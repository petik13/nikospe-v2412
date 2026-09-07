/*---------------------------------------------------------------------------*\
Class
    Foam::functionObjects::nearFieldForm

Description
    Mean second-order wave loads by pressure integration over the hull -- the
    NEARFIELD formulation, following Chen (2022) "Nearfield, midfield and
    farfield formulations of added resistances", JH2022, eq (25) and (27).

    This is an independent check on middleFieldForm.  Both of the objects in
    use so far (middleFieldForm and the original meanWaveLoads) are midfield,
    so their agreement establishes consistency, not correctness; this one
    takes a completely different route to the same quantity.

    Chen orients n INTO the fluid and writes F = -int P n ds.  OpenFOAM's
    patch normal points out of the domain, i.e. out of the fluid, so with
    nu = -n every term below carries the opposite sign to the paper and
    F = +int P nu ds.

    Keeping the classical component (Chen eq 27) and dropping psibar and the
    O(Fr^2) steady-flow blocks F_N^H and F_N^W (eq 30, 31):

      F = < int_H [ -rho/2 |u|^2 nu                         quadratic
                    + rho (S.grad p) nu                     transfer
                    + rho (p - g S_z) (theta ^ nu) ] dS     rotation
        + rho g/2 < oint_G (zeta - S_z)^2 nu / cos(a) dl >  waterline

    with u = grad(phi_T) the first-order velocity, p the first-order KINEMATIC
    pressure the solver already forms, S = xi + theta ^ r the hull
    displacement, S_z its vertical part, zeta the wave elevation at the
    waterline, and cos(a) = sqrt(1 - nu_z^2) the hull inclination there
    (Chen eq 13).

    The moment carries the same integrand crossed with the lever arm, plus the
    term from the lever arm itself moving,

      M = < int_H r ^ f dS > + rho < int_H p (S ^ nu) dS > + waterline

    Everything is written instantaneously, so the mean is taken afterwards by
    the same harmonic fit at 2*omega_e that meanLoads.py applies to the
    midfield output.

    NOTE  The body displacement and rotation are read from the registry, where
    linBodyMotion publishes them as bodyDisp and bodyRot.  The rotation centre
    is CofR, which must match the motion solver's reference point (xG, 0, 0).

Usage
    \table
      Property   | Description                              | Reqd | Default
      patches    | Hull patches                             | yes  |
      freeSurfacePatch | Patch at z = 0, for the waterline  | yes  |
      CofR       | Moment reference and rotation centre     | yes  |
      rhoInf     | Reference density                        | yes  |
      gMag       | Gravity magnitude                        | no   | 9.81
      U          | First-order velocity                     | no   | U
      p          | First-order kinematic pressure           | no   | p
      zeta       | Wave elevation                           | no   | zeta
    \endtable

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_nearFieldForm_H
#define Foam_functionObjects_nearFieldForm_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFields.H"
#include "point.H"
#include "OFstream.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class nearFieldForm Declaration
\*---------------------------------------------------------------------------*/

class nearFieldForm
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        word UName_;
        word pName_;
        word zetaName_;
        word PhiName_;

        //- Take the TANGENTIAL hull velocity from the surface gradient of the
        //  potential instead of from the boundary value of U.
        //
        //  On the hull u.nu is prescribed exactly by the body condition (zero
        //  for a restrained body), so |u|^2 there is almost entirely
        //  tangential -- and U's tangential boundary components are the
        //  zero-order extrapolation of the cell value, i.e. first order.  The
        //  tangential velocity on a surface IS the surface gradient of the
        //  potential, so differencing phi along the patch gives it without any
        //  extrapolation.  The normal part is kept from U, where it is exact.
        //
        //  MEASURED WORSE, DEFAULT OFF.  On a snapped hex mesh PhiD's hull
        //  boundary VALUE is phi_cell + g*d, and g = -u_I.nu carries the
        //  faceting of the snapped normals.  Differencing that along the
        //  surface injects a spurious tangential gradient of order u_I,
        //  i.e. O(1): |u_t| came out ~1.6x too large and the quadratic term
        //  ~2.6x.  The same faceting defeated a curvature-based m-term.
        Switch surfaceTangential_;

        //- Edge-connected neighbours per hull patch face, and the inverse
        //  least-squares moment matrix.  Geometry only, so built once.
        mutable List<labelList> lsNbr_;
        mutable List<tensor> lsInvA_;
        mutable List<label> lsPatchStart_;
        mutable bool lsValid_;

        labelList hullPatchIDs_;

        word freeSurfacePatchName_;
        label freeSurfacePatchID_;

        //- Which of Chen's two algebraically equivalent forms to evaluate:
        //
        //  eq27  int[ -rho/2|u|^2 nu + rho(S.grad p)nu
        //             + rho(p - g S_z)(theta^nu) ] + rho g/2 oint (zeta-S_z)^2 nu
        //
        //  eq29  int[ -rho/2|u|^2 nu + rho(S.nu) grad(p) ]
        //             + rho g/2 oint zeta^2 nu
        //
        //  They are identical on paper (Chen 2007, eq 57), so any difference
        //  between them is pure discretisation error -- a self-check needing
        //  no external reference.  eq29 should be much better conditioned
        //  here: it contracts the displacement with the NORMAL only, which is
        //  the component the body condition prescribes, whereas eq27 needs
        //  S.grad(p) whose tangential part comes from the extrapolated rows of
        //  the boundary gradient.  It also drops the rotation term and needs
        //  no relative elevation at the waterline.
        //
        //  For the MOMENT, eq27 carries the lever-arm term rho p (S ^ nu);
        //  Chen gives the force only, so the eq29 moment is taken as r ^ f,
        //  the natural extension rather than a derived result.
        word formulation_;

        point CofR_;
        scalar rhoRef_;
        scalar gMag_;

        //- Contributions, kept apart: quadratic velocity, the terms linear in
        //  the motion, and the waterline
        vector quadForce_,  quadMoment_;
        vector motForce_,   motMoment_;
        vector wlForce_,    wlMoment_;

        autoPtr<OFstream> forceFilePtr_;
        autoPtr<OFstream> momentFilePtr_;


    // Protected Member Functions

        //- bodyDisp / bodyRot from the registry, zero if absent
        vector bodyVector(const word&) const;

        //- Build the surface least-squares stencil on the hull patches
        void buildSurfaceStencil() const;

        void reset();

        //- Hull integral: quadratic and motion terms
        void hullIntegral();

        //- Waterline integral around the hull
        void waterlineIntegral();

        void createFiles();
        void writeFiles();


public:

    TypeName("nearFieldForm");

    nearFieldForm(const word&, const Time&, const dictionary&);

    virtual ~nearFieldForm() = default;


    // Member Functions

        virtual bool read(const dictionary&);
        virtual bool execute();
        virtual bool write();

        vector force() const;
        vector moment() const;
};


} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
