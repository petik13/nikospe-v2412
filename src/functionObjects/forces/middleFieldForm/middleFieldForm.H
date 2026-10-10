/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::functionObjects::middleFieldForm

Group
    grpForcesFunctionObjects

Description
    Mean wave loads by the midfield formulation: momentum flux through a
    control surface drawn around the body, plus the waterline contribution of
    the strip between z = 0 and the instantaneous free surface.

    Velocity convention is u = +grad(Phi).

    With n the OUTWARD normal of the control volume (pointing away from the
    body) and W the steady basis flow,

        F = rho * int_Sc [ 1/2 |u|^2 n - (u.n) u ] dS          (momentum flux)
          - 1/2 rho g * int_C <zeta^2> n dl                    (wave elevation)
          - rho * int_C <zeta [ (u.n) W + u (W.n) ]> dl        (strip flux, O(U))

    which is Chen (2022) Eq.(43) with his inward normal flipped.  The O(Fr^2)
    waterline group is omitted.  The second-order steady potential psi-bar is
    omitted unless psiBar is on (below).

    psiBar (optional, default off).  The psi-bar terms are the cross momentum
    flux between W and grad(psibar),

        F_psi = rho * int_Sc [ (W.grad psibar) n - (W.n) grad psibar
                               - (grad psibar.n) W ] dS,     M_psi = r x (same)

    (Chen's Eq.(43) lacks the (W.n) grad psibar term; it is needed for the
    flux to be conserved).  Without a net source this gives no horizontal force
    (Grue & Palm 1993, eq.15), but it does give a yaw moment, O(U): a Munk-type
    moment of the mean second-order flow in the steady stream.  psibar is solved
    at leading order in U from running means of its boundary forcing,

        d(psibar)/dz = -<zeta d2phi/dz2>              on z = 0   (G&P 10/47)
        d(psibar)/dn = -<[(X.grad)grad(phi)].n>
                       - <(grad(phi) - Xdot).(theta x n)>         on the hull (G&P 12)

    with X = xi + theta x (r - rRef) from the linBodyMotion body state
    (bodyDisp, bodyRot, bodyVel, bodyOmega), psibar = 0 on the far field and
    zero gradient elsewhere (G&P 46-49).  On a distant control surface the
    moment is the psi(2) term of G&P (19); the near field misses only the
    smaller part (92) of it, so SWAN-type pressure integration without
    psibar is closer than the midfield without it.  The result goes to
    psiBar.dat (force.dat and
    moment.dat are unchanged).  Derivation: researchQuestions/psibarMidfield.ipynb.
    Read the yaw moment only, as for the other terms.

    psiBarForcing picks the free-surface forcing: dynamic (above, the
    default) or kinematic,

        d(psibar)/dz = div_h <zeta grad_h(phi)>                on z = 0

    the mean of the second-order kinematic condition without the mean
    elevation terms.  At zero speed the two are equal (<grad zeta.grad phi>
    = -d/dt <|grad phi|^2>/2g = 0); at forward speed they differ by a
    divergence of O(U).  The kinematic form conserves mass by construction:
    its integral over a free-surface region is the Stokes transport through
    the edge of that region.  It is computed as the conservative surface
    divergence of the running mean of zeta u_h on the free-surface faces.

    Mass-balance diagnostics (psiBar.dat, the columns after Q_fs).  With
    S_C = int_C <zeta u.n> dl the Stokes transport out through the control
    surface (its strip), S_WL the same out of the free surface into the hull
    across the waterline and Phi_C = int_Sc d(psibar)/dn dS,

        Phi_C + S_C = 0       no net mean mass flux out of the control volume
        Q_hull + S_WL = 0     no net mean mass flux into the hull

    With the kinematic forcing the first reduces to the second.  The
    Stokes-transport part of the strip term, Fstokes = -rho int_C <zeta u.n> W dl,
    cancels the matching psi-bar part, Fmass = -rho int_Sc (grad psibar.n) W dS,
    when the first balance holds; otherwise their sum is ~ -rho W (Phi_C + S_C).

    psiBarHullForcing picks the form of the hull condition: gradient (G&P 12
    as above, the default) or divergence,

        d(psibar)/dn = div_S <X_n v> - <X.(n x omega)>,   v = grad(phi) - Xdot

    with X_n = X.n, v tangential (first-order body condition) and omega the
    body angular velocity.  The two are equal pointwise for rigid-body motion
    and an exact first-order solution (Laplace and the first-order body
    condition remove the second derivatives; checked on a translating and
    rotating sphere), but the divergence form needs no second derivatives at
    the wall.  It is computed as the conservative edge divergence over the
    hull faces, so its net flux is exactly that through the waterline plus
    the rotation term (Q_hull_waterline, Q_hull_rotation in psiBar.dat).

    The time average is not taken here: instantaneous values are written every
    step and averaged in post-processing, which lets the mean be extracted by a
    harmonic fit rather than by averaging between zero crossings.

    The waterline is walked from the control-volume faces rather than from
    free-surface edges.  Each control-volume face takes the segment it shares
    with the free surface with weight 1/2 from the cell on either side.  Every
    cell belongs to exactly one rank, so an internal face collects 1/2 + 1/2
    locally and a face on a processor patch collects 1/2 from each rank: exact
    for any decomposition, with no communication and no double counting.  It
    also makes zeta a centred estimate at the segment rather than a one-sided
    one, and accumulates <zeta^2> rather than <zeta>^2.

Usage
    \table
      Property         | Description                        | Reqd | Default
      faceZone         | Control-surface face zone          | yes  |
      cvPoint          | A point inside the control volume  | yes  |
      freeSurfacePatch | Free-surface patch name            | yes  |
      U                | Total unsteady velocity field      | no   | U
      Us               | Steady basis velocity field        | no   | Us
      zeta             | Total free-surface elevation       | no   | zeta
      rhoInf           | Reference density                  | yes  |
      CofR             | Moment reference point             | yes  |
      psiBar           | Add the psi-bar terms (psiBar.dat) | no   | false
      psiBarHullPatches | Hull patches                      | no   | (sphere)
      psiBarFarFieldPatches | psibar = 0 here               | no   | (inlet outlet front back)
      psiBarPeriods    | Average over the last N encounter periods before endTime | no | 6
      psiBarStartTime  | Start of the average (overrides psiBarPeriods) | no |
      psiBarRadius     | Free-surface forcing only within this horizontal distance of CofR (<= 0: all) | no | -1
      psiBarRotationCentre | Rotation centre of the body motion | no | (xG 0 0) of bodyMotionProperties
      psiBarSolver     | fvSolution solver entry            | no   | PhiS
      psiBarCorrectors | Non-orthogonal correctors          | no   | 5
      psiBarSplit      | Also solve free-surface and hull forcing apart | no | true
      psiBarForcing    | Free-surface forcing: dynamic or kinematic | no | dynamic
      psiBarHullForcing | Hull forcing: gradient or divergence | no | gradient
    \endtable

    psibar is solved at every field-write time and at the end of the run.

SourceFiles
    middleFieldForm.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_middleFieldForm_H
#define Foam_functionObjects_middleFieldForm_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "wordRes.H"
#include "surfaceFieldsFwd.H"
#include "point.H"
#include "OFstream.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class middleFieldForm Declaration
\*---------------------------------------------------------------------------*/

class middleFieldForm
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Field names
        word UName_;
        word UsName_;
        word zetaName_;

        //- Disturbance velocity and its potential.  The momentum flux is
        //  dominated by (u.n)u, and on an internal face the normal derivative
        //  of the DISTURBANCE potential is far better resolved by the compact
        //  two-point snGrad -- which carries the non-orthogonal correction --
        //  than by interpolating a cell-centred gradient to the face, which
        //  stacks two uncorrected discretisations.  The incident part is
        //  analytic and is left alone.
        word UDName_;
        word PhiDName_;

        //- Take the disturbance normal component from snGrad(PhiD) rather
        //  than from the interpolated cell gradient
        Switch snGradNormal_;

        //- Control surface and the point used to orient its normals outward
        word faceZoneName_;
        label faceZoneID_;
        point cvPoint_;

        //- Free surface
        word freeSurfacePatchName_;
        label freeSurfacePatchID_;

        //- Reference density and gravity
        scalar rhoRef_;
        scalar gMag_;

        //- Moment reference point
        point CofR_;

        //- Contributions, kept apart so the balance can be inspected: the
        //  total is a difference of larger, individually CV-dependent terms.
        //
        //  WARNING: only the YAW component of the moment is complete.  The
        //  angular-momentum balance carries a waterline term
        //
        //      k x ( 1/2 <zeta^2> q_0 )
        //
        //  from the position vector varying up the strip, r = r_h + z k.  It
        //  is purely horizontal, so it drops out of yaw but is a real part of
        //  roll and pitch -- as is the buoyancy moment, which likewise only
        //  disappeared from yaw because it has no vertical component.  Neither
        //  is computed here.  Read moment().z() only.
        vector surfaceForce_,   surfaceMoment_;
        vector elevationForce_, elevationMoment_;
        vector stripForce_,     stripMoment_;

        //- Output files
        autoPtr<OFstream> forceFilePtr_;
        autoPtr<OFstream> momentFilePtr_;

        //- psi-bar terms (psiBar on); see the class description
        Switch psiBar_;
        wordRes psiBarHullPatches_;
        wordRes psiBarFarFieldPatches_;
        labelList psiBarHullIDs_;
        labelList psiBarFarFieldIDs_;
        word psiBarSolver_;
        label psiBarNCorr_;
        scalar psiBarStartTime_;
        scalar psiBarRadius_;
        Switch psiBarSplit_;
        point psiBarRotRef_;
        scalar omegaE_;

        //- Running time integrals of the boundary forcing: hull faces (all
        //  hull patches in order) and free-surface faces
        scalarField hullForcing_;
        scalarField fsForcing_;
        scalar psiTAcc_;
        scalar psiTPrev_;
        bool psiStarted_;
        scalar psiLastSolved_;

        //- Free-surface forcing: dynamic (-<zeta phi_zz>, G&P 47) or
        //  kinematic (div_h <zeta grad_h phi>, mass conserving)
        word psiBarForcing_;

        //- Mass-balance diagnostics: running time integral of the Stokes
        //  transport zeta u_h on each free-surface face (all faces, no radius
        //  cut), and the Stokes transport through the control surface (the
        //  strip) with its momentum -rho zeta (u.n) W: instantaneous (this
        //  step, reduced) and running time integral
        vectorField fsQ_;
        scalar stokesInst_;
        vector stokesForceInst_;
        vector stokesMomentInst_;
        scalar stokesSum_;
        vector stokesForceSum_;
        vector stokesMomentSum_;

        //- Hull forcing: gradient (G&P 12 as written) or divergence
        //  (div_S <X_n v> - <X.(n x omega)>, the same condition without
        //  second derivatives)
        word psiBarHullForcing_;

        //- Running time integrals for the divergence form, per hull face (all
        //  hull patches in order): X_n v_t and -X.(n x omega)
        vectorField hullP_;
        scalarField hullR_;

        autoPtr<volScalarField> psiBarPtr_;
        autoPtr<OFstream> psiFilePtr_;


    // Protected Member Functions

        //- Zero the accumulators
        void reset();

        //- Momentum flux through the control surface
        void surfaceIntegral
        (
            const surfaceVectorField& Uf,
            const surfaceVectorField& UDf,
            const surfaceScalarField& snGradPhiD
        );

        //- Wave-elevation and strip terms along the waterline
        void waterlineIntegral();

        //- File handling
        void createFiles();
        void writeFiles();

        //- psi-bar: read the options, add this step to the running means,
        //  solve psibar (free-surface and/or hull forcing) and integrate the
        //  cross flux over the control surface, write psiBar.dat
        void readPsiBar(const dictionary& dict);
        void accumulatePsiBar();

        //- Kinematic free-surface forcing div_h <zeta u_h> per free-surface
        //  face (radius cut applied); sWL returns the Stokes transport across
        //  the waterline into the hull (reduced)
        tmp<scalarField> fsKinematicForcing(scalar& sWL) const;

        //- Divergence-form hull forcing per hull face (all hull patches in
        //  order); qWL, qRot return its flux out through the waterline and
        //  the integral of the rotation term (reduced)
        tmp<scalarField> hullDivergenceForcing(scalar& qWL, scalar& qRot) const;

        //- fsGrad, hullGrad: d(psibar)/dn per free-surface / hull face (used
        //  if useFs / useHull).  phiC, Fmass, Mmass: flux of psibar out
        //  through the control surface and its momentum -rho (grad psibar.n) W
        void solvePsiBar
        (
            const scalarField& fsGrad,
            const scalarField& hullGrad,
            const bool useFs,
            const bool useHull,
            vector& F,
            vector& M,
            scalar& qHull,
            scalar& qFs,
            scalar& phiC,
            vector& Fmass,
            vector& Mmass
        );
        void psiBarOutput();


public:

    //- Runtime type information
    TypeName("middleFieldForm");


    // Constructors

        middleFieldForm
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        middleFieldForm(const middleFieldForm&) = delete;

        //- No copy assignment
        void operator=(const middleFieldForm&) = delete;


    //- Destructor
    virtual ~middleFieldForm() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary&);

        //- Evaluate the loads
        virtual bool execute();

        //- Write the loads
        virtual bool write();

        //- End of the run: last psi-bar solve if the final time was not a
        //  write time
        virtual bool end();

        //- Total force and moment
        virtual vector force() const;
        virtual vector moment() const;
};


} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
