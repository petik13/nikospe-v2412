/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::functionObjects::middleFieldFormRot

Group
    grpForcesFunctionObjects

Description
    Mean wave loads by the midfield formulation for a control volume that
    ROTATES with the ship (turning ship, ship-fixed frame).  Extension of
    middleFieldForm; derivation in researchQuestions/
    longWavelengthDiscrepancies.ipynb, "Midfield mean wave loads in a
    rotating ship frame".

    Velocity convention is u = +grad(Phi).  W = Us is the steady flow
    RELATIVE to the ship (with rotation it is not irrotational).

    Balance of the RELATIVE momentum P = int_V rho w dV in the frame rotating
    with Omega (constant, along z), second-order mean:

        F = F_Chen[W] - dP/dt - 2 Omega x P - Omega x int_GH m_c V_S dl
        Mz = Mz_Chen[W] - dH/dt - 2 Omega Q - Omega int_GH m_c r_h.V_S dl

    F_Chen is middleFieldForm unchanged (Chen 2022 eq 43, W from Us), and

        P = rho int_AF zeta grad_h(phi) dA  +  rho int_SB (xi.n) grad(phi) dS
          + int_GH m_c W_h dl
        H = rho int_AF zeta (r x grad(phi))_z dA + rho int_SB (xi.n)(r x grad(phi))_z dS
          + int_GH m_c (r x W)_z dl
        Q = rho int_AF zeta r_h.grad_h(phi) dA + rho int_SB (xi.n) r_h.grad_h(phi) dS
          + int_GH m_c r_h.W_h dl

    AF: free-surface faces inside the control volume; SB: hull, n out of the
    fluid (into the body), xi = bodyDisp + bodyRot x r the hull displacement
    about CofR; r = x - CofR.  CofR MUST be the rotation centre of the frame.

    Waterline corner (cornerTerm, default on): the fluid between z = 0 and
    the free surface that the hull leaves or takes, mass per unit waterline
    length m_c = rho zeta (xi.n_h), the hull taken as wall-sided at the
    waterline (n_h the horizontal unit normal).  GH is the waterline: the
    edges shared by a hull face and a free-surface face of the same cell.
    zeta is the adjacent free-surface face value; W and V_S are averaged
    over the hull and free-surface faces at the edge.  m_c enters P, H, Q
    (and so the Coriolis and storage terms) and the centripetal term
    -Omega x int m_c V_S dl, the corner part of
    rho oint 1/2|V_S|^2 n dS = -rho Omega x int_D V_S dV.

    For Omega = 0 in a periodic state the two extra terms average to zero
    and the result is middleFieldForm.  For the incident wave alone they
    cancel the strip term exactly (the strip then contains +Omega x M_I).

    Everything is written instantaneously; the storage term is the backward
    difference -(P^n - P^{n-1})/dt, which telescopes, so the mean of the
    written total over one encounter period is the corrected mean load.

    Neglected (as in middleFieldForm, plus): the mean second-order potential
    psibar, the mean second-order elevation zeta2bar, the O(xi^2) hull-layer
    terms (xi2bar, 1/2 (xi.n)^2 dW/dn, curvature) in P and in the
    centripetal term, and dOmega/dt, dV0/dt (zero for a circular motion
    test).  Derivation: researchQuestions/rotatingMidfield.ipynb.

    Omega is taken from the registered uniformDimensionedVectorField
    "frameOmega" (published by manFlowPrescribed), otherwise from the entry
    Omega (default zero).

    The free-surface faces inside the control volume are those whose centres
    lie inside the horizontal extent of the faceZone (a box control volume).

Usage
    \table
      Property         | Description                        | Reqd | Default
      faceZone         | Control-surface face zone          | yes  |
      cvPoint          | A point inside the control volume  | yes  |
      freeSurfacePatch | Free-surface patch name            | yes  |
      hullPatch        | Hull patch name                    | no   | sphere
      U                | Total unsteady velocity field      | no   | U
      Us               | Steady relative basis velocity     | no   | Us
      VS               | Earth velocity of ship points V_S  | no   | VS
      cornerTerm       | Include the waterline corner       | no   | true
      zeta             | Total free-surface elevation       | no   | zeta
      UI               | Incident velocity (diagnostics)    | no   | UI
      zetaI            | Incident elevation (diagnostics)   | no   | zetaI
      Omega            | Frame angular velocity if no frameOmega | no | (0 0 0)
      rhoInf           | Reference density                  | yes  |
      CofR             | Moment reference = rotation centre | yes  |
    \endtable

    Output (postProcessing/<name>/<time>/):
        force.dat   Time total chen surface elevation strip coriolis storage
                    centripetal
        moment.dat  same; only the z components are complete
        momentum.dat Time P P_incidentSelf Hz Q Omega_z Pc Hc Qc
                     (P from the free surface, the hull and the corner;
                     P_incidentSelf = rho int_AF zetaI uI_h dA; Pc, Hc, Qc
                     the corner parts, already included in P, Hz, Q)

        total = chen + coriolis + storage + centripetal.  The corner enters
        coriolis and storage through P, H, Q; centripetal is corner only.

SourceFiles
    middleFieldFormRot.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_middleFieldFormRot_H
#define Foam_functionObjects_middleFieldFormRot_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "point.H"
#include "boundBox.H"
#include "OFstream.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class middleFieldFormRot Declaration
\*---------------------------------------------------------------------------*/

class middleFieldFormRot
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Field names
        word UName_;
        word UsName_;
        word VSName_;
        word zetaName_;
        word UDName_;
        word PhiDName_;
        word UIName_;
        word zetaIName_;

        //- Take the disturbance normal component from snGrad(PhiD)
        Switch snGradNormal_;

        //- Include the waterline corner (in P, H, Q and the centripetal term)
        Switch cornerTerm_;

        //- Waterline found and reported
        bool cornerReported_;

        //- Control surface, the point used to orient its normals outward
        word faceZoneName_;
        label faceZoneID_;
        point cvPoint_;

        //- Free surface and hull
        word freeSurfacePatchName_;
        label freeSurfacePatchID_;
        word hullPatchName_;
        label hullPatchID_;

        //- Reference density and gravity
        scalar rhoRef_;
        scalar gMag_;

        //- Moment reference point = rotation centre of the frame
        point CofR_;

        //- Frame angular velocity when no frameOmega is registered
        vector OmegaDict_;

        //- Horizontal extent of the control volume, and the free-surface
        //  faces inside it (local patch indices), built once
        boundBox cvBox_;
        labelList fsInside_;
        bool fsInsideBuilt_;

        //- Chen terms, as in middleFieldForm (only yaw of the moments valid)
        vector surfaceForce_,   surfaceMoment_;
        vector elevationForce_, elevationMoment_;
        vector stripForce_,     stripMoment_;

        //- Rotating-frame terms
        vector coriolisForce_,  coriolisMoment_;
        vector storageForce_,   storageMoment_;
        vector centripetalForce_, centripetalMoment_;

        //- Relative momentum in the control volume (second order):
        //  P, its incident self-part, the yaw angular momentum H and Q
        vector P_, PI_;
        scalar Hz_, Q_;

        //- Waterline-corner parts of P, H, Q (included in P_, Hz_, Q_)
        vector Pc_;
        scalar Hc_, Qc_;

        //- Previous step, for the storage term
        vector Pold_;
        scalar HzOld_;
        scalar tOld_;
        label lastTimeIndex_;
        bool haveOld_;

        //- Output files
        autoPtr<OFstream> forceFilePtr_;
        autoPtr<OFstream> momentFilePtr_;
        autoPtr<OFstream> momentumFilePtr_;


    // Protected Member Functions

        void reset();

        //- Frame angular velocity (frameOmega if registered, else OmegaDict_)
        vector Omega() const;

        //- Horizontal extent of the faceZone and the free-surface faces
        //  inside it
        void buildInside();

        //- Momentum flux through the control surface (as middleFieldForm)
        void surfaceIntegral
        (
            const surfaceVectorField& Uf,
            const surfaceVectorField& UDf,
            const surfaceScalarField& snGradPhiD
        );

        //- Wave-elevation and strip terms along the waterline (as
        //  middleFieldForm)
        void waterlineIntegral();

        //- Relative momentum, angular momentum and Q in the control volume
        void momentumIntegrals();

        //- Waterline corner: its parts of P, H, Q and the centripetal term
        //  (local sums, before the parallel reduction)
        void cornerIntegrals(const vector& Om);

        void createFiles();
        void writeFiles();


public:

    //- Runtime type information
    TypeName("middleFieldFormRot");


    // Constructors

        middleFieldFormRot
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        middleFieldFormRot(const middleFieldFormRot&) = delete;
        void operator=(const middleFieldFormRot&) = delete;


    //- Destructor
    virtual ~middleFieldFormRot() = default;


    // Member Functions

        virtual bool read(const dictionary&);
        virtual bool execute();
        virtual bool write();

        //- Chen (middleFieldForm) part
        vector chenForce() const;
        vector chenMoment() const;

        //- Corrected total
        vector force() const;
        vector moment() const;
};


} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
