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
    Foam::linBodyMotionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    Body boundary condition for the disturbance potential of a linearised,
    freely floating body, coupled to a 6-DOF Newmark integrator.

    Velocity convention is u = +grad(Phi).

    Expanding the exact no-flux condition (grad(Phi) - V_b).N = 0 about the
    mean wetted surface, with S the displacement and theta the rotation,

        grad(phi).n = Sdot.n - n.(S.grad)W - W.(theta x n)

    so that, splitting phi = phi_I + phi_D and solving for the disturbance,

        dphi_D/dn = V_b.n - dphi_I/dn - S.[(n.grad)W] - W.(theta x n)
                    ^^^^^   ^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    motion  incident    m-terms: interaction of the body
                                        motion with the steady basis flow

    The two m-terms vanish for a restrained body and at zero forward speed.
    W.(theta x n) is the change of the normal direction under rotation; it is
    the part of m4..m6 that a translation-only form omits.

    The 6-DOF system

        (M + A) qddot + C qdot + K q = F

    is advanced with Newmark, the loads taken from the pressure on the mean
    wetted surface.  A is an artificial added mass stabilising the explicit
    fluid-structure coupling: see addedMass_ below.

    The force is an integral over the instantaneous surface too, so it needs
    the same transfer onto the mean hull that produced the m-terms.  Expanding
    it leaves three contributions linear in the motion, beyond the unsteady
    pressure over the mean hull,

        (S.grad)p_S n      transfer of the steady pressure onto the
                           displaced hull
        p_S (theta x n)    rotation of the normal in the steady pressure
        p_S (S x n)        movement of the lever arm, moment only

    The hydrostatic part of these is the usual C33..C55 restoring matrix.  The
    part carried by the steady DYNAMIC pressure p_S has no counterpart there,
    and is assembled here as a constant 6x6 added into K: see
    steadyRestoringMatrix().  It contains the Munk moment, so it carries
    negative diagonal entries in pitch and yaw.

    Note that p_S vanishes identically for a uniform basis flow, so these
    terms are zero in a Neumann-Kelvin treatment and only appear once W is a
    genuine double-body flow.

Usage
    \table
      Property   | Description                             | Reqd | Default
      PhiI       | Incident potential field                | no   | PhiI
      Us         | Steady basis velocity field             | no   | Us
      gradUs     | Steady basis velocity gradient field    | no   | gradUs
      pS         | Steady dynamic pressure field           | no   | pS
      p          | Pressure field used for the loads       | no   | p
    \endtable

    Body data are read from constant/bodyMotionProperties:
    \table
      Property        | Description                        | Reqd | Default
      addedMass       | Artificial added mass, per DOF, as a fraction of the
                        diagonal inertia; a single scalar applies to all six
                                                           | no   | 0
      addedMassFromTimeStep | Drive A from the previous time step (stable) or
                        the previous corrector (consistent) | no  | true
      accelRelaxation | Under-relaxation on the new acceleration | no | 1
      accelDamping    | Decay on the new acceleration      | no   | 1
      steadyRestoring | Include the forward-speed restoring carried by the
                        steady dynamic pressure (Munk moment and friends).
                        Set false to recover the hydrostatic-only behaviour
                                                           | no   | true
    \endtable

    Body data are read from constant/bodyMotionProperties, wave data from
    constant/waveConditions.

SourceFiles
    linBodyMotionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_linBodyMotionFvPatchScalarField_H
#define Foam_linBodyMotionFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "FixedList.H"
#include "OFstream.H"
#include "Switch.H"
#include "scalarList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class linBodyMotionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class linBodyMotionFvPatchScalarField
:
    public fixedGradientFvPatchField<scalar>
{
    // Private Types

        //- Generalised displacement/velocity/acceleration, ordered
        //  (surge, sway, heave, roll, pitch, yaw)
        typedef FixedList<scalar, 6> vector6;
        typedef FixedList<vector6, 6> tensor6;


    // Private Data

        //- Field names
        word PhiIName_;
        word UsName_;
        word gradUsName_;
        word pSName_;
        word pName_;

        //- Hold the body at its mean position: the condition then reduces to
        //  the pure diffraction problem
        Switch fixBody_;

        //- Mass, restoring and damping, in body axes
        scalar mass_, xG_, zG_;
        scalar I44_, I55_, I66_, I46_;
        scalar C33_, C44_, C55_, C35_;
        scalar rollDampingFraction_;

        //- Soft mooring: horizontal springs restraining the DOF that have no
        //  hydrostatic restoring (surge, sway, yaw).  Chosen so the natural
        //  period of planar motion is well above the highest wave period, so
        //  they hold the mean position without touching the wave-frequency
        //  response.  Zero leaves the DOF free.
        scalar mooringK11_, mooringK22_, mooringK66_;

        //- Damping on the moored DOF, as a fraction of critical.  Bleeds off
        //  the slow drift oscillation the springs would otherwise sustain.
        scalar mooringDampingFraction_;
        scalar heading_;

        //- Artificial added mass for the explicit fluid-structure coupling,
        //  per DOF, as a fraction of the corresponding diagonal inertia.
        //
        //  The body condition sets dphi_D/dn = V_b.n, so the pressure
        //  -dphi/dt responds to the body ACCELERATION: the coupling is an
        //  added-mass one, evaluated lagged.  It goes unstable once the fluid
        //  added mass exceeds the body mass, which in heave it does.  Putting
        //  A on the left and A*a_old on the right suppresses that.
        FixedList<scalar, 6> addedMass_;

        //- Which acceleration drives the added-mass term.  Driving it from the
        //  previous TIME STEP is what stabilises the coupling, at the price of
        //  an O(dt) consistency error.  Driving it from the previous CORRECTOR
        //  is consistent -- A cancels once the correctors converge -- but then
        //  it provides no stabilisation, so it is only safe with enough
        //  correctors for the coupling to actually converge.
        Switch addedMassFromTimeStep_;

        //- Under-relaxation and decay applied to the new acceleration
        scalar accelRelaxation_;
        scalar accelDamping_;

        //- Newmark parameters
        scalar gamma_, beta_;

        //- Incident wave, from constant/waveConditions
        scalar waveAmp_, waveNumber_, omega_, omegaE_, depth_, rampTime_;
        vector Uinf_;

        //- Reference density for the loads
        scalar rho_;

        //- Include the restoring carried by the steady dynamic pressure.
        //  The hydrostatic restoring above is only half of the first-order
        //  pressure transfer; this is the other half, and it is what makes
        //  the forward-speed formulation consistent with the m-terms in the
        //  boundary condition.  Zero for a uniform basis flow either way.
        Switch steadyRestoring_;

        //- That restoring, as a 6x6 added into K.  W, p_S and the mean hull
        //  are all steady, so it is constant: built once on the first step,
        //  after the solver has finished the basis flow.
        tensor6 Ks_;
        bool KsValid_;

        //- State: displacement, velocity, acceleration at n and n+1
        vector6 q_, qDot_, qDdot_;
        vector6 qNew_, qDotNew_, qDdotNew_;
        vector6 qDdotPrevIter_;

        //- Time index of the last state shift
        label lastTimeIndex_;

        //- Motion history file
        mutable autoPtr<OFstream> motionFilePtr_;


    // Private Member Functions

        //- Read constant/bodyMotionProperties and constant/waveConditions
        void readProperties();

        //- Ramp applied to the incident wave
        scalar ramp(const scalar t) const;

        //- Incident velocity grad(phi_I) at the given points
        tmp<vectorField> incidentVelocity
        (
            const vectorField& pts,
            const scalar t
        ) const;

        //- Pressure load on the mean wetted surface, as (force, moment)
        void bodyLoads(vector& force, vector& moment) const;

        //- Restoring carried by the steady dynamic pressure, in body axes.
        //  Column k is the generalised load from a unit motion in mode k,
        //  negated so that it adds directly into K.
        void steadyRestoringMatrix(tensor6& Ks) const;

        //- Assemble mass, damping and restoring in body axes
        void assembleSystem(tensor6& M, tensor6& C, tensor6& K) const;

        //- Advance the 6-DOF system one Newmark step
        void solveMotion();

        //- Rotate a (linear, angular) pair between body and global axes
        void rotate(vector& lin, vector& ang, const scalar angle) const;

        //- Append the current state to the motion file
        void writeMotion();


public:

    //- Runtime type information
    TypeName("linBodyMotion");


    // Constructors

        linBodyMotionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        linBodyMotionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        linBodyMotionFvPatchScalarField
        (
            const linBodyMotionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        linBodyMotionFvPatchScalarField
        (
            const linBodyMotionFvPatchScalarField&
        );

        linBodyMotionFvPatchScalarField
        (
            const linBodyMotionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new linBodyMotionFvPatchScalarField(*this)
            );
        }

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new linBodyMotionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Advance the body and impose the resulting normal gradient
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif

// ************************************************************************* //
