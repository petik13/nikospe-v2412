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
    Foam::prescribedShipMotion

Description
    Prescribed planar ship motion and the incident wave seen from the ship.
    Header only; shared by manFlowPrescribed, potRotatingFrameBC and
    linBodyMotionRot so that all three evaluate exactly the same incident
    field.

    Velocity convention is u = +grad(Phi) throughout.

    Axes
    ----
    Mesh (ship-fixed) axes: x aft, y starboard, z up, with the bow at -x, as
    in the existing KVLCC2/SOBC cases at heading 0.  The mesh never moves; the
    ship's motion enters through the frame velocity

        V_S(x) = (-u, v, 0) + Omega x (x - x_c),   Omega = (0, 0, -r)

    which is the velocity of the ship-fixed point x relative to the earth,
    expressed in mesh axes.  u, v, r are MMG quantities at the rotation centre
    x_c (x forward, y starboard, z down; r > 0 turns to starboard).

    Earth axes: X north, Y east (MMG / manModel convention), heading psi
    measured clockwise from north.  The ship reference point moves with

        dX/dt = u cos(psi) - v sin(psi),   dY/dt = u sin(psi) + v cos(psi),
        psi(t) = psi0 + r t

    Incident wave
    -------------
    Linear, finite depth, propagating along the earth direction chi
    (waveDirection, the manModel convention: 0 = north, clockwise; chi =
    psi0 + 180 deg is head seas).  In mesh axes its direction is

        e(t) = (-cos(chi - psi), sin(chi - psi), 0)

    and its phase at the mesh point x is

        theta = k [ e_w.X_s(t) + e(t).(x - x_c) ] + k x_c,x - omega t

    (the constant k x_c,x makes the straight head-sea case reproduce the
    phase k x - omega_e t of manFlow exactly).  With

        phi_I = A g/omega C(z) sin(theta),   zeta_I = A cos(theta),
        C = cosh k(z+h)/cosh kh

    the time derivative at a FIXED MESH POINT, which is what the linearised
    Bernoulli equation in the ship frame needs, is

        dphi_I/dt|ship = A g/omega C cos(theta) dtheta/dt,
        dtheta/dt      = -omega + k e.V_S(x)

    i.e. the local encounter frequency.  The wave amplitude is ramped in over
    rampPeriods absolute wave periods, as in manFlow.

    Reads
    -----
        constant/waveConditions : steepness, waveLength, waterDepth,
                                  rampPeriods
        constant/prescribedMotion:
            type            constant;   // circular motion test (only option)
            u               0.318;      // [m/s]
            v               -0.097;     // [m/s], at the rotation centre
            r               0.08;       // [rad/s]
            psi0            0;          // [deg]
            waveDirection   180;        // [deg], propagation, earth
            rotationCentre  (22.4 0 0); // mesh coordinates

\*---------------------------------------------------------------------------*/

#ifndef Foam_prescribedShipMotion_H
#define Foam_prescribedShipMotion_H

#include "IOdictionary.H"
#include "Time.H"
#include "vector.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"
#include "gravityMeshObject.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class prescribedShipMotion Declaration
\*---------------------------------------------------------------------------*/

class prescribedShipMotion
{
    // Private Data

        //- Wave
        scalar g_, amp_, k_, omega_, depth_, rampTime_;

        //- Wave propagation direction, earth axes [rad]
        scalar chi_;

        //- Motion at the rotation centre, MMG axes
        scalar u_, v_, r_, psi0_;

        //- Rotation centre, mesh axes
        vector xc_;


    // Private Member Functions

        static IOdictionary readDict
        (
            const objectRegistry& db,
            const word& name
        )
        {
            return IOdictionary
            (
                IOobject
                (
                    name,
                    db.time().constant(),
                    db,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );
        }


public:

    // Constructors

        explicit prescribedShipMotion(const objectRegistry& db)
        {
            const IOdictionary waveDict(readDict(db, "waveConditions"));
            const IOdictionary motionDict(readDict(db, "prescribedMotion"));

            const word type(motionDict.getOrDefault<word>("type", "constant"));
            if (type != "constant")
            {
                FatalIOErrorInFunction(motionDict)
                    << "Only type constant (circular motion test) is"
                       " implemented, got " << type << exit(FatalIOError);
            }

            g_ = mag(meshObjects::gravity::New(db.time()).value());

            const scalar steepness  = waveDict.get<scalar>("steepness");
            const scalar waveLength = waveDict.get<scalar>("waveLength");
            depth_ = waveDict.get<scalar>("waterDepth");
            const scalar rampPeriods = waveDict.get<scalar>("rampPeriods");

            amp_   = 0.5*steepness*waveLength;
            k_     = constant::mathematical::twoPi/waveLength;
            omega_ = Foam::sqrt(g_*k_*Foam::tanh(k_*depth_));
            rampTime_ = rampPeriods*constant::mathematical::twoPi/omega_;

            u_    = motionDict.get<scalar>("u");
            v_    = motionDict.get<scalar>("v");
            r_    = motionDict.get<scalar>("r");
            psi0_ = degToRad(motionDict.get<scalar>("psi0"));
            chi_  = degToRad(motionDict.get<scalar>("waveDirection"));
            xc_   = motionDict.get<vector>("rotationCentre");
        }


    // Member Functions

        // Access

            scalar u() const { return u_; }
            scalar v() const { return v_; }
            scalar r() const { return r_; }
            scalar omega() const { return omega_; }
            scalar waveNumber() const { return k_; }
            scalar amplitude() const { return amp_; }
            const vector& rotationCentre() const { return xc_; }

            //- Angular velocity of the ship, mesh axes
            vector Omega() const { return vector(0, 0, -r_); }

            //- Translational velocity of the rotation centre, mesh axes
            vector V0() const { return vector(-u_, v_, 0); }


        // Kinematics

            //- Heading [rad]
            scalar psi(const scalar t) const { return psi0_ + r_*t; }

            //- Earth position (X north, Y east) of the rotation centre,
            //  relative to its position at t = 0
            void position(const scalar t, scalar& X, scalar& Y) const
            {
                const scalar p = psi(t);

                if (mag(r_*t) > 1e-9)
                {
                    X = ( u_*(Foam::sin(p) - Foam::sin(psi0_))
                        + v_*(Foam::cos(p) - Foam::cos(psi0_)))/r_;
                    Y = (-u_*(Foam::cos(p) - Foam::cos(psi0_))
                        + v_*(Foam::sin(p) - Foam::sin(psi0_)))/r_;
                }
                else
                {
                    X = (u_*Foam::cos(psi0_) - v_*Foam::sin(psi0_))*t;
                    Y = (u_*Foam::sin(psi0_) + v_*Foam::cos(psi0_))*t;
                }
            }

            //- Velocity of the ship-fixed point x relative to the earth,
            //  mesh axes
            vector VS(const point& x) const
            {
                const vector d(x - xc_);
                // Omega x d with Omega = (0, 0, -r)
                return V0() + vector(r_*d.y(), -r_*d.x(), 0);
            }

            //- Wave propagation direction in mesh axes
            vector waveDir(const scalar t) const
            {
                const scalar a = chi_ - psi(t);
                return vector(-Foam::cos(a), Foam::sin(a), 0);
            }

            //- Encounter angle, manModel convention [deg]: 0 head sea,
            //  positive waves from starboard, in [-180, 180]
            scalar encounterAngle(const scalar t) const
            {
                scalar mu =
                    std::fmod
                    (
                        radToDeg(chi_ - psi(t)) + 180.0 + 3600.0,
                        360.0
                    );
                if (mu > 180.0) mu -= 360.0;
                return mu;
            }


        // Incident wave

            scalar ramp(const scalar t) const
            {
                if (rampTime_ <= SMALL || t >= rampTime_) return 1.0;
                return
                    0.5*(1.0 - Foam::cos(constant::mathematical::pi*t/rampTime_));
            }

            scalar rampDot(const scalar t) const
            {
                if (rampTime_ <= SMALL || t >= rampTime_) return 0.0;
                return
                    0.5*constant::mathematical::pi/rampTime_
                   *Foam::sin(constant::mathematical::pi*t/rampTime_);
            }

            //- Incident potential, velocity, time derivative at a fixed mesh
            //  point, elevation and horizontal elevation gradient at x
            void incident
            (
                const point& x,
                const scalar t,
                scalar& phiI,
                vector& uI,
                scalar& dphiIdt,
                scalar& zetaI,
                vector& gradZetaI
            ) const
            {
                const vector e(waveDir(t));

                scalar X, Y;
                position(t, X, Y);
                const scalar eX = Foam::cos(chi_)*X + Foam::sin(chi_)*Y;

                const scalar theta =
                    k_*(eX + (e & (x - xc_)) + xc_.x()) - omega_*t;
                const scalar dThetaDt = -omega_ + k_*(e & VS(x));

                const scalar cth = Foam::cos(theta);
                const scalar sth = Foam::sin(theta);

                const scalar chDen = Foam::cosh(k_*depth_);
                const scalar C = Foam::cosh(k_*(depth_ + x.z()))/chDen;
                const scalar S = Foam::sinh(k_*(depth_ + x.z()))/chDen;

                const scalar rmp = ramp(t);
                const scalar P0 = amp_*g_/omega_;       // unramped phi_I amplitude
                const scalar Pa = P0*rmp;

                phiI = Pa*C*sth;
                uI = Pa*k_*vector(C*cth*e.x(), C*cth*e.y(), S*sth);
                dphiIdt = Pa*C*cth*dThetaDt + P0*rampDot(t)*C*sth;
                zetaI = amp_*rmp*cth;
                gradZetaI = -amp_*rmp*k_*sth*e;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
