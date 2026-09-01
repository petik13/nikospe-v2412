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

    which is Chen (2022) Eq.(43) with his inward normal flipped.  The
    second-order steady potential psi-bar and the O(Fr^2) waterline group are
    deliberately omitted.

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
    \endtable

SourceFiles
    middleFieldForm.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_middleFieldForm_H
#define Foam_functionObjects_middleFieldForm_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
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
        //  total is a difference of larger, individually CV-dependent terms
        vector surfaceForce_,   surfaceMoment_;
        vector elevationForce_, elevationMoment_;
        vector stripForce_,     stripMoment_;

        //- Discrete waterline closure, reported once.  The strip term's two
        //  halves cancel for a plane wave only if the contour satisfies
        //      int n dl = 0        (closure -- all the FORCE needs)
        //      int r_x n_x dl = int r_y n_y dl = A_cv   (first moments)
        //      int r_x n_y dl = int r_y n_x dl = 0
        //  so a contour that closes can still give a wrong MOMENT.  A
        //  non-zero (rxnx - ryny) leaves a spurious yaw moment that scales
        //  with omega, W and the CV plan area, and so does not decay in long
        //  waves the way a physical drift moment must.
        mutable bool closureReported_;
        vector clNet_;
        scalar clRxnx_, clRyny_, clRxny_, clRynx_, clLen_;

        //- Output files
        autoPtr<OFstream> forceFilePtr_;
        autoPtr<OFstream> momentFilePtr_;


    // Protected Member Functions

        //- Zero the accumulators
        void reset();

        //- Momentum flux through the control surface
        void surfaceIntegral(const surfaceVectorField& Uf);

        //- Wave-elevation and strip terms along the waterline
        void waterlineIntegral();

        //- File handling
        void createFiles();
        void writeFiles();


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

        //- Total force and moment
        virtual vector force() const;
        virtual vector moment() const;
};


} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
