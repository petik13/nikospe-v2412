/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::meanWaveLoadsClaude

Group
    grpForcesFunctionObjects

Description
    Computes meanWaveLoadsClaude and moments over a given list of patches by integrating
    pressure and viscous meanWaveLoadsClaude and moments, and optionally resistance meanWaveLoadsClaude
    and moments from porous zones.

    Forces and moments are output in their total and constituent components:
    - total meanWaveLoadsClaude and moments
    - pressure contributions
    - viscous contributions
    - porous resistance contributions (optional)

    Forces and moments can be computed and output in:
    - the global Cartesian coordinate system (default)
    - a user-defined Cartesian coordinate system

    Operands:
    \table
      Operand       | Type           | Location
      input         | -              | -
      output file   | dat            | postProcessing/\<FO\>/\<time\>/\<file\>s
      output field  | volVectorField | \<time\>/\<outField\>s
    \endtable

    where \c \<file\>s:
    \verbatim
      force.dat        | Forces
      moment.dat       | Moments
    \endverbatim

    where \c \<outField\>s:
    \verbatim
      <namePrefix>:force   | Force field
      <namePrefix>:moment  | Moment field
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <namePrefix>
    {
        // Mandatory entries
        type                meanWaveLoadsClaude;
        libs                (meanWaveLoadsClaude);
        patches             (<wordRes>);

        // Optional entries
        directForceDensity  <bool>;
        porosity            <bool>;
        writeFields         <bool>;
        useNamePrefix       <bool>;

        // Conditional mandatory entries

            // if directForceDensity == true
            fD              <word>;


            // Cartesian coordinate system specification when
            // evaluating meanWaveLoadsClaude and moments, either of the below

            // Define the centre of rotation
            // with implicit directions e1=(1 0 0) and e3=(0 0 1)
            CofR                (0 0 0); // Centre of rotation

            // Define local coordinate system by origin + axes
            origin              (0 0 0);
            e1                  (1 0 0);
            e3                  (0 0 1); // (e1, e2) or (e2, e3) or (e3, e1)

            // General coordinate system specification (always cartesian)
            coordinateSystem
            {
                origin          (0 0 0);
                rotation
                {
                    type        axes;
                    e3          (0 0 1);
                    e1          (1 0 0); // (e1, e2) or (e2, e3) or (e3, e1)
                }
            }

        // Conditional optional entries

            // if directForceDensity == false
            p               <word>;
            U               <word>;
            rho             <word>;
            rhoInf          <scalar>; // enabled if rho=rhoInf
            pRef            <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description               | Type | Reqd    | Deflt
      type       | Type name: meanWaveLoadsClaude         | word | yes     | -
      libs       | Library name: meanWaveLoadsClaude      | word | yes     | -
      patches    | Names of operand patches  | wordRes | yes  | -
      directForceDensity | Flag to directly supply force density <!--
                 -->                         | bool | no      | false
      porosity   | Flag to include porosity contributions | bool | no | false
      writeFields | Flag to write force and moment fields | bool | no | false
      useNamePrefix | Flag to include prefix for field names | bool | no | false
      coordinateSystem | Coordinate system specifier | dictionary | cndtnl | -
      CofR    | Centre of rotation          | vector | cndtnl   | -
      origin  | Origin of coordinate system | vector | cndtnl   | -
      e3      | e3 coordinate axis          | vector | cndtnl   | -
      e1      | e1 coordinate axis          | vector | cndtnl   | -
      fD      | Name of force density field | word   | cndtnl   | -
      p       | Name of pressure field      | word   | cndtnl   | p
      U       | Name of velocity field      | word   | cndtnl   | U
      rho     | Name of density field       | word   | cndtnl   | rho
      rhoInf  | Value of reference density  | scalar | cndtnl   | -
      pRef    | Value of reference pressure | scalar | cndtnl   | 0
    \endtable

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link writeFile.H \endlink
      - \link coordinateSystem.H \endlink

Note
  - For incompressible cases, set \c rho to \c rhoInf.
    You will then be required to provide a \c rhoInf
    value corresponding to the constant freestream density.
  - \c writeControl and \c writeInterval entries of function
    object do control when to output force and moment files and fields.
  - If a \c coordinateSystem entry exists, it is taken in favour of \c CofR.

SourceFiles
    meanWaveLoadsClaude.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_meanWaveLoadsClaude_H
#define Foam_functionObjects_meanWaveLoadsClaude_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "coordinateSystem.H"
#include "Switch.H"
#include "OFstream.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class meanWaveLoadsClaude Declaration
\*---------------------------------------------------------------------------*/

class meanWaveLoadsClaude
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        // Fields

            //- Sum of patch pressure forces
            vector sumPatchForcesP_;

            //- Sum of patch pressure moments
            vector sumPatchMomentsP_;

            //- Sum of internal forces
            vector sumInternalForces_;

            //- Sum of internal moments
            vector sumInternalMoments_;

            vector sumSurfaceForces_;
            vector sumLineForces_;
            vector sumSurfaceMoments_;
            vector sumLineMoments_;

        // File streams

            //- File stream for forces
            autoPtr<OFstream> forceFilePtr_;

            //- File stream for moments
            autoPtr<OFstream> momentFilePtr_;


        // Read from dictionary

            //- Coordinate system used when evaluating meanWaveLoadsClaude and moments
            autoPtr<coordinateSystem> coordSysPtr_;

            //- Selected operand patches
            labelList patchIDs_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Name of Phi field
            word PhiName_;

            //- Name of density field
            word rhoName_;

            //- Flag to write force and moment fields
            bool writeFields_;

            //- Flag of initialisation (internal)
            bool initialised_;

            word faceZoneName_;
            label faceZoneID_;

            // Free-surface elevation and patch
            word zetaName_;
            word freeSurfacePatchName_;
            label freeSurfacePatchID_;
            scalar gMag_;
            point cvPoint_;

            //- Evaluate the incident wave analytically at CV face centres and
            //  at waterline segment midpoints, instead of inheriting it from
            //  the interpolated/face-centre discrete fields.  The incident
            //  field is closed form, and it is the part of the integrand that
            //  must cancel between the surface and waterline terms, so
            //  weighting it inconsistently between the two directly spoils
            //  that cancellation.
            Switch analyticIncident_;

            //- Chen (2022) Eq.(46): the O(Fr^2) waterline group  Fr^2 * F^W_M
            //  of the midfield formulation.  It involves the steady elevation
            //  eta_0 of the double-body flow, which is nil only when the
            //  control surface is at infinity (Chen, section V).  On a CV at a
            //  finite distance from the hull it is not, and the whole group is
            //  O(eps^2 U^2) -- even in U, so it shifts head and following seas
            //  the same way.  Identically zero at U = 0.
            Switch steadyElevationTerms_;

            //- Near-field (direct pressure integration) drift force, Chen
            //  (2022) Eq.(29).  For a FIXED, wall-sided body with psi_bar
            //  dropped it reduces to
            //
            //      F_N = rho * int_hull  1/2 <|u|^2> n dS
            //          - 1/2 rho g * int_WL <zeta^2> n dl        (n out of body)
            //
            //  Note it carries NO explicit forward-speed terms for a fixed
            //  body: all U-dependence enters through the solution.  The
            //  midfield form does carry an explicit O(U) strip term, so
            //  comparing the two isolates whether that term is right.
            Switch nearField_;
            word bodyPatchName_;
            label bodyPatchID_;
            vector sumNearFieldForces_;
            vector sumNearFieldMoments_;
            vector sumNearHullForces_;      // kinetic-energy hull part alone
            vector sumNearWLForces_;        // waterline part alone

            //- Midfield line term split into its U-independent and O(U) parts:
            //    zeta : -1/2 rho g <zeta^2> n        (present at every speed)
            //    strip: -rho <zeta>[(u.n)W + u(W.n)] (Chen Eq.43 "L2", O(U))
            //  The strip term is the ONLY explicit forward-speed term the
            //  midfield has and the near-field does not, so its size bounds how
            //  much of a near-vs-mid gap it could possibly explain.
            vector sumLineZetaForces_;
            vector sumLineStripForces_;
            mutable autoPtr<OFstream> midFilePtr_;
            void writeMidBreakdown();
            mutable autoPtr<OFstream> nearFilePtr_;

            void calcNearField(const volVectorField& U, const scalar tNow);
            void writeNearField();

            //- Diffracted (scattered) fields, used when analyticIncident_
            word UDName_;
            word zetaDName_;

            //- Prescribed incident field, used only for the consistency check
            word UIName_;

            //- Incident wave parameters, read from constant/waveCurConditions
            scalar waveAmp_;
            scalar waveK_;
            scalar waveW_;
            scalar waveOmegaE_;
            scalar waveDepth_;
            scalar waveRampTime_;
            scalar waveU0_;      // |U_inf|, for the O(Fr^2) group
            bool   waveParamsOK_;
            bool   incidentChecked_;

            //- Chen Eq.(46) contribution, accumulated separately so its size
            //  can be seen in the log, then folded into the line term.
            vector sumSteadyElevForces_;
            vector sumSteadyElevMoments_;

            //- Analytic incident wave kinematics (U = -grad(Phi) convention)
            vector incidentU(const point& p, const scalar t) const;
            scalar incidentZeta(const point& p, const scalar t) const;

            //- Read constant/waveCurConditions
            void readWaveConditions();

            //- Compare the analytic incident field against the solver's own
            //  UI field once, so a drift in parameters or formulae is loud
            //  rather than silent.
            void checkIncident();

    // Protected Member Functions

        //- Set the co-ordinate system from dictionary and axes names
        void setCoordinateSystem
        (
            const dictionary& dict,
            const word& e3Name = word::null,
            const word& e1Name = word::null
        );

        //- Return access to the force field
        volVectorField& force();

        //- Return access to the moment field
        volVectorField& moment();

        //- Initialise containers and fields
        void initialise();

        //- Reset containers and fields
        void reset();


    // Evaluation


        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rho if specified otherwise rhoRef for patch
        tmp<scalarField> rho(const label patchi) const;

        //- Return rhoRef if the pressure field is
        //- dynamic (i.e. p/rho), otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- Add patch contributions to force and moment fields
        void addToPatchFields
        (
            const label patchi,
            const vectorField& Md,
            const vectorField& fP
            );

        //- Add cell contributions to force and
        //- moment fields, and include porosity effects
        void addToInternalField
        (
            const labelList& cellIDs,
            const vectorField& Md,
            const vectorField& f
        );


    // I-O

        //- Create the integrated-data files
        void createIntegratedDataFiles();

        //- Write header for an integrated-data file
        void writeIntegratedDataFileHeader
        (
            const word& header,
            OFstream& os
        ) const;

        //- Write integrated data to files
        void writeIntegratedDataFiles();

        //- Write integrated data to a file
        void writeIntegratedDataFile
        (
            const vector& total,
            const vector& surface,
            const vector& line,
            OFstream& os
        ) const;

        //- Write integrated data to stream
        void logIntegratedData
        (
            const string& descriptor,
            const vector& pres,
            const vector& internal
        ) const;


public:

    //- Runtime type information
    TypeName("meanWaveLoadsClaude");


    // Constructors

        //- Construct from Time and dictionary
        meanWaveLoadsClaude
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        meanWaveLoadsClaude
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );

        //- No copy construct
        meanWaveLoadsClaude(const meanWaveLoadsClaude&) = delete;

        //- No copy assignment
        void operator=(const meanWaveLoadsClaude&) = delete;


    //- Destructor
    virtual ~meanWaveLoadsClaude() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate meanWaveLoadsClaude and moments
        virtual void calcForcesMoments();

        //- Return the total force
        virtual vector forceEff() const;

        //- Return the total moment
        virtual vector momentEff() const;

        //- Execute the function object
        virtual bool execute();

        //- Write to data files/fields and to streams
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
