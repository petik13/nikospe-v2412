/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    Foam::functionObjects::meanWaveLoads

Description
    Variant of the user's meanWaveLoads function object that reconstructs the
    control-volume side surface every time-step from a runtime box indicator
    instead of a faceZone. The box is defined by boxMin and boxMax. The side
    surface is the set of faces across which the indicator changes from inside
    to outside, excluding selected physical patches and the free-surface patch.

    The free-surface line integral is then built from those selected CV faces,
    which is generally more robust in parallel than reconstructing the contour
    from fsPatch.meshEdges() plus a faceZone lookup.

SourceFiles
    meanWaveLoads.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_meanWaveLoads_H
#define Foam_functionObjects_meanWaveLoads_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"
#include "wordRes.H"
#include "HashSet.H"

namespace Foam
{
namespace functionObjects
{

class meanWaveLoads
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Integrated totals
    vector sumPatchForcesP_;
    vector sumPatchMomentsP_;
    vector sumInternalForces_;
    vector sumInternalMoments_;

    // File streams
    autoPtr<OFstream> forceFilePtr_;
    autoPtr<OFstream> momentFilePtr_;

    // Dictionary / configuration data
    autoPtr<coordinateSystem> coordSysPtr_;
    labelList patchIDs_;
    scalar rhoRef_;

    word PhiName_;
    word UName_;
    word UcurName_;
    word rhoName_;
    word zetaName_;
    word freeSurfacePatchName_;

    label freeSurfacePatchID_;
    scalar gMag_;
    bool writeFields_;
    bool initialised_;

    // Runtime box indicator definition
    point boxMin_;
    point boxMax_;

    // Additional physical patches to exclude from the side-surface integral
    wordRes excludePatches_;
    HashSet<label> excludePatchIDs_;

    // Protected Member Functions
    void setCoordinateSystem
    (
        const dictionary& dict,
        const word& e3Name = word::null,
        const word& e1Name = word::null
    );

    volVectorField& force();
    volVectorField& moment();

    void initialise();
    void reset();

    tmp<volScalarField> rho() const;
    tmp<scalarField> rho(const label patchi) const;
    scalar rho(const volScalarField& p) const;

    void addToPatchFields
    (
        const label patchi,
        const vectorField& Md,
        const vectorField& fP
    );

    void addToInternalField
    (
        const labelList& cellIDs,
        const vectorField& Md,
        const vectorField& f
    );

    void createIntegratedDataFiles();

    void writeIntegratedDataFileHeader
    (
        const word& header,
        OFstream& os
    ) const;

    void writeIntegratedDataFiles();

    void writeIntegratedDataFile
    (
        const vector& pres,
        const vector& internal,
        OFstream& os
    ) const;

    void logIntegratedData
    (
        const string& descriptor,
        const vector& pres,
        const vector& internal
    ) const;

    inline bool insideBox(const point& p) const
    {
        return
            p.x() >= boxMin_.x() && p.x() <= boxMax_.x()
         && p.y() >= boxMin_.y() && p.y() <= boxMax_.y()
         && p.z() >= boxMin_.z() && p.z() <= boxMax_.z();
    }

    inline bool excludedPatch(const label patchi) const
    {
        return excludePatchIDs_.found(patchi);
    }

public:

    TypeName("meanWaveLoads");

    meanWaveLoads
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict,
        const bool readFields = true
    );

    meanWaveLoads
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        const bool readFields = true
    );

    meanWaveLoads(const meanWaveLoads&) = delete;
    void operator=(const meanWaveLoads&) = delete;

    virtual ~meanWaveLoads() = default;

    virtual bool read(const dictionary& dict);
    virtual void calcForcesMoments();
    virtual vector forceEff() const;
    virtual vector momentEff() const;
    virtual bool execute();
    virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
