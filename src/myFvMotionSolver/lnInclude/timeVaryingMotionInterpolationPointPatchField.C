/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "timeVaryingMotionInterpolationPointPatchField.H"
#include "Time.H"
#include "rawIOField.H"
#include "matchPoints.H"
#include "pointToPointPlanarInterpolation.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingMotionInterpolationPointPatchField<Type>::
timeVaryingMotionInterpolationPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    fieldTableName_(iF.name()),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    endSampleTime_(-1),
    endSampledValues_(0),
    useCustomFolder_(false),
    intOutsideBounds_(true),
    inverseDistRadius_(0),
    startSamplePoints_(Zero),
    startSampleData_(Zero),
    endSamplePoints_(Zero),
    endSampleData_(Zero),
    domainRefPt_(Zero),
    domainMatDm_(Zero),
    domainVoxSz_(Zero)
{}


template<class Type>
Foam::timeVaryingMotionInterpolationPointPatchField<Type>::
timeVaryingMotionInterpolationPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, false),
    fieldTableName_(iF.name()),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    endSampleTime_(-1),
    endSampledValues_(0),
    inputType_
    (
        dict.getOrDefault<word>
        (
            "inputType",
            "unstructured"
        )
    ),
    interpolationType_
    (
        dict.getOrDefault<word>
        (
            "interpolationType",
            "nearest"
        )
    ),
    inputFolderName_(p.name()),
    useCustomFolder_(false),
    intOutsideBounds_(dict.getOrDefault("intOutsideBounds",true)),
    inverseDistRadius_(dict.getOrDefault("inverseDistRadius",-1.0)),
    startSamplePoints_(Zero),
    startSampleData_(Zero),
    endSamplePoints_(Zero),
    endSampleData_(Zero),
    domainRefPt_(Zero),
    domainMatDm_(Zero),
    domainVoxSz_(Zero)
{
    if
    (
        inputType_ != "unstructured"
     && inputType_ != "structured"
    )
    {
        FatalIOErrorInFunction(dict)
            << "inputType should be one of 'unstructured'"
            << ", 'structured', 'struc_mask'" << exit(FatalIOError);
    }

    if
    (
        interpolationType_ != "nearest"
     && interpolationType_ != "trilinear"
     && interpolationType_ != "inverseDist"
    )
    {
        FatalIOErrorInFunction(dict)
            << "interpolationType should be one of 'nearest'"
            << ", 'trilinear', 'inverseDist'" << exit(FatalIOError);
    }

    if
    (
        interpolationType_ == "inverseDist"
     && inverseDistRadius_ <= 0
    )
    {
        FatalIOErrorInFunction(dict)
            << "for interpolationType 'inverseDist'"
            << ", 'inverseDistRadius'>=0 should be defined for B.C."
            << exit(FatalIOError);
    }

    if
    (
        interpolationType_ == "trilinear"
     && inputType_ == "unstructured"
    )
    {
        FatalIOErrorInFunction(dict)
            << "interpolationType 'trilinear' can only be used with "
            << "structured input data (" << inputType_
            << " provided)" << exit(FatalIOError);
    }

    dict.readIfPresent("fieldTableName", fieldTableName_);
    dict.readIfPresent("inputFolderName", inputFolderName_);
    if (inputFolderName_ != p.name())
    {
        useCustomFolder_ = true;
    }

    if (dict.found("value"))
    {
        fixedValuePointPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        // Note: use evaluate to do updateCoeffs followed by a reset
        //       of the pointPatchField::updated_ flag. This is
        //       so if first use is in the next time step it retriggers
        //       a new update.
        pointPatchField<Type>::evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::timeVaryingMotionInterpolationPointPatchField<Type>::
timeVaryingMotionInterpolationPointPatchField
(
    const timeVaryingMotionInterpolationPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(ptf, p, iF, mapper),
    fieldTableName_(ptf.fieldTableName_),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    endSampleTime_(-1),
    endSampledValues_(0),
    inputType_(ptf.inputType_),
    interpolationType_(ptf.interpolationType_),
    inputFolderName_(p.name()),
    useCustomFolder_(ptf.useCustomFolder_),
    intOutsideBounds_(ptf.intOutsideBounds_),
    inverseDistRadius_(ptf.inverseDistRadius_),
    startSamplePoints_(Zero),
    startSampleData_(Zero),
    endSamplePoints_(Zero),
    endSampleData_(Zero),
    domainRefPt_(Zero),
    domainMatDm_(Zero),
    domainVoxSz_(Zero)
{
    if (ptf.useCustomFolder_)
    {
        inputFolderName_ = ptf.inputFolderName_;
    }
}


template<class Type>
Foam::timeVaryingMotionInterpolationPointPatchField<Type>::
timeVaryingMotionInterpolationPointPatchField
(
    const timeVaryingMotionInterpolationPointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf),
    fieldTableName_(ptf.fieldTableName_),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    inputType_(ptf.inputType_),
    interpolationType_(ptf.interpolationType_),
    inputFolderName_(ptf.inputFolderName_),
    useCustomFolder_(ptf.useCustomFolder_),
    intOutsideBounds_(ptf.intOutsideBounds_),
    inverseDistRadius_(ptf.inverseDistRadius_),
    startSamplePoints_(ptf.startSamplePoints_),
    startSampleData_(ptf.startSampleData_),
    endSamplePoints_(ptf.endSamplePoints_),
    endSampleData_(ptf.endSampleData_),
    domainRefPt_(ptf.domainRefPt_),
    domainMatDm_(ptf.domainMatDm_),
    domainVoxSz_(ptf.domainVoxSz_)
{}


template<class Type>
Foam::timeVaryingMotionInterpolationPointPatchField<Type>::
timeVaryingMotionInterpolationPointPatchField
(
    const timeVaryingMotionInterpolationPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    fieldTableName_(ptf.fieldTableName_),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    inputType_(ptf.inputType_),
    interpolationType_(ptf.interpolationType_),
    inputFolderName_(ptf.inputFolderName_),
    useCustomFolder_(ptf.useCustomFolder_),
    intOutsideBounds_(ptf.intOutsideBounds_),
    inverseDistRadius_(ptf.inverseDistRadius_),
    startSamplePoints_(ptf.startSamplePoints_),
    startSampleData_(ptf.startSampleData_),
    endSamplePoints_(ptf.endSamplePoints_),
    endSampleData_(ptf.endSampleData_),
    domainRefPt_(ptf.domainRefPt_),
    domainMatDm_(ptf.domainMatDm_),
    domainVoxSz_(ptf.domainVoxSz_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingMotionInterpolationPointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<Type>::autoMap(m);
    if (startSampledValues_.size())
    {
        startSampledValues_.autoMap(m);
        endSampledValues_.autoMap(m);
    }
    // Clear times list indices
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::timeVaryingMotionInterpolationPointPatchField<Type>::rmap
(
    const pointPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchField<Type>::rmap(ptf, addr);

    const timeVaryingMotionInterpolationPointPatchField<Type>& tiptf =
        refCast<const timeVaryingMotionInterpolationPointPatchField<Type>>(ptf);

    startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    endSampledValues_.rmap(tiptf.endSampledValues_, addr);

    // Clear times list indices
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::timeVaryingMotionInterpolationPointPatchField<Type>::checkTable()
{
    const Time& time = this->db().time();

    const polyMesh& pMesh = this->patch().boundaryMesh().mesh()();

    // Read the initial point position
    pointField meshPts;

    if (pMesh.pointsInstance() == pMesh.facesInstance())
    {
        meshPts = pointField(pMesh.points(), this->patch().meshPoints());
    }
    else
    {
        // Load points from facesInstance
        if (debug)
        {
            Info<< "Reloading points0 from " << pMesh.facesInstance()
                << endl;
        }

        pointIOField points0
        (
            IOobject
            (
                "points",
                pMesh.facesInstance(),
                polyMesh::meshSubDir,
                pMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        meshPts = pointField(points0, this->patch().meshPoints());
    }

    // Initialise
    if (startSampleTime_ == -1 && endSampleTime_ == -1)
    {
        // Structured domain properties file path
        // (Defined here because it's path is used anyways)
        const fileName domainInfoFile
        (
            time.path()
           /time.caseConstant()
           /"boundaryData"
           /inputFolderName_
           /"domainMatrixInfo"
        );

        // Read the times for which data is available

        const fileName samplePointsDir = domainInfoFile.path();
        sampleTimes_ = Time::findTimes(samplePointsDir);

        if (debug)
        {
            Info<< "timeVaryingMotionInterpolationPointPatchField : In directory "
                << samplePointsDir << " found times "
                << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
                << endl;
        }

        // Read structed data domain properties
        if (inputType_ == "structured")
        {
            IOobject ioDomainInfoFile
            (
                domainInfoFile,   // absolute path
                time,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false,              // no need to register
                true                // is global object (currently not used)
            );
            const rawIOField<point> domainInformationData(ioDomainInfoFile, false);
            if (domainInformationData.size() != 3 && domainInformationData.size() != 6)
            {
                FatalErrorInFunction
                    << "Length of file 'domainMatrixInfo' (" << domainInformationData.size()
                    << ") differs from allowed values (3 and 6) in file "
                    << domainInfoFile << exit(FatalError);
            }
            domainRefPt_ = domainInformationData[0];
            domainVoxSz_ = domainInformationData[1];
            domainMatDm_ = domainInformationData[2];
            // Create a variable with the structured point coordinates
            // THIS IS NOT EFFICIENT, BUT SHOULD NOT BE FORBIDDEN!!!
            if (interpolationType_ != "trilinear")
            {
                int m;
                int matSize = domainMatDm_.x()*domainMatDm_.y()*domainMatDm_.z();
                tmp<pointField> tstrPts(new pointField(matSize));
                pointField& strPts = tstrPts.ref();
                for(int ii = 0; ii<domainMatDm_.x(); ii++)
                {
                    for(int jj = 0; jj<domainMatDm_.y(); jj++)
                    {
                        for(int kk = 0; kk<domainMatDm_.z(); kk++)
                        {
                            m = (kk)*(domainMatDm_.x()*domainMatDm_.y()) +
                                (jj)*(domainMatDm_.x()) +
                                (ii);
                            strPts[m] = point
                            (
                                domainRefPt_.x()+ii*domainVoxSz_.x(),
                                domainRefPt_.y()+jj*domainVoxSz_.y(),
                                domainRefPt_.z()+kk*domainVoxSz_.z()
                            );
                        }
                    }
                }
                startSamplePoints_ = tstrPts;
            }
        }
    }

    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;

    bool foundTime = pointToPointPlanarInterpolation::findTime
    (
        sampleTimes_,
        startSampleTime_,
        time.value(),
        lo,
        hi
    );

    if (!foundTime)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << time.value() << nl
            << "Have sampling values for times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory "
            <<  time.path()/time.caseConstant()/"boundaryData"/inputFolderName_
            << "\n    on patch " << this->patch().name()
            << " of field " << fieldTableName_
            << exit(FatalError);
    }
    Info << "loooo" << startSampleTime_ << " " << endSampleTime_ << endl;

    // Update START sampled data fields.
    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;

        if (startSampleTime_ == endSampleTime_)
        {
            // No need to reread since are end values
            if (debug)
            {
                Pout<< "checkTable : Setting startValues to (already read) "
                    <<   "boundaryData"
                        /inputFolderName_
                        /sampleTimes_[startSampleTime_].name()
                    << endl;
            }
            startSampleData_ = endSampleData_;
            if (inputType_ == "unstructured")
            {
                startSamplePoints_ = endSamplePoints_;
            }
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading startValues from "
                    <<   "boundaryData"
                        /inputFolderName_
                        /sampleTimes_[lo].name()
                    << endl;
            }

            // Reread field values at points - path to pointMotionU or pointDisplacement
            const fileName valsFile
            (
                time.path()
               /time.caseConstant()
               /"boundaryData"
               /inputFolderName_
               /sampleTimes_[startSampleTime_].name()
               /fieldTableName_
            );
            IOobject ioField
            (
                valsFile,           // absolute path
                time,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false,              // no need to register
                true                // is global object (currently not used)
            );
            startSampleData_ = rawIOField<Type>(ioField, false);
            // Reread field points coordinates for unstructured
            if (inputType_ == "unstructured")
            {
                // Reread mask data
                const fileName pointsFile
                (
                    time.path()
                   /time.caseConstant()
                   /"boundaryData"
                   /inputFolderName_
                   /sampleTimes_[startSampleTime_].name()
                   /"points"
                );
                Info << "Reading startSamplePoints_ from file: "
                     << pointsFile << endl;
                IOobject ioPoints
                (
                    pointsFile,         // absolute path
                    time,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false,              // no need to register
                    true                // is global object (currently not used)
                );
                startSamplePoints_ = rawIOField<point>(ioPoints, false);

                if (startSampleData_.size() != startSamplePoints_.size())
                {
                    FatalErrorInFunction
                        << "Number of values (" << startSampleData_.size()
                        << ") differs from the number of points ("
                        <<  startSamplePoints_.size()
                        << ") in files " << valsFile
                        << " and " << pointsFile << exit(FatalError);
                }
            }
            else
            {
                // Check size of startSampleData_ for structured data
                int matSize = domainMatDm_.x()*domainMatDm_.y()*domainMatDm_.z();
                if (startSampleData_.size() != matSize)
                {
                    FatalErrorInFunction
                        << "Number of values (" << startSampleData_.size()
                        << ") differs from the number of cells in structured matrix ("
                        <<  matSize
                        << ") in files " << valsFile << exit(FatalError);
                }
            }
        }
    }

    // Update END sampled data fields.
    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;

        if (endSampleTime_ == -1)
        {
            // endTime no longer valid. Might as well clear endValues.
            if (debug)
            {
                Pout<< "checkTable : Clearing endValues" << endl;
            }
            endSampledValues_.clear();
            endSampleData_.clear();
            if (inputType_ == "unstructured")
            {
                endSamplePoints_.clear();
            }
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading endValues from "
                    <<   "boundaryData"
                        /inputFolderName_
                        /sampleTimes_[endSampleTime_].name()
                    << endl;
            }

            // Reread field values at points
            const fileName valsFile
            (
                time.path()
               /time.caseConstant()
               /"boundaryData"
               /inputFolderName_
               /sampleTimes_[endSampleTime_].name()
               /fieldTableName_
            );

            IOobject ioField
            (
                valsFile,           // absolute path
                time,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false,              // no need to register
                true                // is global object (currently not used)
            );
            endSampleData_ = rawIOField<Type>(ioField, false);


            // Reread field points coordinates for unstructured
            if (inputType_ == "unstructured")
            {
                // Reread mask data
                const fileName pointsFile
                (
                    time.path()
                   /time.caseConstant()
                   /"boundaryData"
                   /inputFolderName_
                   /sampleTimes_[endSampleTime_].name()
                   /"points"
                );
                IOobject ioPoints
                (
                    pointsFile,         // absolute path
                    time,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false,              // no need to register
                    true                // is global object (currently not used)
                );
                endSamplePoints_ = rawIOField<point>(ioPoints, false);

                if (endSampleData_.size() != endSamplePoints_.size())
                {
                    FatalErrorInFunction
                        << "Number of values (" << endSampleData_.size()
                        << ") differs from the number of points ("
                        <<  endSamplePoints_.size()
                        << ") in files " << valsFile
                        << " and " << pointsFile << exit(FatalError);
                }
            }
            else
            {
                // Check size of startSampleData_ for structured data
                int matSize = domainMatDm_.x()*domainMatDm_.y()*domainMatDm_.z();
                if (endSampleData_.size() != matSize)
                {
                    FatalErrorInFunction
                        << "Number of values (" << endSampleData_.size()
                        << ") differs from the number of cells in structured matrix ("
                        <<  matSize
                        << ") in files " << valsFile << exit(FatalError);
                }
            }
        }
    }

    // Do the interpolation of the data to mesh coordinates
    scalar deltaTime = time.value() - sampleTimes_[startSampleTime_].value();
    bool interpolateEnd = (endSampleTime_ != -1  && deltaTime > SMALL);

    if (interpolationType_ == "trilinear")
    {
        applyTrilinearInterpolation
        (
            meshPts,
            interpolateEnd
        );
    }
    else if (interpolationType_ == "nearest")
    {
        applyNearestValues
        (
            meshPts,
            interpolateEnd
        );
    }
    else if (interpolationType_ == "inverseDist")
    {
        applyInverseDistanceInterpolation
        (
            meshPts,
            interpolateEnd
        );
    }
    else
    {
        FatalErrorInFunction
            << "Illegal interpolation option."
            << abort(FatalError);
    }
}

// Linear interpolation
template<class Type>
Type Foam::timeVaryingMotionInterpolationPointPatchField<Type>::linearInterpolation
(
    const Type& edgeVal0,
    const Type& edgeVal1,
    const scalar& lf
)
{
    return edgeVal0*(1-lf)+edgeVal1*lf;
}

// Bilinear interpolation
template<class Type>
Type Foam::timeVaryingMotionInterpolationPointPatchField<Type>::bilinearInterpolation
(
    const Type& edgeVal00,
    const Type& edgeVal10,
    const Type& edgeVal01,
    const Type& edgeVal11,
    const scalar& lf1,
    const scalar& lf2
)
{
    Type linA = linearInterpolation(edgeVal00,edgeVal10,lf1);
    Type linB = linearInterpolation(edgeVal01,edgeVal11,lf1);
    return linearInterpolation(linA,linB,lf2);
}

// Trilinear interpolation
template<class Type>
Type Foam::timeVaryingMotionInterpolationPointPatchField<Type>::trilinearInterpolation
(
    const Type& edgeVal000,
    const Type& edgeVal100,
    const Type& edgeVal010,
    const Type& edgeVal110,
    const Type& edgeVal001,
    const Type& edgeVal101,
    const Type& edgeVal011,
    const Type& edgeVal111,
    const scalar& lf1,
    const scalar& lf2,
    const scalar& lf3
)
{
    // First bilinear interpolation
    Type biLinA = bilinearInterpolation
    (
        edgeVal000,
        edgeVal100,
        edgeVal010,
        edgeVal110,
        lf1,
        lf2
    );
    // Second bilinear interpolation
    Type biLinB = bilinearInterpolation
    (
        edgeVal000,
        edgeVal100,
        edgeVal010,
        edgeVal110,
        lf1,
        lf2
    );
    // Interpolate linearly between the two previous
    return linearInterpolation(biLinA,biLinB,lf3);
}


template<class Type>
void Foam::timeVaryingMotionInterpolationPointPatchField<Type>::applyTrilinearInterpolation
(
    const pointField& meshPts,
    const bool& interpolateEnd
)
{
    // Create temp field that will be assigned to start and end values
    // Start
    tmp<Field<Type>> tfldSta(new Field<Type>(meshPts.size()));
    Field<Type>& fldSta = tfldSta.ref();
    // End
    tmp<Field<Type>> tfldEnd;
    if (interpolateEnd)
    {
        tfldEnd = new Field<Type>(meshPts.size());
    }
    else
    {
        // This is ugly, I know.
        tfldEnd = new Field<Type>(0);
    }
    Field<Type>& fldEnd = tfldEnd.ref();

    // Pre-allocated variables
    int iX, iY, iZ;
    int iXd, iYd, iZd;
    int m000, m100, m010, m110, m001, m101, m011, m111;
    scalar l1, l2, l3;
    point invVoxX(1/domainVoxSz_.x(),1/domainVoxSz_.y(),1/domainVoxSz_.z());

    // Check that data matrix is not a plane or line, in which case the
    // interpolation WILL BE PERFORMED as bilinear or linear
    bool inputMat3d = ( domainMatDm_.x()>1 &&
                        domainMatDm_.y()>1 &&
                        domainMatDm_.z()>1);

    // Iterate
    forAll(meshPts, pp)
    {
        // Indices in 3d matrix for point 000 (preliminary)
        iX = (meshPts[pp].x()-domainRefPt_.x())*invVoxX.x();
        iY = (meshPts[pp].y()-domainRefPt_.y())*invVoxX.y();
        iZ = (meshPts[pp].z()-domainRefPt_.z())*invVoxX.z();

        // Check these are within the matrix bounds (false if inputMat3d=false)
        bool withinBounds = (iX >=0 && iX<(domainMatDm_.x()-1) &&
                             iY >=0 && iY<(domainMatDm_.y()-1) &&
                             iZ >=0 && iZ<(domainMatDm_.z()-1));
        // (Case A) Inside the structured grid
        if (withinBounds)
        {
            // Calculate interpolation terms
            l1 = invVoxX.x()*(meshPts[pp].x()-(domainRefPt_.x()+iX*domainVoxSz_.x()));
            l2 = invVoxX.y()*(meshPts[pp].y()-(domainRefPt_.y()+iY*domainVoxSz_.y()));
            l3 = invVoxX.z()*(meshPts[pp].z()-(domainRefPt_.z()+iZ*domainVoxSz_.z()));
            // Calculate matrix indices for cuboid edges
            m000 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX);
            m100 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX+1);
            m010 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY+1)*(domainMatDm_.x()) +
                    (iX);
            m110 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY+1)*(domainMatDm_.x()) +
                    (iX+1);
            m001 =  (iZ+1)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX);
            m101 =  (iZ+1)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX+1);
            m011 =  (iZ+1)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY+1)*(domainMatDm_.x()) +
                    (iX);
            m111 =  (iZ+1)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY+1)*(domainMatDm_.x()) +
                    (iX+1);
            // Do the trilinear interpolation
            fldSta[pp] = trilinearInterpolation
            (
                startSampleData_[m000],
                startSampleData_[m100],
                startSampleData_[m010],
                startSampleData_[m110],
                startSampleData_[m001],
                startSampleData_[m101],
                startSampleData_[m011],
                startSampleData_[m111],
                l1,
                l2,
                l3
            );
            if (interpolateEnd)
            {
                fldEnd[pp] = trilinearInterpolation
                (
                    endSampleData_[m000],
                    endSampleData_[m100],
                    endSampleData_[m010],
                    endSampleData_[m110],
                    endSampleData_[m001],
                    endSampleData_[m101],
                    endSampleData_[m011],
                    endSampleData_[m111],
                    l1,
                    l2,
                    l3
                );
            }
            continue;
        }
        // Case E (no interpolation!)
        else if (inputMat3d && !intOutsideBounds_)
        {
            fldSta[pp] *= 0;
            if (interpolateEnd)
            {
                fldEnd[pp] *= 0;
            }
            continue;
        }

        // Check where the point is when withinBounds=false
        // X - Update indices
        if (iX<0)
        {
            iX = 0;
            iXd = 0;
        }
        else if (iX>=(domainMatDm_.x()-1))
        {
            iX = domainMatDm_.x()-1;
            iXd = 0;
        }
        else
        {
            iXd = 1;
        }
        // Y - Update indices
        if (iY<0)
        {
            iY = 0;
            iYd = 0;
        }
        else if (iY>=(domainMatDm_.y()-1))
        {
            iY = domainMatDm_.y()-1;
            iYd = 0;
        }
        else
        {
            iYd = 1;
        }
        // Z - Update indices
        if (iZ<0)
        {
            iZ = 0;
            iZd = 0;
        }
        else if (iZ>=(domainMatDm_.z()-1))
        {
            iZ = domainMatDm_.z()-1;
            iZd = 0;
        }
        else
        {
            iZd = 1;
        }
        // Vector index of point 000 (first value) used in all cases below
        m000 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                (iY)*(domainMatDm_.x()) +
                (iX);
        // The type of interpolation is based on iXd, iYd, iZd
        // (Case B.1) On one of 2 faces perpendicular to Z
        if (iXd==1 && iYd==1 && iZd==0)
        {
            // Interpolate in X and Y (bilinear)
            l1 = invVoxX.x()*(meshPts[pp].x()-(domainRefPt_.x()+iX*domainVoxSz_.x()));
            l2 = invVoxX.y()*(meshPts[pp].y()-(domainRefPt_.y()+iY*domainVoxSz_.y()));
            m100 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX+1);
            m010 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY+1)*(domainMatDm_.x()) +
                    (iX);
            m110 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY+1)*(domainMatDm_.x()) +
                    (iX+1);
            fldSta[pp] = bilinearInterpolation
            (
                startSampleData_[m000],
                startSampleData_[m100],
                startSampleData_[m010],
                startSampleData_[m110],
                l1,
                l2
            );
            if (interpolateEnd)
            {
                fldEnd[pp] = bilinearInterpolation
                (
                    endSampleData_[m000],
                    endSampleData_[m100],
                    endSampleData_[m010],
                    endSampleData_[m110],
                    l1,
                    l2
                );
            }
        }
        // (Case B.2) On one of 2 faces perpendicular to Y
        else if (iXd==1 && iYd==0 && iZd==1)
        {
            // Interpolate in X and Z (bilinear)
            l1 = invVoxX.x()*(meshPts[pp].x()-(domainRefPt_.x()+iX*domainVoxSz_.x()));
            l3 = invVoxX.z()*(meshPts[pp].z()-(domainRefPt_.z()+iZ*domainVoxSz_.z()));
            m100 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX+1);
            m001 =  (iZ+1)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX);
            m101 =  (iZ+1)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX+1);
            fldSta[pp] = bilinearInterpolation
            (
                startSampleData_[m000],
                startSampleData_[m100],
                startSampleData_[m001],
                startSampleData_[m101],
                l1,
                l3
            );
            if (interpolateEnd)
            {
                fldEnd[pp] = bilinearInterpolation
                (
                    endSampleData_[m000],
                    endSampleData_[m100],
                    endSampleData_[m001],
                    endSampleData_[m101],
                    l1,
                    l3
                );
            }
        }
        // (Case B.3) On one of 2 faces perpendicular to X
        else if (iXd==0 && iYd==1 && iZd==1)
        {
            // Interpolate in Y and Z (bilinear)
            l2 = invVoxX.y()*(meshPts[pp].y()-(domainRefPt_.y()+iY*domainVoxSz_.y()));
            l3 = invVoxX.z()*(meshPts[pp].z()-(domainRefPt_.z()+iZ*domainVoxSz_.z()));
            m010 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY+1)*(domainMatDm_.x()) +
                    (iX);
            m001 =  (iZ+1)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX);
            m011 =  (iZ+1)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY+1)*(domainMatDm_.x()) +
                    (iX);
            fldSta[pp] = bilinearInterpolation
            (
                startSampleData_[m000],
                startSampleData_[m010],
                startSampleData_[m001],
                startSampleData_[m011],
                l2,
                l3
            );
            if (interpolateEnd)
            {
                fldEnd[pp] = bilinearInterpolation
                (
                    endSampleData_[m000],
                    endSampleData_[m010],
                    endSampleData_[m001],
                    endSampleData_[m011],
                    l2,
                    l3
                );
            }
        }
        // (Case C.1) On one of 4 edges aligned with X
        else if (iXd==1 && iYd==0 && iZd==0)
        {
            // Only interpolate in X (linear)
            l1 = invVoxX.x()*(meshPts[pp].x()-(domainRefPt_.x()+iX*domainVoxSz_.x()));
            m100 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX+1);
            fldSta[pp] = linearInterpolation
            (
                startSampleData_[m000],
                startSampleData_[m100],
                l1
            );
            if (interpolateEnd)
            {
                fldEnd[pp] = linearInterpolation
                (
                    endSampleData_[m000],
                    endSampleData_[m100],
                    l1
                );
            }
        }
        // (Case C.2) On one of 4 edges aligned with Y
        else if (iXd==0 && iYd==1 && iZd==0)
        {
            // Only interpolate in Y (linear)
            l2 = invVoxX.y()*(meshPts[pp].y()-(domainRefPt_.y()+iY*domainVoxSz_.y()));
            m010 =  (iZ)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY+1)*(domainMatDm_.x()) +
                    (iX);
            fldSta[pp] = linearInterpolation
            (
                startSampleData_[m000],
                startSampleData_[m010],
                l2
            );
            if (interpolateEnd)
            {
                fldEnd[pp] = linearInterpolation
                (
                    endSampleData_[m000],
                    endSampleData_[m010],
                    l2
                );
            }
        }
        // (Case C.3) On one of 4 edges aligned with Z
        else if (iXd==0 && iYd==0 && iZd==1)
        {
            // Only interpolate in Z (linear)
            l3 = invVoxX.z()*(meshPts[pp].z()-(domainRefPt_.z()+iZ*domainVoxSz_.z()));
            m001 =  (iZ+1)*(domainMatDm_.x()*domainMatDm_.y()) +
                    (iY)*(domainMatDm_.x()) +
                    (iX);
            fldSta[pp] = linearInterpolation
            (
                startSampleData_[m000],
                startSampleData_[m001],
                l3
            );
            if (interpolateEnd)
            {
                fldEnd[pp] = linearInterpolation
                (
                    endSampleData_[m000],
                    endSampleData_[m001],
                    l3
                );
            }
        }
        // (Case D) On one of the 8 regions closer to the corners
        // (iXd==0 && iYd==0 && iZd==0)
        else
        {
            // Assign the nearest point value to fields
            fldSta[pp] = startSampleData_[m000];
            if (interpolateEnd)
            {
                fldEnd[pp] = endSampleData_[m000];
            }
        }
    } // forAll(meshPts, pp)
    // Assign values from temporaty fields to variables
    startSampledValues_ = tfldSta;
    if (interpolateEnd)
    {
        endSampledValues_ = tfldEnd;
    }
}


template<class Type>
void Foam::timeVaryingMotionInterpolationPointPatchField<Type>::applyNearestValues
(
    const pointField& meshPts,
    const bool& interpolateEnd
)
{
    // Create a generic interpolator pointer
    autoPtr<pointToPointPlanarInterpolation> interpPtr;

    // Always interpolate start values
    interpPtr.reset
    (
        new pointToPointPlanarInterpolation
        (
        startSamplePoints_, // sourcePoints
        meshPts,            // destPoints
        0.00,                  // perturb (not used)
        true                // nearestOnly
        )
    );
    startSampledValues_ = interpPtr().interpolate(startSampleData_);

    // Interpolate end
    if (interpolateEnd)
    {
        if (inputType_ == "unstructured")
        {
            interpPtr.reset
            (
                new pointToPointPlanarInterpolation
                (
                endSamplePoints_,   // sourcePoints
                meshPts,            // destPoints
                0.00,                  // perturb (not used)
                true                // nearestOnly
                )
            );
        }

        endSampledValues_ = interpPtr().interpolate(endSampleData_);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::timeVaryingMotionInterpolationPointPatchField<Type>::inverseDistance
(
    const pointField& meshPts,
    const pointField samplePoints,
    const Field<Type> sampleData
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(meshPts.size()));
    Field<Type>& fld = tfld.ref();

    // A sortable list for distances from mesh to point data
    SortableList<scalar> distSorted(samplePoints.size());

    forAll(meshPts, pp)
    {
        // Assign new values to SortableList and sort it
        distSorted = mag(meshPts[pp]-samplePoints);
        distSorted.sort();

        Type interpNum(Zero);
        scalar interpDen(Zero);

        forAll(distSorted,jj)
        {
            scalar pointD = distSorted[jj];
            if (pointD>inverseDistRadius_)
            {
                break;
            }

            label pointI = distSorted.indices()[jj];
            if (pointD < SMALL)
            {
                interpNum = sampleData[pointI];
                interpDen = 1.0;
                break;
            }

            scalar pointInvD = 1/(pointD);

            interpNum += pointInvD * sampleData[pointI];
            interpDen += pointInvD;
        }

        if ( interpDen == 0)
        {
            fld[pp] = Type(Zero);
        }
        else
        {
            fld[pp] = interpNum/interpDen;
        }
    }

    return tfld;
}

template<class Type>
void Foam::timeVaryingMotionInterpolationPointPatchField<Type>::inverseDistance
(
    const pointField& meshPts
)
{
    tmp<Field<Type>> tfldSta(new Field<Type>(meshPts.size()));
    Field<Type>& fldSta = tfldSta.ref();

    tmp<Field<Type>> tfldEnd(new Field<Type>(meshPts.size()));
    Field<Type>& fldEnd = tfldEnd.ref();

    // A sortable list for distances from mesh to point data
    SortableList<scalar> distSorted(startSamplePoints_.size());

    forAll(meshPts, pp)
    {
        // Assign new values to SortableList and sort it
        distSorted = mag(meshPts[pp]-startSamplePoints_);
        distSorted.sort();

        Type interpNumSta(Zero);
        Type interpNumEnd(Zero);
        scalar interpDen(Zero);

        forAll(distSorted,jj)
        {
            scalar pointD = distSorted[jj];
            if (pointD>inverseDistRadius_)
            {
                break;
            }

            label pointI = distSorted.indices()[jj];
            if (pointD < SMALL)
            {
                interpNumSta = startSampleData_[pointI];
                interpNumEnd = endSampleData_[pointI];
                interpDen = 1.0;
                break;
            }

            scalar pointInvD = 1/(pointD);

            interpNumSta += pointInvD * startSampleData_[pointI];
            interpNumEnd += pointInvD * endSampleData_[pointI];
            interpDen += pointInvD;
        }

        if ( interpDen == 0)
        {
            fldSta[pp] = Type(Zero);
            fldEnd[pp] = Type(Zero);
        }
        else
        {
            fldSta[pp] = interpNumSta/interpDen;
            fldEnd[pp] = interpNumEnd/interpDen;
        }
    }
    startSampledValues_ = tfldSta;
    endSampledValues_ = tfldEnd;

}


template<class Type>
void Foam::timeVaryingMotionInterpolationPointPatchField<Type>::applyInverseDistanceInterpolation
(
    const pointField& meshPts,
    const bool& interpolateEnd
)
{
    if (inputType_ == "unstructured")
    {
        startSampledValues_ = inverseDistance
        (
            meshPts,
            startSamplePoints_,
            startSampleData_
        );

        if (interpolateEnd)
        {
            endSampledValues_ = inverseDistance
            (
                meshPts,
                endSamplePoints_,
                endSampleData_
            );
        }
    }
    else
    {
        if (interpolateEnd)
        {
            inverseDistance(meshPts);
        }
        else
        {
            startSampledValues_ = inverseDistance
            (
                meshPts,
                startSamplePoints_,
                startSampleData_
            );
        }
    }
}

template<class Type>
void Foam::timeVaryingMotionInterpolationPointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    Info<< "updateCoeffs: Current time = " << this->db().time().value() << endl;
    checkTable();
    // Interpolate between the sampled data
    scalar deltaTime = this->db().time().value() - sampleTimes_[startSampleTime_].value();
    // Info<< "startSampleTime_ = " << startSampleTime_ << ", endSampleTime_ = " << endSampleTime_ << endl;
    // Info<< "deltaTime = " << deltaTime << ", SMALL = " << SMALL << endl;
    endSampleTime_ = -1;//change - COMMENTED OUT to allow temporal interpolation
    if (endSampleTime_ == -1 || deltaTime < SMALL)
    {
        // only start value
        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, non-interpolated values"
                << " from start time:"
                << sampleTimes_[startSampleTime_].name() << nl;
        }

        this->operator==(startSampledValues_);

        // label nPrint = min(10, startSampledValues_.size());
        // for (label i=0; i<nPrint; ++i)
        // {
        //     Info<< "startSampledValues_[" << i << "] = " << startSampledValues_[i] << nl;
        // }
        // Info<< endl;
    }
    else
    {
        scalar start = sampleTimes_[startSampleTime_].value();
        scalar end = sampleTimes_[endSampleTime_].value();

        scalar s = (this->db().time().value()-start)/(end-start);
        

        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, interpolated values"
                << " between start time:"
                << sampleTimes_[startSampleTime_].name()
                << " and end time:" << sampleTimes_[endSampleTime_].name()
                << " with weight:" << s << endl;
        }
        
        this->operator==((1-s)*startSampledValues_ + s*endSampledValues_);

        // Field<Type> interpVals = (1-s)*startSampledValues_ + s*endSampledValues_;

        // for (label i=0; i<min(50, interpVals.size()); ++i)
        // {
        //     Info<< "interpVals[" << i << "] = " << interpVals[i] << nl;
        // }
    }

    if (debug)
    {
        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this) << endl;
    }

    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingMotionInterpolationPointPatchField<Type>::write
(
    Ostream& os
) const
{
    fixedValuePointPatchField<Type>::write(os);

    os.writeEntryIfDifferent
    (
        "fieldTable",
        this->internalField().name(),
        fieldTableName_
    );

    os.writeEntryIfDifferent<word>
    (
        "inputType",
        "unstructured",
        inputType_
    );

    os.writeEntryIfDifferent<word>
    (
        "interpolationType",
        "nearest",
        interpolationType_
    );

    os.writeEntryIfDifferent<word>
    (
        "inputFolderName",
        this->patch().name(),
        inputFolderName_
    );

    os.writeEntryIfDifferent<Switch>
    (
        "useCustomFolder",
        false,
        useCustomFolder_
    );

    os.writeEntryIfDifferent<Switch>
    (
        "outsideBounds",
        true,
        intOutsideBounds_
    );

    os.writeEntryIfDifferent<scalar>
    (
        "inverseDistRadius",
        -1,
        inverseDistRadius_
    );
}


// ************************************************************************* //
