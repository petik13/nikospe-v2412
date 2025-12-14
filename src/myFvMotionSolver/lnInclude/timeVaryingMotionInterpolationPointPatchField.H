/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::timeVaryingMotionInterpolationPointPatchField

Description
    A time-varying boundary condition for interpolation of motion values from
    provided data fields into the mesh.

See also
    Foam::timeVaryingMotionInterpolationFvPatchField

SourceFiles
    timeVaryingMotionInterpolationPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMotionInterpolationPointPatchField_H
#define timeVaryingMotionInterpolationPointPatchField_H

#include "fixedValuePointPatchField.H"
#include "instantList.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class timeVaryingMotionInterpolationPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingMotionInterpolationPointPatchField
:
    public fixedValuePointPatchField<Type>
{
    // Private data

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- List of boundaryData time directories
        instantList sampleTimes_;

        //- Current starting index in sampleTimes
        label startSampleTime_;

        //- Interpolated values from startSampleTime
        Field<Type> startSampledValues_;

        //- Current end index in sampleTimes
        label endSampleTime_;

        //- Interpolated values from endSampleTime
        Field<Type> endSampledValues_;

        //- Input data type
        word inputType_;

        //- Interpolation type
        word interpolationType_;

        //- Custom boundary subfolder
        word inputFolderName_;

        //- Use custom boundary folder?
        Switch useCustomFolder_;

        //- Outside bounds action
        Switch intOutsideBounds_;

        //- Inverse distance search radius
        scalar inverseDistRadius_;

        //- Sample data coordinates at start
        pointField startSamplePoints_;

        //- Sample data values at start
        Field<Type> startSampleData_;

        //- Sample data coordinates at end
        pointField endSamplePoints_;

        //- Sample data values at end
        Field<Type> endSampleData_;

        //- Reference point for sample data coordinates
        vector domainRefPt_;

        //- Dimensions (number of voxels in x, y, z) from sample data
        vector domainMatDm_;

        //- Voxel sizing for sample data voxels (x, y, z)
        vector domainVoxSz_;


public:

    //- Runtime type information
    TypeName("timeVaryingMotionInterpolation");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMotionInterpolationPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMotionInterpolationPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        timeVaryingMotionInterpolationPointPatchField
        (
            const timeVaryingMotionInterpolationPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMotionInterpolationPointPatchField
        (
            const timeVaryingMotionInterpolationPointPatchField<Type>&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new timeVaryingMotionInterpolationPointPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingMotionInterpolationPointPatchField
        (
            const timeVaryingMotionInterpolationPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new timeVaryingMotionInterpolationPointPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Utility functions

            //- Find boundary data inbetween current time and interpolate
            void checkTable();

            //- Apply trilinear interpolation on a structured grid
            void  applyTrilinearInterpolation
            (
                const pointField& meshPts,
                const bool& interpolateEnd
            );

            //- Apply nearest point values
            void  applyNearestValues
            (
                const pointField& meshPts,
                const bool& interpolateEnd
            );

            //- Inverse distance interpolation
            void applyInverseDistanceInterpolation
            (
                const pointField& meshPts,
                const bool& interpolateEnd
            );

            //- Trililear interpolation
            Type trilinearInterpolation
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
            );

            //- Bilinear interpolation
            // Assumes rectangular grid
            Type bilinearInterpolation
            (
                const Type& edgeVal00,
                const Type& edgeVal10,
                const Type& edgeVal01,
                const Type& edgeVal11,
                const scalar& lf1,
                const scalar& lf2
            );

            //- Linear interpolation
            Type linearInterpolation
            (
                const Type& edgeVal1,
                const Type& edgeVal2,
                const scalar& lf
            );

            //- Perform inverse distance interpolation
            tmp<Field<Type>> inverseDistance
            (
                const pointField& meshPts,
                const pointField samplePoints,
                const Field<Type> sampleData
            ) const;

            //- Perform inverse distance on structured data
            // on both start and end fields
            void inverseDistance
            (
                const pointField& meshPts
            );


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given PointPatchField onto
            //  this PointPatchField
            virtual void rmap
            (
                const pointPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeVaryingMotionInterpolationPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
