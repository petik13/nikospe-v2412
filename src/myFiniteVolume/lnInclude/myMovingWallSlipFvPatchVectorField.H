/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::myMovingWallSlipFvPatchVectorField

Description
    Moving wall velocity with slip tangentially:
    - normal component matches wall motion,
    - tangential component from internal field (slip, zero shear).

Usage
    Example patch in 0/U:
    \verbatim
    movingWall
    {
        type            myMovingWallSlip;
        U               U;
        value           uniform (0 0 0);
    }
    \endverbatim
\*---------------------------------------------------------------------------*/

#ifndef myMovingWallSlipFvPatchVectorField_H
#define myMovingWallSlipFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

class myMovingWallSlipFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
    word UName_;   // Name of velocity field

public:

    //- Runtime type information
    TypeName("myMovingWallSlip");

    // Constructors

        myMovingWallSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        myMovingWallSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        myMovingWallSlipFvPatchVectorField
        (
            const myMovingWallSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        myMovingWallSlipFvPatchVectorField
        (
            const myMovingWallSlipFvPatchVectorField&
        );

        myMovingWallSlipFvPatchVectorField
        (
            const myMovingWallSlipFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new myMovingWallSlipFvPatchVectorField(*this)
            );
        }

        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new myMovingWallSlipFvPatchVectorField(*this, iF)
            );
        }

    // Member functions
        virtual void updateCoeffs();
        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
// ************************************************************************* //
