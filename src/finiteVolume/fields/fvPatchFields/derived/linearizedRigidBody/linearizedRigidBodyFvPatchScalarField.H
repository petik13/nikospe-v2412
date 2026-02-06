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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::linearizedRigidBodyFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a free-stream condition for pressure.

    It is an outlet-inlet condition that uses the velocity orientation to
    continuously blend between zero gradient for normal inlet and fixed value
    for normal outlet flow.

Usage
    \table
        Property        | Description             | Required    | Default value
        U               | velocity field name     | no          | U
        freestreamValue | freestream pressure     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            freestreamPressure;
        freestreamValue uniform 1e5;
    }
    \endverbatim

Note
    This condition is designed to operate with a freestreamVelocity condition

See also
    Foam::mixedFvPatchField
    Foam::freestreamFvPatchField

SourceFiles
    linearizedRigidBodyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef linearizedRigidBodyFvPatchScalarField_H
#define linearizedRigidBodyFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchField.H"
#include "spatialVector.H"
#include "vector.H"
#include "autoPtr.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class linearizedRigidBodyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class linearizedRigidBodyFvPatchScalarField
:
    public fixedGradientFvPatchField<scalar>
{
    // Private data

        //- Name of the velocity field
        word UName_;
	

	mutable autoPtr<Foam::OFstream> motionFilePtr_;
	void createMotionFile();
    void writeMotionHeader(OFstream& os) const;
    void writeMotion(const spatialVector& W);	
	
	// in the .H file
	label lastUpdateTimeIndex_;

	// 6DOF states (angular, linear)
    spatialVector a_old_,  a_new_,  aPrevIter_;
    spatialVector Xb_old_, Xb_new_;
    spatialVector Ub_old_, Ub_new_;

    // Rigid-body parameters
    scalar mass_;
    scalar xG_, zG_;

    scalar I4_, I5_, I6_, I46_;

    scalar C33_, C44_, C55_, C35_;

	scalar heading_;
	
public:

    //- Runtime type information
    TypeName("linearizedRigidBody");


    // Constructors

        //- Construct from patch and internal field
        linearizedRigidBodyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        linearizedRigidBodyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given linearizedRigidBodyFvPatchScalarField onto
        //  a new patch
        linearizedRigidBodyFvPatchScalarField
        (
            const linearizedRigidBodyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        linearizedRigidBodyFvPatchScalarField
        (
            const linearizedRigidBodyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new linearizedRigidBodyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        linearizedRigidBodyFvPatchScalarField
        (
            const linearizedRigidBodyFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new linearizedRigidBodyFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions




        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
		
		spatialVector computeForce();
		
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
