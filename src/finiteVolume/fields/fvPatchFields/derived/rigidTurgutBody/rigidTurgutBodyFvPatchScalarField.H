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
    Foam::rigidTurgutBodyFvPatchScalarField

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
    rigidTurgutBodyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef rigidTurgutBodyFvPatchScalarField_H
#define rigidTurgutBodyFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchField.H"

#include "vector.H"
#include "autoPtr.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class rigidTurgutBodyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class rigidTurgutBodyFvPatchScalarField
:
    public fixedGradientFvPatchField<scalar>
{
    // Private data

        //- Name of the velocity field
        word UName_;
	

	mutable Foam::autoPtr<Foam::OFstream> forceFilePtr_;
	void createForceFile();
    void writeForceHeader(OFstream& os) const;
    void writeForceEntry(const vector& F);	
	
	// in the .H file
	label lastUpdateTimeIndex_;
	vector a_old_;
	vector Xb_old_;
	vector Ub_old_;
	vector a_new_;
	vector Xb_new_;
	vector Ub_new_;
	vector aPrevIter_;
	
	scalar mass_;
	scalar C33_;
	bool masscomputed_=false;
	bool C33computed_=false;
	
	
public:

    //- Runtime type information
    TypeName("rigidTurgutBody");


    // Constructors

        //- Construct from patch and internal field
        rigidTurgutBodyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rigidTurgutBodyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rigidTurgutBodyFvPatchScalarField onto
        //  a new patch
        rigidTurgutBodyFvPatchScalarField
        (
            const rigidTurgutBodyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rigidTurgutBodyFvPatchScalarField
        (
            const rigidTurgutBodyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new rigidTurgutBodyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rigidTurgutBodyFvPatchScalarField
        (
            const rigidTurgutBodyFvPatchScalarField&,
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
                new rigidTurgutBodyFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions




        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
		
		vector computeForce();
		
		virtual void computeMass();
		
		virtual void computeC33();
		
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
