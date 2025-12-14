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
\*---------------------------------------------------------------------------*/

#include "myMovingWallSlipFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "tensorField.H"
#include "fvcMeshPhi.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

myMovingWallSlipFvPatchVectorField::myMovingWallSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    UName_("U")
{}


myMovingWallSlipFvPatchVectorField::myMovingWallSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


myMovingWallSlipFvPatchVectorField::myMovingWallSlipFvPatchVectorField
(
    const myMovingWallSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    UName_(ptf.UName_)
{}


myMovingWallSlipFvPatchVectorField::myMovingWallSlipFvPatchVectorField
(
    const myMovingWallSlipFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf),
    UName_(mwvpvf.UName_)
{}


myMovingWallSlipFvPatchVectorField::myMovingWallSlipFvPatchVectorField
(
    const myMovingWallSlipFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF),
    UName_(mwvpvf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void myMovingWallSlipFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    if (mesh.moving())
    {
        const fvPatch& p = patch();
        const polyPatch& pp = p.patch();
        const pointField& oldPoints = mesh.oldPoints();

        vectorField oldFc(pp.size());
        forAll(oldFc, i)
        {
            oldFc[i] = pp[i].centre(oldPoints);
        }

        const scalar deltaT = mesh.time().deltaTValue();

        // Wall velocity from motion
        const vectorField Up((pp.faceCentres() - oldFc)/deltaT);

        // Lookup velocity field
        const volVectorField& U =
            db().lookupObject<volVectorField>(UName_);

        // Patch flux
        const scalarField phip
        (
            p.patchField<surfaceScalarField, scalar>(fvc::flux(U))
        );

        const vectorField n(p.nf());
        const scalarField& magSf = p.magSf();
        tmp<scalarField> Un = phip / (magSf + VSMALL);
        // scalarField Un = tUn();   // make a full copy (safe)
        const vectorField nHat(p.nf());
//////////////////////////////////////////////////////////////////////////////
        // vectorField newValues = //change
        // (
        //     nHat*(nHat & (Up + n*(Un - (n & Up))))
        // + ((I - sqr(nHat)) & this->patchInternalField())
        // );

        // vectorField newValues = //change
        // (
        //     n*(n & (Up + n*(Un - (n & Up))))
        // + ((I - sqr(n)) & this->patchInternalField())
        // );

        // // Assign it
        // vectorField::operator=(newValues); 

        // --- DEBUG PRINTS ---
        // Print only for one patch (e.g. named "movingWall")
        // if (p.name() == "hull")
        // {
        //     Info << "Debug info for patch: " << p.name() << nl;

        //     Info << "Up (first 10):" << nl;
        //     for (label i = 0; i < min(label(10), Up.size()); i++)
        //         Info << "  Up[" << i << "] = " << Up[i] << nl;

        //     Info << "Un (first 10):" << nl;
        //     for (label i = 0; i < min(label(10), Un.size()); i++)
        //         Info << "  Un[" << i << "] = " << Un[i] << nl;

        //     Info << "n (first 10):" << nl;
        //     for (label i = 0; i < min(label(10), n.size()); i++)
        //         Info << "  n[" << i << "] = " << n[i] << nl;

        //     Info << "newValues (first 10):" << nl;
        //     for (label i = 0; i < min(label(10), newValues.size()); i++)
        //         Info << "  newValues[" << i << "] = " << newValues[i] << nl;

        //     Info << endl;
        // }

//////////////////////////////////////////////////////////////////////////////////
        // Normal velocity = wall motion; tangential = internal field
        vectorField::operator=
        (
            nHat*(nHat & (Up + n*(Un - (n & Up))))
          + ((I - sqr(nHat)) & this->patchInternalField())
        );

    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void myMovingWallSlipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    myMovingWallSlipFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
