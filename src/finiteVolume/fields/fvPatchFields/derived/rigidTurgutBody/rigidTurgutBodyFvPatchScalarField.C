/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "vector.H"
#include "rigidTurgutBodyFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidTurgutBodyFvPatchScalarField::
rigidTurgutBodyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
     fixedGradientFvPatchField<scalar>(p, iF),
    UName_("U"),
	lastUpdateTimeIndex_(0),
    a_new_(vector::zero),
    Xb_new_(vector::zero),
    Ub_new_(vector::zero)
{}


Foam::rigidTurgutBodyFvPatchScalarField::
rigidTurgutBodyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<scalar>(p, iF, dict, IOobjectOption::NO_READ),  // ✅ correct
    UName_(dict.getOrDefault<word>("U", "U")),
    lastUpdateTimeIndex_(0),
    a_new_(vector::zero),
    Xb_new_(vector::zero),
    Ub_new_(vector::zero)
{
    gradient() = Zero;  // safe initialization
}


Foam::rigidTurgutBodyFvPatchScalarField::
rigidTurgutBodyFvPatchScalarField
(
    const rigidTurgutBodyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
     fixedGradientFvPatchField<scalar>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    lastUpdateTimeIndex_(ptf.lastUpdateTimeIndex_),
    a_new_(ptf.a_new_),
    Xb_new_(ptf.Xb_new_),
    Ub_new_(ptf.Ub_new_)
{}


Foam::rigidTurgutBodyFvPatchScalarField::
rigidTurgutBodyFvPatchScalarField
(
    const rigidTurgutBodyFvPatchScalarField& wbppsf
)
:
     fixedGradientFvPatchField<scalar>(wbppsf),
    UName_(wbppsf.UName_),
    lastUpdateTimeIndex_(wbppsf.lastUpdateTimeIndex_),
    a_new_(wbppsf.a_new_),
    Xb_new_(wbppsf.Xb_new_),
    Ub_new_(wbppsf.Ub_new_)
{}


Foam::rigidTurgutBodyFvPatchScalarField::
rigidTurgutBodyFvPatchScalarField
(
    const rigidTurgutBodyFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
     fixedGradientFvPatchField<scalar>(wbppsf, iF),
    UName_(wbppsf.UName_),
    lastUpdateTimeIndex_(wbppsf.lastUpdateTimeIndex_),
    a_new_(wbppsf.a_new_),
    Xb_new_(wbppsf.Xb_new_),
    Ub_new_(wbppsf.Ub_new_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rigidTurgutBodyFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const Time& runTime = db().time();
    const scalar dt = runTime.deltaTValue();
	const scalar aRelax_=0.3;
	const scalar aDamp_  = 1.0;  // like accelerationDamping (<1 adds damping)

	// -- Load Ucur 
    const volVectorField& Ucur = db().lookupObject<volVectorField>("Ucur");

    // current patch index
    const label patchi = patch().index();

    // patch unit normals
    const vectorField np(patch().nf());

    // Gradient of U on the patch as a tensorField (per face)
    tmp<volTensorField> tGradU = fvc::grad(Ucur);   // tmp<volTensorField>
    const volTensorField& gradU = tGradU();         // unwrap tmp -> volTensorField

    const fvPatchTensorField& gradUpP = gradU.boundaryField()[patchi];

    // hard checks
    if (gradUpP.size() != patch().size() || np.size() != patch().size())
    {
        FatalErrorInFunction
            << "Patch size mismatch on patch " << patch().name()
            << " patch.size=" << patch().size()
            << " gradUpP.size=" << gradUpP.size()
            << " np.size=" << np.size()
            << abort(FatalError);
    }

    // make plain Fields (owned storage, no patch refs)
    tensorField gradUp(gradUpP);          // copy
    tensorField gradUpT(gradUp.T());      // transpose into a real field
    vectorField nGradU(patch().size());   // allocate

    forAll(nGradU, i)
    {
        nGradU[i] = gradUpT[i] & np[i];   // Tensor & Vector -> Vector (no temps)
    }


	if (!masscomputed_) computeMass();
	if (!C33computed_)  computeC33();
	if (runTime.timeIndex() == 1)
    {
		a_old_=vector::zero;
		Xb_old_=vector::zero;
		Ub_old_=vector::zero;
		aPrevIter_ = vector::zero;
	}
    
	// Parameters for Newmark
    const scalar gamma = 0.5;
    const scalar beta  = 0.25;
	
	const vector restoring = vector(0.0,0.0,C33_);
	//burayi degistiriyorum
	
	
	if (runTime.timeIndex() != lastUpdateTimeIndex_)
    {
		Xb_old_=Xb_new_;
		a_old_=a_new_;
		Ub_old_=Ub_new_;
		writeForceEntry(Xb_old_);   // write force
	}
	
	
	
	vector Fnew = computeForce();
	vector Fh = restoring * Xb_new_.z();
	
	a_new_ = (Fnew - Fh ) / mass_;
	
	a_new_ = aDamp_ * (aRelax_ * a_new_ + (1.0 - aRelax_) * aPrevIter_);
    aPrevIter_ = a_new_;
	
	
	
	// --- Update velocity
	Ub_new_ = Ub_old_ + dt * ((1 - gamma)*a_old_ + gamma*a_new_);
	
	Xb_new_ =  Xb_old_ + dt*Ub_old_ + dt*dt*((0.5 - beta)*a_old_ + beta*a_new_);
	
	
	//burayi degistiriyorum
	Info << "BC iteration " << runTime.timeIndex() << " aDamp_ " << aDamp_ 
     << " | Xb_new_ = " << Xb_new_.z()
     << " | Fnew = " << Fnew.z() << endl;

	
    // mj term: (n · ∇)U & Xb_new_ on patch faces
    scalarField mj = nGradU & Xb_new_;

	// --- Apply Neumann BC: ∂Φ/∂n = n·Ub
    tmp<vectorField> n = patch().nf();
	scalarField rhs = -(n() & Ub_new_) + mj;
	
	
	this->gradient() = rhs;
  
	lastUpdateTimeIndex_ = runTime.timeIndex();
	
	
    fixedGradientFvPatchField<scalar>::updateCoeffs();  // ✅ call parent
	
}


/*---------------------------------------------------------------------------*\
  if (runTime.timeIndex() != lastUpdateTimeIndex_)
    {
		// --- Compute current total force
		vector F = computeForce();
		
		vector Fh = restoring * Xb_old_.z();

		// --- Compute new acceleration
		vector aNew = (F - Fh ) / mass_;  //burda arti yaptik

		// --- Update velocity
		vector Ub = Ub_old_ + dt * ((1 - gamma)*a_old_ + gamma*aNew);

		// --- Update position
		vector Xb = Xb_old_ + dt*Ub_old_ + dt*dt*((0.5 - beta)*a_old_ + beta*aNew);

		// --- Store new states for next step
		Ub_old_ = Ub;
		a_old_  = aNew;
		Xb_old_ = Xb;
		
		writeForceEntry(Xb);   // write force
	}

\*---------------------------------------------------------------------------*/





// * * * * * * * Helper Functions  * * * * * * * * * * * * * //



namespace Foam {  

void Foam::rigidTurgutBodyFvPatchScalarField::createForceFile()
{
    
	if (!Pstream::master()) return;  // Only master writes

	if (forceFilePtr_.valid()) return;
	
    const Time& runTime = db().time();

    //fileName dir(runTime.path()/"postProcessing"/"rigidBodyMotion");
    fileName dir("./postProcessing/rigidBodyMotion");  
	mkDir(dir);

    fileName filePath(dir/"motionTurgut.dat");
    mkDir(filePath.path());

    Info<< "Creating force output file: " << filePath << endl;

    forceFilePtr_.reset(new OFstream(filePath));
    writeForceHeader(forceFilePtr_());
}


void Foam::rigidTurgutBodyFvPatchScalarField::writeForceHeader(OFstream& os) const
{
    os << "# Motion of Center of Gravity: " << patch().name() << endl;
	os << "# Developed by  Turgut" << endl;	
	os << "#  					 " << endl;	
	os << "# Time   eta1   eta2   eta3" << endl;
    
}


void Foam::rigidTurgutBodyFvPatchScalarField::writeForceEntry(const vector& F)
{
    if (!Pstream::master()) return;  // Only master writes

	
	if (!forceFilePtr_.valid())
        createForceFile();

    const Time& runTime = db().time();

    (*forceFilePtr_)
        << runTime.timeOutputValue() << '\t'
        << F.x() << '\t'
        << F.y() << '\t'
        << F.z() << endl;
}



vector Foam::rigidTurgutBodyFvPatchScalarField::computeForce()
{
    // --- Geometry info
    
	const vectorField& Sfb = patch().Sf();
	
    // --- Look up the pressure field from the object registry
    //const volScalarField& p = db().lookupObject<volScalarField>("p");
        
	const volScalarField& p = db().lookupObject<volScalarField>("p2");
	const scalarField& pPatch = p.boundaryField()[patch().index()];

	
	

    // --- Initialize force vector
    vector F = vector::zero;
	forAll(Sfb, i) F += pPatch[i] * Sfb[i];
	
	reduce(F, sumOp<vector>());
	
	
	//Parallellestirme yapman lazim
	
	//if (this->db().time().timeIndex() == 1)
    //{
      //  Info<< "Skipping computeForce() at first time step (pressure not initialized)" << endl;
       // return vector::zero;
    //}
	//else {
    //forAll(Sfb, i) F += pPatch[i] * Sfb[i];
	
	//}
	
	
	
    return F;
}


void Foam::rigidTurgutBodyFvPatchScalarField::computeMass()
{
    // --- Geometry info
    
	const vectorField& Sfb = patch().Sf();
	const vectorField& Cf = patch().Cf();
    const scalarField zComponents = Cf.component(vector::Z);

	vector mass=vector::zero;

    
	forAll(Sfb, i) mass += - zComponents[i] * Sfb[i];
	
	reduce(mass, sumOp<vector>());
	
	
	mass_ = mass.component(vector::Z);
	
	Info << " mass_ is " << mass_ << endl;
	
	masscomputed_ = true;
}

void Foam::rigidTurgutBodyFvPatchScalarField::computeC33()
{
    // --- Geometry info
    
	const vectorField& Sfb = patch().Sf();
	const scalarField n3Components = Sfb.component(vector::Z);

	
	scalar C33 = 0.0;
    
	forAll(Sfb, i) C33 += 9.81*n3Components[i];
	
	reduce(C33, sumOp<scalar>());
	
	
	C33_ = C33;
	Info << " C33 is " << C33_ << endl;
	
	C33computed_ = true;
}




}




void Foam::rigidTurgutBodyFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        rigidTurgutBodyFvPatchScalarField
    );
}

// ************************************************************************* //
