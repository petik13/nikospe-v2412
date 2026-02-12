/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "fvCFD.H"
#include "fvMesh.H"
#include "waveCurrentPotential3DFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "gravityMeshObject.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
//#include "PrimitivePatchInterpolation.H"
#include "OFstream.H"
//#include "SVD.H"
#include "processorPolyPatch.H"
#include "waveCurPar3DPotUPFD5InlineHelpersInt.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::waveCurrentPotential3DFvPatchScalarField::ddtSchemeType
>
Foam::waveCurrentPotential3DFvPatchScalarField::ddtSchemeTypeNames_
({
    {
        ddtSchemeType::tsEuler,
        fv::EulerDdtScheme<scalar>::typeName_()
    },
    {
        ddtSchemeType::tsCrankNicolson,
        fv::CrankNicolsonDdtScheme<scalar>::typeName_()
    },
    
	
	
	{
        ddtSchemeType::tsBackward,
        fv::backwardDdtScheme<scalar>::typeName_()
    }
	
	
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveCurrentPotential3DFvPatchScalarField::
waveCurrentPotential3DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_("phi"),
    zetaName_("zeta"),
    rhoName_("rho"),
	PHIName_("Phi"),
	UName_("U"),
	shape_(100.0)
	
	
{}


Foam::waveCurrentPotential3DFvPatchScalarField::
waveCurrentPotential3DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    zetaName_(dict.getOrDefault<word>("zeta", "zeta")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
	PHIName_(dict.getOrDefault<word>("Phi", "Phi")),
	UName_(dict.getOrDefault<word>("U", "U")),
	shape_(dict.getOrDefault<scalar>("shape", 100.0)) 
	
	
{{ensureDictLoaded();          // create IOdictionary object
    readParamsFrom(*waveCurDictPtr_); // parse once now
    Info<< "waveCurConditions loaded in constructor." << nl;}}


Foam::waveCurrentPotential3DFvPatchScalarField::
waveCurrentPotential3DFvPatchScalarField
(
    const waveCurrentPotential3DFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    zetaName_(ptf.zetaName_),
    rhoName_(ptf.rhoName_),
	PHIName_(ptf.PHIName_),
	UName_(ptf.UName_),
	shape_(ptf.shape_) 
	
{}


Foam::waveCurrentPotential3DFvPatchScalarField::
waveCurrentPotential3DFvPatchScalarField
(
    const waveCurrentPotential3DFvPatchScalarField& wspsf
)
:
    fixedValueFvPatchScalarField(wspsf),
    phiName_(wspsf.phiName_),
    zetaName_(wspsf.zetaName_),
    rhoName_(wspsf.rhoName_),
	PHIName_(wspsf.PHIName_),
	UName_(wspsf.UName_),
	shape_(wspsf.shape_) 
{}


Foam::waveCurrentPotential3DFvPatchScalarField::
waveCurrentPotential3DFvPatchScalarField
(
    const waveCurrentPotential3DFvPatchScalarField& wspsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wspsf, iF),
    phiName_(wspsf.phiName_),
    zetaName_(wspsf.zetaName_),
    rhoName_(wspsf.rhoName_),
	PHIName_(wspsf.PHIName_),
	UName_(wspsf.UName_),
	shape_(wspsf.shape_) 
	 
	
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::waveCurrentPotential3DFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
	
	if ( db().time().timeIndex() == lastUpdateTimeIndex ){
		Info << "waveCurrentPotential3DFvPatchScalarField--saved BC applied at : " << lastUpdateTimeIndex << endl;
		operator==(data_);
		return;	
	}	

	// -- At top of updateCoeffs() after the early returns:
	ensureDictLoaded();
	if (waveCurDictPtr_->readIfModified() || !params_.init)
	{
		readParamsFrom(*waveCurDictPtr_);
		Info<< "waveCurConditions (re)loaded." << nl;
	}
	
	// Use locals for readability & no overhead in hot loop
	const scalar U0         = params_.U0;
	const scalar steepness  = params_.steepness;
	const scalar wavelength = params_.wavelength;
	const scalar hdepth     = params_.hdepth;
	const scalar rampperiod = params_.rampperiod;

	const scalar v0         = params_.v0;
	const scalar xdamp      = params_.xdamp;
	const scalar Lxdamp      = params_.Lxdamp;
	const scalar ydamp      = params_.ydamp;
	const scalar Lydamp      = params_.Lydamp;

	const scalar xsponge    = params_.xsponge;
	const scalar Lsponge    = params_.Lsponge;

	
		

	Info << "Update status: " << updated() << endl; 
    const label patchi = patch().index();
	const label nFaces = patch().size();
	
	const vectorField& patchFaceCenters = patch().Cf();
	scalarField xComponents = patchFaceCenters.component(vector::X);
	scalarField yComponents = patchFaceCenters.component(vector::Y);
	scalarField zComponents = patchFaceCenters.component(vector::Z);
	
    const scalar dt = db().time().deltaTValue();
	const scalar tt=db().time().value();
	
	Info << "Time in BC : " << tt << endl;
	Info << "Iteration: " << db().time().timeIndex() << endl;
	
    // Retrieve non-const access to zeta field from the database
    volVectorField& zeta = db().lookupObjectRef<volVectorField>(zetaName_);
    vectorField& zetap = zeta.boundaryFieldRef()[patchi]; // current value on patch

	//Turgut
	volScalarField& Phi = db().lookupObjectRef<volScalarField>(PHIName_);
	const volScalarField& Phi0 = Phi.oldTime();
	//scalarField& PhiPatch = Phi.boundaryFieldRef()[patchi];
	const scalarField& Phi0Patch = Phi0.boundaryField()[patchi]; // Previous value on patch
	//Turgut
	
	
		// Lookup d/dt scheme from database for zeta
		const word ddtSchemeName(zeta.mesh().ddtScheme(zeta.name()));
		ddtSchemeType ddtScheme(ddtSchemeTypeNames_[ddtSchemeName]);

		
		// Cache the patch face-normal vectors
		//tmp<vectorField> nf(patch().nf());
		//const vectorField& nfRef = nf.ref();  // Store the reference safely
		const vectorField nfRef(nFaces, vector(0.0,0.0,1.0));
		
		
		const volVectorField& Uvel = db().lookupObjectRef<volVectorField>(UName_);
		vectorField Wn(nfRef* (Uvel.boundaryField()[patchi] & nfRef)); // vertical velocity component (vectorized)
		
		// - Incident wave vertical velocity
		const scalar amp(0.5*steepness*wavelength);
		const scalar wavenumber (2.0*Foam::constant::mathematical::pi/wavelength);
		const scalar w(sqrt(9.81 * wavenumber * tanh(wavenumber*hdepth)));
		const scalar celerity(w/wavenumber);
		const scalar T(2.0*Foam::constant::mathematical::pi/w);
		const scalar ramp_time(3.0 * T);
		const scalar ramp_factor = 0.5 * (1 - cos(Foam::constant::mathematical::pi * min(1.0, tt / ramp_time)));

		Info<< "wavenumber k =" << wavenumber << nl;
		Info<< "angular frequency w =" << w << nl;
		Info<< "Current speed U0 =" << U0 << nl;

		// - Incident wave vertical velocity
		vectorField Wn0
		(
			amp * wavenumber * 9.81 / w 
		 * (Foam::sinh(wavenumber *(hdepth+zComponents))/Foam::cosh(wavenumber*hdepth))
		 * Foam::sin(wavenumber * xComponents- w*tt) * ramp_factor * nfRef
		);

		// - damping factor
		scalarField dampingterm=
		    // x-side damping active for x > xdamp
		    pos(xComponents - xdamp) * v0 * ((xComponents - xdamp) / (Lxdamp)) * ((xComponents - xdamp) / (Lxdamp)) * (nfRef & Wn);
		    // y-side damping active only when x > 5
		//    + pos(xComponents-xsponge)*pos(xdamp-xComponents)*pos(yComponents - ydamp) * v0 * ((yComponents - ydamp) / (Lydamp)) * ((yComponents - ydamp) / (Lydamp)) * (nfRef & Wn)
		//    + pos(xComponents-xsponge)*pos(xdamp-xComponents)*pos(-yComponents - ydamp) * v0 * ((-yComponents - ydamp) / (Lydamp)) * ((-yComponents - ydamp) / (Lydamp)) * (nfRef & Wn);
		    // inlet reflection damping: active for x in [0,5], stronger near x=0
		//   + pos(xsponge - xComponents) * v0 * ((xsponge - xComponents) / (Lsponge)) * ((xsponge - xComponents) / (Lsponge)) * (nfRef & (Wn - Wn0));
		


		// const scalarField wMeas(nfRef & Wn);
		// const scalarField wRef (nfRef & Wn0);
		// const scalarField wErr (wMeas - wRef);
		// const scalarField m(pos(xsponge - xComponents));

		// const scalar nMasked = gSum(m);

		// if (nMasked > SMALL && db().time().timeIndex() % 1 == 0)
		// {
		// 	const scalar maxAbsMeas = gMax(mag(m*wMeas));
		// 	const scalar maxAbsRef  = gMax(mag(m*wRef));
		// 	const scalar maxAbsErr  = gMax(mag(m*wErr));

		// 	const scalar rmsMeas = Foam::sqrt(gSum(m*sqr(wMeas))/nMasked);
		// 	const scalar rmsRef  = Foam::sqrt(gSum(m*sqr(wRef ))/nMasked);
		// 	const scalar rmsErr  = Foam::sqrt(gSum(m*sqr(wErr ))/nMasked);
		// 	const scalar meanErr = gSum(m*wErr)/nMasked;
		// 	const scalar meanAbsErr = gSum(m*mag(wErr))/nMasked;

		// 	Info<< "Sponge diag t=" << tt
		// 		<< "  max|meas|=" << maxAbsMeas
		// 		<< "  max|ref|="  << maxAbsRef
		// 		<< "  max|err|="  << maxAbsErr
		// 		<< "  rms(meas/ref/err)=" << rmsMeas << " " << rmsRef << " " << rmsErr
		// 		<< "  err/ref(max)=" << (maxAbsErr/(maxAbsRef + VSMALL))
		// 		<< "Mean abs error:" << meanAbsErr
		// 		<< nl;
		// }

		// const scalar corr =
		// gSum(m*wMeas*wRef) /
		// (Foam::sqrt(gSum(m*sqr(wMeas))*gSum(m*sqr(wRef))) + VSMALL);

		// Info<< "corr(meas,ref)=" << corr << nl;




		//From SnGrad BUNU DENICEZ SONRA
		//const auto& U2 = db().lookupObject<surfaceScalarField>("U2");
		//vectorField Wn(nf()*U2.boundaryField()[patchi]);
		
		// Retrieve the flux field from the database
        //const auto& phi = db().lookupObject<surfaceScalarField>(phiName_);
		//vectorField dZetap(dt*nf()*phi.boundaryField()[patchi]/patch().magSf());
		
		const volVectorField& zeta0 = zeta.oldTime(); 
		const vectorField& zeta0p = zeta0.boundaryField()[patchi]; // Previous value on patch
		
		scalarField turgut(Phi0Patch.size(), 0.0);
		vectorField UetaDx(Phi0Patch.size(), vector::zero);
		vectorField Wcurdz_zeta0p(Phi0Patch.size(), vector::zero);
		
		if (std::fabs(U0) > SMALL) // if current speed is non-zero
		{
			calcNeigboursV3(); // calc neighboors if not yet calculated (PrePar)
			if (neigboursCalculated_) // set in calcNeigboursV3
			{
				findUpwindDownwindNodesV2();
				
				detectFDSchemes();
				
				
				zetaDx_.setSize(nFaces, vector::zero);
				PhiDx_.setSize(nFaces, Zero);PhiDy_.setSize(nFaces, Zero);
				
				UPFDV2();
				
					
			}
			
			
			
			// Lookup Ucur and compute UetaDx
			const volVectorField& Ucur = db().lookupObjectRef<volVectorField>("Ucur");
			const vectorField& Ucur_p = Ucur.boundaryField()[patchi];
			UetaDx = nfRef * (Ucur_p & zetaDx_);

			// Lookup PhiCurDz2 and compute Wcurdz_zeta0p
			const volScalarField& PhiCurDz2 = db().lookupObjectRef<volScalarField>("PhiCurDz2");
			const scalarField& PhiCurDz2_p = PhiCurDz2.boundaryField()[patchi];
			Wcurdz_zeta0p = PhiCurDz2_p * zeta0p;

			
			
			forAll(turgut, i)
			{
				turgut[i] = Ucur_p[i][vector::X] * PhiDx_[i]+Ucur_p[i][vector::Y] * PhiDy_[i];// + 0.5*((Ucur_p[i] & Ucur_p[i]) - U0*U0) ;   //BUNU modify ettim
			}
		}

		
	

		
		//BUNLAR PARALLELDE CALISMIYOR
		//Info << "zetaDx_p [0] :"  << zetaDx_p[0] << endl;
		//Info << "UetaDX [0] :"  << UetaDx[0] << endl;
		//Info << "PhiDxPatch [0] :"  << PhiDxPatch[0] << endl;
		//Info << "zeta0p[1]:" << zeta0p[1] << endl;
	    //Info << "zetap[1] before the update:" << zetap[1] << endl;
		
		
		// Retriving g
		const uniformDimensionedVectorField& g = meshObjects::gravity::New(db().time());
		const vector gVal = -g.value();
		
		scalarField phiCalc(nFaces, 0.0);
		
	   
		
		
		switch (ddtScheme)
		{
			case tsEuler:
			case tsCrankNicolson:
			{
				
				Info << "Time Scheme = " << ddtSchemeTypeNames_[ddtScheme] << endl;
				
				if ( std::fabs(U0) > 0.0 ){
				
				
					Info << "Current speed is nonzero " << endl;
					
					
					zetap = zeta0p + dt*(Wn+Wcurdz_zeta0p-UetaDx);
					Info << "Turgut ayri applied" << endl;
					phiCalc = ((gVal & zeta0p)+turgut + dampingterm)*dt + Phi0Patch;
					//phiCalc = ((gVal & zetap))*dt + Phi0Patch;
					
				}
				else{
					
				Info << "Current speed is zero " << endl;
				
				zetap = zeta0p + dt*Wn;
				phiCalc = ((gVal & zetap)+ dampingterm)*dt + Phi0Patch;
				
				}
				
				
				
				break;
			}
			case tsBackward:
			{
				
				Info << "Time Scheme = " << ddtSchemeTypeNames_[ddtScheme] << endl;
				
				if ( db().time().timeIndex() == 1 )	 
				{
					if ( std::fabs(U0) > 0.0 )
						{			
						zetap = zeta0p + dt*(Wn+Wcurdz_zeta0p-UetaDx);
						Info << "Turgut ayri applied" << endl;
						phiCalc = ((gVal & zeta0p)+turgut)*dt + Phi0Patch;
						
						WnOld_=(Wn+Wcurdz_zeta0p-UetaDx);
						DPhiold_=((gVal & zeta0p)+turgut);
						}
						else{
							
						Info << "First time step Euler (U0=0) implemented " << endl;
						zetap = zeta0p + dt*Wn;
						phiCalc = (gVal & (zetap))*dt + Phi0Patch;	
							
						WnOld_=Wn;
						DPhiold_=((gVal & zeta0p));
						}
					
					
					
				}	 
				else if ( db().time().timeIndex() == 2 )	 	 
				{
					 if ( std::fabs(U0) > 0.0 ){
									
					//zetaInt = zeta0p+0.5*dt*((Wn-UetaDx)+WnOld_);
					//phiInt = (1.5*((gVal & zeta0p)+(0.5)*turgut)-0.5*DPhiold_)*dt + Phi0Patch;
					
					Info << "Current speed is nonzero " << endl;
					zetap=zeta0p+dt*(1.5*(Wn+Wcurdz_zeta0p-UetaDx)-0.5*WnOld_);
					phiCalc = (1.5*((gVal & zeta0p)+turgut + dampingterm)-0.5*DPhiold_)*dt + Phi0Patch;
					
					
					WnOld2_=WnOld_;
					DPhiold2_=DPhiold_;
					
					WnOld_=(Wn+Wcurdz_zeta0p-UetaDx);
					DPhiold_=((gVal & zeta0p)+turgut + dampingterm);
					
					}
					else{
					
					
					Info << "Current speed U0 is zero " << endl;
					zetap=zeta0p+dt*(1.5*Wn-0.5*WnOld_);
					phiCalc = (1.5*((gVal & zeta0p) + dampingterm)-0.5*DPhiold_)*dt + Phi0Patch;
					
					WnOld2_=WnOld_;
					DPhiold2_=DPhiold_;
					
					WnOld_=Wn;
					DPhiold_=((gVal & zeta0p) + dampingterm);
					
					}
									
					 
					 
				 }
				else
				{
					if ( std::fabs(U0) > 0.0 ){
									
					//zetaInt = zeta0p+0.5*dt*((Wn-UetaDx)+WnOld_);
					//phiInt = (1.5*((gVal & zeta0p)+(0.5)*turgut)-0.5*DPhiold_)*dt + Phi0Patch;
					
					Info << "Current speed is nonzero " << endl;
					zetap=zeta0p+dt*((23.0/12.0)*(Wn+Wcurdz_zeta0p-UetaDx)-(16.0/12.0)*WnOld_+(5.0/12.0)*WnOld2_);
					phiCalc = ((23.0/12.0)*((gVal & zeta0p)+turgut + dampingterm)-(16.0/12.0)*DPhiold_+(5.0/12.0)*DPhiold2_)*dt + Phi0Patch;
					
					
					WnOld2_=WnOld_;
					DPhiold2_=DPhiold_;
					
					
					WnOld_=(Wn+Wcurdz_zeta0p-UetaDx);
					DPhiold_=((gVal & zeta0p)+turgut + dampingterm);
					
					}
					else{
					
					
					Info << "Current speed U0 is zero " << endl;
					zetap=zeta0p+dt*((23.0/12.0)*Wn-(16.0/12.0)*WnOld_+(5.0/12.0)*WnOld2_);
					phiCalc = ((23.0/12.0)*((gVal & zeta0p) + dampingterm)-(16.0/12.0)*DPhiold_+(5.0/12.0)*DPhiold2_)*dt + Phi0Patch;
					
					WnOld2_=WnOld_;
					DPhiold2_=DPhiold_;
					WnOld_=Wn;
					DPhiold_=((gVal & zeta0p) + dampingterm);
					
					} 
					 
					 
					 
				 }		
				
					 
					
				
				break;
			}
			default:
			{
				FatalErrorInFunction
					<< ddtSchemeName << nl
					<< "    on patch " << this->patch().name()
					<< " of field " << this->internalField().name()
					<< " in file " << this->internalField().objectPath()
					<< abort(FatalError);
			}
		}

		

		//Info<< "min/max turgut zetap = " << gMin(zetap & nf()) << ", " << gMax(zetap & nf()) << endl;
			

	Info<< "min/max turgut zetap = " << gMin(zetap & nfRef) << ", "
			<< gMax(zetap & nfRef) << endl;

    
	
	lastUpdateTimeIndex=db().time().timeIndex();
	Info << "lastUpdateTimeIndex:" << lastUpdateTimeIndex << endl;
	
	//BUNLAR PARALLELDE CALISMIYOR
	//Info << "zetap[1]:" << zetap[1] << endl;
	if (gMax(zetap & nfRef) > 100.0)
  {
      FatalErrorInFunction
          << "Surface elevation (zetap) exceeded limit of 100 at time = "
          << db().time().value() << nl
          << "  max(zetap) = " << gMax(zetap & nfRef) << nl
          << abort(FatalError);
  }
	
	// Now assign the computed values using operator==
	data_=phiCalc;
	operator==(phiCalc);
	
	
	
	
	
	
	fixedValueFvPatchScalarField::updateCoeffs();
}


#include "InterpolationsHelpers.H"
#include "2nd_UpwindV6_MQLEAST.H" //numerical scheme
#include "PreParV16D.H"  // neighbours upwind down wind , scheme detection





void Foam::waveCurrentPotential3DFvPatchScalarField::ensureDictLoaded() const
{
    if (!waveCurDictPtr_.valid())
    {
        waveCurDictPtr_.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    "waveCurConditions",
                    db().time().constant(),
                    db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}

void Foam::waveCurrentPotential3DFvPatchScalarField::readParamsFrom(const dictionary& d)
{

    params_.U0         = readScalar(d.lookup("currentspeed"));
    params_.steepness  = readScalar(d.lookup("steepness"));
    params_.wavelength = readScalar(d.lookup("wavelength"));
    params_.hdepth     = readScalar(d.lookup("waterdepth"));
    params_.rampperiod = readScalar(d.lookup("rampperiod"));
    params_.v0         = d.lookupOrDefault<scalar>("v0", 0.0);

	// damping settings
	params_.xdamp      = readScalar(d.lookup("xdamp"));
	params_.Lxdamp      = readScalar(d.lookup("Lxdamp"));
	params_.ydamp      = readScalar(d.lookup("ydamp"));
	params_.Lydamp      = readScalar(d.lookup("Lydamp"));

	params_.xsponge    = readScalar(d.lookup("xsponge"));
	params_.Lsponge    = readScalar(d.lookup("Lsponge"));


    // quick sanity checks (keep minimal)
    if (params_.wavelength <= 0 || params_.hdepth <= 0 || params_.rampperiod <= 0)
        FatalIOErrorInFunction(d) << "Invalid waveCurConditions values." << exit(FatalIOError);

    params_.init = true;
}


void Foam::waveCurrentPotential3DFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("zeta", "zeta", zetaName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
	os.writeEntryIfDifferent<scalar>("shape", 100.0, shape_);
    fvPatchField<scalar>::writeValueEntry(os);
}





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        waveCurrentPotential3DFvPatchScalarField
    );
}

// ************************************************************************* //
