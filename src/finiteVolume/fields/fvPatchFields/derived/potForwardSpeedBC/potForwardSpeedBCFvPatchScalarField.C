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

\*---------------------------------------------------------------------------*/

#include "potForwardSpeedBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "fixedGradientFvPatchFields.H"
#include "gravityMeshObject.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "processorFvPatch.H"
#include "mathematicalConstants.H"

// Free helper functions used by the stencil machinery
#include "waveCurPar3DPotUPFD5InlineHelpersInt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::potForwardSpeedBCFvPatchScalarField::ddtSchemeType
>
Foam::potForwardSpeedBCFvPatchScalarField::ddtSchemeTypeNames_
({
    { ddtSchemeType::tsEuler,         fv::EulerDdtScheme<scalar>::typeName_()         },
    { ddtSchemeType::tsCrankNicolson, fv::CrankNicolsonDdtScheme<scalar>::typeName_() },
    { ddtSchemeType::tsBackward,      fv::backwardDdtScheme<scalar>::typeName_()      }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::potForwardSpeedBCFvPatchScalarField::
potForwardSpeedBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    PhiDName_("PhiD"),
    zetaDName_("zetaD"),
    bodyPatchName_("sphere"),
    shape_(100.0),
    lastTimeIndex_(-1)
{}


Foam::potForwardSpeedBCFvPatchScalarField::
potForwardSpeedBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    PhiDName_(dict.getOrDefault<word>("PhiD", "PhiD")),
    zetaDName_(dict.getOrDefault<word>("zetaD", "zetaD")),
    bodyPatchName_(dict.getOrDefault<word>("bodyPatch", "sphere")),
    shape_(dict.getOrDefault<scalar>("shape", 100.0)),
    lastTimeIndex_(-1)
{
    ensureDictLoaded();
    readParamsFrom(*waveDictPtr_);
}


Foam::potForwardSpeedBCFvPatchScalarField::
potForwardSpeedBCFvPatchScalarField
(
    const potForwardSpeedBCFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    PhiDName_(ptf.PhiDName_),
    zetaDName_(ptf.zetaDName_),
    bodyPatchName_(ptf.bodyPatchName_),
    shape_(ptf.shape_),
    lastTimeIndex_(-1)
{}


Foam::potForwardSpeedBCFvPatchScalarField::
potForwardSpeedBCFvPatchScalarField
(
    const potForwardSpeedBCFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    PhiDName_(ptf.PhiDName_),
    zetaDName_(ptf.zetaDName_),
    bodyPatchName_(ptf.bodyPatchName_),
    shape_(ptf.shape_),
    lastTimeIndex_(-1)
{}


Foam::potForwardSpeedBCFvPatchScalarField::
potForwardSpeedBCFvPatchScalarField
(
    const potForwardSpeedBCFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    PhiDName_(ptf.PhiDName_),
    zetaDName_(ptf.zetaDName_),
    bodyPatchName_(ptf.bodyPatchName_),
    shape_(ptf.shape_),
    lastTimeIndex_(-1)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::potForwardSpeedBCFvPatchScalarField::ensureDictLoaded() const
{
    if (!waveDictPtr_)
    {
        waveDictPtr_.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    "waveConditions",
                    db().time().constant(),
                    db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


void Foam::potForwardSpeedBCFvPatchScalarField::readParamsFrom
(
    const dictionary& dict
)
{
    params_.U0         = dict.get<scalar>("currentSpeed");
    params_.head_ang   = dict.get<scalar>("headingAngle");
    params_.steepness  = dict.get<scalar>("steepness");
    params_.wavelength = dict.get<scalar>("waveLength");
    params_.hdepth     = dict.get<scalar>("waterDepth");
    params_.rampperiod = dict.get<scalar>("rampPeriods");

    params_.v0      = dict.getOrDefault<scalar>("spongeStrength", 0.0);
    params_.xsponge = dict.getOrDefault<scalar>("xInlet",  0.0);
    params_.Lsponge = dict.getOrDefault<scalar>("LInlet",  1.0);
    params_.xdamp   = dict.getOrDefault<scalar>("xOutlet", GREAT);
    params_.Lxdamp  = dict.getOrDefault<scalar>("LOutlet", 1.0);
    params_.ydamp   = dict.getOrDefault<scalar>("ySide",   GREAT);
    params_.Lydamp  = dict.getOrDefault<scalar>("LSide",   1.0);

    if (params_.wavelength <= SMALL || params_.hdepth <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "waveLength and waterDepth must be positive"
            << exit(FatalIOError);
    }

    // Derived: deep/finite-depth dispersion, then the Doppler shift into the
    // ship-fixed frame.  Uinf is the steady stream seen by the body.
    const scalar g = mag(meshObjects::gravity::New(db().time()).value());

    params_.amp    = 0.5*params_.steepness*params_.wavelength;
    params_.k      = constant::mathematical::twoPi/params_.wavelength;
    params_.omega  = Foam::sqrt(g*params_.k*Foam::tanh(params_.k*params_.hdepth));
    params_.Uinf   = params_.U0
                    *vector(Foam::cos(params_.head_ang), Foam::sin(params_.head_ang), 0);
    params_.omegaE = params_.omega + params_.k*params_.Uinf.x();

    params_.rampTime = params_.rampperiod*constant::mathematical::twoPi/params_.omega;
    params_.init = true;
}


Foam::scalar Foam::potForwardSpeedBCFvPatchScalarField::ramp
(
    const scalar t
) const
{
    if (params_.rampTime <= SMALL || t >= params_.rampTime)
    {
        return 1.0;
    }

    return 0.5*(1.0 - Foam::cos(constant::mathematical::pi*t/params_.rampTime));
}


Foam::tmp<Foam::scalarField>
Foam::potForwardSpeedBCFvPatchScalarField::spongeCoeff() const
{
    const scalarField x(patch().Cf().component(vector::X));
    const scalarField y(patch().Cf().component(vector::Y));

    auto tnu = tmp<scalarField>::New(patch().size(), Zero);
    scalarField& nu = tnu.ref();

    if (params_.v0 <= SMALL)
    {
        return tnu;
    }

    // Quadratic ramps: inlet reflections upstream, radiated waves downstream
    // and to the sides.  The side zone is switched off inside the outlet zone
    // so the two do not compound in the corners.
    forAll(nu, i)
    {
        if (x[i] < params_.xsponge)
        {
            nu[i] += params_.v0*sqr((params_.xsponge - x[i])/params_.Lsponge);
        }
        else if (x[i] > params_.xdamp)
        {
            nu[i] += params_.v0*sqr((x[i] - params_.xdamp)/params_.Lxdamp);
        }
        else if (mag(y[i]) > params_.ydamp)
        {
            nu[i] += params_.v0*sqr((mag(y[i]) - params_.ydamp)/params_.Lydamp);
        }
    }

    return tnu;
}


void Foam::potForwardSpeedBCFvPatchScalarField::freeSurfaceRates
(
    const scalar t,
    const scalarField& zetaD,
    scalarField& dZetaDdt,
    scalarField& dPhiDdt
)
{
    const label patchi = patch().index();
    const label nFaces = patch().size();

    const scalar g = mag(meshObjects::gravity::New(db().time()).value());
    const scalar r = ramp(t);

    const scalarField x(patch().Cf().component(vector::X));

    // Vertical velocity of the disturbance.  The patch normal is +z, so with
    // u = +grad(Phi) this is simply the surface-normal gradient.
    const volScalarField& PhiD = db().lookupObject<volScalarField>(PhiDName_);
    const scalarField wD(PhiD.boundaryField()[patchi].snGrad());

    // Steady basis flow and the vertical straining dWz/dz it imposes on the
    // free surface.  Wz vanishes on z = 0 for a double-body flow, so only the
    // horizontal components of W enter the convective terms.
    
    const volVectorField& Us = db().lookupObject<volVectorField>("Us");
    const vectorField& W = Us.boundaryField()[patchi];

    const volScalarField& dUsdzF = db().lookupObject<volScalarField>("dUsdz");
    const scalarField& dWzdz = dUsdzF.boundaryField()[patchi];

    dZetaDdt.setSize(nFaces);
    dPhiDdt.setSize(nFaces);

    const scalar A = params_.amp;
    const scalar k = params_.k;
    const scalar w0 = params_.omega;
    const scalar we = params_.omegaE;
    const scalar dUx = -params_.Uinf.x();   // (W - Uinf)_x is added per face
    const scalar dUy = -params_.Uinf.y();

    forAll(dZetaDdt, i)
    {
        // Incident wave at z = 0 (cosh/cosh = 1 there)
        const scalar theta = k*x[i] - we*t;
        const scalar cth = Foam::cos(theta);
        const scalar sth = Foam::sin(theta);

        const scalar zetaI    =  A*cth*r;          // elevation
        const scalar dZetaIdx = -A*k*sth*r;        // d(zeta_I)/dx
        const scalar uIx      =  A*k*g/w0*cth*r;   // d(phi_I)/dx
        const scalar uIy      =  0.0;              // wave runs along +x

        // Double-body disturbance velocity, W - Uinf: zero in the far field,
        // O(U) near the body.  This is what the incident wave does not see.
        const scalar dWx = W[i].x() + dUx;
        const scalar dWy = W[i].y() + dUy;

        // Kinematic condition
        dZetaDdt[i] =
            wD[i]
          + (zetaD[i] + zetaI)*dWzdz[i]
          - (W[i].x()*zetaDx_[i].x() + W[i].y()*zetaDx_[i].y())
          - (dWx*dZetaIdx);

        // Dynamic condition
        dPhiDdt[i] =
          -  g*zetaD[i]
          - (W[i].x()*PhiDx_[i] + W[i].y()*PhiDy_[i])
          - (dWx*uIx + dWy*uIy);
    }

    // Sponge: -nu*phi_D,z damps the disturbance without touching the incident
    // wave, which is prescribed analytically.
    const tmp<scalarField> tnu = spongeCoeff();
    dPhiDdt -= tnu()*wD;
}


// Mark free-surface faces sharing a vertex or an edge with the body.  The
// finite-difference stencils fall back to a lower order on these, where the
// regular neighbour ring is incomplete.
void Foam::potForwardSpeedBCFvPatchScalarField::findSphereEdgeVertexFaces()
{

	if (sphereEdgeVertexFacesCalculated_) return;
	sphereEdgeVertexFacesCalculated_ = true;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const label nFaces   = patch().size();
    const label startFS  = patch().start();
    ownerHasBodyEdge_.setSize(nFaces, false);
    ownerHasBodyVertex_.setSize(nFaces, false);

    const label spherePatchID = mesh.boundaryMesh().findPatchID(bodyPatchName_);
    if (spherePatchID == -1)
    {
        FatalErrorInFunction << "Patch 'sphere' not found." << abort(FatalError);
    }


    const polyPatch& spherePatch = mesh.boundaryMesh()[spherePatchID];

    // 1) Collect all point labels used by local 'sphere' faces on this proc
    labelHashSet spherePts;
    spherePts.reserve(4*spherePatch.size()); // heuristic

    forAll(spherePatch, sfI)
    {
        const label faceID = spherePatch.start() + sfI;
        const face& f = mesh.faces()[faceID];
        forAll(f, vI)
        {
            spherePts.insert(f[vI]);
        }
    }

    // 2) For each free-surface face, count shared point labels with spherePts
    for (label i = 0; i < nFaces; ++i)
    {
        const label faceID = startFS + i;
        const face& f = mesh.faces()[faceID];
        label shared = 0;
        forAll(f, vI)
        {
            if (spherePts.found(f[vI]))
            {
                ++shared;
                if (shared >= 2) break; // early exit: edge-touch satisfied
            }
        }

        if (shared >= 1) ownerHasBodyVertex_[i] = true;
        if (shared >= 2) ownerHasBodyEdge_[i]   = true;
    }

    // Optional: merge into existing mask if you want
    // for (label i=0;i<nFaces;++i) ownerHasBodyFace_[i] = ownerHasBodyFace_[i] || ownerHasBodyEdge_[i];
}


void Foam::potForwardSpeedBCFvPatchScalarField::applyKreissOligerFilter
(
    scalarField& zetaD,
    const scalarField& zetaDOld
)
{
    // Fourth-difference hyperviscosity, damping only the grid-scale mode that
    // the 4th-order central stencil admits.  Applied where that stencil is in
    // use (scheme code 14); elsewhere the scheme is already dissipative.
    const scalar epsilon = 0.0005;

    auto valueAt = [&](const label gid) -> scalar
    {
        if (globalIdToPackedIdx_.found(gid))
        {
            return zetaDOld[globalIdToPackedIdx_[gid]];
        }
        if (globalIdToWp_.found(gid))
        {
            return globalIdToWp_[gid];
        }

        FatalErrorInFunction
            << "No zeta value for global face " << gid
            << ": the stencil was not exchanged." << abort(FatalError);
        return 0;
    };

    auto secondNeighbour = [&]
    (
        const label first,
        const List<List<label>>& localNodes,
        const HashTable<labelList, label>& remoteNodes,
        const label selfFace
    ) -> label
    {
        if (globalIdToPackedIdx_.found(first))
        {
            return localNodes[globalIdToPackedIdx_[first]][0];
        }

        const auto it = remoteNodes.find(globalIdOfLocalFace_[selfFace]);
        if (it != remoteNodes.end())
        {
            return it()[0];
        }

        FatalErrorInFunction
            << "Missing second neighbour for global face "
            << globalIdOfLocalFace_[selfFace] << abort(FatalError);
        return -1;
    };

    forAll(zetaD, i)
    {
        for (label dir = 0; dir < 2; ++dir)
        {
            const bool xDir = (dir == 0);

            if ((xDir ? schemeCodeX_[i] : schemeCodeY_[i]) != 14)
            {
                continue;
            }

            const List<List<label>>& up   = xDir ? upwindNodesX_   : upwindNodesY_;
            const List<List<label>>& down = xDir ? downwindNodesX_ : downwindNodesY_;
            const auto& remUp   = xDir ? remoteIdToUpwindNodesX_   : remoteIdToUpwindNodesY_;
            const auto& remDown = xDir ? remoteIdToDownwindNodesX_ : remoteIdToDownwindNodesY_;

            if (up[i].empty() || down[i].empty())
            {
                continue;
            }

            const label u1 = up[i][0];
            const label d1 = down[i][0];
            const label u2 = secondNeighbour(u1, up, remUp, i);
            const label d2 = secondNeighbour(d1, down, remDown, i);

            zetaD[i] -= epsilon
               *(
                    valueAt(u2) - 4*valueAt(u1) + 6*zetaDOld[i]
                  - 4*valueAt(d1) + valueAt(d2)
                );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::potForwardSpeedBCFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // updateCoeffs() is called once per non-orthogonal corrector.  Only the
    // first call in a time step advances the free surface; the rest replay it.
    if (db().time().timeIndex() == lastTimeIndex_)
    {
        operator==(cachedPhiD_);
        fixedValueFvPatchScalarField::updateCoeffs();
        return;
    }

    ensureDictLoaded();
    if (waveDictPtr_->readIfModified() || !params_.init)
    {
        readParamsFrom(*waveDictPtr_);
    }

    const label patchi = patch().index();
    const scalar dt = db().time().deltaTValue();

    // The right-hand side is explicit, so every term is taken at the old time
    // level, including the analytic incident wave.
    const scalar tOld = db().time().value() - dt;

    volScalarField& PhiD = db().lookupObjectRef<volScalarField>(PhiDName_);
    volScalarField& zetaD = db().lookupObjectRef<volScalarField>(zetaDName_);

    scalarField& zetaDp = zetaD.boundaryFieldRef()[patchi];
    const scalarField PhiDOld(PhiD.oldTime().boundaryField()[patchi]);
    const scalarField zetaDOld(zetaD.oldTime().boundaryField()[patchi]);

    // Horizontal derivatives of phi_D and zeta_D on the patch.  The stencils,
    // their parallel halo and the near-body fallbacks are built by the
    // finite-difference machinery below; only the first call does the setup.
    zetaDx_.setSize(patch().size(), Zero);
    PhiDx_.setSize(patch().size(), Zero);
    PhiDy_.setSize(patch().size(), Zero);

    findSphereEdgeVertexFaces();
    calcNeigboursV3();
    buildEdgeNeighbours();

    if (neigboursCalculated_)
    {
        findUpwindDownwindNodesV2();
        detectFDSchemes();
        UPFDV2();
    }

    // Rates at the old level
    scalarField dZetaDdt, dPhiDdt;
    freeSurfaceRates(tOld, zetaDOld, dZetaDdt, dPhiDdt);

    // --- Time integration -------------------------------------------------
    const word ddtName(zetaD.mesh().ddtScheme(zetaD.name()));
    const ddtSchemeType scheme = ddtSchemeTypeNames_[ddtName];

    scalarField phiCalc(patch().size(), Zero);

    if (scheme == tsBackward && dZetaDdt2_.size() == dZetaDdt.size())
    {
        // Adams-Bashforth 3
        zetaDp = zetaDOld
               + dt*((23.0/12.0)*dZetaDdt - (16.0/12.0)*dZetaDdt1_ + (5.0/12.0)*dZetaDdt2_);
        phiCalc = PhiDOld
               + dt*((23.0/12.0)*dPhiDdt - (16.0/12.0)*dPhiDdt1_ + (5.0/12.0)*dPhiDdt2_);
    }
    else if (scheme == tsBackward && dZetaDdt1_.size() == dZetaDdt.size())
    {
        // Adams-Bashforth 2 while the third level is still being filled
        zetaDp = zetaDOld + dt*(1.5*dZetaDdt - 0.5*dZetaDdt1_);
        phiCalc = PhiDOld + dt*(1.5*dPhiDdt - 0.5*dPhiDdt1_);
    }
    else
    {
        zetaDp = zetaDOld + dt*dZetaDdt;
        phiCalc = PhiDOld + dt*dPhiDdt;
    }

    // Shift the history
    dZetaDdt2_ = dZetaDdt1_;  dZetaDdt1_ = dZetaDdt;
    dPhiDdt2_  = dPhiDdt1_;   dPhiDdt1_  = dPhiDdt;

    // --- Filter -----------------------------------------------------------
    if (neigboursCalculated_)
    {
        const scalarField zetaBeforeFilter(zetaDp);
        exchangeWpGhostField(zetaBeforeFilter, phiCalc);
        applyKreissOligerFilter(zetaDp, zetaBeforeFilter);
    }

    if (gMax(zetaDp) > 100.0)
    {
        FatalErrorInFunction
            << "Free-surface elevation exceeded 100 m at t = "
            << db().time().value() << ": max(zeta_D) = " << gMax(zetaDp)
            << abort(FatalError);
    }

    Info<< "    zeta_D min/max = " << gMin(zetaDp) << ", " << gMax(zetaDp) << endl;

    lastTimeIndex_ = db().time().timeIndex();
    cachedPhiD_ = phiCalc;

    operator==(phiCalc);
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::potForwardSpeedBCFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("PhiD", "PhiD", PhiDName_);
    os.writeEntryIfDifferent<word>("zetaD", "zetaD", zetaDName_);
    os.writeEntryIfDifferent<scalar>("shape", 100.0, shape_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// Finite-difference machinery: stencil construction, parallel halo exchange,
// scheme detection and the upwind/central derivative evaluation.  Carried over
// unchanged apart from the field renames and the u = +grad(Phi) sign.
// These files use unqualified OpenFOAM type names at file scope.
using namespace Foam;

#include "InterpolationsHelpers.H"
#include "2nd_UpwindV6_MQLEAST.H"
#include "PreParV18D.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        potForwardSpeedBCFvPatchScalarField
    );
}

// ************************************************************************* //
