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

Class
    Foam::potForwardSpeedBCFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    Linearised free-surface condition for the disturbance potential of a body
    advancing with forward speed in waves, on the double-body basis flow.

    Velocity convention is  u = +grad(Phi)  throughout.

    The unsteady potential and the elevation are split into a prescribed
    incident part and a computed disturbance part,

        phi = phi_I + phi_D ,        zeta = zeta_I + zeta_D

    with W = grad(PhiS) the steady double-body flow, tending to Uinf far away.
    The incident wave satisfies the free-surface conditions built on the
    UNIFORM stream Uinf rather than on the local W, so subtracting them leaves
    a forcing proportional to (W - Uinf) that vanishes in the far field:

        zeta_D,t = phi_D,z + (zeta_D + zeta_I) dWz/dz
                 - W_h.grad_h(zeta_D) - (W_h - Uinf).grad_h(zeta_I)

        phi_D,t  = -g zeta_D - W.grad(phi_D) - (W - Uinf).grad(phi_I)
                 + nu phi_D,z                                     (sponge)

    Both are advanced explicitly and phi_D is imposed as a fixed value.
    Horizontal derivatives on the patch come from the finite-difference
    machinery in the included headers, which carries its own parallel halo.

Usage
    \table
      Property | Description                        | Reqd | Default
      PhiD     | Disturbance potential field        | no   | PhiD
      zetaD    | Disturbance elevation field        | no   | zetaD
      shape    | Stencil shape parameter            | no   | 100
    \endtable

    Wave and sponge parameters are read from constant/waveConditions and
    re-read if that file changes.

SourceFiles
    potForwardSpeedBCFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_potForwardSpeedBCFvPatchScalarField_H
#define Foam_potForwardSpeedBCFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Enum.H"
#include "IOdictionary.H"
#include "globalIndex.H"
#include "HashTable.H"
#include "Map.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class potForwardSpeedBCFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class potForwardSpeedBCFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
public:

    //- Time-integration variants, selected from the ddtScheme of zetaD
    enum ddtSchemeType
    {
        tsEuler,
        tsCrankNicolson,
        tsBackward
    };


private:

    // Private Data - names

        word PhiDName_;
        word zetaDName_;
        word bodyPatchName_;
        scalar shape_;

        static const Enum<ddtSchemeType> ddtSchemeTypeNames_;


    // Private Data - explicit time-integration state

        //- Patch value already built at this time index.  updateCoeffs() is
        //  called once per non-orthogonal corrector; only the first call does
        //  work, the rest replay this.
        scalarField cachedPhiD_;
        label lastTimeIndex_;

        //- Right-hand sides at t-dt and t-2dt for Adams-Bashforth
        scalarField dPhiDdt1_, dPhiDdt2_;
        scalarField dZetaDdt1_, dZetaDdt2_;


		
		
		
		
		
		
		
		
		vectorField zetaDx_;
		
		scalarField PhiDx_;
		scalarField PhiDy_;
		
		//List<labelList> neighbours_;
		
		
		List<List<label>> upwindNodesX_;   ///< Stores multiple upwind nodes per face
		List<List<label>> downwindNodesX_; ///< Stores multiple downwind nodes per face
		
		List<List<label>> upwindNodesY_;   ///< Stores multiple upwind nodes per face
		List<List<label>> downwindNodesY_; ///< Stores multiple downwind nodes per face
		
		// Two flags per face: [0] = multiple direct neighbors, [1] = multiple second neighbors
		List<FixedList<bool,2>> interpNeededX_;
		List<FixedList<bool,2>> interpNeededY_;

        List<FixedList<bool,4>> interpNeeded4thX_;
        List<FixedList<bool,4>> interpNeeded4thY_;


		//- Guard so we only build the above once
		boolList procBoundaryFace_; 
		List<labelList> neighbProcNos_; //modified
		Foam::boolList ownerHasInletFace_;
		Foam::boolList ownerHasBodyFace_;
		Foam::boolList ownerHasFirstBodyFace_;
		Foam::Map<Foam::label> ringNodeX_;
		
		//Foam::boolList interpolationNeeded_;
		Map<List<label>> secondRingToFirstRing_;  // key: 2nd-ring local face idx, value: list of first-ring face local idx
		Map<label> ring2NodeX_;  // second-ring face → globalID of its upwind first-ring face
	
		
		Foam::autoPtr<Foam::globalIndex> gIndexPatch_;
		Foam::labelList globalIdOfLocalFace_;     // size = patch().size()
		//Foam::List<Foam::labelList> neighboursGlobal_;  // same shape as neighbours_, but GLOBAL ids
		Foam::Map<label> globalIdToPackedIdx_;
		Foam::Map<Foam::vector> globalIdToCf_;
		
		//Foam::HashTable<Foam::List<Foam::label>, Foam::label> remoteIdToUpwindNodesX_;
		//Foam::HashTable<Foam::List<Foam::label>, Foam::label> remoteIdToDownwindNodesX_;
		//Foam::HashTable<List<label>, label> remoteIdToUpwindNodesY_;     //BUNLAR YENI EKLENDI
		//Foam::HashTable<List<label>, label> remoteIdToDownwindNodesY_;   //BUNLAR YENI EKLENDI
		//HashTable<List<label>, label> remoteNeighbours_;
		
		Foam::HashTable<labelList, label> remoteIdToUpwindNodesX_;
		Foam::HashTable<labelList, label> remoteIdToDownwindNodesX_;
		Foam::HashTable<labelList, label> remoteIdToUpwindNodesY_;
		Foam::HashTable<labelList, label> remoteIdToDownwindNodesY_;
		Foam::HashTable<labelList, label> remoteNeighbours_;
		Foam::HashTable<labelList, label> remoteNeighbours2_;
		
		Foam::HashTable<labelList, label> remIdToUpwindNewX_;
		Foam::HashTable<labelList, label> remIdToDownwindNewX_;
		Foam::HashTable<labelList, label> remIdToUpwindNewY_;
		Foam::HashTable<labelList, label> remIdToDownwindNewY_;
		
		
	    labelHashSet sendTo; // 16D	
		Foam::List<labelList> neighbours_;
		Foam::List<labelList> neighboursGlobal_;
		
		List<List<label>> lsMergedCandidates_;
		//List<List<label>> lsMergedCandidatesY_; 
		
		
		Foam::labelHashSet neededGIDs_;
		Foam::List<Foam::labelHashSet> neededByProc_;
		
		
		Foam::HashTable<scalar, label> globalIdToWp_;
		Foam::HashTable<scalar, label> globalIdToPhi_;
		scalarField schemeCodeX_; 
		scalarField schemeCodeY_; 
		
		
		
		// Flag indicating if the matrix has been computed already
		 bool neigboursCalculated_{false};
         bool sphereEdgeVertexFacesCalculated_{false};
		 bool windsCalculated_{false};
		 bool schemesDetected_{false};
		 bool candidatesCalculated_{false};





    // Private Member Functions - free-surface machinery
    // (defined in the included headers, carried over unchanged; they own the
    //  parallel halo exchange for the patch finite-difference stencils)

        void calcNeigboursV3();
        void buildEdgeNeighbours();
        void buildVertexBasedConnections();
        void initCandidates();
        void findInletFaces();
        void findSphereFaces();
        void findSphereEdgeVertexFaces();
        void findUpwindDownwindNodesV2();
        void detectFDSchemes();
        void UPFDV2();
        void outputter();
        void exchangeWpGhostField(const scalarField& Wp, const scalarField& Wp2);

        void computeFirstDerivativeCoeffs5pt
        (
            const scalar xm2, const scalar xm1, const scalar x0,
            const scalar xp1, const scalar xp2,
            scalar& am2, scalar& am1, scalar& a0,
            scalar& ap1, scalar& ap2
        ) const;

        //- Kreiss-Oliger hyperviscosity, damps the grid-scale mode admitted by
        //  the 4th-order central stencil
        void applyKreissOligerFilter
        (
            scalarField& zetaD,
            const scalarField& zetaDOld
        );


    // Private Member Functions - physics

        //- Read constant/waveConditions into params_
        void readParamsFrom(const dictionary& d);
        void ensureDictLoaded() const;

        //- Ramp applied to the incident wave
        scalar ramp(const scalar t) const;

        //- Sponge coefficient nu(x,y), zero outside the damping zones
        tmp<scalarField> spongeCoeff() const;

        //- Assemble the two explicit right-hand sides at the old time level
        void freeSurfaceRates
        (
            const scalar t,
            const scalarField& zetaD,
            scalarField& dZetaDdt,
            scalarField& dPhiDdt
        );


public:

    //- Runtime type information
    TypeName("potForwardSpeedBC");

    //- Wave and sponge parameters, shared with the FD machinery
    mutable autoPtr<IOdictionary> waveDictPtr_;

    struct Params
    {
        scalar U0{0}, head_ang{0}, steepness{0}, wavelength{0}, hdepth{0},
               v0{0}, xdamp{0}, Lxdamp{0}, ydamp{0}, Lydamp{0},
               xsponge{0}, Lsponge{0}, rampperiod{0};

        // derived
        scalar amp{0}, k{0}, omega{0}, omegaE{0}, rampTime{0};
        vector Uinf{Zero};

        bool init{false};
    } mutable params_;


    // Constructors

        potForwardSpeedBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        potForwardSpeedBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        potForwardSpeedBCFvPatchScalarField
        (
            const potForwardSpeedBCFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        potForwardSpeedBCFvPatchScalarField
        (
            const potForwardSpeedBCFvPatchScalarField&
        );

        potForwardSpeedBCFvPatchScalarField
        (
            const potForwardSpeedBCFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new potForwardSpeedBCFvPatchScalarField(*this)
            );
        }

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new potForwardSpeedBCFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Advance the free surface and impose the new phi_D
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Access
    // Getters and friend declarations required by the inline stencil helpers,
    // which are free template functions and reach the state through these.

        inline const labelList& globalIdOfLocalFace() const
        { return globalIdOfLocalFace_; }

        inline const Map<label>& globalIdToPackedIdx() const
        { return globalIdToPackedIdx_; }

        inline const Map<vector>& globalIdToCf() const
        { return globalIdToCf_; }

        List<FixedList<bool,2>>& interpNeededX() { return interpNeededX_; }
        const List<FixedList<bool,2>>& interpNeededX() const { return interpNeededX_; }

        List<FixedList<bool,2>>& interpNeededY() { return interpNeededY_; }
        const List<FixedList<bool,2>>& interpNeededY() const { return interpNeededY_; }

        List<FixedList<bool,4>>& interpNeeded4thX() { return interpNeeded4thX_; }
        const List<FixedList<bool,4>>& interpNeeded4thX() const { return interpNeeded4thX_; }

        List<FixedList<bool,4>>& interpNeeded4thY() { return interpNeeded4thY_; }
        const List<FixedList<bool,4>>& interpNeeded4thY() const { return interpNeeded4thY_; }

        template<class PatchField>
        friend bool tryUpwind2nd(const PatchField&, label, label&, const vectorField&, const scalarField&, const scalarField&, const List<List<label>>&, const HashTable<List<label>, label>&);

        template<class PatchField>
        friend bool tryCentered2nd(const PatchField&, label, label&, const vectorField&, const scalarField&, const scalarField&, const List<List<label>>&, const List<List<label>>&);

        template<class PatchField>
        friend bool tryDownwind2nd(const PatchField&, label, label&, const vectorField&, const scalarField&, const scalarField&, const List<List<label>>&, const HashTable<List<label>, label>&);

        template<class PatchField>
        friend bool tryFourthOrderCentral(const PatchField&, label, label&, const vectorField&, const scalarField&, const scalarField&, const List<List<label>>&, const List<List<label>>&,
            const HashTable<List<label>, label>&, const HashTable<List<label>, label>&, const word&);

        //- Previous-step derivatives, used by the near-body fallback
        mutable scalarField DzetaZdxOld_, DzetaZdyOld_;
        mutable scalarField DPhidxOld_,   DPhidyOld_;
        mutable bool gradOldInit_ = false;

        //- Body-adjacency masks and edge stencils
        boolList ownerHasBodyEdge_;
        boolList ownerHasBodyVertex_;
        List<labelList> edgeNeighbours_;
        List<List<vector>> edgeNormals2D_;
};


} // End namespace Foam

#endif

// ************************************************************************* //
