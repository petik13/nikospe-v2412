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
    Foam::waveCurMqLeastAB3DPotUPFD5FvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This is a pressure boundary condition, whose value is calculated as
    the hydrostatic pressure based on a given displacement:

        \f[
            p = -\rho*g*\zeta
        \f]

    \vartable
        \rho  | density [kg/m3]
        g     | acceleration due to gravity [m/s2]
        \zeta | wave amplitude [m]
    \endvartable

    The wave amplitude is updated as part of the calculation, derived from the
    local volumetric flux.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        zeta         | wave amplitude field name | no        | zeta
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            waveSurfacePressure;
        phi             phi;
        rho             rho;
        zeta            zeta;
        value           uniform 0;  // place holder
    }
    \endverbatim

    The density field is only required if the flux is mass-based as opposed to
    volumetric-based.

See also
    Foam::fixedValueFvPatchField

SourceFiles
    waveCurMqLeastAB3DPotUPFD5FvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef waveCurMqLeastAB3DPotUPFD5FvPatchScalarField_H
#define waveCurMqLeastAB3DPotUPFD5FvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Enum.H"
// Enum is a wrapper around a list of names/values that represent particular
    // enumeration (or int) values.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class waveCurMqLeastAB3DPotUPFD5FvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class waveCurMqLeastAB3DPotUPFD5FvPatchScalarField
:
    public fixedValueFvPatchScalarField // inherits from fixedValueFvPatchScalarField
{
public:

    // Public data

        //- Enumeration defining the available ddt schemes
        enum ddtSchemeType
        {
            tsEuler,
            tsCrankNicolson,
            tsBackward
        };


private:

    // Private data
		
        //- Flux field name
        word phiName_;

        //- Wave height field name
        word zetaName_;

        //- Density field for mass-based flux evaluations
        word rhoName_;
		
		// Private data
		word PHIName_;
		
		word UName_;
		
		//newly added
		//static constexpr scalar defaultShapeValue = 100.0;
		scalar shape_;
		

        //- Time scheme type names
        static const Enum<ddtSchemeType> ddtSchemeTypeNames_;
		

		scalarField data_;
		
		scalarField DPhiold_;
		vectorField WnOld_;
		
		scalarField DPhiold2_;
		vectorField WnOld2_;
		
		
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
		
		int lastUpdateTimeIndex = -1;
		
		// Flag indicating if the matrix has been computed already
		 bool neigboursCalculated_{false};
		 bool windsCalculated_{false};
		 bool schemesDetected_{false};
		 bool candidatesCalculated_{false};

public:

    //- Runtime type information
    TypeName("waveCurMqLeastAB3DPotUPFD5");

	
	mutable autoPtr<IOdictionary> waveCurDictPtr_;
		struct Params
		{
			scalar U0{0}, steepness{0}, wavelength{0}, hdepth{0}, v0{0}, Ls{1}, xl{100}, rampperiod{0};
			bool init{false};
		} mutable params_;

		void readParamsFrom(const dictionary& d);   // parses d -> params_
		void ensureDictLoaded() const;              // constructs waveCurDictPtr_ if needed
    // Constructors

        //- Construct from patch and internal field
        waveCurMqLeastAB3DPotUPFD5FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveCurMqLeastAB3DPotUPFD5FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given waveCurMqLeastAB3DPotUPFD5FvPatchScalarField
        //  onto a new patch
        waveCurMqLeastAB3DPotUPFD5FvPatchScalarField
        (
            const waveCurMqLeastAB3DPotUPFD5FvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        waveCurMqLeastAB3DPotUPFD5FvPatchScalarField
        (
            const waveCurMqLeastAB3DPotUPFD5FvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new waveCurMqLeastAB3DPotUPFD5FvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveCurMqLeastAB3DPotUPFD5FvPatchScalarField
        (
            const waveCurMqLeastAB3DPotUPFD5FvPatchScalarField&,
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
                new waveCurMqLeastAB3DPotUPFD5FvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
			virtual void calcNeigboursV3();
			
			
			virtual void UPFDV2();
			void exchangeWpGhostField(const scalarField& Wp,const scalarField& Wp2);
			
			virtual void findUpwindDownwindNodesV2();
			void detectFDSchemes(); 
			
			void findInletFaces();
			void findSphereFaces();
			void initCandidates();
			void buildVertexBasedConnections();
			
			void outputter();
			//- Write
			virtual void write(Ostream&) const;

	
	
		//- Safe getters for accessing internal fields
		inline const labelList& globalIdOfLocalFace() const{return globalIdOfLocalFace_;}

		inline const Map<label>& globalIdToPackedIdx() const{return globalIdToPackedIdx_;}

		inline const Map<vector>& globalIdToCf() const{return globalIdToCf_;}
		

		// Accessors for interpolation flags
		List<FixedList<bool,2>>& interpNeededX() { return interpNeededX_; }
		const List<FixedList<bool,2>>& interpNeededX() const { return interpNeededX_; }

		List<FixedList<bool,2>>& interpNeededY() { return interpNeededY_; }
		const List<FixedList<bool,2>>& interpNeededY() const { return interpNeededY_; }
		
		// Allow inline template helpers access to private members
		template<class PatchField>
		friend bool tryUpwind2nd(const PatchField&, label, label&, const vectorField&, const scalarField&, const scalarField&, const List<List<label>>&, const HashTable<List<label>, label>&);

		template<class PatchField>
		friend bool tryCentered2nd(const PatchField&, label, label&, const vectorField&, const scalarField&, const scalarField&, const List<List<label>>&, const List<List<label>>&);

		template<class PatchField>
		friend bool tryDownwind2nd(const PatchField&, label, label&, const vectorField&, const scalarField&, const scalarField&, const List<List<label>>&, const HashTable<List<label>, label>&);

		

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
