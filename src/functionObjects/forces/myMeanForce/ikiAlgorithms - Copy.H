
void Foam::functionObjects::myMeanForce::findFreeSurfaceFaces()
{
    
	if (FSpatchesfound_) return;
	
	
	const fvMesh& mesh = mesh_;
    const labelList& faceOwners = mesh.faceOwner();

    // One list per patch (like before)
    bodyToFreeSurfaceFaces_.setSize(patchIDs_.size());

    forAll(patchIDs_, patchIndex)
    {
        const label patchi = patchIDs_[patchIndex];
        const fvPatch& patch = mesh.boundary()[patchi];
        const label startFace = patch.start();
        const label nFaces = patch.size();

        DynamicList<labelPair> foundPairs; // (localBodyFaceI, globalFreeSurfaceFaceID)

        for (label i = 0; i < nFaces; ++i)
        {
            const label globalBodyFace = startFace + i;
            const label ownerCell = faceOwners[globalBodyFace];
            const cell& ownerCellFaces = mesh.cells()[ownerCell];

            forAll(ownerCellFaces, fJ)
            {
                const label neighFaceID = ownerCellFaces[fJ];
                const label neighPatchI = mesh.boundaryMesh().whichPatch(neighFaceID);

                if (neighPatchI != -1)
                {
                    const polyPatch& neighPatch = mesh.boundaryMesh()[neighPatchI];
                    const word& neighPatchName = neighPatch.name();

                    if (neighPatchName == "freeSurface" || neighPatchName == "upper")
                    {
                        // store (local index, global index)
                        foundPairs.append(labelPair(i, neighFaceID));
                        break;
                    }
                }
            }
        }

         bodyToFreeSurfaceFaces_[patchIndex].transfer(foundPairs);
		if ( Pstream::myProcNo() == 0) Pout << "Proc " << Pstream::myProcNo() << " bodyToFreeSurfaceFaces_ " << bodyToFreeSurfaceFaces_[patchIndex] << endl; 
       
    }
	
	FSpatchesfound_ = true;
	
	
}

	
	/*---------------------------------------------------------------------------*\
	if (debug)
        {
            const word& patchName = mesh.boundary()[patchi].name();
            Pout << "Proc " << Pstream::myProcNo()
                 << " patch " << patchName
                 << " found " << bodyToFreeSurfaceFaces_[patchIndex].size()
                 << " faces sharing owner with freeSurface/upper patch." << endl;

            forAll(bodyToFreeSurfaceFaces_[patchIndex], k)
            {
                const labelPair& pair = bodyToFreeSurfaceFaces_[patchIndex][k];
                Pout << "   bodyLocalI = " << pair.first()
                     << "   freeSurfaceGlobalID = " << pair.second()
                     << endl;
            }
        }
	
	Pout << "Proc " << Pstream::myProcNo() << " bodyToFreeSurfaceFaces_ " << bodyToFreeSurfaceFaces_[patchIndex] << endl; 
	 
	const auto& Cb = mesh_.C().boundaryField();

	forAll(patchIDs_, patchIndex)
	{
		const label patchi = patchIDs_[patchIndex];
		const vectorField& Cpatch = Cb[patchi];  // face centres on this patch

		forAll(ownerHasFreeSurfaceFace_[patchIndex], i)
		{
			if (ownerHasFreeSurfaceFace_[patchIndex][i])
			{
				const vector& Cf = Cpatch[i];  // specific face centre

				Pout << "Proc " << Pstream::myProcNo() << " | Patch: " << mesh_.boundary()[patchi].name()<< " | Face: " << i << " | Centre: (" << Cf.x() << ", " << Cf.y() << ", " << Cf.z() << ")"<< " touches freeSurface/upper" << endl;
			}
		}
	}
	\*---------------------------------------------------------------------------*/

void Foam::functionObjects::myMeanForce::readZetaValues()
{
    const fvMesh& mesh = mesh_;

    // --- Safety check ---
    if (bodyToFreeSurfaceFaces_.empty())
    {
        FatalErrorInFunction
            << "bodyToFreeSurfaceFaces_ is empty. Run findFreeSurfaceFaces() first."
            << exit(FatalError);
    }

	if (!obr_.foundObject<volVectorField>("zeta"))
	{
		FatalErrorInFunction
			<< "Field 'zeta' not found as volVectorField in registry '"
			<< obr_.name() << "'. Available objects: " << nl
			<< obr_.sortedNames() << nl
			<< exit(FatalError);
	}


	// Lookup the volume vector field
	const volVectorField& zeta = obr_.lookupObject<volVectorField>("zeta");

	// Access its boundary field
	const volVectorField::Boundary& zetaB = zeta.boundaryField();



	
   

    // Prepare storage
    zetaAtBodyFaces_.setSize(patchIDs_.size());

    forAll(patchIDs_, patchIndex)
    {
       // const label patchi = patchIDs_[patchIndex];
        //const fvPatch& patch = mesh.boundary()[patchi];

        const List<labelPair>& facePairs = bodyToFreeSurfaceFaces_[patchIndex];
        List<scalar> zetaVals(facePairs.size(), Zero);

        forAll(facePairs, k)
        {
            const label freeSurfFaceID = facePairs[k].second();
			const label fsPatchI = mesh.boundaryMesh().whichPatch(freeSurfFaceID);
			if (fsPatchI == -1)
			{
				WarningInFunction
					<< "Free-surface face " << freeSurfFaceID
					<< " not on a boundary patch. Skipping." << endl;
				continue;
			}

			const fvPatch& fsPatch = mesh.boundary()[fsPatchI];
			const label localFaceI = freeSurfFaceID - fsPatch.start();

			if (localFaceI < 0 || localFaceI >= mesh.boundary()[fsPatchI].size())
			{
				WarningInFunction
					<< "Invalid localFaceI " << localFaceI << " for patch " << fsPatchI << endl;
				continue;
			}

			// Get zeta.z at the free-surface face
			const vector& zetaVec = zetaB[fsPatchI][localFaceI];
			zetaVals[k] = zetaVec.z();

			// --- NEW: Get owner cell center of the body face
			//const label bodyPatchFaceI = facePairs[k].first();
			//const label patchi = patchIDs_[patchIndex];
			//const label globalBodyFaceID = mesh.boundary()[patchi].start() + bodyPatchFaceI;
			//const label ownerCell = mesh.faceOwner()[globalBodyFaceID];
			//const vector& cellCenter = mesh.C()[ownerCell];

			// Debug print
			//if ( Pstream::myProcNo() == 0) Pout << "k = " << k  << " | freeSurfaceFaceID = " << freeSurfFaceID  << " | zeta.z = " << zetaVec.z() << " | bodyFace center = " << mesh.Cf()[globalBodyFaceID] << " | ownerCell center = " << cellCenter << endl;
				
				// --- NEW: Get owner cell center of the free-surface face
			const label fsOwnerCell = mesh.faceOwner()[freeSurfFaceID];
			const vector& fsOwnerCenter = mesh.C()[fsOwnerCell];

			// Debug print
			if (Pstream::myProcNo() == 0)
			{
				Pout << "k = " << k  << " | freeSurfaceFaceID = " << freeSurfFaceID  << " | zeta.z = " << zetaVec.z()  << " | freeSurfaceFace center = " << mesh.Cf()[freeSurfFaceID] << " | ownerCell center (of freeSurface face) = " << fsOwnerCenter << endl;
					
	 
			}
		
				 
				 
				 

        }

        // Store per patch
        zetaAtBodyFaces_[patchIndex].transfer(zetaVals);
		if ( Pstream::myProcNo() == 0) Pout << "Proc " << Pstream::myProcNo() << " zetaAtBodyFaces_ " << zetaAtBodyFaces_[patchIndex] << endl; 
     
    }
}


