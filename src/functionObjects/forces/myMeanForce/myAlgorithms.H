
void Foam::functionObjects::myMeanForce::findFreeSurfaceFaces()
{
    
	if (FSpatchesfound_) return;
	
	
	const fvMesh& mesh = mesh_;
    const labelList& faceOwners = mesh.faceOwner();

    // One list per patch (like before)
	bodyToFreeSurfaceFaces_.setSize(mesh.boundary().size());
	bodyToEdgeLengths_.setSize(mesh.boundary().size());
	
	for (const label patchi : patchIDs_)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        const label startFace = patch.start();
        const label nFaces = patch.size();

        DynamicList<labelPair> foundPairs; // (localBodyFaceI, globalFreeSurfaceFaceID)
		DynamicList<scalar> sharedEdgeLengths; // <-- NEW edge length
		
		
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
                        
						const face& bodyFace = mesh.faces()[globalBodyFace];
                        const face& fsFace   = mesh.faces()[neighFaceID];
                        const pointField& points = mesh.points();

                        // Collect common vertices
                        DynamicList<vector> commonVerts;
						forAll(bodyFace, vi)
						{
							const label vBody = bodyFace[vi]; // vertex index in mesh.points()
							forAll(fsFace, vj)
							{
								const label vFS = fsFace[vj];
								if (vBody == vFS) // same vertex index
								{
									commonVerts.append(points[vBody]);
								}
							}
						}

						Info << "commonVerts" << commonVerts << endl;
                        scalar edgeLength = 0.0;
						if (commonVerts.size() == 2)
						{
							edgeLength = mag(commonVerts[0] - commonVerts[1]);
						}
						else if (commonVerts.size() > 2)
						{
							// degenerate case (shared polygon line) – take first two
							edgeLength = mag(commonVerts[0] - commonVerts[1]);
							WarningInFunction << "Faces share more than one edge ("
											  << commonVerts.size() << " vertices in common)"
											  << " at patch " << patch.name()
											  << ", face " << i << endl;
						}
						else
						{
							WarningInFunction << "No valid shared edge between body face "
											  << globalBodyFace << " and free-surface face "
											  << neighFaceID << endl;
						}
						
						
						
						
						
						
						// store (local index, global index)
                        foundPairs.append(labelPair(i, neighFaceID));
						sharedEdgeLengths.append(edgeLength);
                        break;
                    }
                }
            }
        }
		
		// --- Also check the neighbouring cells of each owner cell ---
		for (label i = 0; i < nFaces; ++i)
		{
			const label globalBodyFace = startFace + i;
			const label ownerCell = faceOwners[globalBodyFace];

			// Get the neighboring cells (cells that share faces with the owner)
			const labelList& neighCells = mesh.cellCells()[ownerCell];

			forAll(neighCells, ncI)
			{
				const label neighCell = neighCells[ncI];
				const cell& neighCellFaces = mesh.cells()[neighCell];

				forAll(neighCellFaces, fJ)
				{
					const label neighFaceID = neighCellFaces[fJ];
					const label neighPatchI = mesh.boundaryMesh().whichPatch(neighFaceID);

					if (neighPatchI != -1)
					{
						const polyPatch& neighPatch = mesh.boundaryMesh()[neighPatchI];
						const word& neighPatchName = neighPatch.name();

						if (neighPatchName == "freeSurface" || neighPatchName == "upper")
						{
							const face& bodyFace = mesh.faces()[globalBodyFace];
							const face& fsFace   = mesh.faces()[neighFaceID];
							const pointField& points = mesh.points();

							// --- Find common vertices (shared edge) ---
							DynamicList<vector> commonVerts;
							forAll(bodyFace, vi)
							{
								const label vBody = bodyFace[vi];
								forAll(fsFace, vj)
								{
									if (vBody == fsFace[vj])
									{
										commonVerts.append(points[vBody]);
									}
								}
							}

							scalar edgeLength = 0.0;
							if (commonVerts.size() == 2)
							{
								// --- prevent duplicates ---
								bool alreadyAdded = false;
								forAll(foundPairs, fpI)
								{
									if (foundPairs[fpI].second() == neighFaceID)
									{
										alreadyAdded = true;
										break;
									}
								}

								if (!alreadyAdded)
								{
									edgeLength = mag(commonVerts[0] - commonVerts[1]);
									foundPairs.append(labelPair(i, neighFaceID));
									sharedEdgeLengths.append(edgeLength);
									break;
							
								}

								
							}
							
							
						}
					}
				}
			}
		}

		
		
         bodyToFreeSurfaceFaces_[patchi].transfer(foundPairs);
		 bodyToEdgeLengths_[patchi].transfer(sharedEdgeLengths);

		//Pout << "Proc " << Pstream::myProcNo() << " bodyToFreeSurfaceFaces_ " << bodyToFreeSurfaceFaces_[patchi] << endl; 
		//Pout << "Proc " << Pstream::myProcNo() << " patch " << patch.name()  << " found " << foundPairs.size() << " pairs" << " | edge lengths: " << bodyToEdgeLengths_ << nl;
      
    }
	
	FSpatchesfound_ = true;
	//Pout << "Proc " << Pstream::myProcNo() << " bodyToFreeSurfaceFaces_ " << bodyToFreeSurfaceFaces_  <<  nl;
	
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
    zetaAtBodyFaces_.setSize(mesh.boundary().size());
	
    for (const label patchi : patchIDs_)
    {
       // const label patchi = patchIDs_[patchi];
        //const fvPatch& patch = mesh.boundary()[patchi];

        const List<labelPair>& facePairs = bodyToFreeSurfaceFaces_[patchi];
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
			//const label patchi = patchIDs_[patchi];
			//const label globalBodyFaceID = mesh.boundary()[patchi].start() + bodyPatchFaceI;
			//const label ownerCell = mesh.faceOwner()[globalBodyFaceID];
			//const vector& cellCenter = mesh.C()[ownerCell];

			// Debug print
			//if ( Pstream::myProcNo() == 0) Pout << "k = " << k  << " | freeSurfaceFaceID = " << freeSurfFaceID  << " | zeta.z = " << zetaVec.z() << " | bodyFace center = " << mesh.Cf()[globalBodyFaceID] << " | ownerCell center = " << cellCenter << endl;
				
				// --- NEW: Get owner cell center of the free-surface face
			//const label fsOwnerCell = mesh.faceOwner()[freeSurfFaceID];
			//const vector& fsOwnerCenter = mesh.C()[fsOwnerCell];

			// Debug print
			//if (Pstream::myProcNo() == 0)
			//{
				//Pout << "k = " << k  << " | freeSurfaceFaceID = " << freeSurfFaceID  << " | zeta.z = " << zetaVec.z()  << " | freeSurfaceFace center = " << mesh.Cf()[freeSurfFaceID] << " | ownerCell center (of freeSurface face) = " << fsOwnerCenter << endl;
					
	 
			//}
		
				 
				 
				 

        }

        // Store per patch
        zetaAtBodyFaces_[patchi].transfer(zetaVals);
		//if ( Pstream::myProcNo() == 0) Pout << "Proc " << Pstream::myProcNo() << " zetaAtBodyFaces_ " << zetaAtBodyFaces_[patchi] << endl; 
     
    }
}


