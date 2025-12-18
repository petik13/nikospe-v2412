/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
neigbouring and scheme FD bu TURGUT

\*---------------------------------------------------------------------------*/


void Foam::waveCurrentPotential3DFvPatchScalarField::calcNeigboursV3()
{
	if (neigboursCalculated_) return;
	
	struct RemoteInfo
	{
		label  globalId;
		vector remoteCc;
		vector localCf;
	};

	HashTable<List<RemoteInfo>, label> sendBuffer_;
	HashTable<List<RemoteInfo>, label> recvBuffer_;	
	
	
	const fvMesh& mesh = this->patch().boundaryMesh().mesh();

	const label startFace = patch().start();
	const label nFaces    = patch().size();
	const label endFace   = startFace + nFaces;
	const labelList& faceOwners = mesh.faceOwner();

	const label nLocal = patch().size();
	globalIdToCf_.clear();
	
	// Automatic global indexing without manual listCombineGather
	gIndexPatch_.reset(new Foam::globalIndex(nLocal));

	globalIdOfLocalFace_.setSize(nLocal);
	for (label i = 0; i < nLocal; ++i)
	{
		globalIdOfLocalFace_[i] = gIndexPatch_->toGlobal(Pstream::myProcNo(), i);
	}
	
	const vectorField Cf = patch().Cf();  
	
	
	// --- Build local arrays for fast access
	List<vector> localCellCentres(nFaces);
	for (label i = 0; i < nFaces; ++i)
	{
		label ownerCell = faceOwners[startFace + i];
		localCellCentres[i] = mesh.C()[ownerCell];
	}


		
	
	globalIdToPackedIdx_.clear();
	for (label i = 0; i < nLocal; ++i)
	{
		globalIdToPackedIdx_.insert(globalIdOfLocalFace_[i], i);
	}	
	
	neighboursGlobal_.setSize(nLocal);
	neighbours_.setSize(nFaces);
	procBoundaryFace_.setSize(nFaces, false);
	neighbProcNos_.setSize(nFaces);
	
	
	// One hash per local face
	List<HashTable<label, label>> neighboursPerProc(nFaces);

	// no need to setSize() again — constructor does it
	forAll(neighboursPerProc, i)
	{
		neighboursPerProc[i].clear();
	}



	
	forAll(neighbProcNos_, i) neighbProcNos_[i] = labelList(0);  // empty list
	

	// --- Loop over patch faces
	for (label i = 0; i < nFaces; ++i)
	{
		const label currentFace = startFace + i;
		const label ownerCell   = faceOwners[currentFace];
		
		// ---- (A) Flag processor-boundary using ONLY the owner cell faces
		const face& fCurrent = mesh.faces()[currentFace];
		// build a lookup set of current face’s vertices
		labelHashSet currentVerts;
		forAll(fCurrent, vI)
		{
			currentVerts.insert(fCurrent[vI]);  // store point label, not coordinate
		}


		const cell& ownerCellFaces = mesh.cells()[ownerCell];
		// Go through all other faces of the owner cell
		forAll(ownerCellFaces, fJ)
		{
			const label faceID = ownerCellFaces[fJ];
			if (faceID == currentFace) continue; // Skip self
			
			

			const label patchi = mesh.boundaryMesh().whichPatch(faceID);

			if (patchi != -1 && isA<processorPolyPatch>(mesh.boundaryMesh()[patchi]))
			{
				const face& fOther = mesh.faces()[faceID];

				// count how many vertices are shared by label
				label shared = 0;
				forAll(fOther, vJ)
				{
					if (currentVerts.found(fOther[vJ]))
						++shared;
				}

				// decide adjacency criterion
				if (shared < 1) continue;  // edge-based adjacency (safer than just 1)

				// ---- at this point, faceID is a relevant proc neighbour
				const processorPolyPatch& pp =
					refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);

				label np = pp.neighbProcNo();
				procBoundaryFace_[i] = true;

				if (!neighbProcNos_[i].found(np))
					neighbProcNos_[i].append(np);


				// local face index on this processor patch
				label localFaceI = faceID - pp.start();
				
				
				// 🔑 Store the mapping
				if (!neighboursPerProc[i].found(np))
				{
					neighboursPerProc[i].insert(np, 1);  // first time we see this proc
				}
				else
				{
					 neighboursPerProc[i][np] += 1;  // increment count
				}
				//label localCell  = pp.faceCells()[localFaceI];
				//const vector& nbrCf = pp.neighbFaceCentres()[localFaceI];
				const vector& nbrCc = pp.neighbFaceCellCentres()[localFaceI];
				
				label gid = globalIdOfLocalFace_[i];
				RemoteInfo ri;
				ri.globalId = gid;
				ri.remoteCc = nbrCc;
				ri.localCf  = Cf[i];

				if (!sendBuffer_.found(np))
				{
					sendBuffer_.insert(np, List<RemoteInfo>());
				}
				sendBuffer_[np].append(ri);
			}
		}

		

		// ---- (B) Build local neighbour list using neighbour cells
		labelHashSet uniqNei; // avoid duplicates

		const labelList& nCells = mesh.cellCells()[ownerCell];
		forAll(nCells, ncI)
		{
			const label nbrCell = nCells[ncI];
			const cell& nbrCellFaces = mesh.cells()[nbrCell];

			forAll(nbrCellFaces, fJ)
			{
				const label faceID = nbrCellFaces[fJ];

				// Only faces on THIS patch (excluding self)
				if (faceID >= startFace && faceID < endFace && faceID != currentFace) //if (mesh.boundaryMesh().whichPatch(faceID) == patchID) why nnot
				{
					const label localNeighbour = faceID - startFace;
					uniqNei.insert(localNeighbour);
				}
			}
		}

		// Write back unique neighbours (local and global)
		const labelList uniq = uniqNei.toc();
		neighbours_[i] = uniq;

		neighboursGlobal_[i].setSize(uniq.size());
		forAll(uniq, jj)
		{
			neighboursGlobal_[i][jj] = globalIdOfLocalFace_[uniq[jj]];
		}
	}
	
	
/*---------------------------------------------------------------------------*\
  using RemoteInfoTable = HashTable<List<RemoteInfo>, label>;

	Pout << "proc " << Pstream::myProcNo() << " sendBuffer_:" << endl;

	forAllConstIter(RemoteInfoTable, sendBuffer_, it)
	{
		label neighbProc = it.key();
		const List<RemoteInfo>& infos = it();

		Pout << "  -> neighbourProc " << neighbProc
			 << " has " << infos.size() << " entries:" << endl;

		forAll(infos, j)
		{
			const RemoteInfo& ri = infos[j];
			Pout << "     gid=" << ri.globalId
				 << " remoteCc=" << ri.remoteCc << endl;
		}
	}
	if (Pstream::myProcNo() ==3) 
		// Loop over my local faces
			for (label i = 0; i < nFaces; ++i)
			{
				Pout << "proc " << Pstream::myProcNo() << " neighboursPerProc_ at:" << globalIdOfLocalFace_[i] << " " << neighboursPerProc[i] << endl;	
			}
	
\*---------------------------------------------------------------------------*/
	
	
	//paralell part
		if (Pstream::parRun())
	{
		// --- Prepare communication
		PstreamBuffers pBuf(Pstream::commsTypes::nonBlocking);
		using RemoteInfoTable = HashTable<List<RemoteInfo>, label>;
		// --- Send phase
		forAllConstIter(RemoteInfoTable, sendBuffer_, it)
		{
			const label neighbProc = it.key();
			const List<RemoteInfo>& infos = it();

			UOPstream os(neighbProc, pBuf);

			// Send number of entries
			os << infos.size();

			// Send each RemoteInfo
			forAll(infos, j)
			{
				 os << infos[j].globalId
					<< infos[j].localCf   // send my Cf
					<< infos[j].remoteCc; // send neighbour cell centre
			}
		}

		// --- Finish sending
		pBuf.finishedSends();


		// --- Receive phase
		forAllConstIter(RemoteInfoTable, sendBuffer_, it)
		{
			const label neighbProc = it.key();
			UIPstream is(neighbProc, pBuf);

			label nRecv; 
			is >> nRecv;

			List<RemoteInfo> rInfos(nRecv);
			for (label j = 0; j < nRecv; ++j)
			{
				is >> rInfos[j].globalId
					>> rInfos[j].localCf
					>> rInfos[j].remoteCc;
			}

			// Store in recvBuffer_
			recvBuffer_.insert(neighbProc, rInfos);
		}

		/*---------------------------------------------------------------------------*\
		// --- Debug output
		Pout << "proc " << Pstream::myProcNo() << " received:" << endl;
		forAllConstIter(RemoteInfoTable, recvBuffer_, it)
		{
			label srcProc = it.key();
			const List<RemoteInfo>& rInfos = it();

			forAll(rInfos, j)
			{
				const RemoteInfo& ri = rInfos[j];
				Pout << "   from proc " << srcProc
					 << " gid=" << ri.globalId
					 << " remoteCc=" << ri.remoteCc << endl;
			}
		}
		
		
		
		\*---------------------------------------------------------------------------*/	
		
		
		forAllConstIter(RemoteInfoTable, recvBuffer_, it)
		{
			label srcProc = it.key();
			const List<RemoteInfo>& rInfos = it();

			// Extract arrays for easier nearest-neighbour scan
			List<vector> remoteCcs(rInfos.size());
			List<label>  remoteGids(rInfos.size());
			forAll(rInfos, j)
			{
				remoteCcs[j]  = rInfos[j].remoteCc;
				remoteGids[j] = rInfos[j].globalId;
				globalIdToCf_.set(remoteGids[j], rInfos[j].localCf);

			}

			// Loop over my local faces
			for (label i = 0; i < nFaces; ++i)
			{
				if (!procBoundaryFace_[i] || !neighbProcNos_[i].found(srcProc)) continue;
				const vector& myCc = localCellCentres[i];

				// track already matched remote indices
				labelHashSet matched;

				// Get number of expected matches with this srcProc
				label N_neighbor = 0;
				if (neighboursPerProc[i].found(srcProc))
					N_neighbor = neighboursPerProc[i][srcProc];
				
				for (label kk = 0; kk < N_neighbor; ++kk)
				{
					scalar bestDist = GREAT;
					label  best = -1;

					// Brute-force nearest neighbour search
					for (label j = 0; j < remoteCcs.size(); ++j)
					{
						if (matched.found(j)) continue; // skip already matched

						const vector& rc = remoteCcs[j];
						scalar dx = rc.x() - myCc.x();
						scalar dy = rc.y() - myCc.y();
						scalar dz = rc.z() - myCc.z();
						scalar dist = sqrt(dx*dx + dy*dy + dz*dz);

						if (dist < bestDist)
						{
							bestDist = dist;
							best = j;
						}
					}

					if (best != -1)
					{
						matched.insert(best); // mark remote index as used

						label gRemote = remoteGids[best];
						if (!neighboursGlobal_[i].found(gRemote))
						{
							neighboursGlobal_[i].append(gRemote);

							//Pout << "proc " << Pstream::myProcNo() << " matched localGID=" << globalIdOfLocalFace_[i] << " with remoteGID=" << gRemote << " from proc " << srcProc << " at distance=" << bestDist << endl;
				
						}
					}
				}
			}

		
		
		}	

		
		
	}

	
	
/*---------------------------------------------------------------------------*\
	Pout<< "proc " << Pstream::myProcNo()
    << " neighboursGlobal_: " << neighboursGlobal_ << nl << endl;
	
	
	Pout << endl;
\*---------------------------------------------------------------------------*/		
	//if (Pstream::myProcNo()==2) Pout<< "proc " << Pstream::myProcNo() << " globalIdToCf_: " << globalIdToCf_ << nl << endl;
	      
	findInletFaces();
	findSphereFaces();
    neigboursCalculated_ = true;
}

void Foam::waveCurrentPotential3DFvPatchScalarField::findInletFaces()
{
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const label startFace = patch().start();
    const label nFaces    = patch().size();
    const labelList& faceOwners = mesh.faceOwner();

	bool hasInletPatchGlobal =
        returnReduce(mesh.boundary().findPatchID("inlet") != -1, orOp<bool>());

    if (!hasInletPatchGlobal)
    {
        FatalErrorInFunction
            << "Global 'inlet' patch not found on any processor." << abort(FatalError);
    }


    ownerHasInletFace_.setSize(nFaces, false);

    for (label i = 0; i < nFaces; ++i)
    {
        const label globalFace = startFace + i;
        const label ownerCell = faceOwners[globalFace];
        const cell& ownerCellFaces = mesh.cells()[ownerCell];

        forAll(ownerCellFaces, fJ)
        {
            const label faceID = ownerCellFaces[fJ];
            const label patchi = mesh.boundaryMesh().whichPatch(faceID);

            if (patchi != -1)
            {
                const polyPatch& bPatch = mesh.boundaryMesh()[patchi];
                const word& patchName = bPatch.name();

                if (patchName == "inlet")
                {
                    ownerHasInletFace_[i] = true;
                    break; // No need to check more faces of this cell
                }
            }
        }
    }
	// for (label i = 0; i < nFaces; ++i)
   // {
	// if (ownerHasInletFace_[i]) Pout<< "proc " << Pstream::myProcNo() << " ownerHasInletFace_ at: " << globalIdOfLocalFace_[i] << nl << endl;
	//}
}



void Foam::waveCurrentPotential3DFvPatchScalarField::findSphereFaces()
{
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
	const word& patchName = patch().name();
	
	
    const label startFace = patch().start();
    const label nFaces    = patch().size();
    const labelList& faceOwners = mesh.faceOwner();
    const labelList& faceNeighbours = mesh.faceNeighbour(); // may be smaller than all faces

    ownerHasBodyFace_.setSize(nFaces, false);
    
	ownerHasFirstBodyFace_.setSize(nFaces, false);
    
	bool bodyfoundLocal = false;

    // --- check globally that sphere exists ---
    bool hasSpherePatchGlobal =
        returnReduce(mesh.boundary().findPatchID("sphere") != -1, orOp<bool>());

    if (!hasSpherePatchGlobal)
    {
        FatalErrorInFunction
            << "Global 'sphere' patch not found on any processor." << abort(FatalError);
    }

    // --- function to test if a cell has a face on 'sphere' ---
    auto touchesSphere = [&](label cellI) -> bool
    {
        if (cellI < 0 || cellI >= mesh.cells().size()) return false;
        const cell& cFaces = mesh.cells()[cellI];
        forAll(cFaces, fI)
        {
            const label fid = cFaces[fI];
            const label patchi = mesh.boundaryMesh().whichPatch(fid);
            if (patchi == -1) continue;
            if (mesh.boundaryMesh()[patchi].name() == "sphere")
                return true;
        }
        return false;
    };


	auto touchesPatch = [&](label cellI, const word& patchNameToCheck) -> bool
	{
		if (cellI < 0 || cellI >= mesh.cells().size()) return false;
		const cell& cFaces = mesh.cells()[cellI];
		forAll(cFaces, fI)
		{
			const label fid = cFaces[fI];
			const label patchi = mesh.boundaryMesh().whichPatch(fid);
			if (patchi == -1) continue;
			if (mesh.boundaryMesh()[patchi].name() == patchNameToCheck)
				return true;
		}
		return false;
	};




    // --- loop faces of this patch ---
    for (label i = 0; i < nFaces; ++i)
    {
        const label globalFace = startFace + i;
        const label owner = faceOwners[globalFace];
        const label neighbour =
            (globalFace < faceNeighbours.size() ? faceNeighbours[globalFace] : -1);

        // direct check on owner/neighbour
        bool found = touchesSphere(owner) || touchesSphere(neighbour);
		bool found2 = touchesSphere(owner);
        
		
		// NEW: also mark as first ring if neighbour touches sphere but not freeSurface for trapezoid
		if (!found2)
		{
			const face& freeSurfFace = mesh.faces()[startFace + i];  // <<< FIX
			const cell& cOwner = mesh.cells()[owner];

			// loop through owner's internal neighbor cells
			forAll(cOwner, fJ)
			{
				const label faceJ = cOwner[fJ];
				if (faceJ >= faceNeighbours.size()) continue;

				label nCell = (faceOwners[faceJ] == owner)
							  ? faceNeighbours[faceJ]
							  : faceOwners[faceJ];

				// Only consider neighbour cells that touch the sphere
				if (!touchesSphere(nCell) || touchesPatch(nCell, patchName))
					continue;

				// Loop over all sphere faces of that neighbour cell
				const cell& nC = mesh.cells()[nCell];
				forAll(nC, fN)
				{
					const label fidN = nC[fN];
					const label patchi = mesh.boundaryMesh().whichPatch(fidN);
					if (patchi == -1) continue;
					if (mesh.boundaryMesh()[patchi].name() != "sphere") continue;

					const face& sphereFace = mesh.faces()[fidN];

					
					/*---------------------------------------------------------------------------*\
					// Debug info
					Pout << "Checking FS face " << globalIdOfLocalFace_[i]
						 << " (owner=" << owner << ") vs sphere face " << fidN
						 << " (cell=" << nCell << ")" << endl;

					Pout << "  FS vertex coords: ";
					forAll(freeSurfFace, vA)
					{
						const point& pA = mesh.points()[freeSurfFace[vA]];
						Pout << "(" << freeSurfFace[vA] << ": "
							 << pA.x() << " " << pA.y() << " " << pA.z() << ") ";
					}
					Pout << endl;

					Pout << "  Sphere vertex coords: ";
					forAll(sphereFace, vB)
					{
						const point& pB = mesh.points()[sphereFace[vB]];
						Pout << "(" << sphereFace[vB] << ": "
							 << pB.x() << " " << pB.y() << " " << pB.z() << ") ";
					}
					Pout << endl;

					\*---------------------------------------------------------------------------*/	
					// Compare vertex coordinates with tolerance
					int sharedCount = 0;
					const scalar tol = 1e-6;  // geometric tolerance

					forAll(freeSurfFace, vA)
					{
						const point& pA = mesh.points()[freeSurfFace[vA]];
						forAll(sphereFace, vB)
						{
							const point& pB = mesh.points()[sphereFace[vB]];
							if (mag(pA - pB) < tol)
							{
								sharedCount++;
								if (sharedCount >= 2) break;
							}
						}
						if (sharedCount >= 2) break;
					}

					if (sharedCount > 0)
						Pout << "  Shared vertices (by coord): " << sharedCount << endl;

					if (sharedCount >= 2)
					{
						//Pout << "proc " << Pstream::myProcNo()
							// << " promoted face " << globalIdOfLocalFace_[i]
							 //<< " to first ring (shares edge geometrically with sphere face)"
							 //<< nl;
						found2 = true;
						break;
					}
				}

				if (found2) break; // stop if found
			}
		}


			
		
		
		// if still not found, check neighbours of those cells
        if (!found)
        {
            const cell& cOwner = mesh.cells()[owner];
            forAll(cOwner, fJ)
            {
                const label faceJ = cOwner[fJ];
                if (faceJ < faceNeighbours.size())
                {
                    label nCell = (faceOwners[faceJ] == owner)
                                  ? faceNeighbours[faceJ]
                                  : faceOwners[faceJ];
                    if (touchesSphere(nCell))
                    {
                        found = true;
                        break;
                    }
                }
            }
        }

        if (found)
        {
            ownerHasBodyFace_[i] = true;
            bodyfoundLocal = true;
        }
		
		if (found2)
        {
			ownerHasFirstBodyFace_[i] = true;
		}
		
	
	}

    bool bodyfoundGlobal = returnReduce(bodyfoundLocal, orOp<bool>());
    if (!bodyfoundGlobal)
    {
        FatalErrorInFunction << "No faces near 'sphere' found globally." << abort(FatalError);
    }
	buildVertexBasedConnections();
    // --- debug output ---
	//for (label i = 0; i < nFaces; ++i)
    //{
		//if (ownerHasFirstBodyFace_[i]) Pout << "proc " << Pstream::myProcNo()  << " ownerHasBodyFace_ at global ID " << globalIdOfLocalFace_[i] <<  endl;
	//}

   
}

void Foam::waveCurrentPotential3DFvPatchScalarField::buildVertexBasedConnections()
{
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    const label startFace = patch().start();
    const label nFaces = patch().size();
	
	
	List<List<label>> firstRingConnections_;
    firstRingConnections_.setSize(nFaces);
    forAll(firstRingConnections_, i)
        firstRingConnections_[i].clear();

    // --- Step 1: Collect all vertices used by first-ring faces ---
    // Build a map from vertexID -> list of faces using that vertex
    Map<DynamicList<label>> vertexToFaces;

    for (label i = 0; i < nFaces; ++i)
    {
        if (!ownerHasFirstBodyFace_[i]) continue;

        const label faceI = startFace + i;
        const face& f = mesh.faces()[faceI];

        forAll(f, vI)
        {
            label vertexID = f[vI];
            if (!vertexToFaces.found(vertexID))
                vertexToFaces.insert(vertexID, DynamicList<label>());
            vertexToFaces[vertexID].append(i); // store patch-local index
        }
    }

    // --- Step 2: For each vertex, connect all faces sharing it ---
    forAllIter(Map<DynamicList<label>>, vertexToFaces, it)
    {
        const DynamicList<label>& facesAtVertex = it.val();

        forAll(facesAtVertex, a)
        {
            label fa = facesAtVertex[a];
            for (label b = a + 1; b < facesAtVertex.size(); ++b)
            {
                label fb = facesAtVertex[b];

                if (fa == fb) continue;

                // connect fa <-> fb
                if (!firstRingConnections_[fa].found(fb)) firstRingConnections_[fa].append(fb);
				if (!firstRingConnections_[fb].found(fa)) firstRingConnections_[fb].append(fa);
                
            
			}
        }
    }
	const vectorField Cf = patch().Cf();  
     // --- Step 3: Debug print ---
    for (label i = 0; i < nFaces; ++i)
    {
        if (!ownerHasFirstBodyFace_[i]) continue;
        const List<label>& conn = firstRingConnections_[i];
        if (conn.empty()) continue;

             
        forAll(conn, j)
		{
			if (Cf[conn[j]].x() <  Cf[i].x())
			{
				if (!ringNodeX_.found(i))
					ringNodeX_.insert(i, globalIdOfLocalFace_[conn[j]]);
				else
					ringNodeX_[i] = globalIdOfLocalFace_[conn[j]];  // overwrite if already assigned
			
			//label node = ringNodeX_[i];
			//Pout << "proc " << Pstream::myProcNo() << " vertex-connected face " << globalIdOfLocalFace_[i] << " linked to: " << node << endl;;
        
			}
			
			
		}     
		
    }
	
	/*---------------------------------------------------------------------------*\
	for (label i = 0; i < nFaces; ++i)
    {
        if (!ownerHasFirstBodyFace_[i]) continue;
        const List<label>& conn = firstRingConnections_[i];
        if (conn.empty()) continue;

        Pout << "proc " << Pstream::myProcNo() << " vertex-connected face " << globalIdOfLocalFace_[i] << " linked to: ";
             
             

        forAll(conn, j)
           Pout << globalIdOfLocalFace_[conn[j]] << " ";
        Pout << endl;
    }


	\*---------------------------------------------------------------------------*/	
	Map<List<label>> firstRingToSecondRing_;
	Map<DynamicList<label>> vertexToFirstRingFaces;

	forAll(ownerHasFirstBodyFace_, i)
	{
		if (!ownerHasFirstBodyFace_[i]) continue;
		const face& f = mesh.faces()[startFace + i];
		forAll(f, v)
		{
			label vid = f[v];
			if (!vertexToFirstRingFaces.found(vid))
				vertexToFirstRingFaces.insert(vid, DynamicList<label>());
			vertexToFirstRingFaces[vid].append(i);
		}
	}


	

	for (label i = 0; i < nFaces; ++i)
	{
		if (ownerHasFirstBodyFace_[i]) continue;  // skip ring0 and ring1
		
		
		if (!ownerHasBodyFace_[i]) continue;  // bunu yapmassak vertex paylasanlari da buluyor
		
		
		const face& f = mesh.faces()[startFace + i];
		DynamicList<label> linkedFirstRingFaces;

		forAll(f, v)
		{
			label vid = f[v];
			if (vertexToFirstRingFaces.found(vid))
			{
				const DynamicList<label>& fr = vertexToFirstRingFaces[vid];
				forAll(fr, k)
				{
					if (!linkedFirstRingFaces.found(fr[k]))
						linkedFirstRingFaces.append(fr[k]);
				}
			}
		}

		if (linkedFirstRingFaces.size() > 0)
		{
			secondRingToFirstRing_.insert(i, linkedFirstRingFaces);
			// optional reverse link:
			forAll(linkedFirstRingFaces, k)
			{
				label frI = linkedFirstRingFaces[k];
				if (!firstRingToSecondRing_.found(frI))
					firstRingToSecondRing_.insert(frI, DynamicList<label>());
				firstRingToSecondRing_[frI].append(i);
			}
		}
	}
	/*---------------------------------------------------------------------------*\
	forAllConstIter(Map<List<label>>, secondRingToFirstRing_, it)
	{
		label sFace = it.key();
		const List<label>& frList = it.val();
		Pout << "proc " << Pstream::myProcNo()
			 << " 2nd-ring face " << globalIdOfLocalFace_[sFace]
			 << " connected to first-ring faces: ";
		forAll(frList, j)
			Pout << globalIdOfLocalFace_[frList[j]] << " ";
		Pout << endl;
	}
	\*---------------------------------------------------------------------------*/	
	
	 for (label i = 0; i < nFaces; ++i)
    {
        if (ownerHasFirstBodyFace_[i]) continue;
		if (!ownerHasBodyFace_[i]) continue;  
		
        const List<label>& conn = secondRingToFirstRing_[i];
		if (conn.empty()) continue;

		// Pick the first-ring face with smaller x (upwind)
		const vectorField& Cf = patch().Cf();
		scalar minDx = GREAT;
		label bestConn = -1;

		forAll(conn, j)
		{
			label frIdx = conn[j]; // first-ring local index
			scalar dx = Cf[i].x() - Cf[frIdx].x();

			// look for upstream connection (smaller x)
			if (dx > VSMALL && dx < minDx)
			{
				minDx = dx;
				bestConn = frIdx;
			}
		}

		if (bestConn != -1)
		{
			ring2NodeX_.insert(i, globalIdOfLocalFace_[bestConn]);
			//Pout << "proc " << Pstream::myProcNo()  << " second-ring face " << globalIdOfLocalFace_[i] << " linked (−x) to first-ring face " << globalIdOfLocalFace_[bestConn] << nl;
	 
				 
		} 

		if (bestConn == -1)
		{
			forAll(conn, j)
			{
				label frIdx = conn[j]; // first-ring local index
				//scalar dx = Cf[i].x() - Cf[frIdx].x();
				scalar dx = Cf[frIdx].x()-Cf[i].x();

				// look for upstream connection (smaller x)
				if (dx < minDx)
				{
					minDx = dx;
					bestConn = frIdx;
				}
			}

			
			
			ring2NodeX_.insert(i, globalIdOfLocalFace_[bestConn]);
			//Pout << "proc " << Pstream::myProcNo()  << " second-ring face " << globalIdOfLocalFace_[i] << " linked (−x) to first-ring face " << globalIdOfLocalFace_[bestConn] << nl;
	 
				 
		} 
			
			
    }
	
	
	
	    // ===============================================================
    // Find global maximum Cf.x among ownerHasFirstBodyFace_ == true
    // ===============================================================

   // ===============================================================
// Global maximum of face center X (robust for missing faces)
// ===============================================================

	scalar localMaxX = -GREAT;

	// --- compute local maximum only for valid faces ---
	forAll(Cf, i)
	{
		if (ownerHasFirstBodyFace_[i])
		{
			localMaxX = Foam::max(localMaxX, Cf[i].x());
		}
	}

	// --- ensure all processors have a defined value ---
	scalar globalMaxX = localMaxX;

	// --- if a proc has no faces, its localMaxX = -GREAT, so max reduction ignores it ---
	Foam::reduce(globalMaxX, Foam::maxOp<scalar>());

	// --- optional info ---
	if (Pstream::master())
	{
		Info << "Global maximum Cf.x = " << globalMaxX << nl;
	}


    
    // ===============================================================
    // 4. Debug output 
    // ===============================================================
    //Pout << "Processor " << Pstream::myProcNo() << " knows that global maximum Cf.x = " << globalMaxX << nl;
       
         
    const scalar limitingX = 0.05;     
	
	 forAll(Cf, i)
    {
        
		//if (ownerHasFirstBodyFace_[i])
		if (ownerHasBodyFace_[i]) //V13B
		{	scalar x00 = Cf[i].x();  // <-- or .y() or .z() depending on what you want
			if (x00 > (globalMaxX-limitingX))  //burasi 0.01 idi
			{
				Pout << globalIdOfLocalFace_[i] << " stay outside" << endl;
				ownerHasBodyFace_[i]= false;
			}
		}	
    }
	
	
}


#include "Outputter.H" 

void Foam::waveCurrentPotential3DFvPatchScalarField::findUpwindDownwindNodesV2()
{
    if (windsCalculated_)
    {
        return;
    }

    Info << "Finding upwind and downwind nodes in X and Y..." << endl;

    const vectorField& patchFaceCenters = patch().Cf();
    scalarField xComponents = patchFaceCenters.component(vector::X);
    scalarField yComponents = patchFaceCenters.component(vector::Y);

    const label nFaces = patch().size();
	const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    upwindNodesX_.setSize(nFaces, List<label>());
    downwindNodesX_.setSize(nFaces, List<label>());
    upwindNodesY_.setSize(nFaces, List<label>());
    downwindNodesY_.setSize(nFaces, List<label>());
	const label startFace = patch().start();

    for (label i = 0; i < nFaces; i++)
    {
        scalar x_current = xComponents[i];
        scalar y_current = yComponents[i];

        List<label> upwindX, downwindX;
        List<label> upwindY, downwindY;
		 
		//const face& f = mesh.faces()[globalIdOfLocalFace_[i]];
		const face& f = mesh.faces()[startFace + i];
		scalar minX = GREAT, maxX = -GREAT;
		scalar minY = GREAT, maxY = -GREAT;

		forAll(f, vI)
		{
		const point& pt = mesh.points()[f[vI]];
		minX = min(minX, pt.x());
		maxX = max(maxX, pt.x());
		minY = min(minY, pt.y());
		maxY = max(maxY, pt.y());
						}

		
        for (label j = 0; j < neighbours_[i].size(); j++)
        {
            
			bool foundneighbor = false; 
			
			label neighborLocal = neighbours_[i][j];
            scalar x_neighbor = xComponents[neighborLocal];
            scalar y_neighbor = yComponents[neighborLocal];

            scalar dx = x_current - x_neighbor;
            scalar dy = y_current - y_neighbor;

		
            bool centerInX = (x_neighbor >= minX && x_neighbor <= maxX);
			bool centerInY = (y_neighbor >= minY && y_neighbor <= maxY);
			
			//if (globalIdOfLocalFace_[i]==983) Pout << "Lokal neighbour :" << globalIdOfLocalFace_[neighborLocal] << " centerInX :" << centerInX  << " centerInY :" << centerInY << endl;
				
			
			
			
			if (dx > VSMALL && centerInY) {upwindX.append(globalIdOfLocalFace_[neighborLocal]);foundneighbor = true;}
			if (dx < -VSMALL && centerInY) {downwindX.append(globalIdOfLocalFace_[neighborLocal]);foundneighbor = true;}
			if (dy > VSMALL && centerInX) {upwindY.append(globalIdOfLocalFace_[neighborLocal]);foundneighbor = true;}
			if (dy < -VSMALL && centerInX) {downwindY.append(globalIdOfLocalFace_[neighborLocal]);foundneighbor = true;}
			
			if (!foundneighbor){
				const face& fj = mesh.faces()[startFace + neighborLocal];
				scalar minXj = GREAT, maxXj = -GREAT;
				scalar minYj = GREAT, maxYj = -GREAT;
				forAll(fj, vi)
				{
					const point& pt = mesh.points()[fj[vi]];
					minXj = min(minXj, pt.x());
					maxXj = max(maxXj, pt.x());
					minYj = min(minYj, pt.y());
					maxYj = max(maxYj, pt.y());
				}
				bool centerOfIinXj = (x_current >= minXj && x_current <= maxXj);
				bool centerOfIinYj = (y_current >= minYj && y_current <= maxYj);
				
				if (dx > VSMALL && centerOfIinYj)
					upwindX.append(globalIdOfLocalFace_[neighborLocal]);
				if (dx < -VSMALL && centerOfIinYj)
					downwindX.append(globalIdOfLocalFace_[neighborLocal]);

				if (dy > VSMALL && centerOfIinXj)
					upwindY.append(globalIdOfLocalFace_[neighborLocal]);
				if (dy < -VSMALL && centerOfIinXj)
					downwindY.append(globalIdOfLocalFace_[neighborLocal]);
				
				
			}
	 
	   }
		

        for (label j = 0; j < neighboursGlobal_[i].size(); j++)
        {
            const label globalNei = neighboursGlobal_[i][j];

            if (!globalIdToCf_.found(globalNei)) continue;
			
			bool foundneighbor = false; 
			
            const vector& CfNei = globalIdToCf_[globalNei];
            scalar x_neighbor = CfNei.x();
            scalar y_neighbor = CfNei.y();

            scalar dx = x_current - x_neighbor;
            scalar dy = y_current - y_neighbor;
			bool centerInX = (x_neighbor >= minX && x_neighbor <= maxX);
			bool centerInY = (y_neighbor >= minY && y_neighbor <= maxY);

            if (dx > VSMALL && centerInY)
                if (!upwindX.found(globalNei)) {upwindX.append(globalNei);foundneighbor = true;}

            if (dx < -VSMALL && centerInY)
				if (!downwindX.found(globalNei)) {downwindX.append(globalNei);foundneighbor = true;}

			if (dy > VSMALL && centerInX)
                if (!upwindY.found(globalNei)) {upwindY.append(globalNei);foundneighbor = true;}

			if (dy < -VSMALL && centerInX)
				if (!downwindY.found(globalNei)) {downwindY.append(globalNei);foundneighbor = true;}
        
			
			//minX maxX eksik o yuzden eski usul (can fail when surface face is not close to square)
			if (!foundneighbor)
			{	if (dx > VSMALL && mag(dx) > mag(dy))
					if (!upwindX.found(globalNei)) upwindX.append(globalNei);

				if (dx < -VSMALL && mag(dx) > mag(dy))
					if (!downwindX.found(globalNei)) downwindX.append(globalNei);

				if (dy > VSMALL && mag(dy) > mag(dx))
					if (!upwindY.found(globalNei)) upwindY.append(globalNei);

				if (dy < -VSMALL && mag(dy) > mag(dx))
					if (!downwindY.found(globalNei)) downwindY.append(globalNei);
			}
			//minX maxX eksik o yuzden eski usul
		
		}

        upwindNodesX_[i]   = upwindX;
        downwindNodesX_[i] = downwindX;
        upwindNodesY_[i]   = upwindY;
        downwindNodesY_[i] = downwindY;
    }

    //BURAYI EKLIYORUK VERSION 5
	/*---------------------------------------------------------------------------*\
	
	

	
	\*---------------------------------------------------------------------------*/	


	//if (Pstream::myProcNo()==2){Pout << "Before proc  " << Pstream::myProcNo() << " globalIdToCf_: " << globalIdToCf_ << nl << endl;}
	
    if (Pstream::parRun())
    {
        labelHashSet sendTo;
		for (label i = 0; i < patch().size(); ++i)
		{
			if (!procBoundaryFace_[i]) continue;

			forAll(neighbProcNos_[i], j)
			{
				label proc = neighbProcNos_[i][j];
				if (proc >= 0)
				{
					sendTo.insert(proc);
				}
			}
		}


        List<labelList> sendToLists(Pstream::nProcs());
        sendToLists[Pstream::myProcNo()] = sendTo.toc();
        Pstream::gatherList(sendToLists);
        Pstream::scatterList(sendToLists);

        labelHashSet recvFrom;
        for (label src = 0; src < Pstream::nProcs(); ++src)
        {
            if (src == Pstream::myProcNo()) continue;
            const labelList& lst = sendToLists[src];
            if (lst.found(Pstream::myProcNo()))
                recvFrom.insert(src);
        }

        PstreamBuffers pBuf(Pstream::commsTypes::nonBlocking);

        forAllConstIter(labelHashSet, sendTo, it)
        {
            const label tgt = it.key();
            UOPstream os(tgt, pBuf);

            label n = 0;
            for (label i = 0; i < patch().size(); ++i)
			{
				if (procBoundaryFace_[i] && neighbProcNos_[i].found(tgt)) ++n;
			}
            os << n;

            for (label i = 0; i < patch().size(); ++i)
			{
				if (procBoundaryFace_[i] && neighbProcNos_[i].found(tgt))
				{
					label gid = globalIdOfLocalFace_[i];
					const List<label>& upListX = upwindNodesX_[i];
					const List<label>& downListX = downwindNodesX_[i];
					const List<label>& upListY = upwindNodesY_[i];
					const List<label>& downListY = downwindNodesY_[i];

					os << gid << upListX << downListX << upListY << downListY;

					os << patch().Cf()[i];

					//You’ll also have to add a send for the Y upwind/downwind coordinates if needed — this is optional for correctness if not used.
					forAll(upListX, k)
					{
						const label gU = upListX[k];
						if (globalIdToPackedIdx_.found(gU))
						{
							os << patchFaceCenters[globalIdToPackedIdx_[gU]];
						}
						else
							os << vector::zero;
					}
					forAll(downListX, k)
					{
						const label gD = downListX[k];
						if (globalIdToPackedIdx_.found(gD))
						{
							os << patchFaceCenters[globalIdToPackedIdx_[gD]];
						}
						else
							os << vector::zero;
					}
					forAll(upListY, k)
					{
						const label gU = upListY[k];
						if (globalIdToPackedIdx_.found(gU))
						{
							os << patchFaceCenters[globalIdToPackedIdx_[gU]];
						}
						else
							os << vector::zero;
					}
					forAll(downListY, k)
					{
						const label gD = downListY[k];
						if (globalIdToPackedIdx_.found(gD))
						{
							os << patchFaceCenters[globalIdToPackedIdx_[gD]];
						}
						else
							os << vector::zero;
					}
					
					
				}
			}

        }

        pBuf.finishedSends();

        forAllConstIter(labelHashSet, recvFrom, it)
        {
            const label src = it.key();
            UIPstream is(src, pBuf);

            label nRecv; is >> nRecv;
            for (label j = 0; j < nRecv; ++j)
            {
                label gid;
                List<label> upListX, downListX, upListY, downListY;
				vector CfGid;

				is >> gid >> upListX >> downListX >> upListY >> downListY >> CfGid;

                remIdToUpwindNewX_.insert(gid, upListX);
				remIdToDownwindNewX_.insert(gid, downListX);
				remIdToUpwindNewY_.insert(gid, upListY);
				remIdToDownwindNewY_.insert(gid, downListY);
				
				
                
				forAll(upwindNodesX_, i)
				{
					const List<label>& myUpList = upwindNodesX_[i];
					if (myUpList.found(gid))
						remoteIdToUpwindNodesX_.insert(globalIdOfLocalFace_[i], upListX);
						
				}

				forAll(downwindNodesX_, i)
				{
					const List<label>& myDownList = downwindNodesX_[i];
					if (myDownList.found(gid))
						remoteIdToDownwindNodesX_.insert(globalIdOfLocalFace_[i], downListX);
				}

				forAll(upwindNodesY_, i)
				{
					const List<label>& myUpY = upwindNodesY_[i];
					if (myUpY.found(gid))
						remoteIdToUpwindNodesY_.insert(globalIdOfLocalFace_[i], upListY);
				}
								
				forAll(downwindNodesY_, i)
				{
					const List<label>& myDownY = downwindNodesY_[i];
					if (myDownY.found(gid))
						remoteIdToDownwindNodesY_.insert(globalIdOfLocalFace_[i], downListY);
				}
								
				
				if (!remoteNeighbours_.found(gid))
				{
					remoteNeighbours_.insert(gid, List<label>());
				}

				remoteNeighbours_[gid].append(upListX);
				remoteNeighbours_[gid].append(downListX);
				remoteNeighbours_[gid].append(upListY);
				remoteNeighbours_[gid].append(downListY);
				
				
                globalIdToCf_.insert(gid, CfGid);
				

				//Matching
				
                forAll(upListX, k)
                {
                    vector CfUp;
                    is >> CfUp;
                    if (!globalIdToCf_.found(upListX[k]) && !globalIdToPackedIdx_.found(upListX[k]))
                        globalIdToCf_.insert(upListX[k], CfUp);
                }

                forAll(downListX, k)
                {
                    vector CfDown;
                    is >> CfDown;
                    if (!globalIdToCf_.found(downListX[k]) && !globalIdToPackedIdx_.found(downListX[k]))
                        globalIdToCf_.insert(downListX[k], CfDown);
                }
				
				 forAll(upListY, k)
                {
                    vector CfUp;
                    is >> CfUp;
                    if (!globalIdToCf_.found(upListY[k]) && !globalIdToPackedIdx_.found(upListY[k]))
                        globalIdToCf_.insert(upListY[k], CfUp);
                }

                forAll(downListY, k)
                {
                    vector CfDown;
                    is >> CfDown;
                    if (!globalIdToCf_.found(downListY[k]) && !globalIdToPackedIdx_.found(downListY[k]))
                        globalIdToCf_.insert(downListY[k], CfDown);
                }
				
				
				
				
            }
        }
    

	
	
	}
	
	
	
	remoteNeighbours2_.clear();
	
	//VERSION 14 IMPLEMENTATION (EXTENSION OF remoteNeighbours2_)
	
	
	if (Pstream::parRun())
	{
		// ----------------------
		// Build sendTo / recvFrom
		// ----------------------
		labelHashSet sendTo;
		for (label i = 0; i < patch().size(); ++i)
		{
			if (!procBoundaryFace_[i]) continue;

			forAll(neighbProcNos_[i], j)
			{
				label proc = neighbProcNos_[i][j];
				if (proc >= 0) sendTo.insert(proc);
			}
		}

		List<labelList> sendToLists(Pstream::nProcs());
		sendToLists[Pstream::myProcNo()] = sendTo.toc();
		Pstream::gatherList(sendToLists);
		Pstream::scatterList(sendToLists);

		labelHashSet recvFrom;
		for (label src = 0; src < Pstream::nProcs(); ++src)
		{
			if (src == Pstream::myProcNo()) continue;
			const labelList& lst = sendToLists[src];
			if (lst.found(Pstream::myProcNo())) recvFrom.insert(src);
		}

		// ----------------------
		// Send phase (2nd ring)
		// ----------------------
		PstreamBuffers pBuf(Pstream::commsTypes::nonBlocking);

		forAllConstIter(labelHashSet, sendTo, it)
		{
			const label tgt = it.key();
			UOPstream os(tgt, pBuf);

			// count how many faces we’ll send
			label n = 0;
			forAll(procBoundaryFace_, i)
				if (procBoundaryFace_[i] && neighbProcNos_[i].found(tgt)) ++n;
			os << n;

			forAll(procBoundaryFace_, i)
			{
				if (procBoundaryFace_[i] && neighbProcNos_[i].found(tgt))
				{
					label gid = globalIdOfLocalFace_[i];   // my global face ID
					const List<label>& firstRingLocal = neighbours_[i];

					// Send gid + number of first-ring neighbours
					os << gid << firstRingLocal.size();

					forAll(firstRingLocal, jj)
					{
                    label locId1 = firstRingLocal[jj];
                    label gid1   = globalIdOfLocalFace_[locId1];

                    const List<label>& upX   = upwindNodesX_[locId1];
                    const List<label>& downX = downwindNodesX_[locId1];
                    const List<label>& upY   = upwindNodesY_[locId1];
                    const List<label>& downY = downwindNodesY_[locId1];

                    // center of gid1 on this proc (we know locId1 is local)
                    const vector& Cf_gid1 = patchFaceCenters[locId1];

                    // Send basic lists + Cf(gid1)
                    os << gid1 << upX << downX << upY << downY << Cf_gid1;

                    // Now send Cf for each candidate in upX/downX/upY/downY,
                    // if we know it; otherwise send vector::zero as placeholder.
                    forAll(upX, k)
						{
                        const label cand = upX[k];
                        if (globalIdToPackedIdx_.found(cand))
                        {
                            os << patchFaceCenters[globalIdToPackedIdx_[cand]];
                        }
                        else if (globalIdToCf_.found(cand))
                        {
                            os << globalIdToCf_[cand];
                        }
                        else
                        {
                            os << vector::zero;
                        }
						}

                    forAll(downX, k)
						{
                        const label cand = downX[k];
                        if (globalIdToPackedIdx_.found(cand))
                        {
                            os << patchFaceCenters[globalIdToPackedIdx_[cand]];
                        }
                        else if (globalIdToCf_.found(cand))
                        {
                            os << globalIdToCf_[cand];
                        }
                        else
                        {
                            os << vector::zero;
                        }
						}

                    forAll(upY, k)
						{
                        const label cand = upY[k];
                        if (globalIdToPackedIdx_.found(cand))
                        {
                            os << patchFaceCenters[globalIdToPackedIdx_[cand]];
                        }
                        else if (globalIdToCf_.found(cand))
                        {
                            os << globalIdToCf_[cand];
                        }
                        else
                        {
                            os << vector::zero;
                        }
						}

                    forAll(downY, k)
						{
                        const label cand = downY[k];
                        if (globalIdToPackedIdx_.found(cand))
                        {
                            os << patchFaceCenters[globalIdToPackedIdx_[cand]];
                        }
                        else if (globalIdToCf_.found(cand))
                        {
                            os << globalIdToCf_[cand];
                        }
                        else
                        {
                            os << vector::zero;
                        }
						}
					}

				}
			}
			
			
			
			
			
			
		}


		pBuf.finishedSends();

		// ----------------------
		// Receive phase (2nd ring)
		// ----------------------
		forAllConstIter(labelHashSet, recvFrom, it)
		{
			const label src = it.key();
			UIPstream is(src, pBuf);

			label nRecv; 
			is >> nRecv;

			for (label j = 0; j < nRecv; ++j)
			{
				label gid; 
				is >> gid;

				label nFirst; 
				is >> nFirst;

				for (label jj = 0; jj < nFirst; ++jj)
				{
					label gid1;
					List<label> upX, downX, upY, downY;
					vector Cf_gid1;

					is >> gid1 >> upX >> downX >> upY >> downY >> Cf_gid1;

					// Make sure we store Cf for this first-ring neighbour (gid1)
					if (!globalIdToCf_.found(gid1))
					{
						globalIdToCf_.insert(gid1, Cf_gid1);
					}

					// Read Cf for all candidates in upX/downX/upY/downY
					// in the same order they were sent:

					// upX
					List<vector> Cf_upX(upX.size());
					forAll(upX, k)
					{
						is >> Cf_upX[k];
						const label cand = upX[k];
						if (!globalIdToCf_.found(cand) && (mag(Cf_upX[k]) > VSMALL) && !globalIdToPackedIdx_.found(cand))
						{
							globalIdToCf_.insert(cand, Cf_upX[k]);
						}
					}

					// downX
					List<vector> Cf_downX(downX.size());
					forAll(downX, k)
					{
						is >> Cf_downX[k];
						const label cand = downX[k];
						if (!globalIdToCf_.found(cand) && (mag(Cf_downX[k]) > VSMALL) && !globalIdToPackedIdx_.found(cand))
						{
							globalIdToCf_.insert(cand, Cf_downX[k]);
						}
					}

					// upY
					List<vector> Cf_upY(upY.size());
					forAll(upY, k)
					{
						is >> Cf_upY[k];
						const label cand = upY[k];
						if (!globalIdToCf_.found(cand) && (mag(Cf_upY[k]) > VSMALL) && !globalIdToPackedIdx_.found(cand))
						{
							globalIdToCf_.insert(cand, Cf_upY[k]);
						}
					}

					// downY
					List<vector> Cf_downY(downY.size());
					forAll(downY, k)
					{
						is >> Cf_downY[k];
						const label cand = downY[k];
						if (!globalIdToCf_.found(cand) && (mag(Cf_downY[k]) > VSMALL) && !globalIdToPackedIdx_.found(cand))
						{
							globalIdToCf_.insert(cand, Cf_downY[k]);
						}
					}

					// Now build 2nd-ring neighbour lists, WITHOUT filtering by globalIdToCf_
					if (!remoteNeighbours_.found(gid1))
					{   // store as 2nd-ring neighbours
						if (!remoteNeighbours2_.found(gid1))
						{
							remoteNeighbours2_.insert(gid1, List<label>());
						}

						// upX
						forAll(upX, kk)
						{
							label cand = upX[kk];
							if (!remoteNeighbours2_[gid1].found(cand))
							{
								remoteNeighbours2_[gid1].append(cand);
							}
						}

						// downX
						forAll(downX, kk)
						{
							label cand = downX[kk];
							if (!remoteNeighbours2_[gid1].found(cand))
							{
								remoteNeighbours2_[gid1].append(cand);
							}
						}

						// upY
						forAll(upY, kk)
						{
							label cand = upY[kk];
							if (!remoteNeighbours2_[gid1].found(cand))
							{
								remoteNeighbours2_[gid1].append(cand);
							}
						}

						// downY
						forAll(downY, kk)
						{
							label cand = downY[kk];
							if (!remoteNeighbours2_[gid1].found(cand))
							{
								remoteNeighbours2_[gid1].append(cand);
							}
						}

						// Also store directional lists for gid1
						remIdToUpwindNewX_.insert(gid1, upX);
						remIdToDownwindNewX_.insert(gid1, downX);
						remIdToUpwindNewY_.insert(gid1, upY);
						remIdToDownwindNewY_.insert(gid1, downY);
					}
				}
			}
		}
		
	
	}


    Info << "VELOCITY BC Upwind and downwind nodes in X and Y identified." << endl;

	
	//if (Pstream::myProcNo()==2) Pout << "After " << Pstream::myProcNo() << " globalIdToCf_: " << globalIdToCf_ << nl << endl;
	//for (label i = 0; i < patch().size(); ++i) Pout << "Proc " << Pstream::myProcNo() << " upwindNodesX_ at: " << globalIdOfLocalFace_[i] << " " << upwindNodesX_[i] << nl << endl;
	//for (label i = 0; i < patch().size(); ++i) Pout << "Proc " << Pstream::myProcNo() << " upwindNodesY_ at: " << globalIdOfLocalFace_[i] << " " << upwindNodesY_[i] << nl << endl;
	
	



	if (Pstream::parRun())
	{
		neededGIDs_.clear();

		using TableType = HashTable<labelList, label>;

		// 1) First-ring neighbours (remoteNeighbours_)
		forAllConstIter(TableType, remoteNeighbours_, it)
		{
			const labelList& myGIDs = it();  // Faces that remote proc expects from me
			forAll(myGIDs, i)
			{
				neededGIDs_.insert(myGIDs[i]);
			}
		}

		// 2) Second-ring neighbours (remoteNeighbours2_)
		forAllConstIter(TableType, remoteNeighbours2_, it2)
		{
			const labelList& myGIDs2 = it2();
			forAll(myGIDs2, i)
			{
				// labelHashSet / set semantics -> duplicates automatically ignored
				neededGIDs_.insert(myGIDs2[i]);
			}
		}

		// Communicate to all processors
		neededByProc_.setSize(Pstream::nProcs());
		neededByProc_[Pstream::myProcNo()] = neededGIDs_;

		Pstream::gatherList(neededByProc_);
		Pstream::scatterList(neededByProc_);
	}
		//WRITING pART
	

	/*---------------------------------------------------------------------------*\
	if (Pstream::myProcNo() == 2) 
	{
		forAllConstIter(Foam::Map<Foam::vector>, globalIdToCf_, it)
		{
			const label gid = it.key();
			const vector& Cf = it();
			if (!globalIdToPackedIdx_.found(gid)) Pout << "After proc " << Pstream::myProcNo()<< " globalIdToCf_[" << gid << "]: " << Cf << nl << endl;
		
		}
	}
	
	
	\*---------------------------------------------------------------------------*/	
		
		
	//Pout << "After proc " << Pstream::myProcNo()<< " neededByProc_:" << neededByProc_  << nl << endl;
	
	

	
    windsCalculated_ = true;
}

void Foam::waveCurrentPotential3DFvPatchScalarField::detectFDSchemes()
{
    if (schemesDetected_) return;

    const vectorField& Cf = patch().Cf();
    const scalarField x = Cf.component(vector::X);
    const scalarField y = Cf.component(vector::Y);
    const label nFaces = patch().size();

    schemeCodeX_.setSize(nFaces, -1);
    schemeCodeY_.setSize(nFaces, -1);
    
	interpNeededX_.setSize(nFaces, FixedList<bool,2>(false));
	interpNeededY_.setSize(nFaces, FixedList<bool,2>(false));

	
	const volVectorField& Ucur = db().lookupObjectRef<volVectorField>("Ucur");
	const vectorField& Ucur_p = Ucur.boundaryField()[patch().index()];
	
	//vectorField Ucur_p(nFaces);
	//for (label i = 0; i < nFaces; ++i) Ucur_p[i] = vector(1.0, 1.0, 0);

    for (label i = 0; i < nFaces; ++i)
	{
		label schemeX = -1;
		label schemeY = -1;

		const scalar uX = Ucur_p[i].x();
		const scalar uY = Ucur_p[i].y();

		// === X-direction scheme detection
		{
			const word dir = "x";
			if (
				(uX > SMALL && (
					tryUpwind2nd(*this, i, schemeX, Cf, x, y, upwindNodesX_, remoteIdToUpwindNodesX_, dir) ||
					tryCentered2nd(*this, i, schemeX, Cf, x, y, upwindNodesX_, downwindNodesX_, dir) ||
					tryDownwind2nd(*this, i, schemeX, Cf, x, y, downwindNodesX_, remoteIdToDownwindNodesX_, dir)
				)) ||
				(uX < -SMALL && (
					tryDownwind2nd(*this, i, schemeX, Cf, x, y, downwindNodesX_, remoteIdToDownwindNodesX_, dir) ||
					tryCentered2nd(*this, i, schemeX, Cf, x, y, upwindNodesX_, downwindNodesX_, dir) ||
					tryUpwind2nd(*this, i, schemeX, Cf, x, y, upwindNodesX_, remoteIdToUpwindNodesX_, dir)
				)) ||
				(mag(uX) <= SMALL && (
					tryCentered2nd(*this, i, schemeX, Cf, x, y, upwindNodesX_, downwindNodesX_, dir) ||
					tryUpwind2nd(*this, i, schemeX, Cf, x, y, upwindNodesX_, remoteIdToUpwindNodesX_, dir) ||
					tryDownwind2nd(*this, i, schemeX, Cf, x, y, downwindNodesX_, remoteIdToDownwindNodesX_, dir)
				))
			)
			{
				schemeCodeX_[i] = schemeX;
			}
			else
			{
				const List<label>& L = upwindNodesX_[i];
				if (L.size())
				{
					auto dxdyL = dxdy(Cf, globalIdToPackedIdx_, globalIdToCf_, L[0], x[i], y[i]);
					if (aligned(dxdyL.first, dxdyL.second, "x"))
						schemeCodeX_[i] = 4;
				}

				const List<label>& R = downwindNodesX_[i];
				if (R.size())
				{
					auto dxdyR = dxdy(Cf, globalIdToPackedIdx_, globalIdToCf_, R[0], x[i], y[i]);
					if (aligned(dxdyR.first, dxdyR.second, "x"))
						schemeCodeX_[i] = 5;
				}
			}
		}

		// === Y-direction scheme detection
		{
			const word dir = "y";
			if (
				(uY > SMALL && (
					tryUpwind2nd(*this, i, schemeY, Cf, x, y, upwindNodesY_, remoteIdToUpwindNodesY_, dir) ||
					tryCentered2nd(*this, i, schemeY, Cf, x, y, upwindNodesY_, downwindNodesY_, dir) ||
					tryDownwind2nd(*this, i, schemeY, Cf, x, y, downwindNodesY_, remoteIdToDownwindNodesY_, dir)
				)) ||
				(uY < -SMALL && (
					tryDownwind2nd(*this, i, schemeY, Cf, x, y, downwindNodesY_, remoteIdToDownwindNodesY_, dir) ||
					tryCentered2nd(*this, i, schemeY, Cf, x, y, upwindNodesY_, downwindNodesY_, dir) ||
					tryUpwind2nd(*this, i, schemeY, Cf, x, y, upwindNodesY_, remoteIdToUpwindNodesY_, dir)
				)) ||
				(mag(uY) <= SMALL && (
					tryCentered2nd(*this, i, schemeY, Cf, x, y, upwindNodesY_, downwindNodesY_, dir) ||
					tryUpwind2nd(*this, i, schemeY, Cf, x, y, upwindNodesY_, remoteIdToUpwindNodesY_, dir) ||
					tryDownwind2nd(*this, i, schemeY, Cf, x, y, downwindNodesY_, remoteIdToDownwindNodesY_, dir)
				))
			)
			{
				schemeCodeY_[i] = schemeY;
			}
			else
			{
				const List<label>& L = upwindNodesY_[i];
				if (L.size())
				{
					auto dxdyL = dxdy(Cf, globalIdToPackedIdx_, globalIdToCf_, L[0], x[i], y[i]);
					if (aligned(dxdyL.first, dxdyL.second, "y")){
						schemeCodeY_[i] = 4;
						WarningInFunction << "Face " << globalIdOfLocalFace_[i] << ": First Order Detected !!!!" << endl;
					}
				}

				const List<label>& R = downwindNodesY_[i];
				if (R.size())
				{
					auto dxdyR = dxdy(Cf, globalIdToPackedIdx_, globalIdToCf_, R[0], x[i], y[i]);
					if (aligned(dxdyR.first, dxdyR.second, "y")){
						schemeCodeY_[i] = 5;
						WarningInFunction << "Face " << globalIdOfLocalFace_[i] << ": First Order Detected !!!" << endl;
					}
				}
			}
		}

		// Fallback (optional)
		//if (schemeCodeX_[i] == -1) schemeCodeX_[i] = 0; // Or other fallback
		//if (schemeCodeY_[i] == -1) schemeCodeY_[i] = 0; // Or other fallback
			
	}

	//for (label i = 0; i < nFaces; ++i) if( schemeCodeY_[i] !=6) {Pout << "proc " << Pstream::myProcNo() << " Schemes implemented at " << globalIdOfLocalFace_[i] << " " << schemeCodeY_[i] << nl;}
	//for (label i = 0; i < nFaces; ++i) if( schemeCodeX_[i] !=6) {Pout << "proc " << Pstream::myProcNo() << " Schemes implemented at " << globalIdOfLocalFace_[i] << " " << schemeCodeX_[i] << nl;}
	//for (label i = 0; i < nFaces; ++i)  {Pout << "proc " << Pstream::myProcNo() << " interpNeededX_ at " << globalIdOfLocalFace_[i] << " " << interpNeededX_[i] << nl;}
	//for (label i = 0; i < nFaces; ++i) Pout << "proc " << Pstream::myProcNo() << " interpNeededY_ at " << globalIdOfLocalFace_[i] << " " << interpNeededY_[i] << nl;
    
	Info << "Scheme detection complete." << endl;
	
	initCandidates();
	
	
	schemesDetected_ = true;
}


void Foam::waveCurrentPotential3DFvPatchScalarField::initCandidates()
{
    // If both already built, return
    if (candidatesCalculated_)
    {
        return;
    }

    const label nFaces = patch().size();
	lsMergedCandidates_.setSize(nFaces);
	const vectorField& Cf = patch().Cf();
    const scalarField xComponents = Cf.component(vector::X);
    const scalarField yComponents = Cf.component(vector::Y);
	
	
	
	// ---- Utility: compute distance between two global IDs ----
	auto getDistBetweenGIDs = [&](label gid1, label gid2) -> scalar
	{
		vector C1(Zero), C2(Zero);

		// --- Fetch coordinates for gid1 ---
		if (globalIdToPackedIdx_.found(gid1))
			C1 = Cf[globalIdToPackedIdx_[gid1]];
		else if (globalIdToCf_.found(gid1))
			C1 = globalIdToCf_[gid1];
		else
		{
			FatalErrorInFunction
				<< "getDistBetweenGIDs(): gid1 " << gid1
				<< " not found in local or remote Cf maps." << nl << abort(FatalError);
		}

		// --- Fetch coordinates for gid2 ---
		if (globalIdToPackedIdx_.found(gid2))
			C2 = Cf[globalIdToPackedIdx_[gid2]];
		else if (globalIdToCf_.found(gid2))
			C2 = globalIdToCf_[gid2];
		else
		{
			FatalErrorInFunction
				<< "getDistBetweenGIDs(): gid2 " << gid2
				<< " not found in local or remote Cf maps." << nl << abort(FatalError);
		}

		// --- Return Euclidean distance ---
		return mag(C2 - C1);
	};

	// ---- Utility: remove candidates based on Uy criterion ----
	auto removeCandidatesByUx_13D = [&](const List<label>& candidates,
									scalar UxRef,
									label selfGid) -> List<label>
	{
		DynamicList<label> filtered;

		// Reference coordinate of self face
		vector Corig(Zero);
		if (globalIdToPackedIdx_.found(selfGid))
			Corig = Cf[globalIdToPackedIdx_[selfGid]];
		else if (globalIdToCf_.found(selfGid))
			Corig = globalIdToCf_[selfGid];
		else
		{
			FatalErrorInFunction
				<< "removeCandidatesByUx(): Corig selfGid " << selfGid
				<< " not found in local or remote Cf maps." << endl;
			return candidates;
		}

		// Loop over all candidate faces
		forAll(candidates, j)
		{
			const label gid = candidates[j];

			vector Ccand(Zero);
			if (globalIdToPackedIdx_.found(gid))
				Ccand = Cf[globalIdToPackedIdx_[gid]];
			else if (globalIdToCf_.found(gid))
				Ccand = globalIdToCf_[gid];
			else
				FatalErrorInFunction << "removeCandidatesByUx(): Ccand selfGid " << selfGid << " not found in local or remote Cf maps." << endl; // Skip if Cf not found anywhere
				

			// --- Geometric filtering criterion based on Uy direction ---
			if (UxRef > 0.0)
			{
				// Keep candidates below (smaller x)
				if (Corig.x() >= (Ccand.x()-0.001))
					filtered.append(gid);
			}
			else if (UxRef < 0.0)
			{
				// Keep candidates above (larger x)
				if (Corig.x() <= (Ccand.x()+0.001))
					filtered.append(gid);
			}
			else
			{
				// Ux ≈ 0 → keep all
				filtered.append(gid);
			}
		}

		if (filtered.empty())
		{
			FatalErrorInFunction
				<< "removeCandidatesByUx(): No valid candidates kept for selfGid "
				<< selfGid << " (UxRef=" << UxRef << ")." << endl;
		}

		return List<label>(filtered);
	};


	
	// ---- Utility: merge with uniqueness ----
    auto mergeUniquewoSelf = [](const List<label>& a, const List<label>& b,const label selfGid)
    {
        DynamicList<label> merged;

        forAll(a, i)
        {
            if (a[i] != selfGid && !merged.found(a[i])) 
                merged.append(a[i]);
        }

        forAll(b, j)
        {
            if (b[j] != selfGid && !merged.found(b[j])) //
                merged.append(b[j]);
        }

        return List<label>(merged);
    };
	
	
    // ---- Utility: safe candidate getter ----
    auto getCandidates = [&](label selfGid, label gid, scalar Uy) -> const labelList&
	{
		// static container to safely return references (cannot return local variable)
		static labelList sendt;

		// Clear previous contents to avoid accumulation
		sendt.clear();

		// ========================
		// SchemeCode 1 – General (X/Y combined)
		// ========================
		if (mag(Uy) <= VSMALL)
		{
			
			// --- Start with local neighbours ---
			if (globalIdToPackedIdx_.found(gid))
			{
				const labelList& local = neighboursGlobal_[globalIdToPackedIdx_[gid]];
				sendt.append(local);
			}

			// --- Add 1st-ring remote neighbours if exist ---
			if (remoteNeighbours_.found(gid))
			{
				const labelList& rem1 = remoteNeighbours_[gid];
				sendt.append(rem1);
			}

			// --- Add 2nd-ring remote neighbours if exist ---
			if (remoteNeighbours2_.found(gid))
			{
				const labelList& rem2 = remoteNeighbours2_[gid];
				sendt.append(rem2);
			}


			return sendt;
		}

		// ========================
		// SchemeCode 6 – Upwind Y-direction
		// ========================
		if (Uy > VSMALL)
		{
			

			// ---------- Local branch ----------
			if (globalIdToPackedIdx_.found(gid))
			{
				const labelList& local = upwindNodesY_[globalIdToPackedIdx_[gid]];
				sendt.append(local);

				// --- 2nd ring: try local of the first upwind neighbour
				if (!local.empty())
				{
					// Find nearest candidate among 'local' faces to selfGid
					label nearest = local[0];
					scalar minDist = GREAT;

					forAll(local, j)
					{
						scalar d = getDistBetweenGIDs(gid, local[j]);
						if (d < minDist)
						{
							minDist = d;
							nearest = local[j];
						}
					}

					if (globalIdToPackedIdx_.found(nearest))
					{
						const labelList& local2 = upwindNodesY_[globalIdToPackedIdx_[nearest]];
						sendt.append(local2);
					}
					else if (remIdToUpwindNewY_.found(nearest))
					{
						const labelList& rem = remIdToUpwindNewY_[nearest];
						sendt.append(rem);
					}
				}

				return sendt;
			}

			// ---------- Remote branch ----------
			if (remIdToUpwindNewY_.found(gid))
			{
				const labelList& rem1 = remIdToUpwindNewY_[gid];
				sendt.append(rem1);

				// --- 2nd ring: remote of the first upwind neighbour
				if (!rem1.empty())
				{
					const label first = rem1[0];

					if (remIdToUpwindNewY_.found(first))
					{
						const labelList& rem2 = remIdToUpwindNewY_[first];
						sendt.append(rem2);
					}
					else if (globalIdToPackedIdx_.found(first))
					{
						const labelList& local2 = upwindNodesY_[globalIdToPackedIdx_[first]];
						sendt.append(local2);
					}
				}

				return sendt;
			}

			
		}

		// ========================
		// SchemeCode 7 – Downwind Y-direction
		// ========================
		if (Uy < -VSMALL)
		{
			static labelList sendt;
			sendt.clear();

			// ---------- Local branch ----------
			if (globalIdToPackedIdx_.found(gid))
			{
				const labelList& local = downwindNodesY_[globalIdToPackedIdx_[gid]];
				sendt.append(local);

				// --- 2nd ring: try the first downwind neighbour
				if (!local.empty())
				{
					label nearest = local[0];
					scalar minDist = GREAT;

					forAll(local, j)
					{
						scalar d = getDistBetweenGIDs(gid, local[j]);
						if (d < minDist)
						{
							minDist = d;
							nearest = local[j];
						}
					}

					if (globalIdToPackedIdx_.found(nearest))
					{
						const labelList& local2 = downwindNodesY_[globalIdToPackedIdx_[nearest]];
						sendt.append(local2);
					}
					else if (remIdToDownwindNewY_.found(nearest))
					{
						const labelList& rem = remIdToDownwindNewY_[nearest];
						sendt.append(rem);
					}
				}

				return sendt;
			}

			// ---------- Remote branch ----------
			if (remIdToDownwindNewY_.found(gid))
			{
				const labelList& rem1 = remIdToDownwindNewY_[gid];
				sendt.append(rem1);

				// --- 2nd ring: remote of the first downwind neighbour
				if (!rem1.empty())
				{
					const label first = rem1[0];

					if (remIdToDownwindNewY_.found(first))
					{
						const labelList& rem2 = remIdToDownwindNewY_[first];
						sendt.append(rem2);
					}
					else if (globalIdToPackedIdx_.found(first))
					{
						const labelList& local2 = downwindNodesY_[globalIdToPackedIdx_[first]];
						sendt.append(local2);
					}
				}

				return sendt;
			}

			
		}

		
		FatalErrorInFunction << "should have return something so far for selfGid " << selfGid << abort(FatalError);
				
				
		
		static const labelList dummy;
		return dummy;
	};
    
	const volVectorField& Ucur = db().lookupObjectRef<volVectorField>("Ucur");
	const vectorField& Ucur_p = Ucur.boundaryField()[patch().index()];

                    
	for (label i = 0; i < nFaces; ++i)
	{
		
		if (ownerHasBodyFace_[i]) 
		{		
			const scalar uY = Ucur_p[i].y();
			//sconst scalar uX = Ucur_p[i].x();
			List<label> mergedCands;
			const label selfGid = globalIdOfLocalFace_[i];
			
			if (ownerHasFirstBodyFace_[i])
			{	
				
				
				//if (selfGid == 205) Pout << "selfGid:" << selfGid << " has ringNodeX_[i]" << ringNodeX_[i] << endl;
				if (ringNodeX_.found(i))  
				{	
					
					
					const List<label>& candsSelf = getCandidates(selfGid, selfGid, uY);
					mergedCands = mergeUniquewoSelf(mergedCands, candsSelf, selfGid);label conn = ringNodeX_[i];          // single neighbor
					mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn), selfGid);
					
					
					
					//if (selfGid == 205) Pout << "selfGid:" << selfGid << " has conn" << conn << endl;
				
					
					if (globalIdToPackedIdx_.found(conn))
					{	
						const labelList& ring = neighboursGlobal_[globalIdToPackedIdx_[conn]];
						mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn), selfGid);
						mergedCands = mergeUniquewoSelf(mergedCands, ring , selfGid);
						lsMergedCandidates_[i]=mergedCands;
						//if (selfGid == 205) Pout << "selfGid:" << selfGid << " has mergedCands" << mergedCands << endl;
				
						
					}
					else
					{
						mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn), selfGid);
						lsMergedCandidates_[i]=mergedCands;
					}	
					
					if (ringNodeX_.found(globalIdToPackedIdx_[conn]))
					{
						label conn2 = ringNodeX_[globalIdToPackedIdx_[conn]];
						mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn2), selfGid);
						
						const labelList& ring = neighboursGlobal_[globalIdToPackedIdx_[conn2]];
						
						//const List<label>& candsX2 = getCandidates(selfGid,conn2,uY); 
						mergedCands = mergeUniquewoSelf(mergedCands, ring , selfGid);
						//if (selfGid == 205) Pout << "selfGid:" << selfGid << " has mergedCands" << mergedCands << endl;
						lsMergedCandidates_[i]=mergedCands;
					
					
						if (ringNodeX_.found(globalIdToPackedIdx_[conn2]))
						{
							label conn3 = ringNodeX_[globalIdToPackedIdx_[conn2]];
							mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn3), selfGid);
							
							const labelList& ring = neighboursGlobal_[globalIdToPackedIdx_[conn3]]; //V12
							
							//const List<label>& candsX2 = getCandidates(selfGid,conn3,uY); //V11
							mergedCands = mergeUniquewoSelf(mergedCands, ring , selfGid);
							//if (selfGid == 205) Pout << "selfGid:" << selfGid << " has mergedCands" << mergedCands << endl;
							lsMergedCandidates_[i]=mergedCands;
							
							
							
							//13D burayi ekledik
							if (ringNodeX_.found(globalIdToPackedIdx_[conn3]))
							{
								label conn4 = ringNodeX_[globalIdToPackedIdx_[conn3]];
								mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn4), selfGid);
							
								lsMergedCandidates_[i]=mergedCands;
							
							}
							//13D burayi ekledik
							
						}
					}
					const List<label>& L = upwindNodesX_[i];
					
					//if (selfGid == 191) Pout << "selfGid:" << selfGid << " has L" << L << endl;
					
					if (!L.empty())
					{
						label up1 = L[0];
						const List<label>& candsSelf = getCandidates(selfGid,selfGid,uY); 
						mergedCands = mergeUniquewoSelf(mergedCands,candsSelf,selfGid);
						
						mergedCands = mergeUniquewoSelf(mergedCands,L,selfGid);
						const List<label>& candsX1 = getCandidates(selfGid,up1,uY); 
						mergedCands = mergeUniquewoSelf(mergedCands,candsX1,selfGid);
						
						List<label> up2List;

						if (globalIdToPackedIdx_.found(up1))
						{
						
						up2List = upwindNodesX_[globalIdToPackedIdx_[up1]];
						
						}
						else
						{
							const auto it = remoteIdToUpwindNodesX_.find(selfGid);
							if (it != remoteIdToUpwindNodesX_.end())
								up2List = it();
						}
						mergedCands = mergeUniquewoSelf(mergedCands,up2List,selfGid);
				
						
						lsMergedCandidates_[i]=mergedCands;
						
					}	
			
					lsMergedCandidates_[i]=mergedCands;
						
				}
				
			}
			
			if (ring2NodeX_.found(i))
			{
				label conn = ring2NodeX_[i];          // single neighbor
				mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn), selfGid);
				const List<label>& candsSelf = getCandidates(selfGid, selfGid, uY);
				mergedCands = mergeUniquewoSelf(mergedCands, candsSelf, selfGid);
				
				
				if (selfGid == 707) Pout << "selfGid:" << selfGid << " has conn" << conn << endl;
				
				
				if (globalIdToPackedIdx_.found(conn))
				{	
					mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn), selfGid);
					const List<label>& candsX1 = getCandidates(selfGid,conn,uY); 
					//const labelList& ring = neighboursGlobal_[globalIdToPackedIdx_[conn]];
					
					mergedCands = mergeUniquewoSelf(mergedCands, candsX1 , selfGid);
						
					
					
				}
				
				if (ringNodeX_.found(globalIdToPackedIdx_[conn]))
				{
					label conn2 = ringNodeX_[globalIdToPackedIdx_[conn]];
					mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn2), selfGid);
					const labelList& ring2 = neighboursGlobal_[globalIdToPackedIdx_[conn2]];
					mergedCands = mergeUniquewoSelf(mergedCands, ring2 , selfGid);
					if (ringNodeX_.found(globalIdToPackedIdx_[conn2]))
					{
						label conn3 = ringNodeX_[globalIdToPackedIdx_[conn2]];
						mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn3), selfGid);
						const labelList& ring3 = neighboursGlobal_[globalIdToPackedIdx_[conn3]];
						mergedCands = mergeUniquewoSelf(mergedCands, ring3 , selfGid);
								
					
						//13D burayi ekledik
						if (ringNodeX_.found(globalIdToPackedIdx_[conn3]))
						{
							label conn4 = ringNodeX_[globalIdToPackedIdx_[conn3]];
							mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn4), selfGid);
							
							lsMergedCandidates_[i]=mergedCands;
							
						}
						//13D burayi ekledik
					
					
					
					}
					
					
					
					
				}
				const List<label>& L = upwindNodesX_[i];
				if (selfGid == 707) Pout << "selfGid:" << selfGid << " has mergedCands" << mergedCands << endl;
					
				if (!L.empty())	
				{	label up1 = L[0];
					
					const List<label>& candsSelf = getCandidates(selfGid,selfGid,uY); 
					if (selfGid == 973)  Pout << "selfGid:" << selfGid << " has candsSelf " << candsSelf << endl;
				
					mergedCands = mergeUniquewoSelf(mergedCands,L,selfGid);
					mergedCands = mergeUniquewoSelf(mergedCands,candsSelf,selfGid);
					
					
					
					const List<label>& candsX1 = getCandidates(selfGid,up1,uY); 
					mergedCands = mergeUniquewoSelf(mergedCands,candsX1,selfGid);
				
					List<label> up2List;

					if (globalIdToPackedIdx_.found(up1))
					{
					
					up2List = upwindNodesX_[globalIdToPackedIdx_[up1]];
					
					}
					else
					{
						const auto it = remoteIdToUpwindNodesX_.find(selfGid);
						if (it != remoteIdToUpwindNodesX_.end())
							up2List = it();
					}
					mergedCands = mergeUniquewoSelf(mergedCands,up2List,selfGid);
			
					
					lsMergedCandidates_[i]=mergedCands;
				}
			}
			
			if (selfGid == 678)  Pout << "selfGid:" << selfGid << " has mergedCands " << mergedCands << endl;
			
			
			if (mergedCands.empty())
			{
				const List<label>& L = upwindNodesX_[i];

				if (L.empty())
				{
					WarningInFunction
						<< "Face " << selfGid << " has no upwindNodesX_. Skipping candidate build." << endl;
					continue; // nothing to merge safely
				}

				label up1 = L[0];
				List<label> up2List;

				if (globalIdToPackedIdx_.found(up1))
				{
					up2List = upwindNodesX_[globalIdToPackedIdx_[up1]];
				}
				else
				{
					const auto it = remoteIdToUpwindNodesX_.find(selfGid);
					if (it != remoteIdToUpwindNodesX_.end())
					{
						up2List = it();
					}
					else
					{
						WarningInFunction
							<< "Face " << selfGid << " missing remote second upwind node — skipping second order."
							<< endl;
						up2List.clear(); // fallback to 1st order
					}
				}

				if (up2List.empty())
				{
					WarningInFunction
						<< "Face " << selfGid << " has empty up2List — fallback to first order." << endl;
					lsMergedCandidates_[i] = L;
					continue;
				}

				label up2 = up2List[0];

				const List<label>& candsSelf = getCandidates(selfGid, selfGid, uY);
				mergedCands = mergeUniquewoSelf(candsSelf, L, selfGid);
				mergedCands = mergeUniquewoSelf(mergedCands, up2List, selfGid);

				const List<label>& candsX1 = getCandidates(selfGid, up1, uY);
				mergedCands = mergeUniquewoSelf(mergedCands, candsX1, selfGid);

				const List<label>& candsX2 = getCandidates(selfGid, up2, uY);
				mergedCands = mergeUniquewoSelf(mergedCands, candsX2, selfGid);

				lsMergedCandidates_[i] = mergedCands;
			}

		
			//mergedCands=removeCandidatesByUx(mergedCands,uX,selfGid);
			lsMergedCandidates_[i] = mergedCands;
		}
		
	}	
  
    //outputter();               
    Info << "PreParV14D used " << endl;               
	candidatesCalculated_ = true;
    Info << "Candidates for Mq built on " << patch().name() << endl;
         
}


