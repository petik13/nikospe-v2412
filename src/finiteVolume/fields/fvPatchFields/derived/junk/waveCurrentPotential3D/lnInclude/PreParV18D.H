/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
   Parallel-aware Neighbor Calculation using Global Mesh Files (V18D)
\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "IOobject.H"
#include "PrimitivePatch.H"
#include "SortableList.H"

void Foam::waveCurrentPotential3DFvPatchScalarField::calcNeigboursV3()
{
    if (neigboursCalculated_) return;

    Info << "PreParV18D: Calculating neighbors using Global Mesh Files strategy..." << endl;

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    const word patchName = patch().name();

    // ---------------------------------------------------------
    // 1. Read 'boundary' file to find the patch info (start, nFaces) globally
    // ---------------------------------------------------------
    label globalStart = -1;
    label globalSize = -1;
    {
        IOobject ioBound("boundary", mesh.time().constant(), "polyMesh",
                         mesh.time(), IOobject::MUST_READ, IOobject::NO_WRITE);
        
        if (!ioBound.typeHeaderOk<polyBoundaryMesh>(false))
        {
            FatalErrorInFunction 
                << "Cannot find boundary file in " << ioBound.path() 
                << abort(FatalError);
        }

        // Parse manually
        IFstream is(ioBound.objectPath());
        
        // Skip header if present (FoamFile { ... })
        token t(is);
        if (t.isWord() && t.wordToken() == "FoamFile")
        {
            dictionary header(is);
        }
        else
        {
            is.putBack(t);
        }

        PtrList<entry> patches(is);
        
        forAll(patches, i)
        {
            if (patches[i].keyword() == patchName)
            {
                const dictionary& dict = patches[i].dict();
                globalSize = readLabel(dict.lookup("nFaces"));
                globalStart = readLabel(dict.lookup("startFace"));
                Info << "Found global patch '" << patchName 
                     << "' start=" << globalStart << " size=" << globalSize << endl;
                break;
            }
        }
    }

    if (globalStart == -1 || globalSize == -1)
    {
        FatalErrorInFunction << "Could not find patch '" << patchName 
            << "' in global boundary file." << abort(FatalError);
    }

    // ---------------------------------------------------------
    // 2. Read Global 'faces' (subset for patch)
    // ---------------------------------------------------------
    List<face> globalPatchFaces(globalSize);
    {
        IOobject ioFaces("faces", mesh.time().constant(), "polyMesh",
                         mesh.time(), IOobject::MUST_READ, IOobject::NO_WRITE);
        
        // We read all faces because IOstream seeking for binary/compressed is hard
        // This assumes master or local reads - Concurrent read on FS usually OK
        // Using faceIOList to read standard faces file
        faceIOList allFaces(ioFaces);
        
        if (allFaces.size() < globalStart + globalSize)
        {
             FatalErrorInFunction << "Global faces file too short." << abort(FatalError);
        }

        for(label i=0; i<globalSize; ++i)
        {
            globalPatchFaces[i] = allFaces[globalStart + i];
        }
        Info << "Read " << globalSize << " global patch faces." << endl;
    }

    // ---------------------------------------------------------
    // 3. Read Global 'points'
    // ---------------------------------------------------------
    pointIOField globalPoints
    (
        IOobject("points", mesh.time().constant(), "polyMesh",
                 mesh.time(), IOobject::MUST_READ, IOobject::NO_WRITE)
    );
    Info << "Read " << globalPoints.size() << " global points." << endl;

    Info << "DEBUG: CHECKPOINT A - Constructing Global PrimitivePatch..." << endl;

    // ---------------------------------------------------------
    // 4. Build Global PrimitivePatch
    // ---------------------------------------------------------
    // PrimitivePatch<FaceList, PointField>
    PrimitivePatch<List<face>, pointField> globalPatch
    (
        globalPatchFaces,
        globalPoints
    );

    Info << "DEBUG: CHECKPOINT B - Calculation Global Face Centres..." << endl;
    const vectorField& globalCf = globalPatch.faceCentres();
    
    Info << "DEBUG: CHECKPOINT C - Calculating Global Face-Face Topology..." << endl;
    const labelListList& globalFaceFaces = globalPatch.faceFaces();

    // ---------------------------------------------------------
    // 5. Build Global-to-Global Stencils
    // ---------------------------------------------------------
    Info << "DEBUG: CHECKPOINT D - Calculating Stencils and Rays..." << endl;

    List<labelList> gUpwindX(globalSize), gDownwindX(globalSize);
    List<labelList> gUpwindY(globalSize), gDownwindY(globalSize);
    List<labelList> gRing(globalSize); // Vertex connections

    // Build vertex connectivity (Ring)
    const labelListList& globalPointFaces = globalPatch.pointFaces();

    Info << "DEBUG: CHECKPOINT E - Entering Loop over " << globalSize << " global faces..." << endl;
    for(label i=0; i<globalSize; ++i)
    {
        const vector& C = globalCf[i];
        
        // Neighbor (Face-to-Face) Analysis
        const labelList& nbrs = globalFaceFaces[i];
        forAll(nbrs, k)
        {
            label j = nbrs[k];
            const vector& Cj = globalCf[j];
            vector d = C - Cj;
            
            // Use a robust tolerance (1e-6) to avoid near-coincident faces causing FPE later
            if (d.x() > 1e-6) gUpwindX[i].append(j);
            else if (d.x() < -1e-6) gDownwindX[i].append(j);

            if (d.y() > 1e-6) gUpwindY[i].append(j);
            else if (d.y() < -1e-6) gDownwindY[i].append(j);
        }

        // Ring Analysis
        // Just collecting all point-connected faces
        const face& f = globalPatchFaces[i];
        labelHashSet ringSet;
        forAll(f, v)
        {
             label pGlobal = f[v];
             // whichPoint returns -1 if point not in patch support
             label pLocal = globalPatch.whichPoint(pGlobal);
             
             if (pLocal != -1)
             {
                 if (pLocal >= globalPointFaces.size()) 
                 {
                     FatalErrorInFunction << "Logic Error: pLocal out of bounds" << abort(FatalError);
                 }
                 
                 const labelList& pf = globalPointFaces[pLocal];
                 forAll(pf, m)
                 {
                     label pf_idx = pf[m];
                     if (pf_idx != i) ringSet.insert(pf_idx);
                 }
             }
        }
        gRing[i] = ringSet.toc();
    }
    Info << "DEBUG: CHECKPOINT F - Loop Complete." << endl;


    // ---------------------------------------------------------
    // 6. Map Local Faces to Global Faces
    // ---------------------------------------------------------
    Info << "DEBUG: CHECKPOINT G - Sorting Global X..." << endl;
    
    // SortableList for fast X-lookup
    scalarList globalX(globalSize);
    scalar minGlobalX = GREAT;
    forAll(globalCf, i) {
        globalX[i] = globalCf[i].x();
        if (globalX[i] < minGlobalX) minGlobalX = globalX[i];
    }
    Info << "Global Min X = " << minGlobalX << endl;
    
    SortableList<scalar> sortedGlobalX(globalX); 
    
    const label localSize = patch().size();
    labelList localToGlobal(localSize, -1);
    
    Info << "DEBUG: CHECKPOINT H - Mapping " << localSize << " local faces to global..." << endl;
    
    // Initialize data structures
    Info << "DEBUG: Resizing arrays..." << endl;
    upwindNodesX_.setSize(localSize);
    downwindNodesX_.setSize(localSize);
    upwindNodesY_.setSize(localSize);
    downwindNodesY_.setSize(localSize);
    lsMergedCandidates_.setSize(localSize);
    neighboursGlobal_.setSize(localSize);
    
    // Critical flags initialization for UPFDV2 safety
    ownerHasBodyFace_.setSize(localSize, false);
    ownerHasInletFace_.setSize(localSize, false);
    ownerHasFirstBodyFace_.setSize(localSize, false);
    
    schemeCodeX_.setSize(localSize, 1);
    schemeCodeY_.setSize(localSize, 1);
    interpNeededX_.setSize(localSize); 
    forAll(interpNeededX_, i) { interpNeededX_[i] = false; }
    interpNeededY_.setSize(localSize);
    forAll(interpNeededY_, i) { interpNeededY_[i] = false; }

    globalIdToCf_.clear();
    globalIdOfLocalFace_.setSize(localSize);
    globalIdToPackedIdx_.clear();
    
    remoteIdToUpwindNodesX_.clear();
    remoteIdToDownwindNodesX_.clear();
    remoteIdToUpwindNodesY_.clear();
    remoteIdToDownwindNodesY_.clear();
    remoteNeighbours_.clear();
    
    Info << "DEBUG: Arrays resized. Entering mapping loop..." << endl;

    // ---------------------------------------------------------
    // 6c. Parallel Topology Handshake (Part 1 - Local Mapping)
    // ---------------------------------------------------------
    
    scalar tol = 1e-4; // Tolerance
    
    for(label i=0; i<localSize; ++i)
    {
        if (i < 5 || i % 5000 == 0) Info << "DEBUG: Mapping local face " << i << endl;
        
        const vector& C = patch().Cf()[i];
        
        scalar minDist = GREAT;
        label match = -1;
        
        // Brute force search (optimized by X-check)
        // Note: globalCf is a reference to globalPatch.faceCentres()
        for(label k=0; k<globalSize; ++k)
        {
            if (mag(globalCf[k].x() - C.x()) < tol) // Pre-filter X
            {
                scalar d = mag(globalCf[k] - C);
                if (d < minDist)
                {
                    minDist = d;
                    match = k;
                }
            }
        }
        
        if (match == -1 || minDist > 1e-3)
        {
            FatalErrorInFunction 
                << "Could not map local face " << i << " at " << C 
                << " to any global face. MinDist=" << minDist 
                << abort(FatalError);
        }
        
        localToGlobal[i] = match;
    }

    Info << "DEBUG: Local mapping complete." << endl;

    // Determine GlobalID -> ProcID Ownership
    Map<label> globalToProc;
    if (Pstream::parRun())
    {
        List<labelList> allGlobalIDs(Pstream::nProcs());
        allGlobalIDs[Pstream::myProcNo()] = localToGlobal;
        Pstream::gatherList(allGlobalIDs);
        Pstream::scatterList(allGlobalIDs);

        forAll(allGlobalIDs, procI) 
        {
            const labelList& gIDs = allGlobalIDs[procI];
            forAll(gIDs, k) globalToProc.insert(gIDs[k], procI);
        }
    }
    else
    {
        forAll(localToGlobal, i) globalToProc.insert(localToGlobal[i], 0);
    }
    
    // Track what we need from others
    List<labelHashSet> myRequests(Pstream::nProcs());
    
    auto addRequest = [&](label gID)
    {
        if (gID < 0) return;
        if (globalIdToPackedIdx_.found(gID)) return; // Local
        
        if (globalToProc.found(gID))
        {
            label owner = globalToProc[gID];
            if (owner != Pstream::myProcNo())
            {
                myRequests[owner].insert(gID);
            }
            // Ensure we have geometry for calculation
            if (!globalIdToCf_.found(gID)) globalIdToCf_.insert(gID, globalCf[gID]);
        }
    };

    // ---------------------------------------------------------
    // 7. Transfer Global Stencils to Local Lists & Populating Maps
    // ---------------------------------------------------------
    
    // First Pass: Populate Local Lookup (Needed for isLocal check)
    for(label i=0; i<localSize; ++i)
    {
        label gID = localToGlobal[i];
        globalIdOfLocalFace_[i] = gID;
        globalIdToPackedIdx_.insert(gID, i);
    }

    // Second Pass: Stencils and Dependencies
    for(label i=0; i<localSize; ++i)
    {
        label gID = localToGlobal[i];
        
        // --- X Direction ---
        upwindNodesX_[i] = gUpwindX[gID];
        downwindNodesX_[i] = gDownwindX[gID];
        
        // Handle Remote Stencils for UPFDV2 support
        if (upwindNodesX_[i].size() > 0) 
        {
            label up1 = upwindNodesX_[i][0];
            addRequest(up1);
            if (!globalIdToPackedIdx_.found(up1)) 
            {
                // up1 is remote. Store up2List (stencil of up1) in remoteId map
                remoteIdToUpwindNodesX_.insert(gID, gUpwindX[up1]);
                if (gUpwindX[up1].size() > 0) addRequest(gUpwindX[up1][0]); 
            }
        }
        if (downwindNodesX_[i].size() > 0)
        {
            label dn1 = downwindNodesX_[i][0];
            addRequest(dn1);
            if (!globalIdToPackedIdx_.found(dn1))
            {
                remoteIdToDownwindNodesX_.insert(gID, gDownwindX[dn1]);
                if (gDownwindX[dn1].size() > 0) addRequest(gDownwindX[dn1][0]);
            }
        }

        // --- Y Direction ---
        upwindNodesY_[i] = gUpwindY[gID];
        downwindNodesY_[i] = gDownwindY[gID];
        
        if (upwindNodesY_[i].size() > 0) 
        {
            label up1 = upwindNodesY_[i][0];
            addRequest(up1);
            if (!globalIdToPackedIdx_.found(up1)) 
            {
                remoteIdToUpwindNodesY_.insert(gID, gUpwindY[up1]);
                if (gUpwindY[up1].size() > 0) addRequest(gUpwindY[up1][0]);
            }
        }
        if (downwindNodesY_[i].size() > 0)
        {
            label dn1 = downwindNodesY_[i][0];
            addRequest(dn1);
            if (!globalIdToPackedIdx_.found(dn1))
            {
                remoteIdToDownwindNodesY_.insert(gID, gDownwindY[dn1]);
                if (gDownwindY[dn1].size() > 0) addRequest(gDownwindY[dn1][0]);
            }
        }

        // --- Candidates & Remote Neighbours ---
        labelHashSet candidates;
        
        auto processCandidate = [&](label c) {
            candidates.insert(c);
            addRequest(c);
            
            // If candidate is remote, store ITS candidates for interpolation
            if (!globalIdToPackedIdx_.found(c))
            {
                if (!remoteNeighbours_.found(c))
                {
                    labelHashSet rc;
                    const labelList& r1 = gRing[c];
                    forAll(r1, k) { rc.insert(r1[k]); addRequest(r1[k]); }
                    const labelList& uX = gUpwindX[c];
                    forAll(uX, k) { rc.insert(uX[k]); addRequest(uX[k]); }
                    const labelList& dX = gDownwindX[c];
                    forAll(dX, k) { rc.insert(dX[k]); addRequest(dX[k]); }
                    
                    remoteNeighbours_.insert(c, rc.toc());
                }
            }
        };

        // Aggressive candidate collection (Ring + Upwind + Downwind)
        const labelList& ring1 = gRing[gID];
        forAll(ring1, k) processCandidate(ring1[k]);
        const labelList& upX = gUpwindX[gID];
        forAll(upX, k) processCandidate(upX[k]);
        const labelList& downX = gDownwindX[gID];
        forAll(downX, k) processCandidate(downX[k]);

        candidates.erase(gID);
        lsMergedCandidates_[i] = candidates.toc();
        neighboursGlobal_[i] = lsMergedCandidates_[i]; 
        
        // Ensure Geometry
        forAll(lsMergedCandidates_[i], k)
        {
            label candID = lsMergedCandidates_[i][k];
            if (!globalIdToCf_.found(candID)) globalIdToCf_.insert(candID, globalCf[candID]);
        }
        if (!globalIdToCf_.found(gID)) globalIdToCf_.insert(gID, globalCf[gID]);
        
        // Fallback Schemes if missing neighbors at boundary
        if (upwindNodesX_[i].empty()) schemeCodeX_[i] = 5; 
        if (downwindNodesX_[i].empty()) schemeCodeX_[i] = 4;
        if (upwindNodesX_[i].empty() && downwindNodesX_[i].empty()) schemeCodeX_[i] = 100;

        if (upwindNodesY_[i].empty()) schemeCodeY_[i] = 5;
        if (downwindNodesY_[i].empty()) schemeCodeY_[i] = 4;
        if (upwindNodesY_[i].empty() && downwindNodesY_[i].empty()) schemeCodeY_[i] = 100;
    }
    
    // ---------------------------------------------------------
    // 7b. Safety Check for FPE (Floating Point Exceptions)
    // ---------------------------------------------------------
    // Validate stencils to ensure UPFDV2 won't crash with DivByZero
    for(label i=0; i<localSize; ++i)
    {
        label gID = localToGlobal[i];
        const vector& C = globalCf[gID]; // local Cf might slightly differ? Prefer global for consistency
        
        auto getDx = [&](label nID) -> scalar {
             return mag(C.x() - globalCf[nID].x());
        };
        auto getDy = [&](label nID) -> scalar {
             return mag(C.y() - globalCf[nID].y());
        };

        // Check X Scheme
        label sX = schemeCodeX_[i];
        if (sX == 1) // Central
        {
            if (upwindNodesX_[i].empty() || downwindNodesX_[i].empty()) {
                FatalErrorInFunction << "Face " << i << " assigned Scheme 1 but missing neighbors!" << abort(FatalError);
            }
            scalar h1 = getDx(upwindNodesX_[i][0]);
            scalar h2 = getDx(downwindNodesX_[i][0]);
            if (h1 < VSMALL || h2 < VSMALL) {
                // Determine if we should fallback or abort
                Info << "WARNING: Face " << i << " (Global " << gID << ") has zero X-node distance. h1=" << h1 << " h2=" << h2 << ". Falling back to 100." << endl;
                schemeCodeX_[i] = 100;
            }
        }
        else if (sX == 4) // Backward
        {
             if (upwindNodesX_[i].empty()) {
                 FatalErrorInFunction << "Face " << i << " assigned Scheme 4 but missing upwind!" << abort(FatalError);
             }
             scalar h1 = getDx(upwindNodesX_[i][0]);
             if (h1 < VSMALL) {
                  Info << "WARNING: Face " << i << " has zero Upwind X distance. Fallback to 100." << endl;
                  schemeCodeX_[i] = 100;
             }
        }
        else if (sX == 5) // Forward
        {
             if (downwindNodesX_[i].empty()) {
                  FatalErrorInFunction << "Face " << i << " assigned Scheme 5 but missing downwind!" << abort(FatalError);
             }
             scalar h1 = getDx(downwindNodesX_[i][0]);
             if (h1 < VSMALL) {
                  Info << "WARNING: Face " << i << " has zero Downwind X distance. Fallback to 100." << endl;
                  schemeCodeX_[i] = 100;
             }
        }

        // Check Y Scheme
        label sY = schemeCodeY_[i];
        if (sY == 1) // Central
        {
            if (upwindNodesY_[i].empty() || downwindNodesY_[i].empty()) {
                FatalErrorInFunction << "Face " << i << " assigned Scheme 1 but missing neighbors!" << abort(FatalError);
            }
            scalar h1 = getDy(upwindNodesY_[i][0]);
            scalar h2 = getDy(downwindNodesY_[i][0]);
            if (h1 < VSMALL || h2 < VSMALL) {
                Info << "WARNING: Face " << i << " has zero Y-node distance. h1=" << h1 << " h2=" << h2 << ". Falling back to 100." << endl;
                schemeCodeY_[i] = 100;
            }
        }
    }

    // ---------------------------------------------------------
    // 8. Finalize Communication Lists (neededByProc_)
    // ---------------------------------------------------------
    if (Pstream::parRun())
    {
        // Exchange requests to build neededByProc_
        // 1. Prepare data
        List<List<label>> requestsToSend(Pstream::nProcs());
        labelList myTargets; 
        {
            DynamicList<label> t(Pstream::nProcs());
            for(label p=0; p<Pstream::nProcs(); ++p) {
                requestsToSend[p] = myRequests[p].toc();
                if (requestsToSend[p].size() > 0) t.append(p);
            }
            myTargets.transfer(t);
        }
        
        // 2. Send Data
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
        forAll(myTargets, i)
        {
            label p = myTargets[i];
            UOPstream os(p, pBufs);
            os << requestsToSend[p];
        }
        pBufs.finishedSends();
        
        // 3. Topology Exchange (to know who sent to us)
        List<labelList> allTargets(Pstream::nProcs());
        allTargets[Pstream::myProcNo()] = myTargets;
        Pstream::gatherList(allTargets);
        Pstream::scatterList(allTargets);

        // 4. Receive Data
        neededByProc_.setSize(Pstream::nProcs());
        sendTo.clear();
        
        for(label src=0; src<Pstream::nProcs(); ++src)
        {
             if (src == Pstream::myProcNo()) continue;
             
             // Did src send to me?
             if (allTargets[src].found(Pstream::myProcNo()))
             {
                 UIPstream is(src, pBufs);
                 List<label> reqs;
                 is >> reqs;
                 neededByProc_[src].insert(reqs);
                 sendTo.insert(src);
             }
        }
        
        Info << "V18D: Parallel handshake complete. Neighbors: " << sendTo.size() << endl;
    }

    // Set completion flags
    neigboursCalculated_ = true;
    candidatesCalculated_ = false; 
    schemesDetected_ = false; 
    
    Info << "Global Stencil Calculation Complete. " << endl;
    
    findInletFaces();
    findSphereFaces();
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
            label cellFaceI = ownerCellFaces[fJ];
            label patchID = mesh.boundaryMesh().whichPatch(cellFaceI);

            if (patchID != -1 && mesh.boundaryMesh()[patchID].name() == "inlet")
            {
                ownerHasInletFace_[i] = true;
                break;
            }
        }
    }
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
            const face& freeSurfFace = mesh.faces()[startFace + i];  
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

                    if (sharedCount >= 2)
                    {
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
        Info << "WARNING: No faces near 'sphere' found globally." << endl;
    }
    
    buildVertexBasedConnections();
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
            vertexToFaces[vertexID].append(i); 
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

                if (!firstRingConnections_[fa].found(fb)) firstRingConnections_[fa].append(fb);
				if (!firstRingConnections_[fb].found(fa)) firstRingConnections_[fb].append(fa);
			}
        }
    }
	const vectorField Cf = patch().Cf();  
     // --- Step 3: Populate ringNodeX_ ---
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
					ringNodeX_[i] = globalIdOfLocalFace_[conn[j]]; 
			}
		}     
    }
	
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
		if (ownerHasFirstBodyFace_[i]) continue;  
		if (!ownerHasBodyFace_[i]) continue;  
		
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
		}
	}
	
	 for (label i = 0; i < nFaces; ++i)
    {
        if (ownerHasFirstBodyFace_[i]) continue;
		if (!ownerHasBodyFace_[i]) continue;  
		
		if (!secondRingToFirstRing_.found(i)) continue;
        const List<label>& conn = secondRingToFirstRing_[i];
		if (conn.empty()) continue;

		scalar minDx = GREAT;
		label bestConn = -1;

		forAll(conn, j)
		{
			label frIdx = conn[j]; // first-ring local index
			scalar dx = Cf[i].x() - Cf[frIdx].x();

			// look for upstream connection
			if (dx > VSMALL && dx < minDx)
			{
				minDx = dx;
				bestConn = frIdx;
			}
		}

		if (bestConn != -1)
		{
			ring2NodeX_.insert(i, globalIdOfLocalFace_[bestConn]);
		} 
        else
		{
			forAll(conn, j)
			{
				label frIdx = conn[j]; 
				scalar dx = Cf[frIdx].x()-Cf[i].x(); // check forward ?

				if (dx < minDx)
				{
					minDx = dx;
					bestConn = frIdx;
				}
			}
			if (bestConn != -1) ring2NodeX_.insert(i, globalIdOfLocalFace_[bestConn]);
		} 
    }
	
	scalar localMaxX = -GREAT;
	// --- compute local maximum only for valid faces ---
	forAll(Cf, i)
	{
		if (ownerHasFirstBodyFace_[i])
		{
			localMaxX = Foam::max(localMaxX, Cf[i].x());
		}
	}
	scalar globalMaxX = localMaxX;
	Foam::reduce(globalMaxX, Foam::maxOp<scalar>());

    const scalar limitingX = 0.05;     
	 forAll(Cf, i)
    {
		if (ownerHasBodyFace_[i]) 
		{	scalar x00 = Cf[i].x(); 
			if (x00 > (globalMaxX-limitingX)) 
			{
				// exclude tail?
				ownerHasBodyFace_[i]= false;
			}
		}	
    }
}
void Foam::waveCurrentPotential3DFvPatchScalarField::initCandidates()
{
    if (candidatesCalculated_) return;

    const label nFaces = patch().size();
	lsMergedCandidates_.setSize(nFaces);
	const vectorField& Cf = patch().Cf();
    const scalarField xComponents = Cf.component(vector::X);
    const scalarField yComponents = Cf.component(vector::Y);
	
	auto getDistBetweenGIDs = [&](label gid1, label gid2) -> scalar
	{
		vector C1(Zero), C2(Zero);
		if (globalIdToPackedIdx_.found(gid1)) C1 = Cf[globalIdToPackedIdx_[gid1]];
		else if (globalIdToCf_.found(gid1)) C1 = globalIdToCf_[gid1];
        else return GREAT;

		if (globalIdToPackedIdx_.found(gid2)) C2 = Cf[globalIdToPackedIdx_[gid2]];
		else if (globalIdToCf_.found(gid2)) C2 = globalIdToCf_[gid2];
        else return GREAT;

		return mag(C2 - C1);
	};

    auto mergeUniquewoSelf = [](const List<label>& a, const List<label>& b,const label selfGid)
    {
        DynamicList<label> merged;
        forAll(a, i) if (a[i] != selfGid && !merged.found(a[i])) merged.append(a[i]);
        forAll(b, j) if (b[j] != selfGid && !merged.found(b[j])) merged.append(b[j]);
        return List<label>(merged);
    };
	
    auto getCandidates = [&](label selfGid, label gid, scalar Uy) -> const labelList&
	{
		static labelList sendt;
		sendt.clear();

		// Scheme 1
		if (mag(Uy) <= VSMALL)
		{
			if (globalIdToPackedIdx_.found(gid))
			{
				const labelList& local = neighboursGlobal_[globalIdToPackedIdx_[gid]];
				sendt.append(local);
			}
			if (remoteNeighbours_.found(gid))
			{
				const labelList& rem1 = remoteNeighbours_[gid];
				sendt.append(rem1);
			}
            // remoteNeighbours2_ handling omitted for V18D
			return sendt;
		}

		// Scheme 6 – Upwind Y
		if (Uy > VSMALL)
		{
			if (globalIdToPackedIdx_.found(gid))
			{
				const labelList& local = upwindNodesY_[globalIdToPackedIdx_[gid]];
				sendt.append(local);

				if (!local.empty())
				{
					label nearest = local[0];
					scalar minDist = GREAT;
					forAll(local, j) {
						scalar d = getDistBetweenGIDs(gid, local[j]);
						if (d < minDist) { minDist = d; nearest = local[j]; }
					}

					if (globalIdToPackedIdx_.found(nearest)) {
						const labelList& local2 = upwindNodesY_[globalIdToPackedIdx_[nearest]];
						sendt.append(local2);
					}
					else if (remoteIdToUpwindNodesY_.found(nearest)) {
						const labelList& rem = remoteIdToUpwindNodesY_[nearest];
						sendt.append(rem);
					}
				}
				return sendt;
			}
			if (remoteIdToUpwindNodesY_.found(gid))
			{
				const labelList& rem1 = remoteIdToUpwindNodesY_[gid];
				sendt.append(rem1);

				if (!rem1.empty())
				{
					const label first = rem1[0];
					if (remoteIdToUpwindNodesY_.found(first)) {
						const labelList& rem2 = remoteIdToUpwindNodesY_[first];
						sendt.append(rem2);
					}
					else if (globalIdToPackedIdx_.found(first)) {
						const labelList& local2 = upwindNodesY_[globalIdToPackedIdx_[first]];
						sendt.append(local2);
					}
				}
				return sendt;
			}
		}

		// Scheme 7 – Downwind Y
		if (Uy < -VSMALL)
		{
			if (globalIdToPackedIdx_.found(gid))
			{
				const labelList& local = downwindNodesY_[globalIdToPackedIdx_[gid]];
				sendt.append(local);

				if (!local.empty())
				{
					label nearest = local[0];
					scalar minDist = GREAT;
					forAll(local, j) {
						scalar d = getDistBetweenGIDs(gid, local[j]);
						if (d < minDist) { minDist = d; nearest = local[j]; }
					}

					if (globalIdToPackedIdx_.found(nearest)) {
						const labelList& local2 = downwindNodesY_[globalIdToPackedIdx_[nearest]];
						sendt.append(local2);
					}
					else if (remoteIdToDownwindNodesY_.found(nearest)) {
						const labelList& rem = remoteIdToDownwindNodesY_[nearest];
						sendt.append(rem);
					}
				}
				return sendt;
			}
			if (remoteIdToDownwindNodesY_.found(gid))
			{
				const labelList& rem1 = remoteIdToDownwindNodesY_[gid];
				sendt.append(rem1);

				if (!rem1.empty())
				{
					const label first = rem1[0];
					if (remoteIdToDownwindNodesY_.found(first)) {
						const labelList& rem2 = remoteIdToDownwindNodesY_[first];
						sendt.append(rem2);
					}
					else if (globalIdToPackedIdx_.found(first)) {
						const labelList& local2 = downwindNodesY_[globalIdToPackedIdx_[first]];
						sendt.append(local2);
					}
				}
				return sendt;
			}
		}
		
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
			List<label> mergedCands;
			const label selfGid = globalIdOfLocalFace_[i];
			
			if (ownerHasFirstBodyFace_[i])
			{	
				if (ringNodeX_.found(i))  
				{	
					const List<label>& candsSelf = getCandidates(selfGid, selfGid, uY);
					mergedCands = mergeUniquewoSelf(mergedCands, candsSelf, selfGid);
                    label conn = ringNodeX_[i];
					mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn), selfGid);
					
					if (globalIdToPackedIdx_.found(conn))
					{	
						const labelList& ring = neighboursGlobal_[globalIdToPackedIdx_[conn]];
						mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn), selfGid);
						mergedCands = mergeUniquewoSelf(mergedCands, ring , selfGid);
						lsMergedCandidates_[i]=mergedCands;
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
						mergedCands = mergeUniquewoSelf(mergedCands, ring , selfGid);
						lsMergedCandidates_[i]=mergedCands;
					
						if (ringNodeX_.found(globalIdToPackedIdx_[conn2]))
						{
							label conn3 = ringNodeX_[globalIdToPackedIdx_[conn2]];
							mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn3), selfGid);
							const labelList& ring2 = neighboursGlobal_[globalIdToPackedIdx_[conn3]];
							mergedCands = mergeUniquewoSelf(mergedCands, ring2 , selfGid);
							lsMergedCandidates_[i]=mergedCands;
							
							if (ringNodeX_.found(globalIdToPackedIdx_[conn3]))
							{
								label conn4 = ringNodeX_[globalIdToPackedIdx_[conn3]];
								mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn4), selfGid);
								lsMergedCandidates_[i]=mergedCands;
							}
						}
					}
					const List<label>& L = upwindNodesX_[i];
					if (!L.empty())
					{
						label up1 = L[0];
						const List<label>& candsSelf2 = getCandidates(selfGid,selfGid,uY); 
						mergedCands = mergeUniquewoSelf(mergedCands,candsSelf2,selfGid);
						mergedCands = mergeUniquewoSelf(mergedCands,L,selfGid);
						const List<label>& candsX1 = getCandidates(selfGid,up1,uY); 
						mergedCands = mergeUniquewoSelf(mergedCands,candsX1,selfGid);
						
						List<label> up2List;
						if (globalIdToPackedIdx_.found(up1)) up2List = upwindNodesX_[globalIdToPackedIdx_[up1]];
						else {
							const auto it = remoteIdToUpwindNodesX_.find(selfGid);
							if (it != remoteIdToUpwindNodesX_.end()) up2List = it();
						}
						mergedCands = mergeUniquewoSelf(mergedCands,up2List,selfGid);
						lsMergedCandidates_[i]=mergedCands;
					}	
					lsMergedCandidates_[i]=mergedCands;
				}
			}
			
			if (ring2NodeX_.found(i))
			{
				label conn = ring2NodeX_[i];
				mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn), selfGid);
				const List<label>& candsSelf = getCandidates(selfGid, selfGid, uY);
				mergedCands = mergeUniquewoSelf(mergedCands, candsSelf, selfGid);
				
				if (globalIdToPackedIdx_.found(conn))
				{	
					mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn), selfGid);
					const List<label>& candsX1 = getCandidates(selfGid,conn,uY); 
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
						if (ringNodeX_.found(globalIdToPackedIdx_[conn3]))
						{
							label conn4 = ringNodeX_[globalIdToPackedIdx_[conn3]];
							mergedCands = mergeUniquewoSelf(mergedCands, List<label>(1, conn4), selfGid);
							lsMergedCandidates_[i]=mergedCands;
						}
					}
				}
				const List<label>& L = upwindNodesX_[i];
				if (!L.empty())	
				{	
                    label up1 = L[0];
					const List<label>& candsSelf2 = getCandidates(selfGid,selfGid,uY); 
					mergedCands = mergeUniquewoSelf(mergedCands,L,selfGid);
					mergedCands = mergeUniquewoSelf(mergedCands,candsSelf2,selfGid);
					const List<label>& candsX1 = getCandidates(selfGid,up1,uY); 
					mergedCands = mergeUniquewoSelf(mergedCands,candsX1,selfGid);
				
					List<label> up2List;
					if (globalIdToPackedIdx_.found(up1)) up2List = upwindNodesX_[globalIdToPackedIdx_[up1]];
					else {
						const auto it = remoteIdToUpwindNodesX_.find(selfGid);
						if (it != remoteIdToUpwindNodesX_.end()) up2List = it();
					}
					mergedCands = mergeUniquewoSelf(mergedCands,up2List,selfGid);
					lsMergedCandidates_[i]=mergedCands;
				}
			}
			
			if (mergedCands.empty())
			{
				const List<label>& L = upwindNodesX_[i];
				if (L.empty()) continue; 

				label up1 = L[0];
				List<label> up2List;

				if (globalIdToPackedIdx_.found(up1)) up2List = upwindNodesX_[globalIdToPackedIdx_[up1]];
				else {
					const auto it = remoteIdToUpwindNodesX_.find(selfGid);
					if (it != remoteIdToUpwindNodesX_.end()) up2List = it();
					else up2List.clear(); 
				}

				if (up2List.empty()) {
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
			lsMergedCandidates_[i] = mergedCands;
		}
	}	
	candidatesCalculated_ = true;
    Info << "Candidates for Mq built on " << patch().name() << endl;
}
void Foam::waveCurrentPotential3DFvPatchScalarField::findUpwindDownwindNodesV2() {}
void Foam::waveCurrentPotential3DFvPatchScalarField::detectFDSchemes()
{
    if (schemesDetected_) return;

    const vectorField& Cf = patch().Cf();
    const scalarField x = Cf.component(vector::X);
    const scalarField y = Cf.component(vector::Y);
    const label nFaces = patch().size();

    schemeCodeX_.setSize(nFaces, -1);
    schemeCodeY_.setSize(nFaces, -1);
    
    // interpNeededX_ sized in calcNeigboursV3

	const volVectorField& Ucur = db().lookupObjectRef<volVectorField>("Ucur");
	const vectorField& Ucur_p = Ucur.boundaryField()[patch().index()];
	
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
						//WarningInFunction << "Face " << globalIdOfLocalFace_[i] << ": First Order Detected !!!!" << endl;
					}
				}

				const List<label>& R = downwindNodesY_[i];
				if (R.size())
				{
					auto dxdyR = dxdy(Cf, globalIdToPackedIdx_, globalIdToCf_, R[0], x[i], y[i]);
					if (aligned(dxdyR.first, dxdyR.second, "y")){
						schemeCodeY_[i] = 5;
						//WarningInFunction << "Face " << globalIdOfLocalFace_[i] << ": First Order Detected !!!" << endl;
					}
				}
			}
		}
	}

	Info << "Scheme detection complete." << endl;
	
	initCandidates();
	schemesDetected_ = true;
}
