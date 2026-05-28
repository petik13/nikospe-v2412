/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
PreParV19D — SVD interpolation weight cache builder (Change 1)
Precomputes MQ RBF weight vectors once (mesh is fixed) to avoid repeated
SVDinv() calls inside UPFDV2().
\*---------------------------------------------------------------------------*/

void Foam::waveCurrentPotential3DFvPatchScalarField::buildInterpCaches()
{
    if (interpCachesBuilt_) return;

    const vectorField& Cf = patch().Cf();
    const label nFaces = patch().size();

    interp2ndCacheX_.setSize(nFaces);
    interp2ndCacheY_.setSize(nFaces);
    interp4thCacheX_.setSize(nFaces);
    interp4thCacheY_.setSize(nFaces);

    // ---- helper: get face-centre vector for a global ID ----
    auto getCf = [&](label gid) -> const vector&
    {
        const auto it = globalIdToPackedIdx_.find(gid);
        if (it != globalIdToPackedIdx_.end()) return Cf[*it];
        const auto it2 = globalIdToCf_.find(gid);
        if (it2 != globalIdToCf_.end()) return *it2;
        FatalErrorInFunction << "buildInterpCaches: gid " << gid
            << " not found in Cf maps" << abort(FatalError);
        static const vector dummy(Zero); return dummy;
    };

    // ---- helper: get candidates for a global ID ----
    auto getCandidates = [&](label gid) -> const labelList&
    {
        {
            const auto it = globalIdToPackedIdx_.find(gid);
            if (it != globalIdToPackedIdx_.end())
                return neighboursGlobal_[*it];
        }
        if (remoteNeighbours_.found(gid))
            return remoteNeighbours_[gid];
        if (remoteNeighbours2_.found(gid))
            return remoteNeighbours2_[gid];
        static const labelList dummy;
        return dummy;
    };

    // ---- helper: build one cache entry ----
    // selfGid: the node being interpolated
    // candidates: neighbour GIDs (getCandidates result)
    // x_interp, y_interp: target interpolation point
    // epsilon: MQ shape parameter
    auto buildCache = [&](
        label selfGid,
        const labelList& candidates,
        scalar x_interp,
        scalar y_interp,
        scalar epsilon,
        InterpWeightCache& cache)
    {
        if (candidates.empty()) return; // leave valid=false

        // Collect points: [selfGid, cand0, cand1, ...]
        DynamicList<label> gids;
        DynamicList<vector> points;

        auto addPoint = [&](label gid)
        {
            const vector& c = getCf(gid);
            gids.append(gid);
            points.append(vector(c.x(), c.y(), 0.0));
        };

        addPoint(selfGid);
        forAll(candidates, j) addPoint(candidates[j]);

        const label N = gids.size();
        if (N < 2) return; // not enough for RBF interpolation

        // Build MQ RBF matrix A
        scalarSquareMatrix A(N, 0.0);
        for (label ii = 0; ii < N; ++ii)
            for (label jj = 0; jj < N; ++jj)
            {
                scalar r = mag(points[ii] - points[jj]);
                A(ii, jj) = sqrt(r*r + epsilon*epsilon);
            }

        // Invert A via SVD
        scalarRectangularMatrix Ainv = SVDinv(A, 1e-12);

        // Evaluate RBF at interpolation point
        vector pInterp(x_interp, y_interp, 0.0);
        scalarField phiVec(N, 0.0);
        for (label jj = 0; jj < N; ++jj)
        {
            scalar r = mag(pInterp - points[jj]);
            phiVec[jj] = sqrt(r*r + epsilon*epsilon);
        }

        // w[j] = sum_k phiVec[k] * Ainv(k,j)   (row of phi^T * Ainv)
        cache.gids = List<label>(gids);
        cache.w.setSize(N, 0.0);
        for (label jj = 0; jj < N; ++jj)
            for (label kk = 0; kk < N; ++kk)
                cache.w[jj] += phiVec[kk] * Ainv(kk, jj);

        cache.valid = true;
    };

    // ---- helper: get the second upwind/downwind node list for a face ----
    auto getUp2ListX = [&](label localFaceI, label up1) -> List<label>
    {
        const auto it = globalIdToPackedIdx_.find(up1);
        if (it != globalIdToPackedIdx_.end())
            return upwindNodesX_[*it];
        const auto it2 = remoteIdToUpwindNodesX_.find(globalIdOfLocalFace_[localFaceI]);
        if (it2 != remoteIdToUpwindNodesX_.end())
            return it2();
        return List<label>();
    };

    auto getDown2ListX = [&](label localFaceI, label down1) -> List<label>
    {
        const auto it = globalIdToPackedIdx_.find(down1);
        if (it != globalIdToPackedIdx_.end())
            return downwindNodesX_[*it];
        const auto it2 = remoteIdToDownwindNodesX_.find(globalIdOfLocalFace_[localFaceI]);
        if (it2 != remoteIdToDownwindNodesX_.end())
            return it2();
        return List<label>();
    };

    auto getUp2ListY = [&](label localFaceI, label up1) -> List<label>
    {
        const auto it = globalIdToPackedIdx_.find(up1);
        if (it != globalIdToPackedIdx_.end())
            return upwindNodesY_[*it];
        const auto it2 = remoteIdToUpwindNodesY_.find(globalIdOfLocalFace_[localFaceI]);
        if (it2 != remoteIdToUpwindNodesY_.end())
            return it2();
        return List<label>();
    };

    auto getDown2ListY = [&](label localFaceI, label down1) -> List<label>
    {
        const auto it = globalIdToPackedIdx_.find(down1);
        if (it != globalIdToPackedIdx_.end())
            return downwindNodesY_[*it];
        const auto it2 = remoteIdToDownwindNodesY_.find(globalIdOfLocalFace_[localFaceI]);
        if (it2 != remoteIdToDownwindNodesY_.end())
            return it2();
        return List<label>();
    };

    const scalarField xC = Cf.component(vector::X);
    const scalarField yC = Cf.component(vector::Y);

    derivCoeffsX_.setSize(nFaces);
    derivCoeffsY_.setSize(nFaces);

    for (label i = 0; i < nFaces; ++i)
    {
        if (ownerHasBodyFace_[i]) continue; // body faces handled by GreenGauss

        const scalar epsilon = shape_ / patch().deltaCoeffs()[i];
        const scalar y0 = yC[i];
        const scalar x0 = xC[i];

        // ============================================================
        // X-direction
        // ============================================================
        const label sX = label(schemeCodeX_[i]);
        const List<label>& LX = upwindNodesX_[i];
        const List<label>& RX = downwindNodesX_[i];

        if (sX == 6 && !LX.empty())
        {
            label up1 = LX[0];

            // slot 0: interp at (up1.x, y0)
            if (interpNeededX_[i][0])
            {
                const scalar x_up1 = getCf(up1).x();
                const labelList& cands = getCandidates(up1);
                buildCache(up1, cands, x_up1, y0, epsilon, interp2ndCacheX_[i][0]);
            }

            // slot 1: interp at (up2.x, y0)
            if (interpNeededX_[i][1])
            {
                List<label> up2List = getUp2ListX(i, up1);
                if (!up2List.empty())
                {
                    label up2 = up2List[0];
                    const scalar x_up2 = getCf(up2).x();
                    const labelList& cands = getCandidates(up2);
                    buildCache(up2, cands, x_up2, y0, epsilon, interp2ndCacheX_[i][1]);
                }
            }
        }
        else if (sX == 7 && !RX.empty())
        {
            label down1 = RX[0];

            // slot 0: interp at (down1.x, y0)
            if (interpNeededX_[i][0])
            {
                const scalar x_dn1 = getCf(down1).x();
                const labelList& cands = getCandidates(down1);
                buildCache(down1, cands, x_dn1, y0, epsilon, interp2ndCacheX_[i][0]);
            }

            // slot 1: interp at (down2.x, y0)
            if (interpNeededX_[i][1])
            {
                List<label> dn2List = getDown2ListX(i, down1);
                if (!dn2List.empty())
                {
                    label down2 = dn2List[0];
                    const scalar x_dn2 = getCf(down2).x();
                    const labelList& cands = getCandidates(down2);
                    buildCache(down2, cands, x_dn2, y0, epsilon, interp2ndCacheX_[i][1]);
                }
            }
        }
        else if (sX == 1 && !LX.empty() && !RX.empty())
        {
            // slot 0: upwind node L[0]
            if (interpNeededX_[i][0])
            {
                label up1 = LX[0];
                const scalar x_up1 = getCf(up1).x();
                const labelList& cands = getCandidates(up1);
                buildCache(up1, cands, x_up1, y0, epsilon, interp2ndCacheX_[i][0]);
            }

            // slot 1: downwind node R[0]
            if (interpNeededX_[i][1])
            {
                label dn1 = RX[0];
                const scalar x_dn1 = getCf(dn1).x();
                const labelList& cands = getCandidates(dn1);
                buildCache(dn1, cands, x_dn1, y0, epsilon, interp2ndCacheX_[i][1]);
            }
        }
        else if (sX == 14 && !LX.empty() && !RX.empty())
        {
            label up1   = LX[0];
            label down1 = RX[0];
            List<label> up2List   = getUp2ListX(i, up1);
            List<label> down2List = getDown2ListX(i, down1);

            // slot 0: up1
            if (interpNeeded4thX_[i][0])
            {
                const scalar x_up1 = getCf(up1).x();
                const labelList& cands = getCandidates(up1);
                buildCache(up1, cands, x_up1, y0, epsilon, interp4thCacheX_[i][0]);
            }

            // slot 1: up2
            if (interpNeeded4thX_[i][1] && !up2List.empty())
            {
                label up2 = up2List[0];
                const scalar x_up2 = getCf(up2).x();
                const labelList& cands = getCandidates(up2);
                buildCache(up2, cands, x_up2, y0, epsilon, interp4thCacheX_[i][1]);
            }

            // slot 2: down1
            if (interpNeeded4thX_[i][2])
            {
                const scalar x_dn1 = getCf(down1).x();
                const labelList& cands = getCandidates(down1);
                buildCache(down1, cands, x_dn1, y0, epsilon, interp4thCacheX_[i][2]);
            }

            // slot 3: down2
            if (interpNeeded4thX_[i][3] && !down2List.empty())
            {
                label down2 = down2List[0];
                const scalar x_dn2 = getCf(down2).x();
                const labelList& cands = getCandidates(down2);
                buildCache(down2, cands, x_dn2, y0, epsilon, interp4thCacheX_[i][3]);
            }
        }

        // ============================================================
        // Y-direction
        // ============================================================
        const label sY = label(schemeCodeY_[i]);
        const List<label>& LY = upwindNodesY_[i];
        const List<label>& RY = downwindNodesY_[i];

        if (sY == 6 && !LY.empty())
        {
            label up1 = LY[0];

            // slot 0: interp at (x0, up1.y)
            if (interpNeededY_[i][0])
            {
                const scalar y_up1 = getCf(up1).y();
                const labelList& cands = getCandidates(up1);
                buildCache(up1, cands, x0, y_up1, epsilon, interp2ndCacheY_[i][0]);
            }

            // slot 1: interp at (x0, up2.y)
            if (interpNeededY_[i][1])
            {
                List<label> up2List = getUp2ListY(i, up1);
                if (!up2List.empty())
                {
                    label up2 = up2List[0];
                    const scalar y_up2 = getCf(up2).y();
                    const labelList& cands = getCandidates(up2);
                    buildCache(up2, cands, x0, y_up2, epsilon, interp2ndCacheY_[i][1]);
                }
            }
        }
        else if (sY == 7 && !RY.empty())
        {
            label down1 = RY[0];

            // slot 0: interp at (x0, down1.y)
            if (interpNeededY_[i][0])
            {
                const scalar y_dn1 = getCf(down1).y();
                const labelList& cands = getCandidates(down1);
                buildCache(down1, cands, x0, y_dn1, epsilon, interp2ndCacheY_[i][0]);
            }

            // slot 1: interp at (x0, down2.y)
            if (interpNeededY_[i][1])
            {
                List<label> dn2List = getDown2ListY(i, down1);
                if (!dn2List.empty())
                {
                    label down2 = dn2List[0];
                    const scalar y_dn2 = getCf(down2).y();
                    const labelList& cands = getCandidates(down2);
                    buildCache(down2, cands, x0, y_dn2, epsilon, interp2ndCacheY_[i][1]);
                }
            }
        }
        else if (sY == 1 && !LY.empty() && !RY.empty())
        {
            // slot 0: up1
            if (interpNeededY_[i][0])
            {
                label up1 = LY[0];
                const scalar y_up1 = getCf(up1).y();
                const labelList& cands = getCandidates(up1);
                buildCache(up1, cands, x0, y_up1, epsilon, interp2ndCacheY_[i][0]);
            }

            // slot 1: dn1
            if (interpNeededY_[i][1])
            {
                label dn1 = RY[0];
                const scalar y_dn1 = getCf(dn1).y();
                const labelList& cands = getCandidates(dn1);
                buildCache(dn1, cands, x0, y_dn1, epsilon, interp2ndCacheY_[i][1]);
            }
        }
        else if (sY == 14 && !LY.empty() && !RY.empty())
        {
            label up1   = LY[0];
            label down1 = RY[0];
            List<label> up2List   = getUp2ListY(i, up1);
            List<label> down2List = getDown2ListY(i, down1);

            // slot 0: up1
            if (interpNeeded4thY_[i][0])
            {
                const scalar y_up1 = getCf(up1).y();
                const labelList& cands = getCandidates(up1);
                buildCache(up1, cands, x0, y_up1, epsilon, interp4thCacheY_[i][0]);
            }

            // slot 1: up2
            if (interpNeeded4thY_[i][1] && !up2List.empty())
            {
                label up2 = up2List[0];
                const scalar y_up2 = getCf(up2).y();
                const labelList& cands = getCandidates(up2);
                buildCache(up2, cands, x0, y_up2, epsilon, interp4thCacheY_[i][1]);
            }

            // slot 2: down1
            if (interpNeeded4thY_[i][2])
            {
                const scalar y_dn1 = getCf(down1).y();
                const labelList& cands = getCandidates(down1);
                buildCache(down1, cands, x0, y_dn1, epsilon, interp4thCacheY_[i][2]);
            }

            // slot 3: down2
            if (interpNeeded4thY_[i][3] && !down2List.empty())
            {
                label down2 = down2List[0];
                const scalar y_dn2 = getCf(down2).y();
                const labelList& cands = getCandidates(down2);
                buildCache(down2, cands, x0, y_dn2, epsilon, interp4thCacheY_[i][3]);
            }
        }

        // ---- Change A: cache 5pt coefficients for s==14 faces ----
        if (sX == 14 && !LX.empty() && !RX.empty())
        {
            label up1 = LX[0], dn1 = RX[0];
            List<label> up2L = getUp2ListX(i, up1);
            List<label> dn2L = getDown2ListX(i, dn1);
            if (!up2L.empty() && !dn2L.empty())
            {
                const scalar h1 = mag(getCf(up1).x()     - x0);
                const scalar h2 = mag(getCf(up2L[0]).x() - x0);
                const scalar h3 = mag(getCf(dn1).x()     - x0);
                const scalar h4 = mag(getCf(dn2L[0]).x() - x0);
                scalar am2, am1, a0c, ap1, ap2;
                computeFirstDerivativeCoeffs5pt(-h2, -h1, 0.0, h3, h4, am2, am1, a0c, ap1, ap2);
                derivCoeffsX_[i][0] = am2; derivCoeffsX_[i][1] = am1;
                derivCoeffsX_[i][2] = a0c; derivCoeffsX_[i][3] = ap1; derivCoeffsX_[i][4] = ap2;
            }
        }

        if (sY == 14 && !LY.empty() && !RY.empty())
        {
            label up1 = LY[0], dn1 = RY[0];
            List<label> up2L = getUp2ListY(i, up1);
            List<label> dn2L = getDown2ListY(i, dn1);
            if (!up2L.empty() && !dn2L.empty())
            {
                const scalar h1 = mag(getCf(up1).y()     - y0);
                const scalar h2 = mag(getCf(up2L[0]).y() - y0);
                const scalar h3 = mag(getCf(dn1).y()     - y0);
                const scalar h4 = mag(getCf(dn2L[0]).y() - y0);
                scalar am2, am1, a0c, ap1, ap2;
                computeFirstDerivativeCoeffs5pt(-h2, -h1, 0.0, h3, h4, am2, am1, a0c, ap1, ap2);
                derivCoeffsY_[i][0] = am2; derivCoeffsY_[i][1] = am1;
                derivCoeffsY_[i][2] = a0c; derivCoeffsY_[i][3] = ap1; derivCoeffsY_[i][4] = ap2;
            }
        }
    }

    // ---- Changes B+C: cache comm topology and pre-build send index lists ----
    if (Pstream::parRun() && !commTopologyBuilt_)
    {
        List<labelList> sendToLists(Pstream::nProcs());
        sendToLists[Pstream::myProcNo()] = sendTo.toc();
        Pstream::gatherList(sendToLists);
        Pstream::scatterList(sendToLists);

        for (label src = 0; src < Pstream::nProcs(); ++src)
        {
            if (src == Pstream::myProcNo()) continue;
            if (sendToLists[src].found(Pstream::myProcNo()))
                recvFrom_.insert(src);
        }

        sendGids_.setSize(Pstream::nProcs());
        sendPackedIdxs_.setSize(Pstream::nProcs());

        forAllConstIter(labelHashSet, sendTo, stIt)
        {
            const label tgt = stIt.key();
            const labelHashSet& targetNeeds = neededByProc_[tgt];
            DynamicList<label> gids, pkidxs;
            forAllConstIter(labelHashSet, targetNeeds, gIt)
            {
                const label gID = gIt.key();
                const auto pkIt = globalIdToPackedIdx_.find(gID);
                if (pkIt != globalIdToPackedIdx_.end())
                {
                    gids.append(gID);
                    pkidxs.append(*pkIt);
                }
            }
            sendGids_[tgt]       = labelList(gids);
            sendPackedIdxs_[tgt] = labelList(pkidxs);
        }

        commTopologyBuilt_ = true;
    }

    interpCachesBuilt_ = true;
    Info << "buildInterpCaches: SVD weight caches built for " << nFaces << " faces." << endl;
}
