/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
centered FD bu TURGUT

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
	
	if ( db().time().timeIndex() == 1 )
	{
	  for (label i = 0; i < nFaces; ++i) {
	  //if (Pstream::myProcNo()==2)	
	  if (schemeCode[i]!=6)	  Pout << "proc " << Pstream::myProcNo() << " Schemes implemented[i] : " << schemeCode[i] << "at location:" << x[i] <<nl;
	  }
	   //for (label i = 0; i < nFaces; ++i) {
	   //Pout << "proc " << Pstream::myProcNo() << " PhiDx_[i] : " << PhiDx_[i] << " at location:" << i <<nl;
	   //}
	 }
		if (globalIdToPackedIdx_.found(up1)){
				const vector& cUp1 = Cf[globalIdToPackedIdx_[up1]];}
				else {
				const vector& cUp1 = globalIdToCf_[up1];}	
				
				
	
				
\*---------------------------------------------------------------------------*/
	
void Foam::waveCurrentPotential3DFvPatchScalarField::UPFDV2()
{
    const vectorField& Cf = patch().Cf();
    const scalarField x = Cf.component(vector::X);
    const scalarField y = Cf.component(vector::Y);
    const label nFaces = patch().size();
    const label patchi = patch().index();

    const volVectorField& zeta = db().lookupObjectRef<volVectorField>(zetaName_);
    const scalarField zetaZ = zeta.boundaryField()[patchi].component(vector::Z);

    const volScalarField& Phi = db().lookupObjectRef<volScalarField>(PHIName_);
    const scalarField& PhiPatch_p = Phi.boundaryField()[patchi];

    scalarField PhiPatch(nFaces, 0.0);
    for (label i = 0; i < nFaces; ++i) PhiPatch[i] = PhiPatch_p[i];

    exchangeWpGhostField(zetaZ, PhiPatch);

    scalarField DzetaZdx(nFaces, 0.0), DzetaZdy(nFaces, 0.0);
    scalarField DPhidx(nFaces, 0.0),   DPhidy(nFaces, 0.0);

    auto getWp = [&](label globalId) -> scalar
    {
        if (globalIdToPackedIdx_.found(globalId)) return zetaZ[globalIdToPackedIdx_[globalId]];
        if (globalIdToWp_.found(globalId))        return globalIdToWp_[globalId];
        FatalErrorInFunction << "Missing zetaZ for gid " << globalId << nl << abort(FatalError);
        return 0.0;
    };

    auto getPhip = [&](label globalId) -> scalar
    {
        if (globalIdToPackedIdx_.found(globalId)) return PhiPatch[globalIdToPackedIdx_[globalId]];
        if (globalIdToPhi_.found(globalId))       return globalIdToPhi_[globalId];
        FatalErrorInFunction << "Missing PhiPatch for gid " << globalId << nl << abort(FatalError);
        return 0.0;
    };

    auto getC = [&](label gid) -> const vector&
    {
        if (globalIdToPackedIdx_.found(gid)) return Cf[globalIdToPackedIdx_[gid]];
        if (globalIdToCf_.found(gid))        return globalIdToCf_[gid];
        FatalErrorInFunction << "Missing Cf for gid " << gid << nl << abort(FatalError);
        static const vector dummy(Zero);
        return dummy;
    };

    auto getCandidates = [&](label gid) -> const labelList&
    {
        if (globalIdToPackedIdx_.found(gid)) return neighboursGlobal_[globalIdToPackedIdx_[gid]];
        if (remoteNeighbours_.found(gid))    return remoteNeighbours_[gid];
        if (remoteNeighbours2_.found(gid))   return remoteNeighbours2_[gid];

        FatalErrorInFunction << "No neighbours for gid " << gid << nl << abort(FatalError);
        static const labelList dummy;
        return dummy;
    };

    // ---- Robustness knobs (start here)
    const scalar regCoeff   = 1e-2;  // 1e-3..1e-1 (bigger = more stable, more diffusive)
    const scalar wMaxFactor = 50.0;  // caps max weight relative to typical 1/h^2 scale
    const scalar Cgrad      = 5.0;   // gradient limiter strength (3..10)

    for (label i = 0; i < nFaces; ++i)
    {
        const scalar x0 = x[i];
        const scalar y0 = y[i];
        const scalar f0 = zetaZ[i];
        const scalar g0 = PhiPatch[i];
        const label  gid0 = globalIdOfLocalFace_[i];

        // Prefer your merged LSQ stencil if available
        const labelList* Nptr = nullptr;
        if (i < lsMergedCandidates_.size() && !lsMergedCandidates_[i].empty())
            Nptr = &lsMergedCandidates_[i];
        else
            Nptr = &getCandidates(gid0);

        const labelList& N = *Nptr;

        // local spacing estimate (z=0 patch, so XY spacing is fine)
        const scalar h   = 1.0/max(VSMALL, patch().deltaCoeffs()[i]);
        const scalar eps = sqr(1e-6*h);                 // MUCH larger than before (important)
        const scalar wMax = wMaxFactor/(h*h + VSMALL);  // cap weights

        scalar M00=0.0, M01=0.0, M11=0.0;
        scalar bf0=0.0, bf1=0.0;
        scalar bg0=0.0, bg1=0.0;

        scalar dfMax = 0.0;
        scalar dgMax = 0.0;

        forAll(N, k)
        {
            const label gidj = N[k];
            if (gidj == gid0) continue;

            const vector& Cj = getC(gidj);
            const scalar dx = Cj.x() - x0;
            const scalar dy = Cj.y() - y0;
            const scalar r2 = dx*dx + dy*dy;

            if (r2 < eps) continue;

            // bounded weight
            scalar w = 1.0/(r2 + eps);
            w = min(w, wMax);

            const scalar fj = getWp(gidj);
            const scalar gj = getPhip(gidj);

            const scalar df = fj - f0;
            const scalar dg = gj - g0;

            dfMax = max(dfMax, mag(df));
            dgMax = max(dgMax, mag(dg));

            M00 += w*dx*dx;
            M01 += w*dx*dy;
            M11 += w*dy*dy;

            bf0 += w*dx*df;
            bf1 += w*dy*df;

            bg0 += w*dx*dg;
            bg1 += w*dy*dg;
        }

        // Tikhonov regularization
        const scalar trace = M00 + M11;
        const scalar lambda = regCoeff*trace + VSMALL;

        M00 += lambda;
        M11 += lambda;

        const scalar det = M00*M11 - M01*M01;

        scalar fx=0.0, fy=0.0, gx=0.0, gy=0.0;

        if (mag(det) > VSMALL && trace > VSMALL)
        {
            fx = ( M11*bf0 - M01*bf1)/det;
            fy = (-M01*bf0 + M00*bf1)/det;

            gx = ( M11*bg0 - M01*bg1)/det;
            gy = (-M01*bg0 + M00*bg1)/det;
        }
        else
        {
            // If it’s STILL singular, zero it (better than blow-up)
            fx = fy = gx = gy = 0.0;
        }

        // Gradient limiter: |grad| <= C * dfMax / h
        // (prevents “one bad face” from destroying the run)
        const scalar fLim = (dfMax > 0 ? (Cgrad*dfMax/(h + VSMALL)) : 0.0);
        const scalar gLim = (dgMax > 0 ? (Cgrad*dgMax/(h + VSMALL)) : 0.0);

        const scalar magF = Foam::sqrt(fx*fx + fy*fy);
        if (magF > fLim && magF > VSMALL)
        {
            const scalar s = fLim/magF;
            fx *= s; fy *= s;
        }

        const scalar magG = Foam::sqrt(gx*gx + gy*gy);
        if (magG > gLim && magG > VSMALL)
        {
            const scalar s = gLim/magG;
            gx *= s; gy *= s;
        }

        // Inlet override: keep your exact behavior
        if (ownerHasInletFace_[i])
        {
            const scalar rampperiod = params_.rampperiod;
            const scalar U0         = params_.U0;
            const scalar wavelength = params_.wavelength;
            const scalar steepness  = params_.steepness;
            const scalar hdepth     = params_.hdepth;

            const scalar kWave(2.0*Foam::constant::mathematical::pi / wavelength);
            const scalar amp(0.5 * steepness * wavelength);
            const scalar w(Foam::sqrt(9.81 * kWave * Foam::tanh(kWave*hdepth)));

            scalar ramp_factor = 1.0;
            if (rampperiod > 0)
            {
                const scalar T(2.0*Foam::constant::mathematical::pi/w);
                const scalar ramp_time(rampperiod * T);
                ramp_factor = 0.5 * (1 - cos(Foam::constant::mathematical::pi
                               * min(1.0, db().time().value()/ramp_time)));
            }

            fx = -amp * kWave
               * Foam::sin(kWave*x0 - (w + U0*kWave)*db().time().value())
               * ramp_factor;

            gx =  amp * kWave * (9.81/w)
               * Foam::cos(kWave*x0 - (w + U0*kWave)*db().time().value())
               * ramp_factor;

            fy = 0.0;
            gy = 0.0;
        }

        DzetaZdx[i] = fx;
        DzetaZdy[i] = fy;
        DPhidx[i]   = gx;
        DPhidy[i]   = gy;
    }

    for (label j = 0; j < nFaces; ++j)
    {
        zetaDx_[j][vector::X] = DzetaZdx[j];
        zetaDx_[j][vector::Y] = DzetaZdy[j];

        PhiDx_[j] = -DPhidx[j];
        PhiDy_[j] = -DPhidy[j];
    }
}


void Foam::waveCurrentPotential3DFvPatchScalarField::exchangeWpGhostField(const scalarField& zetaZ,const scalarField& Phi)
{
    if (!Pstream::parRun()) return;

    globalIdToWp_.clear();
	globalIdToPhi_.clear();
    



    // --- Gather communication partners
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
	
	
	
	

	//using TableType = HashTable<List<label>, label>;

	
	const List<labelHashSet>& neededByProc = neededByProc_;
	
	
	// --- Send phase
    PstreamBuffers pBuf(Pstream::commsTypes::nonBlocking);

    forAllConstIter(labelHashSet, sendTo, it)
	{
		const label tgt = it.key();
		UOPstream os(tgt, pBuf);

		Foam::Map<scalar> zetaSendMap;
		Foam::Map<scalar> phiSendMap;

		

		// NEW: Also send what processor tgt requested
		const labelHashSet& targetNeeds = neededByProc[tgt];

		forAllConstIter(labelHashSet, targetNeeds, itG)
		{
			const label gID = itG.key();
			if (globalIdToPackedIdx_.found(gID))
			{
				const label packedIdx = globalIdToPackedIdx_[gID];
				zetaSendMap.insert(gID, zetaZ[packedIdx]);
				phiSendMap.insert(gID, Phi[packedIdx]);
			}
		}

		// Send packed data
		os << zetaSendMap.size();
		forAllConstIter(Foam::Map<scalar>, zetaSendMap, it)
		{
			const label gid = it.key();
			os << gid << it();
			os << phiSendMap[gid];  // Make sure phiSendMap has the same keys
		}
		//if ( Pstream::myProcNo()==0) Pout << "proc " << Pstream::myProcNo() << "target:" << tgt << " zetaSendMap " << zetaSendMap << endl;
	
	}
	
    pBuf.finishedSends();

    // --- Receive phase
    // --- Receive phase
    forAllConstIter(labelHashSet, recvFrom, it)
    {
        const label src = it.key();
        UIPstream is(src, pBuf);

        label nRecv; is >> nRecv;
		for (label j = 0; j < nRecv; ++j)
		{
			label gid;
			scalar zetaVal, phiVal;
			is >> gid >> zetaVal >> phiVal;
			globalIdToWp_.insert(gid, zetaVal);
			globalIdToPhi_.insert(gid, phiVal);
		}
    }
    // Optional debug:
    //if ( Pstream::myProcNo()==1) Pout << "proc " << Pstream::myProcNo() << " globalIdToWp_ " << globalIdToWp_ << endl;
}

