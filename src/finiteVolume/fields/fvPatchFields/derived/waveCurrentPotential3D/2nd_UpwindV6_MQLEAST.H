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
	
	//if ( db().time().timeIndex() == 1 )	Pout << "proc " << Pstream::myProcNo() << " Schemes implemented : " << schemeCode <<nl;
void Foam::waveCurrentPotential3DFvPatchScalarField::UPFDV2()
{
    const vectorField& Cf = patch().Cf();
    const scalarField x = Cf.component(vector::X);
    const scalarField y = Cf.component(vector::Y);
    const label nFaces = patch().size();
    const label patchi = patch().index();

    const volVectorField& zeta = db().lookupObjectRef<volVectorField>(zetaName_);
    const scalarField zetaZ = zeta.boundaryField()[patchi].component(vector::Z); //scalarField zetaZ(nFaces, 0.0); zetaZ[i] = x[i]*0.99;
	
	const volVectorField& Ucur = db().lookupObjectRef<volVectorField>("Ucur");
	const vectorField& Ucur_p = Ucur.boundaryField()[patchi];
	
    const volScalarField& Phi = db().lookupObjectRef<volScalarField>(PHIName_);
    const scalarField& PhiPatch_p = Phi.boundaryField()[patchi];
    scalarField PhiPatch(nFaces, 0.0);
    
	for (label i = 0; i < nFaces; i++)
    {
        PhiPatch[i] = PhiPatch_p[i];// - PhiCur_p[i];
		//zetaZ[i] = sin(2*Foam::constant::mathematical::pi*x[i]/L) + cos(2*Foam::constant::mathematical::pi*y[i]/L);//zetaZ[i] = x[i]*x[i] + y[i]*y[i];
		
		
	}

    exchangeWpGhostField(zetaZ, PhiPatch);

    scalarField DzetaZdx(nFaces, 0.0);
    scalarField DPhidx(nFaces, 0.0);
	scalarField DzetaZdy(nFaces, 0.0);
    scalarField DPhidy(nFaces, 0.0);
	//Pout << "proc " << Pstream::myProcNo() << " globalIdToWp_ " << globalIdToWp_ << nl;


    auto getWp = [&](label globalId) -> scalar
    {
        if (globalIdToPackedIdx_.found(globalId))
            return zetaZ[globalIdToPackedIdx_[globalId]];
        else if (globalIdToWp_.found(globalId))
            return globalIdToWp_[globalId];
        else
        {
			FatalErrorInFunction << "Missing zetaZ value for globalId " << globalId << " (not found in local or ghost maps)." << " This indicates that the stencil was not exchanged properly." << nl << abort(FatalError);
			return -12345; // Never reached, silences compiler warning
		}
    };

    auto getPhip = [&](label globalId) -> scalar
    {
        if (globalIdToPackedIdx_.found(globalId))
            return PhiPatch[globalIdToPackedIdx_[globalId]];
        else if (globalIdToPhi_.found(globalId))
            return globalIdToPhi_[globalId];
        else
        {
            FatalErrorInFunction << "Missing PhiPatch value for globalId " << globalId << " (not found in local or ghost maps)." << " This indicates that the stencil was not exchanged properly." << nl << abort(FatalError);
			return -12345; // Never reached, silences compiler warning
        }
    };

    auto dx = [&](label gid, scalar x0) -> scalar
    {
        const vector& c = globalIdToPackedIdx_.found(gid) ? Cf[globalIdToPackedIdx_[gid]] : globalIdToCf_[gid];
        return x0 - c.x();
    };
	
	auto dy = [&](label gid, scalar y0) -> scalar
    {
        const vector& c = globalIdToPackedIdx_.found(gid)
                            ? Cf[globalIdToPackedIdx_[gid]]
                            : globalIdToCf_[gid];
        return y0 - c.y();
    };
	
	
	auto getCandidates = [&](label gid) -> const labelList&
	{
		if (globalIdToPackedIdx_.found(gid))
		{
			// Local face → safe direct lookup
			return neighboursGlobal_[globalIdToPackedIdx_[gid]];
		}

		if (remoteNeighbours_.found(gid))
		{
			return remoteNeighbours_[gid];
		}

		if (remoteNeighbours2_.found(gid))
		{
			return remoteNeighbours2_[gid];
		}

		// FatalErrorInFunction
		// 	<< "getCandidates: No neighbours for gid "
		// 	<< gid << abort(FatalError); // Turned off... 

		static const labelList dummy;  // safe static fallback
		return dummy;
	};


	scalar shalpeVal = shape_;
	Info << " shalpeVal : " << shalpeVal << endl;
	
	
    for (label i = 0; i < nFaces; ++i)
    {
        if (!ownerHasBodyFace_[i])
		{	
			const scalar dist =  shalpeVal/patch().deltaCoeffs()[i];
			const scalar x0 = x[i];
			const scalar y0 = y[i];
			const scalar f0 = zetaZ[i], g0 = PhiPatch[i];
			const List<label>& L = upwindNodesX_[i];
			const List<label>& R = downwindNodesX_[i];
			label s = schemeCodeX_[i];

			scalar h1 = 0, h2 = 0, f1 = 0, f2 = 0, g1 = 0, g2 = 0;

			if (s == 6)
			{
				
				label up1 = L[0];
				List<label> up2List;

				if (globalIdToPackedIdx_.found(up1))
				{
					
					up2List = upwindNodesX_[globalIdToPackedIdx_[up1]];
				
				}
				else
				{
					const auto it = remoteIdToUpwindNodesX_.find(globalIdOfLocalFace_[i]);
					if (it != remoteIdToUpwindNodesX_.end())
						up2List = it();
					else
						FatalErrorInFunction << "Missing second upwind node for globalId " << globalIdOfLocalFace_[i]
											 << " in scheme 6. This violates detection guarantee." << nl << abort(FatalError);
				}
				////INTERPOLATION CHECK
				 
				// f1 and g1
				if (interpNeededX_[i][0])
				{
					const vector& cUp1 = globalIdToPackedIdx_.found(up1) ? Cf[globalIdToPackedIdx_[up1]] : globalIdToCf_[up1];
					const scalar x_up1 = cUp1.x();
					const List<label>& candsX1 = getCandidates(up1); 
					
					
					interpolatePairMQAturgut(x_up1, y0, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,up1, candsX1, dist,f1, g1);
					//f1 = interpolateMQAturgut(x_up1, y0, Cf, globalIdToPackedIdx_, globalIdToCf_, getWp, up1 ,candsX1,dist);
					//g1 = interpolateMQAturgut(x_up1, y0, Cf, globalIdToPackedIdx_, globalIdToCf_, getPhip,up1 ,candsX1,dist);
				}
				else
				{
					f1 = getWp(up1);
					g1 = getPhip(up1);
				}

				// f2 and g2
				label up2 = up2List[0];
				if (interpNeededX_[i][1])
				{
					const vector& cUp2 = globalIdToPackedIdx_.found(up2)? Cf[globalIdToPackedIdx_[up2]]: globalIdToCf_[up2];
					const scalar x_up2 = cUp2.x();
					const List<label>& candsX1 = getCandidates(up2); 
					
					
					interpolatePairMQAturgut(x_up2, y0, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,up2, candsX1, dist,f2, g2);
					
				}
				else
				{
					f2 = getWp(up2);
					g2 = getPhip(up2);
				}			
				
				////INTERPOLATION CHECK
			   
				h1 = mag(dx(up1, x0)); h2 = mag(dx(up2, x0));
				DzetaZdx[i] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h2 - h1));
				DPhidx[i] = (-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h2 - h1));
				
							 
			}
			else if (s == 1)
			{
				
				////INTERPOLATION CHECK
				 
				// f1 and g1
				if (interpNeededX_[i][0])
				{
					label up1 = L[0];
					const vector& cUp1 = globalIdToPackedIdx_.found(up1) ? Cf[globalIdToPackedIdx_[up1]] : globalIdToCf_[up1];
					const scalar x_up1 = cUp1.x();
					
					const List<label>& candsX1 = getCandidates(up1);

					interpolatePairMQAturgut(x_up1, y0, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,up1, candsX1, dist,f1, g1);				
					
				}
				else
				{
					f1 = getWp(L[0]);
					g1 = getPhip(L[0]);
				}

				// f2 and g2
				
				if (interpNeededX_[i][1])
				{
					label dn1 = R[0];
					const vector& cDn1 = globalIdToPackedIdx_.found(dn1)? Cf[globalIdToPackedIdx_[dn1]]: globalIdToCf_[dn1];
					const scalar x_dn1 = cDn1.x();
					const List<label>& candsX1 = getCandidates(dn1); 
					
					interpolatePairMQAturgut(x_dn1, y0, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,dn1, candsX1, dist,f2, g2);
					
				}
				else
				{
					f2 = getWp(R[0]);
					g2 = getPhip(R[0]);
				}			
				
				////INTERPOLATION CHECK
				
				h1 = mag(dx(L[0], x0));
				h2 = mag(dx(R[0], x0));
			   
				DzetaZdx[i] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h1 + h2));
				DPhidx[i] = (-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h1 + h2));
			}
			else if (s == 7)
			{
				label down1 = R[0];
				List<label> down2List;

				if (globalIdToPackedIdx_.found(down1))
				{
					down2List = downwindNodesX_[globalIdToPackedIdx_[down1]];
				}
				else
				{
					const auto it = remoteIdToDownwindNodesX_.find(globalIdOfLocalFace_[i]);
					if (it != remoteIdToDownwindNodesX_.end())
					{
						down2List = it();
					}
					else
					{
						FatalErrorInFunction
							<< "Missing second downwind node for globalId " << globalIdOfLocalFace_[i]
							<< " in scheme 7. This violates detection guarantee." << nl
							<< abort(FatalError);
					}
				}
				////INTERPOLATION CHECK
				 
				// f1 and g1
				if (interpNeededX_[i][0])
				{
					const vector& cUp1 = globalIdToPackedIdx_.found(down1) ? Cf[globalIdToPackedIdx_[down1]] : globalIdToCf_[down1];
					const scalar x_down1 = cUp1.x();
					const List<label>& candsX1 = getCandidates(down1); 
					
					interpolatePairMQAturgut(x_down1, y0, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,down1, candsX1, dist,f1, g1);
					
				
				}
				else
				{
					f1 = getWp(down1); 
					g1 = getPhip(down1); 
				}

				// f2 and g2
				label down2 = down2List[0];
				if (interpNeededX_[i][1])
				{
					const vector& cUp2 = globalIdToPackedIdx_.found(down2)? Cf[globalIdToPackedIdx_[down2]]: globalIdToCf_[down2];
					const scalar x_down2 = cUp2.x();
					const List<label>& candsX1 = getCandidates(down2); 
					
					interpolatePairMQAturgut(x_down2, y0, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,down2, candsX1, dist,f2, g2);
					
				}
				else
				{
					f2 = getWp(down2);
					g2 = getPhip(down2);
				}			
				
				////INTERPOLATION CHECK
				
				
				h1 = mag(dx(down1, x0)); h2 = mag(dx(down2, x0));
				DzetaZdx[i] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h1 - h2));
				DPhidx[i] = (-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h1 - h2));
				
			}
			else if (s == 4 && L.size())
			{
				h1 = mag(dx(L[0], x0));f1 = getWp(L[0]);g1 = getPhip(L[0]);
				DzetaZdx[i] = (f0 - f1) / h1;DPhidx[i] = (g0 - g1) / h1;
				
			}
			else if (s == 5 && R.size())
			{
				h1 = mag(dx(R[0], x0));f1 = getWp(R[0]);g1 = getPhip(R[0]);
				DzetaZdx[i] = (f1 - f0) / h1;DPhidx[i] = (g1 - g0) / h1;
			}
			else
			{
				 FatalErrorInFunction
						<< "Face " << i << ": No stencil match for schemeCode = " << s
						<< ". This indicates a bug or unsupported stencil. Aborting!" << nl
						<< "Global face ID: " << globalIdOfLocalFace_[i] << nl
						<< "Upwind nodes: " << upwindNodesX_[i] << nl
						<< "Downwind nodes: " << downwindNodesX_[i] << nl
						<< abort(FatalError);
			}
			
			
			if (ownerHasInletFace_[i]){
				
				const scalar rampperiod =  params_.rampperiod;
				const scalar U0          = params_.U0;
				const scalar head_ang   = params_.head_ang;
				const scalar wavelength  = params_.wavelength;
				const scalar steepness  = params_.steepness;
				const scalar hdepth  = params_.hdepth;
				
				const scalar wavenumber(2.0 * Foam::constant::mathematical::pi / wavelength);
				const scalar amp(0.5 * steepness * wavelength);
				const scalar w(Foam::sqrt(9.81 * wavenumber * Foam::tanh(wavenumber*hdepth)));
				const scalar celerity(w/wavenumber); 
				
				//Info << "Input wave given on:" << x0 << endl; 
				scalar ramp_factor = 1.0;
				if (rampperiod > 0){
					const scalar T(2.0*Foam::constant::mathematical::pi/w);
					const scalar ramp_time(rampperiod * T);
					ramp_factor = 0.5 * (1 - cos(Foam::constant::mathematical::pi * min(1.0, db().time().value() / ramp_time)));
				//Info << "Ramp implemented for nPeriod:" << rampperiod << endl;
				
				}
				
				DzetaZdx[i] = -amp * wavenumber * Foam::sin(wavenumber * x0 -(w + U0 * wavenumber * Foam::cos(head_ang)) * db().time().value())*ramp_factor;
				DPhidx[i] =  amp * wavenumber * (9.81 / w) * Foam::cos(wavenumber * x0 - (w + U0 * wavenumber * Foam::cos(head_ang)) * db().time().value())*ramp_factor;
			}

			// if (ownerHasInletFace_[i]){
				
			// 	const scalar rampperiod =  params_.rampperiod;
			// 	const scalar U0          = params_.U0;
			// 	const scalar head_ang   = params_.head_ang;
			// 	const scalar wavelength  = params_.wavelength;
			// 	const scalar steepness  = params_.steepness;
			// 	const scalar hdepth  = params_.hdepth;
				
			// 	const scalar wavenumber(2.0 * Foam::constant::mathematical::pi / wavelength);
			// 	const scalar amp(0.5 * steepness * wavelength);
			// 	const scalar w(Foam::sqrt(9.81 * wavenumber * Foam::tanh(wavenumber*hdepth)));
			// 	const scalar celerity(w/wavenumber); 
				
			// 	//Info << "Input wave given on:" << x0 << endl; 
			// 	scalar ramp_factor = 1.0;
			// 	if (rampperiod > 0){
			// 		const scalar T(2.0*Foam::constant::mathematical::pi/w);
			// 		const scalar ramp_time(rampperiod * T);
			// 		ramp_factor = 0.5 * (1 - cos(Foam::constant::mathematical::pi * min(1.0, db().time().value() / ramp_time)));
			// 	//Info << "Ramp implemented for nPeriod:" << rampperiod << endl;
				
			// 	}
				
				
				
			// 	DzetaZdx[i] = -amp * wavenumber * Foam::sin(wavenumber * x0 -(w + U0 * wavenumber * Foam::cos(head_ang)) * db().time().value())*ramp_factor;
			// 	DPhidx[i] =  amp * wavenumber * (9.81 / w) * Foam::cos(wavenumber * x0 - (w + U0 * wavenumber * Foam::cos(head_ang)) * db().time().value())*ramp_factor;
			// }
		}
	}
	
		
			

	for (label i = 0; i < nFaces; ++i)
    {
		if (!ownerHasBodyFace_[i])
		{	
			const scalar dist =  shalpeVal/patch().deltaCoeffs()[i];
			
			const scalar x0 = x[i];
			const scalar y0 = y[i];
			const scalar f0 = zetaZ[i], g0 = PhiPatch[i];
			const List<label>& L = upwindNodesY_[i];
			const List<label>& R = downwindNodesY_[i];
			label s = schemeCodeY_[i];

			scalar h1 = 0, h2 = 0, f1 = 0, f2 = 0, g1 = 0, g2 = 0;

			if (s == 6 )
			{
				label up1 = L[0];
				List<label> up2List;

				if (globalIdToPackedIdx_.found(up1))
				{
					up2List = upwindNodesY_[globalIdToPackedIdx_[up1]];
				}
				else
				{
					const auto it = remoteIdToUpwindNodesY_.find(globalIdOfLocalFace_[i]);
					if (it != remoteIdToUpwindNodesY_.end())
						up2List = it();
					else
						FatalErrorInFunction << "Missing 2nd upwind node for face " << i << nl << abort(FatalError);
				}

				label up2 = up2List[0];
				
				// f1,g1 (at up1)
				if (interpNeededY_[i][0])
				{
					const vector& cUp1 = globalIdToPackedIdx_.found(up1)
										   ? Cf[globalIdToPackedIdx_[up1]]
										   : globalIdToCf_[up1];
					const scalar y_up1 = cUp1.y();

					const List<label>& candsY1 = getCandidates(up1);
					
					
					interpolatePairMQAturgut(x0, y_up1, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,up1, candsY1, dist,f1, g1);
					
				}
				else
				{
					f1 = getWp(up1);
					g1 = getPhip(up1);
				}

				// f2,g2 (at up2)
				if (interpNeededY_[i][1])
				{
					const vector& cUp2 = globalIdToPackedIdx_.found(up2)
										   ? Cf[globalIdToPackedIdx_[up2]]
										   : globalIdToCf_[up2];
					const scalar y_up2 = cUp2.y();

					const List<label>& candsY2 = getCandidates(up2);

					interpolatePairMQAturgut(x0, y_up2, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,up2, candsY2, dist,f2, g2);
				
				}
				else
				{
					f2 = getWp(up2);
					g2 = getPhip(up2);
				}

				
				h1 = mag(dy(up1, y0)); h2 = mag(dy(up2, y0));
				DzetaZdy[i] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h2 - h1));
				DPhidy[i] = (-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h2 - h1));
			}
			else if (s == 1)
			{
				
					if (interpNeededY_[i][0])
				{
					label up1 = L[0];
					const vector& cUp1 = globalIdToPackedIdx_.found(up1)
										   ? Cf[globalIdToPackedIdx_[up1]]
										   : globalIdToCf_[up1];
					const scalar y_up1 = cUp1.y();

					const List<label>& candsY1 = getCandidates(up1);
					
					interpolatePairMQAturgut(x0, y_up1, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,up1, candsY1, dist,f1, g1);
					
				}
				else
				{
					f1 = getWp(L[0]);
					g1 = getPhip(L[0]);
				}

				// f2,g2
				if (interpNeededY_[i][1])
				{
					label dn1 = R[0];
					const vector& cDn1 = globalIdToPackedIdx_.found(dn1)
										   ? Cf[globalIdToPackedIdx_[dn1]]
										   : globalIdToCf_[dn1];
					const scalar y_dn1 = cDn1.y();

					const List<label>& candsY2 = getCandidates(dn1);
					
					interpolatePairMQAturgut(x0, y_dn1, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,dn1, candsY2, dist,f2, g2);
					
				}
				else
				{
					f2 = getWp(R[0]);
					g2 = getPhip(R[0]);
				}
				h1 = mag(dy(L[0], y0));
				h2 = mag(dy(R[0], y0));
				
				DzetaZdy[i] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h1 + h2));
				DPhidy[i] = (-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h1 + h2));
			}
			else if (s == 7)
			{
				label down1 = R[0];
				List<label> down2List;

				if (globalIdToPackedIdx_.found(down1))
				{
					down2List = downwindNodesY_[globalIdToPackedIdx_[down1]];
				}
				else
				{
					const auto it = remoteIdToDownwindNodesY_.find(globalIdOfLocalFace_[i]);
					if (it != remoteIdToDownwindNodesY_.end())
						down2List = it();
					else
						FatalErrorInFunction << "Missing 2nd downwind node for face " << i << nl << abort(FatalError);
				}

				label down2 = down2List[0];
				
				// f1,g1
				if (interpNeededY_[i][0])
				{
					const vector& cDn1 = globalIdToPackedIdx_.found(down1)
										   ? Cf[globalIdToPackedIdx_[down1]]
										   : globalIdToCf_[down1];
					const scalar y_dn1 = cDn1.y();

					const List<label>& candsY1 = getCandidates(down1);
					
					
					interpolatePairMQAturgut(x0, y_dn1, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,down1, candsY1, dist,f1, g1);
					
				}
				else
				{
					f1 = getWp(down1);
					g1 = getPhip(down1);
				}

				// f2,g2
				if (interpNeededY_[i][1])
				{
					const vector& cDn2 = globalIdToPackedIdx_.found(down2)
										   ? Cf[globalIdToPackedIdx_[down2]]
										   : globalIdToCf_[down2];
					const scalar y_dn2 = cDn2.y();

					const List<label>& candsY2 = getCandidates(down2);
					
					interpolatePairMQAturgut(x0, y_dn2, Cf,globalIdToPackedIdx_, globalIdToCf_,getWp, getPhip,down2, candsY2, dist,f2, g2);
				
				}
				else
				{
					f2 = getWp(down2);
					g2 = getPhip(down2);
				}

				
				
				h1 = mag(dy(down1, y0)); h2 = mag(dy(down2, y0));
				DzetaZdy[i] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h1 - h2));
				DPhidy[i] = (-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h1 - h2));
			}
			else if (s == 4 && L.size())
			{
				h1 = mag(dy(L[0], y0));
				f1 = getWp(L[0]); g1 = getPhip(L[0]);
				DzetaZdy[i] = (f0 - f1) / h1;
				DPhidy[i] = (g0 - g1) / h1;
			}
			else if (s == 5 && R.size())
			{
				h1 = mag(dy(R[0], y0));
				f1 = getWp(R[0]); g1 = getPhip(R[0]);
				DzetaZdy[i] = (f1 - f0) / h1;
				DPhidy[i] = (g1 - g0) / h1;
			}
			else
			{
				FatalErrorInFunction << "Face " << i << ": Invalid schemeCodeY = " << s << nl << abort(FatalError);
			}
			
			if (ownerHasInletFace_[i]){ 
				DzetaZdy[i] = 0.0;
				DPhidy[i] = 0.0;
			}
		
		}
		
    }

	for (label i = 0; i < nFaces; ++i)
    {
		
	    if (ownerHasBodyFace_[i] || ownerHasBodyEdge_[i] || ownerHasBodyVertex_[i])
		{ 
			
			const label selfGid = globalIdOfLocalFace_[i];
			if (lsMergedCandidates_[i].empty()){
				 Pout << "lsMergedCandidates_ should not be empty local i: " 
				 << i << " GlobalID : " << globalIdOfLocalFace_[i] << nl << endl;
				 continue;
			}
			
			const scalar dist =  shalpeVal/patch().deltaCoeffs()[i];
			scalar fx = 0, fy = 0, gx = 0, gy = 0;	
			// DerivativeMqRBFxy(dist,Cf,globalIdToPackedIdx_,globalIdToCf_,getWp,getPhip,selfGid,lsMergedCandidates_[i],fx,fy,gx,gy);
			// DerivativeLSTaylor(Cf,globalIdToPackedIdx_,globalIdToCf_,getWp,getPhip,selfGid,lsMergedCandidates_[i],fx,fy,gx,gy);
			// DerivativeLSRobustTaylor(
			// 	Cf, globalIdToPackedIdx_, globalIdToCf_,
			// 	getWp, getPhip,
			// 	selfGid, lsMergedCandidates_[i],
			// 	fx, fy, gx, gy
			// );

			GreenGaussUpwindXY
			(
				patch(),                 // fvPatch
				i,                       // face index on patch
				selfGid,
				edgeNeighbours_[i],
				edgeNormals2D_[i],
				Ucur_p[i],               // advecting velocity at this patch face
				getWp,
				getPhip,
				fx, fy, gx, gy
			);

			DzetaZdx[i] = fx;
			DPhidx[i] = gx;
			DzetaZdy[i] = fy;
			DPhidy[i] = gy;
		}
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

