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

void Foam::waveCurrentPotential3DFvPatchScalarField::UPFDV2()
{
    const vectorField& Cf = patch().Cf();
    const scalarField x = Cf.component(vector::X);
    const scalarField y = Cf.component(vector::Y);
    const label nFaces = patch().size();
    const label patchi = patch().index();

    const volVectorField& zeta = db().lookupObjectRef<volVectorField>(zetaName_);
    const scalarField zetaZ = zeta.boundaryField()[patchi].component(vector::Z);

    const volVectorField& Ucur = db().lookupObjectRef<volVectorField>("Ucur");
    const vectorField& Ucur_p = Ucur.boundaryField()[patchi];

    const volScalarField& Phi = db().lookupObjectRef<volScalarField>(PHIName_);
    const scalarField& PhiPatch_p = Phi.boundaryField()[patchi];
    scalarField PhiPatch(nFaces, 0.0);

    for (label i = 0; i < nFaces; i++)
    {
        PhiPatch[i] = PhiPatch_p[i];
    }

    exchangeWpGhostField(zetaZ, PhiPatch);

    // ---- Change 3: single find() helper for Cf lookup ----
    auto getCf = [&](label gid) -> const vector&
    {
        const auto it = globalIdToPackedIdx_.find(gid);
        if (it != globalIdToPackedIdx_.end()) return Cf[*it];
        const auto it2 = globalIdToCf_.find(gid);
        if (it2 != globalIdToCf_.end()) return *it2;
        FatalErrorInFunction << "gid " << gid << " not found in Cf maps" << abort(FatalError);
        static const vector dummy(Zero); return dummy;
    };

    // ---- Change 3: getWp/getPhip using find() ----
    auto getWp = [&](label globalId) -> scalar
    {
        const auto it = globalIdToPackedIdx_.find(globalId);
        if (it != globalIdToPackedIdx_.end())
            return zetaZ[*it];
        const auto it2 = globalIdToWp_.find(globalId);
        if (it2 != globalIdToWp_.end())
            return *it2;
        FatalErrorInFunction << "Missing zetaZ value for globalId " << globalId
            << " (not found in local or ghost maps)."
            << " This indicates that the stencil was not exchanged properly."
            << nl << abort(FatalError);
        return -12345;
    };

    auto getPhip = [&](label globalId) -> scalar
    {
        const auto it = globalIdToPackedIdx_.find(globalId);
        if (it != globalIdToPackedIdx_.end())
            return PhiPatch[*it];
        const auto it2 = globalIdToPhi_.find(globalId);
        if (it2 != globalIdToPhi_.end())
            return *it2;
        FatalErrorInFunction << "Missing PhiPatch value for globalId " << globalId
            << " (not found in local or ghost maps)."
            << " This indicates that the stencil was not exchanged properly."
            << nl << abort(FatalError);
        return -12345;
    };

    // ---- Change 5: combined getWpPhip helper ----
    auto getWpPhip = [&](label gid, scalar& f, scalar& g)
    {
        const auto it = globalIdToPackedIdx_.find(gid);
        if (it != globalIdToPackedIdx_.end())
        {
            f = zetaZ[*it];
            g = PhiPatch[*it];
            return;
        }
        const auto it2 = globalIdToWp_.find(gid);
        if (it2 != globalIdToWp_.end()) f = *it2;
        else FatalErrorInFunction << "Missing zetaZ for gid " << gid << abort(FatalError);
        const auto it3 = globalIdToPhi_.find(gid);
        if (it3 != globalIdToPhi_.end()) g = *it3;
        else FatalErrorInFunction << "Missing Phi for gid " << gid << abort(FatalError);
    };

    auto dx = [&](label gid, scalar x0) -> scalar
    {
        return x0 - getCf(gid).x();
    };

    auto dy = [&](label gid, scalar y0) -> scalar
    {
        return y0 - getCf(gid).y();
    };

    auto getCandidates = [&](label gid) -> const labelList&
    {
        const auto it = globalIdToPackedIdx_.find(gid);
        if (it != globalIdToPackedIdx_.end())
            return neighboursGlobal_[*it];
        if (remoteNeighbours_.found(gid))
            return remoteNeighbours_[gid];
        if (remoteNeighbours2_.found(gid))
            return remoteNeighbours2_[gid];
        static const labelList dummy;
        return dummy;
    };

    scalar shalpeVal = shape_;

    // ---- Change 2: merged single loop over faces for X and Y ----
    for (label i = 0; i < nFaces; ++i)
    {
        if (!ownerHasBodyFace_[i])
        {
            const scalar dist = shalpeVal / patch().deltaCoeffs()[i];
            const scalar x0 = x[i];
            const scalar y0 = y[i];
            const scalar f0 = zetaZ[i], g0 = PhiPatch[i];

            // ============================================================
            // X-direction
            // ============================================================
            {
                const List<label>& L = upwindNodesX_[i];
                const List<label>& R = downwindNodesX_[i];
                label s = label(schemeCodeX_[i]);

                scalar h1 = 0, h2 = 0, f1 = 0, f2 = 0, g1 = 0, g2 = 0;
                scalar f3 = 0, g3 = 0, f4 = 0, g4 = 0;

                if (s == 6 || (s == 14 && Ucur_p[i].x() > SMALL
                    && (interpNeeded4thX_[i][0] || interpNeeded4thX_[i][1]
                        || interpNeeded4thX_[i][2] || interpNeeded4thX_[i][3])))
                {
                    label up1 = L[0];

                    // ---- Change 4: pointer instead of copy ----
                    const List<label>* up2Ptr = nullptr;
                    {
                        const auto pkIt = globalIdToPackedIdx_.find(up1);
                        if (pkIt != globalIdToPackedIdx_.end())
                            up2Ptr = &upwindNodesX_[*pkIt];
                        else
                        {
                            const auto it = remoteIdToUpwindNodesX_.find(globalIdOfLocalFace_[i]);
                            if (it != remoteIdToUpwindNodesX_.end())
                                up2Ptr = &it();
                            else
                                FatalErrorInFunction
                                    << "Missing second upwind node for globalId "
                                    << globalIdOfLocalFace_[i]
                                    << " in scheme 6. This violates detection guarantee."
                                    << nl << abort(FatalError);
                        }
                    }

                    // s==14 fallback: use 4th-order cache slots (up1=0, up2=1)
                    const InterpWeightCache& cUp1X = (s == 14)
                        ? interp4thCacheX_[i][0] : interp2ndCacheX_[i][0];
                    const InterpWeightCache& cUp2X = (s == 14)
                        ? interp4thCacheX_[i][1] : interp2ndCacheX_[i][1];

                    // f1, g1
                    if (interpNeededX_[i][0])
                    {
                        if (cUp1X.valid)
                        {
                            f1 = g1 = 0.0;
                            forAll(cUp1X.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cUp1X.gids[j], ww, pp);
                                f1 += cUp1X.w[j]*ww;
                                g1 += cUp1X.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_up1 = getCf(up1).x();
                            const List<label>& candsX1 = getCandidates(up1);
                            interpolatePairMQAturgut(x_up1, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up1, candsX1, dist, f1, g1);
                        }
                    }
                    else
                    {
                        getWpPhip(up1, f1, g1);
                    }

                    // f2, g2
                    label up2 = (*up2Ptr)[0];
                    if (interpNeededX_[i][1])
                    {
                        if (cUp2X.valid)
                        {
                            f2 = g2 = 0.0;
                            forAll(cUp2X.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cUp2X.gids[j], ww, pp);
                                f2 += cUp2X.w[j]*ww;
                                g2 += cUp2X.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_up2 = getCf(up2).x();
                            const List<label>& candsX2 = getCandidates(up2);
                            interpolatePairMQAturgut(x_up2, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up2, candsX2, dist, f2, g2);
                        }
                    }
                    else
                    {
                        getWpPhip(up2, f2, g2);
                    }

                    h1 = mag(dx(up1, x0)); h2 = mag(dx(up2, x0));
                    zetaDx_[i][vector::X] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h2 - h1));
                    PhiDx_[i] = -((-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h2 - h1)));
                }
                else if (s == 1)
                {
                    // f1, g1
                    if (interpNeededX_[i][0])
                    {
                        label up1 = L[0];
                        if (interp2ndCacheX_[i][0].valid)
                        {
                            const auto& cache = interp2ndCacheX_[i][0];
                            f1 = g1 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f1 += cache.w[j]*ww;
                                g1 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_up1 = getCf(up1).x();
                            const List<label>& candsX1 = getCandidates(up1);
                            interpolatePairMQAturgut(x_up1, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up1, candsX1, dist, f1, g1);
                        }
                    }
                    else
                    {
                        getWpPhip(L[0], f1, g1);
                    }

                    // f2, g2
                    if (interpNeededX_[i][1])
                    {
                        label dn1 = R[0];
                        if (interp2ndCacheX_[i][1].valid)
                        {
                            const auto& cache = interp2ndCacheX_[i][1];
                            f2 = g2 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f2 += cache.w[j]*ww;
                                g2 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_dn1 = getCf(dn1).x();
                            const List<label>& candsX1 = getCandidates(dn1);
                            interpolatePairMQAturgut(x_dn1, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, dn1, candsX1, dist, f2, g2);
                        }
                    }
                    else
                    {
                        getWpPhip(R[0], f2, g2);
                    }

                    h1 = mag(dx(L[0], x0));
                    h2 = mag(dx(R[0], x0));
                    zetaDx_[i][vector::X] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h1 + h2));
                    PhiDx_[i] = -((-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h1 + h2)));
                }
                else if (s == 7 || (s == 14 && Ucur_p[i].x() < -SMALL
                    && (interpNeeded4thX_[i][0] || interpNeeded4thX_[i][1]
                        || interpNeeded4thX_[i][2] || interpNeeded4thX_[i][3])))
                {
                    label down1 = R[0];

                    // ---- Change 4: pointer instead of copy ----
                    const List<label>* down2Ptr = nullptr;
                    {
                        const auto pkIt = globalIdToPackedIdx_.find(down1);
                        if (pkIt != globalIdToPackedIdx_.end())
                            down2Ptr = &downwindNodesX_[*pkIt];
                        else
                        {
                            const auto it = remoteIdToDownwindNodesX_.find(globalIdOfLocalFace_[i]);
                            if (it != remoteIdToDownwindNodesX_.end())
                                down2Ptr = &it();
                            else
                                FatalErrorInFunction
                                    << "Missing second downwind node for globalId "
                                    << globalIdOfLocalFace_[i]
                                    << " in scheme 7. This violates detection guarantee."
                                    << nl << abort(FatalError);
                        }
                    }

                    // s==14 fallback: use 4th-order cache slots (dn1=2, dn2=3)
                    const InterpWeightCache& cDn1X = (s == 14)
                        ? interp4thCacheX_[i][2] : interp2ndCacheX_[i][0];
                    const InterpWeightCache& cDn2X = (s == 14)
                        ? interp4thCacheX_[i][3] : interp2ndCacheX_[i][1];

                    // f1, g1
                    if (interpNeededX_[i][0])
                    {
                        if (cDn1X.valid)
                        {
                            f1 = g1 = 0.0;
                            forAll(cDn1X.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cDn1X.gids[j], ww, pp);
                                f1 += cDn1X.w[j]*ww;
                                g1 += cDn1X.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_dn1 = getCf(down1).x();
                            const List<label>& candsX1 = getCandidates(down1);
                            interpolatePairMQAturgut(x_dn1, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, down1, candsX1, dist, f1, g1);
                        }
                    }
                    else
                    {
                        getWpPhip(down1, f1, g1);
                    }

                    // f2, g2
                    label down2 = (*down2Ptr)[0];
                    if (interpNeededX_[i][1])
                    {
                        if (cDn2X.valid)
                        {
                            f2 = g2 = 0.0;
                            forAll(cDn2X.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cDn2X.gids[j], ww, pp);
                                f2 += cDn2X.w[j]*ww;
                                g2 += cDn2X.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_dn2 = getCf(down2).x();
                            const List<label>& candsX2 = getCandidates(down2);
                            interpolatePairMQAturgut(x_dn2, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, down2, candsX2, dist, f2, g2);
                        }
                    }
                    else
                    {
                        getWpPhip(down2, f2, g2);
                    }

                    h1 = mag(dx(down1, x0)); h2 = mag(dx(down2, x0));
                    zetaDx_[i][vector::X] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h1 - h2));
                    PhiDx_[i] = -((-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h1 - h2)));
                }
                else if (s == 14)
                {
                    label up1   = L[0];
                    label down1 = R[0];

                    // ---- Change 4: pointers ----
                    const List<label>* up2Ptr = nullptr;
                    {
                        const auto pkIt = globalIdToPackedIdx_.find(up1);
                        if (pkIt != globalIdToPackedIdx_.end())
                            up2Ptr = &upwindNodesX_[*pkIt];
                        else
                        {
                            const auto it = remoteIdToUpwindNodesX_.find(globalIdOfLocalFace_[i]);
                            if (it != remoteIdToUpwindNodesX_.end())
                                up2Ptr = &it();
                            else
                                FatalErrorInFunction
                                    << "Missing second upwind node for globalId "
                                    << globalIdOfLocalFace_[i]
                                    << " in scheme 14. This violates detection guarantee."
                                    << nl << abort(FatalError);
                        }
                    }
                    const List<label>* down2Ptr = nullptr;
                    {
                        const auto pkIt = globalIdToPackedIdx_.find(down1);
                        if (pkIt != globalIdToPackedIdx_.end())
                            down2Ptr = &downwindNodesX_[*pkIt];
                        else
                        {
                            const auto it = remoteIdToDownwindNodesX_.find(globalIdOfLocalFace_[i]);
                            if (it != remoteIdToDownwindNodesX_.end())
                                down2Ptr = &it();
                            else
                                FatalErrorInFunction
                                    << "Missing second downwind node for globalId "
                                    << globalIdOfLocalFace_[i]
                                    << " in scheme 14. This violates detection guarantee."
                                    << nl << abort(FatalError);
                        }
                    }

                    // f1, g1 (up1)
                    if (interpNeeded4thX_[i][0])
                    {
                        if (interp4thCacheX_[i][0].valid)
                        {
                            const auto& cache = interp4thCacheX_[i][0];
                            f1 = g1 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f1 += cache.w[j]*ww;
                                g1 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_up1 = getCf(up1).x();
                            const List<label>& candsX1 = getCandidates(up1);
                            interpolatePairMQAturgut(x_up1, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up1, candsX1, dist, f1, g1);
                        }
                    }
                    else
                    {
                        getWpPhip(up1, f1, g1);
                    }

                    // f2, g2 (up2)
                    label up2 = (*up2Ptr)[0];
                    if (interpNeeded4thX_[i][1])
                    {
                        if (interp4thCacheX_[i][1].valid)
                        {
                            const auto& cache = interp4thCacheX_[i][1];
                            f2 = g2 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f2 += cache.w[j]*ww;
                                g2 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_up2 = getCf(up2).x();
                            const List<label>& candsX2 = getCandidates(up2);
                            interpolatePairMQAturgut(x_up2, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up2, candsX2, dist, f2, g2);
                        }
                    }
                    else
                    {
                        getWpPhip(up2, f2, g2);
                    }

                    // f3, g3 (down1)
                    if (interpNeeded4thX_[i][2])
                    {
                        if (interp4thCacheX_[i][2].valid)
                        {
                            const auto& cache = interp4thCacheX_[i][2];
                            f3 = g3 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f3 += cache.w[j]*ww;
                                g3 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_dn1 = getCf(down1).x();
                            const List<label>& candsX3 = getCandidates(down1);
                            interpolatePairMQAturgut(x_dn1, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, down1, candsX3, dist, f3, g3);
                        }
                    }
                    else
                    {
                        getWpPhip(down1, f3, g3);
                    }

                    // f4, g4 (down2)
                    label down2 = (*down2Ptr)[0];
                    if (interpNeeded4thX_[i][3])
                    {
                        if (interp4thCacheX_[i][3].valid)
                        {
                            const auto& cache = interp4thCacheX_[i][3];
                            f4 = g4 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f4 += cache.w[j]*ww;
                                g4 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar x_dn2 = getCf(down2).x();
                            const List<label>& candsX4 = getCandidates(down2);
                            interpolatePairMQAturgut(x_dn2, y0, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, down2, candsX4, dist, f4, g4);
                        }
                    }
                    else
                    {
                        getWpPhip(down2, f4, g4);
                    }

                    const scalar am2 = derivCoeffsX_[i][0], am1 = derivCoeffsX_[i][1];
                    const scalar a0  = derivCoeffsX_[i][2];
                    const scalar ap1 = derivCoeffsX_[i][3], ap2 = derivCoeffsX_[i][4];

                    zetaDx_[i][vector::X] = am2*f2 + am1*f1 + a0*f0 + ap1*f3 + ap2*f4;
                    PhiDx_[i]             = -(am2*g2 + am1*g1 + a0*g0 + ap1*g3 + ap2*g4);
                }
                else if (s == 4 && L.size())
                {
                    h1 = mag(dx(L[0], x0));
                    getWpPhip(L[0], f1, g1);
                    zetaDx_[i][vector::X] = (f0 - f1) / h1;
                    PhiDx_[i] = -((g0 - g1) / h1);
                }
                else if (s == 5 && R.size())
                {
                    h1 = mag(dx(R[0], x0));
                    getWpPhip(R[0], f1, g1);
                    zetaDx_[i][vector::X] = (f1 - f0) / h1;
                    PhiDx_[i] = -((g1 - g0) / h1);
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

                // Inlet face override for X
                if (ownerHasInletFace_[i])
                {
                    const scalar rampperiod = params_.rampperiod;
                    const scalar U0         = params_.U0;
                    const scalar head_ang   = params_.head_ang;
                    const scalar wavelength = params_.wavelength;
                    const scalar steepness  = params_.steepness;
                    const scalar hdepth     = params_.hdepth;

                    const scalar wavenumber(2.0*Foam::constant::mathematical::pi/wavelength);
                    const scalar amp(0.5*steepness*wavelength);
                    const scalar w(Foam::sqrt(9.81*wavenumber*Foam::tanh(wavenumber*hdepth)));

                    scalar ramp_factor = 1.0;
                    if (rampperiod > 0)
                    {
                        const scalar T(2.0*Foam::constant::mathematical::pi/w);
                        const scalar ramp_time(rampperiod*T);
                        ramp_factor = 0.5*(1 - cos(Foam::constant::mathematical::pi*min(1.0, db().time().value()/ramp_time)));
                    }

                    zetaDx_[i][vector::X] = -amp*wavenumber*Foam::sin(wavenumber*x0 - (w + U0*wavenumber*Foam::cos(head_ang))*db().time().value())*ramp_factor;
                    PhiDx_[i] = -(amp*wavenumber*(9.81/w)*Foam::cos(wavenumber*x0 - (w + U0*wavenumber*Foam::cos(head_ang))*db().time().value())*ramp_factor);
                }
            }

            // ============================================================
            // Y-direction
            // ============================================================
            {
                const List<label>& L = upwindNodesY_[i];
                const List<label>& R = downwindNodesY_[i];
                label s = label(schemeCodeY_[i]);

                scalar h1 = 0, h2 = 0, f1 = 0, f2 = 0, g1 = 0, g2 = 0;
                scalar f3 = 0, g3 = 0, f4 = 0, g4 = 0;

                if (s == 6 || (s == 14 && Ucur_p[i].y() > SMALL
                    && (interpNeeded4thY_[i][0] || interpNeeded4thY_[i][1]
                        || interpNeeded4thY_[i][2] || interpNeeded4thY_[i][3])))
                {
                    label up1 = L[0];

                    // ---- Change 4: pointer ----
                    const List<label>* up2Ptr = nullptr;
                    {
                        const auto pkIt = globalIdToPackedIdx_.find(up1);
                        if (pkIt != globalIdToPackedIdx_.end())
                            up2Ptr = &upwindNodesY_[*pkIt];
                        else
                        {
                            const auto it = remoteIdToUpwindNodesY_.find(globalIdOfLocalFace_[i]);
                            if (it != remoteIdToUpwindNodesY_.end())
                                up2Ptr = &it();
                            else
                                FatalErrorInFunction
                                    << "Missing 2nd upwind node for face " << i
                                    << nl << abort(FatalError);
                        }
                    }

                    label up2 = (*up2Ptr)[0];

                    // s==14 fallback: use 4th-order cache slots (up1=0, up2=1)
                    const InterpWeightCache& cUp1Y = (s == 14)
                        ? interp4thCacheY_[i][0] : interp2ndCacheY_[i][0];
                    const InterpWeightCache& cUp2Y = (s == 14)
                        ? interp4thCacheY_[i][1] : interp2ndCacheY_[i][1];

                    // f1, g1
                    if (interpNeededY_[i][0])
                    {
                        if (cUp1Y.valid)
                        {
                            f1 = g1 = 0.0;
                            forAll(cUp1Y.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cUp1Y.gids[j], ww, pp);
                                f1 += cUp1Y.w[j]*ww;
                                g1 += cUp1Y.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_up1 = getCf(up1).y();
                            const List<label>& candsY1 = getCandidates(up1);
                            interpolatePairMQAturgut(x0, y_up1, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up1, candsY1, dist, f1, g1);
                        }
                    }
                    else
                    {
                        getWpPhip(up1, f1, g1);
                    }

                    // f2, g2
                    if (interpNeededY_[i][1])
                    {
                        if (cUp2Y.valid)
                        {
                            f2 = g2 = 0.0;
                            forAll(cUp2Y.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cUp2Y.gids[j], ww, pp);
                                f2 += cUp2Y.w[j]*ww;
                                g2 += cUp2Y.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_up2 = getCf(up2).y();
                            const List<label>& candsY2 = getCandidates(up2);
                            interpolatePairMQAturgut(x0, y_up2, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up2, candsY2, dist, f2, g2);
                        }
                    }
                    else
                    {
                        getWpPhip(up2, f2, g2);
                    }

                    h1 = mag(dy(up1, y0)); h2 = mag(dy(up2, y0));
                    zetaDx_[i][vector::Y] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h2 - h1));
                    PhiDy_[i] = -((-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h2 - h1)));
                }
                else if (s == 1)
                {
                    // f1, g1
                    if (interpNeededY_[i][0])
                    {
                        label up1 = L[0];
                        if (interp2ndCacheY_[i][0].valid)
                        {
                            const auto& cache = interp2ndCacheY_[i][0];
                            f1 = g1 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f1 += cache.w[j]*ww;
                                g1 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_up1 = getCf(up1).y();
                            const List<label>& candsY1 = getCandidates(up1);
                            interpolatePairMQAturgut(x0, y_up1, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up1, candsY1, dist, f1, g1);
                        }
                    }
                    else
                    {
                        getWpPhip(L[0], f1, g1);
                    }

                    // f2, g2
                    if (interpNeededY_[i][1])
                    {
                        label dn1 = R[0];
                        if (interp2ndCacheY_[i][1].valid)
                        {
                            const auto& cache = interp2ndCacheY_[i][1];
                            f2 = g2 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f2 += cache.w[j]*ww;
                                g2 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_dn1 = getCf(dn1).y();
                            const List<label>& candsY2 = getCandidates(dn1);
                            interpolatePairMQAturgut(x0, y_dn1, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, dn1, candsY2, dist, f2, g2);
                        }
                    }
                    else
                    {
                        getWpPhip(R[0], f2, g2);
                    }

                    h1 = mag(dy(L[0], y0));
                    h2 = mag(dy(R[0], y0));
                    zetaDx_[i][vector::Y] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h1 + h2));
                    PhiDy_[i] = -((-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h1 + h2)));
                }
                else if (s == 7 || (s == 14 && Ucur_p[i].y() < -SMALL
                    && (interpNeeded4thY_[i][0] || interpNeeded4thY_[i][1]
                        || interpNeeded4thY_[i][2] || interpNeeded4thY_[i][3])))
                {
                    label down1 = R[0];

                    // ---- Change 4: pointer ----
                    const List<label>* down2Ptr = nullptr;
                    {
                        const auto pkIt = globalIdToPackedIdx_.find(down1);
                        if (pkIt != globalIdToPackedIdx_.end())
                            down2Ptr = &downwindNodesY_[*pkIt];
                        else
                        {
                            const auto it = remoteIdToDownwindNodesY_.find(globalIdOfLocalFace_[i]);
                            if (it != remoteIdToDownwindNodesY_.end())
                                down2Ptr = &it();
                            else
                                FatalErrorInFunction
                                    << "Missing 2nd downwind node for face " << i
                                    << nl << abort(FatalError);
                        }
                    }

                    label down2 = (*down2Ptr)[0];

                    // s==14 fallback: use 4th-order cache slots (dn1=2, dn2=3)
                    const InterpWeightCache& cDn1Y = (s == 14)
                        ? interp4thCacheY_[i][2] : interp2ndCacheY_[i][0];
                    const InterpWeightCache& cDn2Y = (s == 14)
                        ? interp4thCacheY_[i][3] : interp2ndCacheY_[i][1];

                    // f1, g1
                    if (interpNeededY_[i][0])
                    {
                        if (cDn1Y.valid)
                        {
                            f1 = g1 = 0.0;
                            forAll(cDn1Y.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cDn1Y.gids[j], ww, pp);
                                f1 += cDn1Y.w[j]*ww;
                                g1 += cDn1Y.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_dn1 = getCf(down1).y();
                            const List<label>& candsY1 = getCandidates(down1);
                            interpolatePairMQAturgut(x0, y_dn1, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, down1, candsY1, dist, f1, g1);
                        }
                    }
                    else
                    {
                        getWpPhip(down1, f1, g1);
                    }

                    // f2, g2
                    if (interpNeededY_[i][1])
                    {
                        if (cDn2Y.valid)
                        {
                            f2 = g2 = 0.0;
                            forAll(cDn2Y.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cDn2Y.gids[j], ww, pp);
                                f2 += cDn2Y.w[j]*ww;
                                g2 += cDn2Y.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_dn2 = getCf(down2).y();
                            const List<label>& candsY2 = getCandidates(down2);
                            interpolatePairMQAturgut(x0, y_dn2, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, down2, candsY2, dist, f2, g2);
                        }
                    }
                    else
                    {
                        getWpPhip(down2, f2, g2);
                    }

                    h1 = mag(dy(down1, y0)); h2 = mag(dy(down2, y0));
                    zetaDx_[i][vector::Y] = (-sqr(h2)*f1 + sqr(h1)*f2 - (sqr(h1)-sqr(h2))*f0) / (h1*h2*(h1 - h2));
                    PhiDy_[i] = -((-sqr(h2)*g1 + sqr(h1)*g2 - (sqr(h1)-sqr(h2))*g0) / (h1*h2*(h1 - h2)));
                }
                else if (s == 14)
                {
                    label up1   = L[0];
                    label down1 = R[0];

                    // ---- Change 4: pointers ----
                    const List<label>* up2Ptr = nullptr;
                    {
                        const auto pkIt = globalIdToPackedIdx_.find(up1);
                        if (pkIt != globalIdToPackedIdx_.end())
                            up2Ptr = &upwindNodesY_[*pkIt];
                        else
                        {
                            const auto it = remoteIdToUpwindNodesY_.find(globalIdOfLocalFace_[i]);
                            if (it != remoteIdToUpwindNodesY_.end())
                                up2Ptr = &it();
                            else
                                FatalErrorInFunction
                                    << "Missing second upwind Y node for globalId "
                                    << globalIdOfLocalFace_[i]
                                    << nl << abort(FatalError);
                        }
                    }
                    const List<label>* down2Ptr = nullptr;
                    {
                        const auto pkIt = globalIdToPackedIdx_.find(down1);
                        if (pkIt != globalIdToPackedIdx_.end())
                            down2Ptr = &downwindNodesY_[*pkIt];
                        else
                        {
                            const auto it = remoteIdToDownwindNodesY_.find(globalIdOfLocalFace_[i]);
                            if (it != remoteIdToDownwindNodesY_.end())
                                down2Ptr = &it();
                            else
                                FatalErrorInFunction
                                    << "Missing second downwind Y node for globalId "
                                    << globalIdOfLocalFace_[i]
                                    << nl << abort(FatalError);
                        }
                    }

                    // f1, g1 (up1)
                    if (interpNeeded4thY_[i][0])
                    {
                        if (interp4thCacheY_[i][0].valid)
                        {
                            const auto& cache = interp4thCacheY_[i][0];
                            f1 = g1 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f1 += cache.w[j]*ww;
                                g1 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_up1 = getCf(up1).y();
                            const List<label>& candsY1 = getCandidates(up1);
                            interpolatePairMQAturgut(x0, y_up1, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up1, candsY1, dist, f1, g1);
                        }
                    }
                    else
                    {
                        getWpPhip(up1, f1, g1);
                    }

                    // f2, g2 (up2)
                    label up2 = (*up2Ptr)[0];
                    if (interpNeeded4thY_[i][1])
                    {
                        if (interp4thCacheY_[i][1].valid)
                        {
                            const auto& cache = interp4thCacheY_[i][1];
                            f2 = g2 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f2 += cache.w[j]*ww;
                                g2 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_up2 = getCf(up2).y();
                            const List<label>& candsY2 = getCandidates(up2);
                            interpolatePairMQAturgut(x0, y_up2, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, up2, candsY2, dist, f2, g2);
                        }
                    }
                    else
                    {
                        getWpPhip(up2, f2, g2);
                    }

                    // f3, g3 (down1)
                    if (interpNeeded4thY_[i][2])
                    {
                        if (interp4thCacheY_[i][2].valid)
                        {
                            const auto& cache = interp4thCacheY_[i][2];
                            f3 = g3 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f3 += cache.w[j]*ww;
                                g3 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_dn1 = getCf(down1).y();
                            const List<label>& candsY3 = getCandidates(down1);
                            interpolatePairMQAturgut(x0, y_dn1, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, down1, candsY3, dist, f3, g3);
                        }
                    }
                    else
                    {
                        getWpPhip(down1, f3, g3);
                    }

                    // f4, g4 (down2)
                    label down2 = (*down2Ptr)[0];
                    if (interpNeeded4thY_[i][3])
                    {
                        if (interp4thCacheY_[i][3].valid)
                        {
                            const auto& cache = interp4thCacheY_[i][3];
                            f4 = g4 = 0.0;
                            forAll(cache.gids, j)
                            {
                                scalar ww, pp;
                                getWpPhip(cache.gids[j], ww, pp);
                                f4 += cache.w[j]*ww;
                                g4 += cache.w[j]*pp;
                            }
                        }
                        else
                        {
                            const scalar y_dn2 = getCf(down2).y();
                            const List<label>& candsY4 = getCandidates(down2);
                            interpolatePairMQAturgut(x0, y_dn2, Cf, globalIdToPackedIdx_,
                                globalIdToCf_, getWp, getPhip, down2, candsY4, dist, f4, g4);
                        }
                    }
                    else
                    {
                        getWpPhip(down2, f4, g4);
                    }

                    const scalar am2 = derivCoeffsY_[i][0], am1 = derivCoeffsY_[i][1];
                    const scalar a0  = derivCoeffsY_[i][2];
                    const scalar ap1 = derivCoeffsY_[i][3], ap2 = derivCoeffsY_[i][4];

                    zetaDx_[i][vector::Y] = am2*f2 + am1*f1 + a0*f0 + ap1*f3 + ap2*f4;
                    PhiDy_[i]             = -(am2*g2 + am1*g1 + a0*g0 + ap1*g3 + ap2*g4);
                }
                else if (s == 4 && L.size())
                {
                    h1 = mag(dy(L[0], y0));
                    getWpPhip(L[0], f1, g1);
                    zetaDx_[i][vector::Y] = (f0 - f1) / h1;
                    PhiDy_[i] = -((g0 - g1) / h1);
                }
                else if (s == 5 && R.size())
                {
                    h1 = mag(dy(R[0], y0));
                    getWpPhip(R[0], f1, g1);
                    zetaDx_[i][vector::Y] = (f1 - f0) / h1;
                    PhiDy_[i] = -((g1 - g0) / h1);
                }
                else
                {
                    FatalErrorInFunction
                        << "Face " << i << ": Invalid schemeCodeY = " << s
                        << nl << abort(FatalError);
                }

                // Inlet face override for Y
                if (ownerHasInletFace_[i])
                {
                    zetaDx_[i][vector::Y] = 0.0;
                    PhiDy_[i] = 0.0;
                }
            }
        }
    }

    // ---- Body-face GreenGauss loop (separate code path) ----
    for (label i = 0; i < nFaces; ++i)
    {
        if (ownerHasBodyFace_[i] || ownerHasBodyEdge_[i] || ownerHasBodyVertex_[i])
        {
            const label selfGid = globalIdOfLocalFace_[i];
            if (lsMergedCandidates_[i].empty())
            {
                Pout << "lsMergedCandidates_ should not be empty local i: "
                     << i << " GlobalID : " << globalIdOfLocalFace_[i] << nl << endl;
                continue;
            }

            scalar fx = 0, fy = 0, gx = 0, gy = 0;

            GreenGaussUpwindXY
            (
                patch(),
                i,
                selfGid,
                edgeNeighbours_[i],
                edgeNormals2D_[i],
                Ucur_p[i],
                getWp,
                getPhip,
                fx, fy, gx, gy
            );

            zetaDx_[i][vector::X] = fx;
            zetaDx_[i][vector::Y] = fy;
            PhiDx_[i] = -gx;
            PhiDy_[i] = -gy;
        }
    }
}


void Foam::waveCurrentPotential3DFvPatchScalarField::exchangeWpGhostField(const scalarField& zetaZ,const scalarField& Phi)
{
    if (!Pstream::parRun()) return;

    globalIdToWp_.clear();
    globalIdToPhi_.clear();

    // --- Send phase: use pre-built index lists (no Map allocation, no collective) ---
    PstreamBuffers pBuf(Pstream::commsTypes::nonBlocking);

    forAllConstIter(labelHashSet, sendTo, it)
    {
        const label tgt = it.key();
        UOPstream os(tgt, pBuf);

        const labelList& gids   = sendGids_[tgt];
        const labelList& pkidxs = sendPackedIdxs_[tgt];
        os << gids.size();
        forAll(gids, j)
        {
            os << gids[j] << zetaZ[pkidxs[j]] << Phi[pkidxs[j]];
        }
    }

    pBuf.finishedSends();

    // --- Receive phase: use cached recvFrom_ ---
    forAllConstIter(labelHashSet, recvFrom_, it)
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
}


void Foam::waveCurrentPotential3DFvPatchScalarField::computeFirstDerivativeCoeffs5pt
(
    const scalar xm2,
    const scalar xm1,
    const scalar x0,
    const scalar xp1,
    const scalar xp2,
    scalar& am2,
    scalar& am1,
    scalar& a0,
    scalar& ap1,
    scalar& ap2
) const
{
    // Offsets relative to stencil center
    scalar h[5] = {xm2, xm1, x0, xp1, xp2};

    // Augmented matrix [A|b]
    scalar A[5][6];

    for (label j = 0; j < 5; ++j)
    {
        A[0][j] = 1.0;
        A[1][j] = h[j];
        A[2][j] = h[j]*h[j];
        A[3][j] = h[j]*h[j]*h[j];
        A[4][j] = h[j]*h[j]*h[j]*h[j];
    }

    A[0][5] = 0.0;
    A[1][5] = 1.0;
    A[2][5] = 0.0;
    A[3][5] = 0.0;
    A[4][5] = 0.0;

    // Gaussian elimination with partial pivoting
    for (label k = 0; k < 5; ++k)
    {
        label pivotRow = k;
        scalar maxVal = Foam::mag(A[k][k]);

        for (label i = k + 1; i < 5; ++i)
        {
            const scalar val = Foam::mag(A[i][k]);
            if (val > maxVal)
            {
                maxVal = val;
                pivotRow = i;
            }
        }

        if (maxVal < SMALL)
        {
            FatalErrorInFunction
                << "Singular 5-point derivative stencil matrix."
                << nl
                << "Offsets = (" << xm2 << ", " << xm1 << ", " << x0
                << ", " << xp1 << ", " << xp2 << ")"
                << nl << abort(FatalError);
        }

        if (pivotRow != k)
        {
            for (label j = k; j < 6; ++j)
            {
                Swap(A[k][j], A[pivotRow][j]);
            }
        }

        const scalar pivot = A[k][k];
        for (label j = k; j < 6; ++j)
        {
            A[k][j] /= pivot;
        }

        for (label i = k + 1; i < 5; ++i)
        {
            const scalar factor = A[i][k];
            for (label j = k; j < 6; ++j)
            {
                A[i][j] -= factor*A[k][j];
            }
        }
    }

    // Back substitution
    scalar a[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

    for (label i = 4; i >= 0; --i)
    {
        scalar sum = A[i][5];
        for (label j = i + 1; j < 5; ++j)
        {
            sum -= A[i][j]*a[j];
        }
        a[i] = sum;
    }

    am2 = a[0];
    am1 = a[1];
    a0  = a[2];
    ap1 = a[3];
    ap2 = a[4];
}
