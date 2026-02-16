/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
InterpolationsHelpers by TURGUT
f1 = interpolateAtXY( x0, y0, Cf,globalIdToPackedIdx_, globalIdToCf_, getWp, upwindNodesX_[i]);
\*---------------------------------------------------------------------------*/
// Interpolates value of a scalar field at (x_interp, y_interp)
// using globalId-based neighbors (local or remote).
#ifndef InterpolationsHelpers_H
#define InterpolationsHelpers_H

#include <cmath>
#include "scalar.H"
#include "vector.H"
#include "label.H"
#include "Map.H"
#include "List.H"

namespace Foam
{
    template<class GetValueFunc>
    inline Foam::scalar interpolateAtXY(
        const Foam::scalar x_interp,
        const Foam::scalar y_interp,
        const Foam::vectorField& Cf,
        const Foam::Map<Foam::label>& globalIdToPackedIdx,
        const Foam::Map<Foam::vector>& globalIdToCf,
        const GetValueFunc& getValue,
        const Foam::List<Foam::label>& candidates
    )
    {
        Foam::scalar totalWeight = 0.0;
        Foam::scalar interpolatedValue = 0.0;

        forAll(candidates, j)
        {
            const Foam::label gid = candidates[j];

            Foam::vector c = Foam::vector::zero;
            if (globalIdToPackedIdx.found(gid))
                c = Cf[globalIdToPackedIdx[gid]];
            else if (globalIdToCf.found(gid))
                c = globalIdToCf[gid];
            else
                FatalErrorInFunction << "interpolateAtXY: gid " << gid << " not found in local or ghost maps." << abort(FatalError);
                    
                    
                    

            Foam::scalar dx = x_interp - c.x();
            Foam::scalar dy = y_interp - c.y();
            Foam::scalar r  = Foam::sqrt(dx*dx + dy*dy) + VSMALL;

            Foam::scalar w   = 1.0 / r;
            Foam::scalar val = getValue(gid);

            interpolatedValue += w * val;
            totalWeight       += w;
        }

        return totalWeight > VSMALL ? interpolatedValue / totalWeight : 0.0;
    }

template<class GetValueFunc>
inline Foam::scalar interpolateMQAturgut(
    const Foam::scalar x_interp,
    const Foam::scalar y_interp,
    const Foam::vectorField& Cf,
    const Foam::Map<Foam::label>& globalIdToPackedIdx,
    const Foam::Map<Foam::vector>& globalIdToCf,
    const GetValueFunc& getValue,
	const Foam::label selfGid,  
    const Foam::List<Foam::label>& candidates,
	const Foam::scalar epsilon
)
{
    
    
	
	if (candidates.empty()) FatalErrorInFunction << "interpolateMQAtXY: No candidates provided" << abort(FatalError);
		
	if (candidates.size() == 1) {
    const label only = candidates[0];
	scalar val = getValue(only);
    Pout << "DEBUG candidate-only branch: selfGid=" << selfGid << " gid=" << only << nl;
    //Pout << "value return as =" << getValue(only) << nl;
	
    if (!globalIdToPackedIdx.found(only) && !globalIdToCf.found(only)) {
        FatalErrorInFunction << "Candidate " << only << " not found in Cf maps" << abort(FatalError);
    }
   
    return val;
	}
	


    // ---- Collect coordinates and values ----
    DynamicList<Foam::vector> points;
    DynamicList<Foam::scalar> values;
	
	// Always add self point first
    {
        Foam::vector c = Foam::vector::zero;
        if (globalIdToPackedIdx.found(selfGid))
            c = Cf[globalIdToPackedIdx[selfGid]];
        else if (globalIdToCf.found(selfGid))
            c = globalIdToCf[selfGid];
        else
            FatalErrorInFunction << "interpolateMQAtXY: selfGid " << selfGid << " not found in local or ghost maps."  << abort(FatalError);

        points.append(Foam::vector(c.x(), c.y(), 0.0));
        values.append(getValue(selfGid));
    }
	
    forAll(candidates, j)
    {
        const Foam::label gid = candidates[j];
		
		//if (gid == 113) Pout << "interpolateMQAtXY: proc " << Pstream::myProcNo() << " selfGid " << selfGid  << " has candidate 113 in the list: " << candidates << nl;
	
        Foam::vector c = Foam::vector::zero;
        if (globalIdToPackedIdx.found(gid))
            c = Cf[globalIdToPackedIdx[gid]];
        else if (globalIdToCf.found(gid))
            c = globalIdToCf[gid];
        else
            FatalErrorInFunction
                << "interpolateMQAtXY: gid " << gid
                << " not found in local or ghost maps."
                << abort(FatalError);
		// ---- Check if still vector::zero ----
		if (c == Foam::vector::zero)
		{
			FatalErrorInFunction << "interpolateMQAtXY: gid " << gid  << " has Cf == vector::zero (uninitialized coordinate)." << abort(FatalError);
	 
		}
        points.append(Foam::vector(c.x(), c.y(), 0.0));
        values.append(getValue(gid));
    }
	 


    const Foam::label N = points.size();
    Foam::scalarSquareMatrix A(N, 0.0);
    Foam::scalarField f(values);

    //const Foam::scalar epsilon = 1.0;
	
	
	
    // ---- Build MQ RBF matrix ----
    for (label i = 0; i < N; ++i)
    {
        for (label j = 0; j < N; ++j)
        {
            Foam::scalar r = mag(points[i] - points[j]);
            A(i, j) = Foam::sqrt(r*r + epsilon*epsilon);
        }
    }

    // ---- Solve A * w = f (via pseudo-inverse) ----
    Foam::scalarRectangularMatrix Ainv = SVDinv(A, 1e-12);
    Foam::scalarField w = Ainv * f;

    // ---- Evaluate interpolant ----
    Foam::scalar interpolatedValue = 0.0;
    Foam::vector pInterp(x_interp, y_interp, 0.0);

    for (label i = 0; i < N; ++i)
    {
        Foam::scalar r = mag(pInterp - points[i]);
        interpolatedValue += w[i] * Foam::sqrt(r*r + epsilon*epsilon);
    }

    return interpolatedValue;
}

template<class GetValueFunc1, class GetValueFunc2>
inline void interpolatePairMQAturgut(
    const scalar x_interp,
    const scalar y_interp,
    const vectorField& Cf,
    const Map<label>& globalIdToPackedIdx,
    const Map<vector>& globalIdToCf,
    const GetValueFunc1& getValue1,
    const GetValueFunc2& getValue2,
    const label selfGid,
    const List<label>& candidates,
    const scalar epsilon,
    scalar& out1,
    scalar& out2
)
{
    if (candidates.empty())
        FatalErrorInFunction << "interpolatePairMQAturgut: No candidates provided" << abort(FatalError);

    // Special case: one candidate → direct evaluation
    if (candidates.size() == 1)
    {
        const label only = candidates[0];

        if (!globalIdToPackedIdx.found(only) && !globalIdToCf.found(only))
        {
            FatalErrorInFunction
                << "Candidate " << only << " not found in Cf maps" << abort(FatalError);
        }

        out1 = getValue1(only);
        out2 = getValue2(only);
        return;
    }

    // ---- Collect coordinates and values ----
    DynamicList<vector> points;
    DynamicList<scalar> vals1, vals2;

    auto collectValue = [&](label gid)
    {
        vector c = vector::zero;
        if (globalIdToPackedIdx.found(gid))
            c = Cf[globalIdToPackedIdx[gid]];
        else if (globalIdToCf.found(gid))
            c = globalIdToCf[gid];
        else
            FatalErrorInFunction
                << "interpolatePairMQAturgut: gid " << gid
                << " not found in local or ghost maps." << abort(FatalError);

        if (c == vector::zero)
            FatalErrorInFunction
                << "interpolatePairMQAturgut: gid " << gid
                << " has Cf == vector::zero (uninitialized coordinate)." << abort(FatalError);

        points.append(vector(c.x(), c.y(), 0.0));
        vals1.append(getValue1(gid));
        vals2.append(getValue2(gid));
    };

    collectValue(selfGid);
    forAll(candidates, j) collectValue(candidates[j]);

    const label N = points.size();
    scalarSquareMatrix A(N, 0.0);

    // ---- Build MQ RBF matrix ----
    for (label i = 0; i < N; ++i)
    {
        for (label j = 0; j < N; ++j)
        {
            scalar r = mag(points[i] - points[j]);
            A(i, j) = sqrt(r*r + epsilon*epsilon);
        }
    }

    // Reuse one SVD inverse
    scalarRectangularMatrix Ainv = SVDinv(A, 1e-12);

    scalarField w1 = Ainv * scalarField(vals1);
    scalarField w2 = Ainv * scalarField(vals2);

    vector pInterp(x_interp, y_interp, 0.0);

    out1 = out2 = 0.0;
    for (label i = 0; i < N; ++i)
    {
        scalar r = mag(pInterp - points[i]);
        scalar mqrbf = sqrt(r*r + epsilon*epsilon);
        out1 += w1[i] * mqrbf;
        out2 += w2[i] * mqrbf;
    }
}

template<class GetValueFunc1, class GetValueFunc2>
inline void DerivativeMqRBFxy(
    const scalar epsilon,
	const vectorField& Cf,
    const Map<label>& globalIdToPackedIdx,
    const Map<vector>& globalIdToCf,
    const GetValueFunc1& getValue1,
    const GetValueFunc2& getValue2,
    const label selfGid,
    const List<label>& candidates,
    scalar& out1x,
	scalar& out1y,
    scalar& out2x,
	scalar& out2y
)
{
    if (candidates.empty())
        FatalErrorInFunction
            << "DerivativeMqRBF: No candidates provided" << abort(FatalError);

    if (candidates.size() < 5)
        FatalErrorInFunction
            << "DerivativeMqRBF: Not enough candidates at " << selfGid
            << " for DerivativeMqRBF candidates: " << candidates
            << abort(FatalError);

// ---- Collect coordinates and values ----
    DynamicList<vector> points;
    DynamicList<scalar> vals1, vals2;

    auto collectValue = [&](label gid)
    {
        vector c = vector::zero;
        if (globalIdToPackedIdx.found(gid))
            c = Cf[globalIdToPackedIdx[gid]];
        else if (globalIdToCf.found(gid))
            c = globalIdToCf[gid];
        else
            FatalErrorInFunction
                << "DerivativeMqRBF: gid " << gid
                << " not found in local or ghost maps." << abort(FatalError);

        if (c == vector::zero)
            FatalErrorInFunction
                << "DerivativeMqRBF: gid " << gid
                << " has Cf == vector::zero (uninitialized coordinate)." << abort(FatalError);

        points.append(vector(c.x(), c.y(), 0.0));
        vals1.append(getValue1(gid));
        vals2.append(getValue2(gid));
    };

    collectValue(selfGid);
    forAll(candidates, j) collectValue(candidates[j]);
	const vector& c0 = Cf[globalIdToPackedIdx[selfGid]];


    const label N = points.size();
    scalarSquareMatrix A(N, 0.0);
	


    // ---- Build MQ RBF matrix ----
    for (label i = 0; i < N; ++i)
    {
        for (label j = 0; j < N; ++j)
        {
            scalar r = mag(points[i] - points[j]);
            A(i, j) = sqrt(r*r + epsilon*epsilon);
        }
    }

    // Reuse one SVD inverse
    scalarRectangularMatrix Ainv = SVDinv(A, 1e-12);

    scalarField w1 = Ainv * scalarField(vals1);
    scalarField w2 = Ainv * scalarField(vals2);

    vector pInterp( c0.x(),  c0.y(), 0.0);

    out1x = out2x = out1y = out2y = 0.0;
    

		for (label i = 0; i < N; ++i)
		{
			scalar r = mag(pInterp - points[i]);
			scalar mqrbfx = (c0.x()-points[i].x())/sqrt(r*r + epsilon*epsilon);
			scalar mqrbfy = (c0.y()-points[i].y())/sqrt(r*r + epsilon*epsilon);
			
			
			out1x += w1[i] * mqrbfx;
			out2x += w2[i] * mqrbfx;
			out1y += w1[i] * mqrbfy;
			out2y += w2[i] * mqrbfy;
		}
		
	
   
}

// New method to handle few neighbor cases with linear least squares
template<class GetValueFunc1, class GetValueFunc2>
inline void DerivativeLSRobustTaylor(
    const vectorField& Cf,
    const Map<label>& globalIdToPackedIdx,
    const Map<vector>& globalIdToCf,
    const GetValueFunc1& getValue1,
    const GetValueFunc2& getValue2,
    const label selfGid,
    const List<label>& candidates,
    scalar& out1x,
    scalar& out1y,
    scalar& out2x,
    scalar& out2y,
    const scalar svdTol = 1e-12,      // SVD tolerance
    const label  minQuadPts = 8,      // require more than 5 for quadratic robustness
    const scalar weightP = 2.0,       // distance weight exponent
    const scalar epsR = 1e-12,        // avoid division by zero in weights
    const scalar maxGradFactor = 5.0  // gradient clamp factor (see below)
)
{
    // ----------------------------
    // Helpers: coordinate lookup
    // ----------------------------
    auto getCoord = [&](label gid) -> vector
    {
        if (globalIdToPackedIdx.found(gid))
        {
            return Cf[globalIdToPackedIdx[gid]];
        }
        if (globalIdToCf.found(gid))
        {
            return globalIdToCf[gid];
        }

        FatalErrorInFunction
            << "DerivativeLSRobustTaylor: gid " << gid
            << " not found in local or ghost Cf maps." << nl
            << abort(FatalError);
        return vector::zero;
    };

    // ----------------------------
    // Gather points/values
    // ----------------------------
    if (candidates.empty())
    {
        FatalErrorInFunction
            << "DerivativeLSRobustTaylor: No candidates provided for selfGid "
            << selfGid << nl << abort(FatalError);
    }

    const vector c0 = getCoord(selfGid);
    const scalar f0 = getValue1(selfGid);
    const scalar g0 = getValue2(selfGid);

    // Remove duplicates and self (defensive)
    DynamicList<label> uniq;
    uniq.reserve(candidates.size());
    forAll(candidates, j)
    {
        const label gid = candidates[j];
        if (gid == selfGid) continue;
        if (!uniq.found(gid)) uniq.append(gid);
    }

    if (uniq.empty())
    {
        // No usable candidates: safest output
        out1x = out1y = out2x = out2y = 0.0;
        return;
    }

    const label N = uniq.size();

    // Compute a local length scale h from candidate distances
    // (robust scaling dramatically improves conditioning)
    scalar sumR = 0.0;
    label  nR   = 0;
    forAll(uniq, j)
    {
        const vector cj = getCoord(uniq[j]);
        const scalar dx = cj.x() - c0.x();
        const scalar dy = cj.y() - c0.y();
        const scalar r  = Foam::sqrt(dx*dx + dy*dy);
        if (r > VSMALL)
        {
            sumR += r;
            ++nR;
        }
    }
    scalar h = (nR > 0 ? sumR / nR : 1.0);
    if (h < VSMALL) h = 1.0;

    // ----------------------------
    // Weighted LS builder/solver
    // ----------------------------
    auto solveLS = [&](const label nUnknowns,
                       const bool quadratic,
                       scalar& fx, scalar& fy,
                       scalar& gx, scalar& gy) -> bool
    {
        // A: N x nUnknowns, rhs: N
        scalarRectangularMatrix A(N, nUnknowns, 0.0);
        scalarField rhs1(N, 0.0), rhs2(N, 0.0);

        for (label i = 0; i < N; ++i)
        {
            const label gid = uniq[i];
            const vector ci = getCoord(gid);

            // scaled coordinates
            const scalar dx = (ci.x() - c0.x())/h;
            const scalar dy = (ci.y() - c0.y())/h;

            const scalar r2 = dx*dx + dy*dy;
            const scalar w  = 1.0 / (Foam::pow(Foam::sqrt(r2) + epsR, weightP)); // 1/(r^p)
            const scalar sw = Foam::sqrt(w);

            // linear terms
            A(i,0) = sw*dx;
            A(i,1) = sw*dy;

            if (quadratic)
            {
                // quadratic Taylor terms (still scaled)
                A(i,2) = sw*(0.5*dx*dx);
                A(i,3) = sw*(dx*dy);
                A(i,4) = sw*(0.5*dy*dy);
            }

            rhs1[i] = sw*(getValue1(gid) - f0);
            rhs2[i] = sw*(getValue2(gid) - g0);
        }

        // Solve directly with SVD on A (NOT normal equations)
        // coeff contains derivatives w.r.t scaled coordinates:
        // fx_scaled = d f / d(dx) = h * df/dx  => df/dx = fx_scaled / h
        scalarRectangularMatrix Ainv = SVDinv(A, svdTol);

        // If A is rank-deficient, SVDinv will still return something,
        // but we detect badness via residual + sanity checks.
        const scalarField c1 = Ainv * rhs1;
        const scalarField c2 = Ainv * rhs2;

        if (c1.size() < 2 || c2.size() < 2)
            return false;

        // unscale gradients back to physical derivatives
        fx = c1[0]/h;
        fy = c1[1]/h;
        gx = c2[0]/h;
        gy = c2[1]/h;

        // Residual check (if fit is nonsense, fall back)
        // Compute ||A c - rhs|| / ||rhs|| for both fields
        auto relResidual = [&](const scalarField& c, const scalarField& rhs) -> scalar
        {
            scalar num = 0.0;
            scalar den = 0.0;
            for (label i = 0; i < N; ++i)
            {
                scalar pred = 0.0;
                for (label k = 0; k < nUnknowns; ++k)
                    pred += A(i,k)*c[k];
                const scalar r = pred - rhs[i];
                num += r*r;
                den += rhs[i]*rhs[i];
            }
            if (den < SMALL) return Foam::sqrt(num); // absolute residual if rhs ~ 0
            return Foam::sqrt(num/den);
        };

        const scalar rr1 = relResidual(c1, rhs1);
        const scalar rr2 = relResidual(c2, rhs2);

        // Hard sanity: NaN/Inf or absurd numbers
        auto badScalar = [&](scalar v) -> bool
        {
            return !std::isfinite(v) || mag(v) > GREAT/1000;
        };

        if (badScalar(fx) || badScalar(fy) || badScalar(gx) || badScalar(gy))
            return false;

        // Residual thresholds: quadratic is more fragile → require better fit
        const scalar rrLimit = quadratic ? 0.5 : 2.0;
        if (rr1 > rrLimit || rr2 > rrLimit)
            return false;

        return true;
    };

    // ----------------------------
    // Clamp helper (prevents blow-ups)
    // ----------------------------
    auto clampGradients = [&](scalar& fx, scalar& fy, scalar& gx, scalar& gy)
    {
        // Estimate a safe gradient magnitude from neighbor differences
        // (based on max finite-difference slopes)
        scalar fxMax = 0.0, fyMax = 0.0, gxMax = 0.0, gyMax = 0.0;

        forAll(uniq, j)
        {
            const label gid = uniq[j];
            const vector cj = getCoord(gid);

            const scalar dx = cj.x() - c0.x();
            const scalar dy = cj.y() - c0.y();

            const scalar df = getValue1(gid) - f0;
            const scalar dg = getValue2(gid) - g0;

            if (mag(dx) > VSMALL)
            {
                fxMax = Foam::max(fxMax, mag(df/dx));
                gxMax = Foam::max(gxMax, mag(dg/dx));
            }
            if (mag(dy) > VSMALL)
            {
                fyMax = Foam::max(fyMax, mag(df/dy));
                gyMax = Foam::max(gyMax, mag(dg/dy));
            }
        }

        // If we couldn't estimate anything (all dx or all dy ~ 0), do not clamp hard
        const scalar fxLim = (fxMax > VSMALL ? maxGradFactor*fxMax : GREAT);
        const scalar fyLim = (fyMax > VSMALL ? maxGradFactor*fyMax : GREAT);
        const scalar gxLim = (gxMax > VSMALL ? maxGradFactor*gxMax : GREAT);
        const scalar gyLim = (gyMax > VSMALL ? maxGradFactor*gyMax : GREAT);

        fx = Foam::max(-fxLim, Foam::min(fxLim, fx));
        fy = Foam::max(-fyLim, Foam::min(fyLim, fy));
        gx = Foam::max(-gxLim, Foam::min(gxLim, gx));
        gy = Foam::max(-gyLim, Foam::min(gyLim, gy));
    };

    // ----------------------------
    // Strategy:
    // 1) Try quadratic WLS only if we have enough points
    // 2) Otherwise do linear WLS
    // 3) If still fails, do a very safe 1D fallback
    // 4) Always clamp
    // ----------------------------
    bool ok = false;
    scalar fx = 0.0, fy = 0.0, gx = 0.0, gy = 0.0;

    if (N >= minQuadPts)
    {
        ok = solveLS(/*nUnknowns*/5, /*quadratic*/true, fx, fy, gx, gy);
    }

    if (!ok)
    {
        ok = solveLS(/*nUnknowns*/2, /*quadratic*/false, fx, fy, gx, gy);
    }

    if (!ok)
    {
        // Very safe last resort:
        // pick the closest candidate and do directional finite difference
        label best = uniq[0];
        scalar bestR2 = GREAT;

        forAll(uniq, j)
        {
            const vector cj = getCoord(uniq[j]);
            const scalar dx = cj.x() - c0.x();
            const scalar dy = cj.y() - c0.y();
            const scalar r2 = dx*dx + dy*dy;
            if (r2 < bestR2)
            {
                bestR2 = r2;
                best = uniq[j];
            }
        }

        const vector cb = getCoord(best);
        const scalar dx = cb.x() - c0.x();
        const scalar dy = cb.y() - c0.y();
        const scalar r2 = dx*dx + dy*dy;

        if (r2 > VSMALL)
        {
            const scalar df = getValue1(best) - f0;
            const scalar dg = getValue2(best) - g0;

            // Project onto x/y with a simple least-norm split
            fx = df * dx / r2;
            fy = df * dy / r2;
            gx = dg * dx / r2;
            gy = dg * dy / r2;
        }
        else
        {
            fx = fy = gx = gy = 0.0;
        }
    }

    // Always clamp to prevent blow-ups from any weird stencil
    clampGradients(fx, fy, gx, gy);

    out1x = fx;
    out1y = fy;
    out2x = gx;
    out2y = gy;
}

inline void diffuseBodyGradientsCandidates
(
    scalarField& G,
    const boolList& ownerHasBodyFace,
    const List<labelList>& lsMergedCandidates,
    const Map<label>& globalIdToPackedIdx,
    const scalar mu,
    const label nSweeps
)
{
    scalarField Gold(G);

    for (label sweep = 0; sweep < nSweeps; ++sweep)
    {
        Gold = G;

        forAll(G, i)
        {
            if (!ownerHasBodyFace[i]) continue;

            const labelList& cands = lsMergedCandidates[i];
            if (cands.empty()) continue;

            scalar sum = 0.0;
            label n = 0;

            forAll(cands, k)
            {
                const label gid = cands[k];

                if (globalIdToPackedIdx.found(gid))
                {
                    const label j = globalIdToPackedIdx[gid];

                    sum += Gold[j];
                    ++n;
                }
            }

            if (n > 0)
            {
                const scalar avg = sum / n;
                G[i] = (1.0 - mu)*Gold[i] + mu*avg;
            }
        }
    }
}

// ----------------------------------------------------------------------------

// --------- Green Gauss gradient reconstruction implementation --------------------

template<class GetW, class GetP>
static inline void GreenGaussUpwindXY
(
    const fvPatch& patch,
    const label faceI,
    const label selfGid,
    const labelList& neiGids,
    const List<vector>& neiEdgeNormals2D,
    const vector& Uadv,
    const GetW& getWp,
    const GetP& getPhip,
    scalar& fx, scalar& fy,
    scalar& gx, scalar& gy
)
{
    const scalar A = mag(patch.Sf()[faceI]);
    if (A < VSMALL)
    {
        fx = fy = gx = gy = 0;
        return;
    }

    const scalar wSelf = getWp(selfGid);
    const scalar pSelf = getPhip(selfGid);

    vector gradW(0,0,0);
    vector gradP(0,0,0);

    bool usedSelfAtLeastOnce = false;

    // ---------- pass 1: upwind ----------
    forAll(neiGids, k)
    {
        const label gidN = neiGids[k];
        const vector& nk = neiEdgeNormals2D[k];

        const scalar flux = (Uadv & nk);

        if (flux < 0)
        {
            usedSelfAtLeastOnce = true;
            gradW += getWp(gidN)   * nk;
            gradP += getPhip(gidN) * nk;
        }
        else
        {
            gradW += wSelf * nk;
            gradP += pSelf * nk;
        }
    }

    // ---------- fallback: if ALL edges were inflow, switch to downwind everywhere ----------
    if (!usedSelfAtLeastOnce)
    {
        gradW = vector::zero;
        gradP = vector::zero;

        forAll(neiGids, k)
        {
            const label gidN = neiGids[k];
            const vector& nk = neiEdgeNormals2D[k];

            // downwind everywhere: always neighbour
            gradW += getWp(gidN)   * nk;
            gradP += getPhip(gidN) * nk;
        }
    }

    gradW /= A;
    gradP /= A;

    fx = gradW.x();
    fy = gradW.y();
    gx = gradP.x();
    gy = gradP.y();
}




// --------- Least squares Taylor gradient reconstruction implementation ----------------



template<class GetValueFunc1, class GetValueFunc2>
inline void DerivativeLSTaylor(
    const vectorField& Cf,
    const Map<label>& globalIdToPackedIdx,
    const Map<vector>& globalIdToCf,
    const GetValueFunc1& getValue1,
    const GetValueFunc2& getValue2,
    const label selfGid,
    const List<label>& candidates,
    scalar& out1x,
	scalar& out1y,
    scalar& out2x,
	scalar& out2y
)
{
    if (candidates.empty())
        FatalErrorInFunction
            << "DerivativeLSTaylor: No candidates provided" << abort(FatalError);

    // if (candidates.size() < 5)
    //     FatalErrorInFunction
    //         << "DerivativeLSTaylor: Not enough candidates at " << selfGid
    //         << " for quadratic LS (need ≥5), candidates: " << candidates
    //         << abort(FatalError);

    // ---- Collect coordinates and values ----
    DynamicList<vector> points;
    DynamicList<scalar> vals1, vals2;

    auto collectValue = [&](label gid)
    {
        vector c = vector::zero;
        if (globalIdToPackedIdx.found(gid)) // 
            c = Cf[globalIdToPackedIdx[gid]];
        else if (globalIdToCf.found(gid))
            c = globalIdToCf[gid];
        else
            FatalErrorInFunction
                << "DerivativeLSTaylor: gid " << gid
                << " not found in local or ghost maps." << abort(FatalError);

        if (c == vector::zero)
            FatalErrorInFunction
                << "DerivativeLSTaylor: gid " << gid
                << " has Cf == vector::zero (uninitialized coordinate)." << abort(FatalError);

        points.append(vector(c.x(), c.y(), 0.0));
        vals1.append(getValue1(gid));
        vals2.append(getValue2(gid));
    };

    forAll(candidates, j) collectValue(candidates[j]);

    const label N = points.size();
    const vector& c0 = Cf[globalIdToPackedIdx[selfGid]];
    const scalar f0 = getValue1(selfGid);
    const scalar g0 = getValue2(selfGid);

    // ---- Build design matrix for quadratic LS ----
    scalarRectangularMatrix A(N, 5, 0.0);
    scalarField rhs1(N, 0.0), rhs2(N, 0.0);

    for (label i = 0; i < N; ++i)
    {
        scalar dx = points[i].x() - c0.x();
        scalar dy = points[i].y() - c0.y();

        A(i,0) = dx;
        A(i,1) = dy;
        A(i,2) = 0.5*dx*dx;
        A(i,3) = dx*dy;
        A(i,4) = 0.5*dy*dy;

        rhs1[i] = vals1[i] - f0;
        rhs2[i] = vals2[i] - g0;
    }

    // ---- Solve least squares via SVD ----
    

   // ---- Compute AᵀA and Aᵀb for both fields ----
	scalarSquareMatrix ATA(5, 0.0);
	scalarField ATb1(5, 0.0);
	scalarField ATb2(5, 0.0);

	for (label i = 0; i < N; ++i)
	{
		for (label j = 0; j < 5; ++j)
		{
			for (label k = 0; k < 5; ++k)
				ATA[j][k] += A(i,j) * A(i,k);

			ATb1[j] += A(i,j) * rhs1[i];
			ATb2[j] += A(i,j) * rhs2[i];
		}
	}

	// ---- Solve (AᵀA) c = Aᵀb via LU decomposition ----
	
	scalarRectangularMatrix Ainv = SVDinv(ATA, 1e-12);
	scalarField coeff1 = Ainv * ATb1;  // [fx, fy, fxx, fxy, fyy]
    scalarField coeff2 = Ainv * ATb2;
	
	//scalarRectangularMatrix Ainv = SVDinv(A, 1e-12);
	//scalarField coeff1 = Ainv * rhs1;  // [fx, fy, fxx, fxy, fyy]
    //scalarField coeff2 = Ainv * rhs2;
	


  
        out1x = coeff1[0];  // fx
        out2x = coeff2[0];  // gx
  
        out1y = coeff1[1];  // fy
        out2y = coeff2[1];  // gy
   


}

}
#endif
