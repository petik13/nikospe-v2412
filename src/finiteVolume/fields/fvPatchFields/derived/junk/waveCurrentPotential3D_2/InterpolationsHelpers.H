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

    if (candidates.size() < 5)
        FatalErrorInFunction
            << "DerivativeLSTaylor: Not enough candidates at " << selfGid
            << " for quadratic LS (need ≥5), candidates: " << candidates
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
