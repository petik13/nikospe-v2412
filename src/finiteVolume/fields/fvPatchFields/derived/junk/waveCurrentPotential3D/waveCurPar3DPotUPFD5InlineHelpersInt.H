#ifndef waveCurPar3DPotUPFD5InlineHelpers_H
#define waveCurPar3DPotUPFD5InlineHelpers_H

#include "vector.H"
#include "scalar.H"
#include "label.H"
#include "Map.H"
#include "List.H"
#include <utility>  // For std::pair

namespace Foam
{

//--- Direction-aware inline function to check alignment
inline bool aligned(scalar dx, scalar dy, const word& direction)
{
    if (direction == "x")
        return mag(dx) > VSMALL && (mag(dy) / mag(dx)) < 0.001; //0.0175
    else if (direction == "y")
        return mag(dy) > VSMALL && (mag(dx) / mag(dy)) < 0.001; // 0.08725
    else
        FatalErrorInFunction << "Unknown direction: " << direction << abort(FatalError);

    return false;
}

//--- Inline function to get dx and dy from global ID
inline std::pair<scalar, scalar> dxdy
(
    const vectorField& Cf,
    const Map<label>& globalIdToPackedIdx_,
    const Map<vector>& globalIdToCf_,
    label gid,
    scalar x0,
    scalar y0
)
{
    
	
	if (!globalIdToPackedIdx_.found(gid) && !globalIdToCf_.found(gid) )
	{
		FatalErrorInFunction
			<< "gid " << gid << " not found in globalIdToPackedIdx_ or globalIdToCf_"
			<< abort(FatalError);
	}
	
	const vector& c = globalIdToPackedIdx_.found(gid)
        ? Cf[globalIdToPackedIdx_[gid]]
        : globalIdToCf_[gid];

    return std::make_pair(x0 - c.x(), y0 - c.y());
	
}

template<class PatchField>
inline bool tryUpwind2nd
(
    PatchField& self,
    label i,
    label& scheme,
    const vectorField& Cf,
    const scalarField& x,
    const scalarField& y,
    const List<List<label>>& upwindNodes,
    const HashTable<List<label>, label>& remoteUpwindNodes,
    const word& direction
)
{
    const bool xDir = direction == "x";
    if (!xDir && direction != "y") FatalErrorInFunction << "Unknown direction: " << direction << abort(FatalError);
   
    auto& interpNeeded = xDir ? self.interpNeededX() : self.interpNeededY();

    const List<label>& L = upwindNodes[i];
    if (L.empty()) return false;

    // Flag [0]: multiple direct upwind neighbors
    if (L.size() > 1)
        interpNeeded[i][0] = true;

    
    label up1 = L[0];
    List<label> up2List;
	
	//if( self.globalIdOfLocalFace()[i] ==375) {Pout << "proc " << Pstream::myProcNo() << " self.globalIdOfLocalFace()[i]" << self.globalIdOfLocalFace()[i] << "  up1 : " << up1 << nl;}
	
	
    if (self.globalIdToPackedIdx().found(up1))
    {
        label packedIdx = self.globalIdToPackedIdx()[up1];
        if (packedIdx < upwindNodes.size())
            up2List = upwindNodes[packedIdx];
    }

    if (up2List.empty() && remoteUpwindNodes.found(self.globalIdOfLocalFace()[i]))
    {
        up2List = remoteUpwindNodes[self.globalIdOfLocalFace()[i]];
    }
	
	if (up2List.size() > 1)
        interpNeeded[i][1] = true;
	
	
	
	
    if (up2List.size())
    {
		
        label up2 = up2List[0];
        scheme = 6;
		auto dxdy1 = dxdy(Cf, self.globalIdToPackedIdx(), self.globalIdToCf(), up1, x[i], y[i]);
        
		const vector& c0 = self.globalIdToPackedIdx().found(up1) ? Cf[self.globalIdToPackedIdx()[up1]]  : self.globalIdToCf()[up1];
    
   
		
		auto dxdy2 = dxdy(Cf, self.globalIdToPackedIdx(), self.globalIdToCf(), up2, c0.x(), c0.y());

        if (!aligned(dxdy1.first, dxdy1.second, direction))
            interpNeeded[i][0] = true;

        if (!aligned(dxdy2.first, dxdy2.second, direction))
            interpNeeded[i][1] = true;
		
		
        return true;
    }
    

    return false;
}



//--- Generalized tryCentered2nd()
template<class PatchField>
inline bool tryCentered2nd
(
    PatchField& self,
    label i,
    label& scheme,
    const vectorField& Cf,
    const scalarField& x,
    const scalarField& y,
    const List<List<label>>& upwindNodes,
    const List<List<label>>& downwindNodes,
    const word& direction
)
{
    const List<label>& L = upwindNodes[i];
    const List<label>& R = downwindNodes[i];

    if (L.empty() || R.empty()) return false;

    const bool xDir = direction == "x";
    if (!xDir && direction != "y")
    {
        FatalErrorInFunction << "Unknown direction: " << direction << abort(FatalError);
    }

    auto& interpNeeded = xDir ? self.interpNeededX() : self.interpNeededY();

    if (L.size() > 1) interpNeeded[i][0] = true;
    if (R.size() > 1) interpNeeded[i][1] = true;

    auto dxdyL = dxdy(Cf, self.globalIdToPackedIdx(), self.globalIdToCf(), L[0], x[i], y[i]);
    auto dxdyR = dxdy(Cf, self.globalIdToPackedIdx(), self.globalIdToCf(), R[0], x[i], y[i]);

    scheme = 1;  // default: centered second-order

    bool alignedL = aligned(dxdyL.first, dxdyL.second, direction);
    bool alignedR = aligned(dxdyR.first, dxdyR.second, direction);

    if (!alignedL) interpNeeded[i][0] = true;
    if (!alignedR) interpNeeded[i][1] = true;

    return true;
}



//--- Generalized tryDownwind2nd() with chain-based logic
template<class PatchField>
inline bool tryDownwind2nd
(
    PatchField& self,
    label i,
    label& scheme,
    const vectorField& Cf,
    const scalarField& x,
    const scalarField& y,
    const List<List<label>>& downwindNodes,
    const HashTable<List<label>, label>& remoteDownwindNodes,
    const word& direction
)
{
    const bool xDir = direction == "x";
    if (!xDir && direction != "y")
    {
        FatalErrorInFunction << "Unknown direction: " << direction << abort(FatalError);
    }

    auto& interpNeeded = xDir ? self.interpNeededX() : self.interpNeededY();

    const List<label>& R = downwindNodes[i];
    if (R.empty()) return false;

    // Flag [0]: multiple direct downwind neighbors
    if (R.size() > 1)
        interpNeeded[i][0] = true;

   
    label dn1 = R[0];
    List<label> dn2List;

    // Local lookup for dn2
    if (self.globalIdToPackedIdx().found(dn1))
    {
        label packedIdx = self.globalIdToPackedIdx()[dn1];
        if (packedIdx < downwindNodes.size())
            dn2List = downwindNodes[packedIdx];
    }

    // Remote lookup for dn2
    if (dn2List.empty() && remoteDownwindNodes.found(self.globalIdOfLocalFace()[i]))
    {
        dn2List = remoteDownwindNodes[self.globalIdOfLocalFace()[i]];
    }
	
	if (dn2List.size() > 1)
        interpNeeded[i][1] = true;

    if (dn2List.size())
    {
        label dn2 = dn2List[0];
        scheme = 7; // example scheme for 2nd-order downwind

        auto dxdy1 = dxdy(Cf, self.globalIdToPackedIdx(), self.globalIdToCf(), dn1, x[i], y[i]);
        
		const vector& c0 = self.globalIdToPackedIdx().found(dn1) ? Cf[self.globalIdToPackedIdx()[dn1]]  : self.globalIdToCf()[dn1];
    
		auto dxdy2 = dxdy(Cf, self.globalIdToPackedIdx(), self.globalIdToCf(), dn2,  c0.x(), c0.y());

        if (!aligned(dxdy1.first, dxdy1.second, direction))
            interpNeeded[i][0] = true;

        if (!aligned(dxdy2.first, dxdy2.second, direction))
            interpNeeded[i][1] = true;

        return true;
    }
    

    return false;
}



} // namespace Foam

#endif
