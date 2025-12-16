
void Foam::functionObjects::myMeanForce::findFreeSurfaceFaces()
{
    const fvMesh& mesh = mesh_;
    const labelList& faceOwners = mesh.faceOwner();

    ownerHasFreeSurfaceFace_.setSize(patchIDs_.size());

    forAll(patchIDs_, patchIndex)
    {
        const label patchi = patchIDs_[patchIndex];
        const fvPatch& patch = mesh.boundary()[patchi];
        const label startFace = patch.start();
        const label nFaces = patch.size();

        ownerHasFreeSurfaceFace_[patchIndex].setSize(nFaces, false);

        for (label i = 0; i < nFaces; ++i)
        {
            const label globalFace = startFace + i;
            const label ownerCell = faceOwners[globalFace];
            const cell& ownerCellFaces = mesh.cells()[ownerCell];

            forAll(ownerCellFaces, fJ)
            {
                const label faceID = ownerCellFaces[fJ];
                const label neighbPatchI = mesh.boundaryMesh().whichPatch(faceID);

                if (neighbPatchI != -1)
                {
                    const polyPatch& bPatch = mesh.boundaryMesh()[neighbPatchI];
                    const word& neighbPatchName = bPatch.name();

                    if (neighbPatchName == "freeSurface" || neighbPatchName == "upper")
                    {
                        ownerHasFreeSurfaceFace_[patchIndex][i] = true;
                        break;
                    }
                }
            }
        }
    }
	
	
	label nTurgut = 0;

	forAll(patchIDs_, patchIndex)
	{
		nTurgut = 0;  // reset per patch (optional, depending on your goal)

		forAll(ownerHasFreeSurfaceFace_[patchIndex], i)
		{
			if (ownerHasFreeSurfaceFace_[patchIndex][i])
			{
				nTurgut++;   // ✅ increment by 1
			}
		}

		Pout << "Proc " << Pstream::myProcNo()
			 << " has " << nTurgut
			 << " freeSurface faces on patch "
			 << mesh_.boundary()[patchIDs_[patchIndex]].name()
			 << endl;
	}
	
	
	/*---------------------------------------------------------------------------*\
	 
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


}
