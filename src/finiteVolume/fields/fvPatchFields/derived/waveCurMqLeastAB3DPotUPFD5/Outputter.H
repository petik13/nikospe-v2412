void Foam::waveCurMqLeastAB3DPotUPFD5FvPatchScalarField::outputter()
{
	const vectorField& patchFaceCenters = patch().Cf();
    scalarField xComponents = patchFaceCenters.component(vector::X);
    scalarField yComponents = patchFaceCenters.component(vector::Y);

    
	const fvMesh& mesh = this->patch().boundaryMesh().mesh();
	const label startFace = patch().start();
	
	fileName faceXYFile = "faceXY_" + Foam::name(Pstream::myProcNo()) + ".dat";
	OFstream faceXYOS(
		faceXYFile,
		Foam::IOstreamOption::ASCII,
		Foam::IOstreamOption::UNCOMPRESSED
	); // 

	if (faceXYOS.opened())
	{
		for (label i = 0; i < patch().size(); ++i)
		{
			faceXYOS 
				<< globalIdOfLocalFace_[i] << "    "  // Global face ID
				<< xComponents[i] << "    "          // X component of face center
				<< yComponents[i] << endl;           // Y component of face center
		}
	}
	
	const polyMesh& pMesh = mesh;
    const faceList& allFaces = pMesh.faces();
    const pointField& points = pMesh.points();

   
    fileName vertFile = "faceVertices_" + Foam::name(Pstream::myProcNo()) + ".dat";
    OFstream vertOS(
        vertFile,
        Foam::IOstreamOption::ASCII,
        Foam::IOstreamOption::UNCOMPRESSED,
        Foam::IOstreamOption::NON_APPEND
    );

    if (vertOS.opened())
    {
        for (label i = 0; i < patch().size(); ++i)
        {
            const face& f = allFaces[startFace + i];

            vertOS << globalIdOfLocalFace_[i] << " [ ";
            forAll(f, j)
            {
                const point& vert = points[f[j]];
                vertOS << "(" << vert.x() << " " << vert.y() << " " << vert.z() << ")";
                if (j < f.size()-1) vertOS << ", ";
            }
            vertOS << " ]" << endl;
        }
    }
	
	fileName upwindFile = "upwindX_" + Foam::name(Pstream::myProcNo()) + ".dat";
	OFstream upwindOS(upwindFile,Foam::IOstreamOption::ASCII, Foam::IOstreamOption::UNCOMPRESSED, Foam::IOstreamOption::NON_APPEND); // , Foam::IOstreamOption::NON_APPEND

	if (upwindOS.opened())
	{
		for (label i = 0; i < upwindNodesX_.size(); i++) {
			upwindOS << globalIdOfLocalFace_[i]  << "    " << upwindNodesX_[i].size() << ", [";
			for (label j = 0; j < upwindNodesX_[i].size(); j++) {
				upwindOS << upwindNodesX_[i][j];
				if (j < upwindNodesX_[i].size() - 1) {
					upwindOS << ", ";
				}
			}
			upwindOS << "];" << endl;
		}
	}
	
	fileName downwindFile = "downwindX_" + Foam::name(Pstream::myProcNo()) + ".dat";
	OFstream downwindOS(downwindFile,Foam::IOstreamOption::ASCII,Foam::IOstreamOption::UNCOMPRESSED, Foam::IOstreamOption::NON_APPEND); // 

	if (downwindOS.opened())
	{
		for (label i = 0; i < downwindNodesX_.size(); i++)
		{
			downwindOS << globalIdOfLocalFace_[i] << "    " << downwindNodesX_[i].size() << ", [";

			for (label j = 0; j < downwindNodesX_[i].size(); j++)
			{
				downwindOS << downwindNodesX_[i][j];
				if (j < downwindNodesX_[i].size() - 1)
				{
					downwindOS << ", ";
				}
			}

			downwindOS << "];" << endl;
		}
	}
	
	
	fileName upwindFileY = "upwindY_" + Foam::name(Pstream::myProcNo()) + ".dat";
	OFstream upwindYOS(upwindFileY,Foam::IOstreamOption::ASCII, Foam::IOstreamOption::UNCOMPRESSED, Foam::IOstreamOption::NON_APPEND); // 

	if (upwindYOS.opened())
	{
		for (label i = 0; i < upwindNodesY_.size(); i++) {
			upwindYOS << globalIdOfLocalFace_[i]  << "    " << upwindNodesY_[i].size() << ", [";
			for (label j = 0; j < upwindNodesY_[i].size(); j++) {
				upwindYOS << upwindNodesY_[i][j];
				if (j < upwindNodesY_[i].size() - 1) {
					upwindYOS << ", ";
				}
			}
			upwindYOS << "];" << endl;
		}
	}
	
	fileName downwindFileY = "downwindY_" + Foam::name(Pstream::myProcNo()) + ".dat";
	OFstream downwindYOS(downwindFileY,Foam::IOstreamOption::ASCII,Foam::IOstreamOption::UNCOMPRESSED, Foam::IOstreamOption::NON_APPEND); 

	if (downwindYOS.opened())
	{
		for (label i = 0; i < downwindNodesY_.size(); i++)
		{
			downwindYOS << globalIdOfLocalFace_[i] << "    " << downwindNodesY_[i].size() << ", [";

			for (label j = 0; j < downwindNodesY_[i].size(); j++)
			{
				downwindYOS << downwindNodesY_[i][j];
				if (j < downwindNodesY_[i].size() - 1)
				{
					downwindYOS << ", ";
				}
			}

			downwindYOS << "];" << endl;
		}
	}
	
	
	fileName neighFile = "neighboursGlobal_" + Foam::name(Pstream::myProcNo()) + ".dat";

	OFstream neighOS(
		neighFile,
		Foam::IOstreamOption::ASCII,
		Foam::IOstreamOption::UNCOMPRESSED,
		Foam::IOstreamOption::NON_APPEND
	); // ,
		// Foam::IOstreamOption::NON_APPEND

	if (neighOS.opened())
	{
		// For safety, ensure we don't go out of bounds
		const label n = Foam::min(patch().size(), neighboursGlobal_.size());

		for (label i = 0; i < n; ++i)
		{
			const labelList& nei = neighboursGlobal_[i];

			neighOS << globalIdOfLocalFace_[i] << "    "  // Global face ID
					<< nei.size() << ", [";

			for (label j = 0; j < nei.size(); ++j)
			{
				neighOS << nei[j];
				if (j + 1 < nei.size())
				{
					neighOS << ", ";
				}
			}

			neighOS << "];" << endl;
		}
	}
	
	fileName remNeighFile = "remoteNeighbours_" + Foam::name(Pstream::myProcNo()) + ".dat";

	OFstream remNeighOS(
		remNeighFile,
		Foam::IOstreamOption::ASCII,
		Foam::IOstreamOption::UNCOMPRESSED,
		Foam::IOstreamOption::NON_APPEND
	); // ,
		// IOstreamOption::NON_APPEND

	if (remNeighOS.opened())
	{
		using TableType = HashTable<List<label>, label>;

		forAllConstIter(TableType, remoteNeighbours_, it)
		{
			const label gid = it.key();
			const List<label>& nei = it();

			remNeighOS << gid << "    " << nei.size() << ", [";

			for (label j = 0; j < nei.size(); ++j)
			{
				remNeighOS << nei[j];
				if (j + 1 < nei.size())
				{
					remNeighOS << ", ";
				}
			}

			remNeighOS << "];" << endl;
		}
	}

	fileName remNeighFile2 = "remoteNeighbours2_" + Foam::name(Pstream::myProcNo()) + ".dat";

	OFstream remNeighOS2(
		remNeighFile2,
		Foam::IOstreamOption::ASCII,
		Foam::IOstreamOption::UNCOMPRESSED,
		Foam::IOstreamOption::NON_APPEND
	); // ,
        // IOstreamOption::NON_APPEND

	if (remNeighOS2.opened())
	{
		using TableType = HashTable<List<label>, label>;

		forAllConstIter(TableType, remoteNeighbours2_, it)
		{
			const label gid = it.key();
			const List<label>& nei = it();

			remNeighOS2 << gid << "    " << nei.size() << ", [";

			for (label j = 0; j < nei.size(); ++j)
			{
				remNeighOS2 << nei[j];
				if (j + 1 < nei.size())
				{
					remNeighOS2 << ", ";
				}
			}

			remNeighOS2 << "];" << endl;
		}
	}


	fileName interpFile =
    "interpFlags_" + Foam::name(Pstream::myProcNo()) + ".dat";

	OFstream interpOS(
		interpFile,
		Foam::IOstreamOption::ASCII,
		Foam::IOstreamOption::UNCOMPRESSED,
		Foam::IOstreamOption::NON_APPEND
	);// ,
        // IOstreamOption::NON_APPEND

	if (interpOS.opened())
	{
		const label n = Foam::min(interpNeededX_.size(), interpNeededY_.size());

		for (label i = 0; i < n; ++i)
		{
			const FixedList<bool,2>& flagsX = interpNeededX_[i];
			const FixedList<bool,2>& flagsY = interpNeededY_[i];

			interpOS << globalIdOfLocalFace_[i]
					 << " interpX=[" << flagsX[0] << "," << flagsX[1] << "]"
					 << " interpY=[" << flagsY[0] << "," << flagsY[1] << "]"
					 << endl;
		}
	}
		
	fileName cfFile = "globalIdToCf_" + Foam::name(Pstream::myProcNo()) + ".dat";

	OFstream cfOS(
		cfFile,
		Foam::IOstreamOption::ASCII,
		Foam::IOstreamOption::UNCOMPRESSED,
		Foam::IOstreamOption::NON_APPEND
	); // ,
		// Foam::IOstreamOption::NON_APPEND

	if (cfOS.opened())
	{
		forAllConstIter(Foam::Map<Foam::vector>, globalIdToCf_, iter)
		{
			const label globalId = iter.key();
			const vector& Cf = iter();

			cfOS << globalId << " "
				 << Cf.x() << " " << Cf.y() << " " << Cf.z() << endl;
		}
	}

	using TableType = HashTable<List<label>, label>;

	fileName RemUpwindFile =
		"RemupwindX_" + Foam::name(Pstream::myProcNo()) + ".dat";

	OFstream RemupwindOS(
		RemUpwindFile,
		IOstreamOption::ASCII,
		IOstreamOption::UNCOMPRESSED,	
		IOstreamOption::NON_APPEND
	); // ,
        // IOstreamOption::NON_APPEND

	if (RemupwindOS.opened())
	{
		forAllConstIter(TableType, remoteIdToUpwindNodesX_, it)
		{
			const label globalId = it.key();
			const List<label>& upList = it();

			RemupwindOS << globalId << "    " << upList.size() << ", [";
			forAll(upList, j)
			{
				RemupwindOS << upList[j];
				if (j + 1 < upList.size()) RemupwindOS << ", ";
			}
			RemupwindOS << "];" << endl;
		}
	}

	
// === RemdownwindX ===
{
    fileName file = "RemdownwindX_" + Foam::name(Pstream::myProcNo()) + ".dat";
    OFstream os(file, IOstreamOption::ASCII, IOstreamOption::UNCOMPRESSED, IOstreamOption::NON_APPEND); // , IOstreamOption::NON_APPEND

    if (os.opened())
    {
        forAllConstIter(TableType, remoteIdToDownwindNodesX_, it)
        {
            const label globalId = it.key();
            const List<label>& lst = it();

            os << globalId << "    " << lst.size() << ", [";
            forAll(lst, j)
            {
                os << lst[j];
                if (j + 1 < lst.size()) os << ", ";
            }
            os << "];" << endl;
        }
    }
}

// === RemupwindY ===
{
    fileName file = "RemupwindY_" + Foam::name(Pstream::myProcNo()) + ".dat";
    OFstream os(file, IOstreamOption::ASCII, IOstreamOption::UNCOMPRESSED, IOstreamOption::NON_APPEND); // , IOstreamOption::NON_APPEND

    if (os.opened())
    {
        forAllConstIter(TableType, remoteIdToUpwindNodesY_, it)
        {
            const label globalId = it.key();
            const List<label>& lst = it();

            os << globalId << "    " << lst.size() << ", [";
            forAll(lst, j)
            {
                os << lst[j];
                if (j + 1 < lst.size()) os << ", ";
            }
            os << "];" << endl;
        }
    }
}

// === RemdownwindY ===
{
    fileName file = "RemdownwindY_" + Foam::name(Pstream::myProcNo()) + ".dat";
    OFstream os(file, IOstreamOption::ASCII, IOstreamOption::UNCOMPRESSED, IOstreamOption::NON_APPEND); // , IOstreamOption::NON_APPEND

    if (os.opened())
    {
        forAllConstIter(TableType, remoteIdToDownwindNodesY_, it)
        {
            const label globalId = it.key();
            const List<label>& lst = it();

            os << globalId << "    " << lst.size() << ", [";
            forAll(lst, j)
            {
                os << lst[j];
                if (j + 1 < lst.size()) os << ", ";
            }
            os << "];" << endl;
        }
    }
}
	
// Write neededByProc_ into a dat file
{
    fileName outFile = "neededByProc_" + Foam::name(Pstream::myProcNo()) + ".dat";
    OFstream os(
        outFile,
        IOstreamOption::ASCII,
        IOstreamOption::UNCOMPRESSED,
		IOstreamOption::NON_APPEND
    ); // ,
        // IOstreamOption::NON_APPEND

    if (os.opened())
    {
        os << "Proc " << Pstream::myProcNo() 
           << " neededByProc_ size: " << neededByProc_.size() << nl;

        forAll(neededByProc_, procI)
        {
            os << "  For proc " << procI << ": [";
            const labelHashSet& gids = neededByProc_[procI];
            label count = 0;
            forAllConstIter(labelHashSet, gids, it)
            {
                os << it.key();
                if (++count < gids.size()) os << ", ";
            }
            os << "]" << nl;
        }
    }
    else
    {
        FatalErrorInFunction << "Could not open file " << outFile << abort(FatalError);
    }
}
{
	fileName candsFile = "candidates_" + Foam::name(Pstream::myProcNo()) + ".dat";

	OFstream candsOS(
		candsFile,
		Foam::IOstreamOption::ASCII,
		Foam::IOstreamOption::UNCOMPRESSED,
		Foam::IOstreamOption::NON_APPEND
	); // ,
		// Foam::IOstreamOption::NON_APPEND

	if (candsOS.opened())
	{
		// For safety, ensure we don't go out of bounds
		const label n = Foam::min(patch().size(), lsMergedCandidates_.size());

		for (label i = 0; i < n; ++i)
		{
			const labelList& nei = lsMergedCandidates_[i];

			candsOS << globalIdOfLocalFace_[i] << "    "  // Global face ID
					<< nei.size() << ", [";

			for (label j = 0; j < nei.size(); ++j)
			{
				candsOS << nei[j];
				if (j + 1 < nei.size())
				{
					candsOS << ", ";
				}
			}

			candsOS << "];" << endl;
		}
	}
}	
	
}


//END