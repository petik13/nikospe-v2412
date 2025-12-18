const dictionary& potentialFlow
(
    mesh.solutionDict().subDict("turgutFlow")
);

const int nNonOrthCorr
(
    potentialFlow.getOrDefault<int>("nNonOrthogonalCorrectors", 0)
);
