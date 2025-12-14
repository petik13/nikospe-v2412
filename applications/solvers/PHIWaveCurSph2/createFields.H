Info<< "Reading velocity field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volVectorField Ucur
(
    IOobject
    (
        "Ucur",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedVector(dimVelocity, Zero)
);

// Initialise the velocity internal field to zero
//U = dimensionedVector(U.dimensions(), Zero);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);

phi = fvc::flux(U);

Info<< "Creating field zeta\n" << endl;
volVectorField zeta
(
    IOobject
    (
        "zeta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimLength, Zero)
);

volVectorField zetaDx
(
    IOobject
    (
        "zetaDx",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimLength, Zero)
);




volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar(dimVelocity*dimVelocity, Zero)
);

volScalarField p2
(
    IOobject
    (
        "p2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar(dimVelocity*dimVelocity, Zero)
);

volScalarField p3
(
    IOobject
    (
        "p3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar(dimVelocity*dimVelocity, Zero)
);


Info<< "Constructing velocity potential field Phi\n" << endl;
volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE,
        IOobject::REGISTER
    ),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, Zero)
    
);

Info<< "Constructing velocity potential field Phi\n" << endl;
volScalarField PhiCur
(
    IOobject
    (
        "PhiCur",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, Zero)
    
);


volScalarField PhiDx
(
    IOobject
    (
        "PhiDx",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimVelocity, Zero)
    
);


volScalarField PhiDy
(
    IOobject
    (
        "PhiDy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimVelocity, Zero)
    
);



volScalarField PhiCurDz2
(
    IOobject
    (
        "PhiCurDz2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
    
);




label PhiRefCell = 0;
scalar PhiRefValue = 0;
setRefCell
(
    Phi,
    turgutFlow.dict(),
    PhiRefCell,
    PhiRefValue
);
mesh.setFluxRequired(Phi.name());

IOdictionary waveCurConditions
(
    IOobject
    (
        "waveCurConditions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


#include "createMRF.H"



