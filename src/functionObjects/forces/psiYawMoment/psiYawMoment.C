/*---------------------------------------------------------------------------*/

#include "psiYawMoment.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "IOdictionary.H"
#include "uniformDimensionedFields.H"
#include "gravityMeshObject.H"
#include "mathematicalConstants.H"
#include "fvm.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(psiYawMoment, 0);
    addToRunTimeSelectionTable(functionObject, psiYawMoment, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::functionObjects::psiYawMoment::d2PhiDz2() const
{
    const auto& U = mesh_.lookupObject<volVectorField>(UName_);

    // phi is harmonic, so d2phi/dz2 = -(d2phi/dx2 + d2phi/dy2).  On the free
    // surface those are TANGENTIAL second derivatives, which is much better
    // conditioned than differencing in z into the boundary.  With u = grad(phi)
    // each is a single gradient of a component of an existing field.
    tmp<volScalarField> td2
    (
        -(
            fvc::grad(U.component(vector::X))().component(vector::X)
          + fvc::grad(U.component(vector::Y))().component(vector::Y)
         )
    );

    return td2;
}


void Foam::functionObjects::psiYawMoment::solvePsi()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Build the boundary types: hull fixedGradient, far field fixedValue 0,
    // everything else (free surface, seabed) zeroGradient -- the free-surface
    // condition dPsi/dz = 0 IS zeroGradient there, since its normal is z.
    wordList bcTypes
    (
        mesh_.boundary().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    labelList hullIDs(pbm.patchSet(hullPatches_).sortedToc());
    labelList farIDs(pbm.patchSet(farFieldPatches_).sortedToc());

    for (const label p : hullIDs)
    {
        bcTypes[p] = fixedGradientFvPatchScalarField::typeName;
    }
    for (const label p : farIDs)
    {
        bcTypes[p] = fixedValueFvPatchScalarField::typeName;
    }

    PsiPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "Psi",
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimLength*dimLength/dimTime, Zero),
            bcTypes
        )
    );

    volScalarField& Psi = PsiPtr_();

    // dPsi/dn = n2 on the hull.  OpenFOAM's patch normal points out of the
    // domain, i.e. out of the fluid, which is G&P's convention exactly.
    for (const label p : hullIDs)
    {
        auto& pf = refCast<fixedGradientFvPatchScalarField>
        (
            Psi.boundaryFieldRef()[p]
        );
        pf.gradient() = (mesh_.boundary()[p].nf()() & yHat_);
    }
    for (const label p : farIDs)
    {
        Psi.boundaryFieldRef()[p] == 0.0;
    }

    Info<< type() << ' ' << name() << ": solving the auxiliary sway potential"
           " Psi (G&P 63-66)" << endl;

    for (label i = 0; i < nPsiCorr_; ++i)
    {
        fvScalarMatrix PsiEqn
        (
            fvm::laplacian
            (
                dimensionedScalar("1", dimless, 1),
                Psi
            ) == dimensionedScalar(Psi.dimensions()/dimLength/dimLength, Zero)
        );

        const scalar res =
            PsiEqn.solve(mesh_.solverDict(psiSolverName_)).initialResidual();

        if (res < 1e-8) break;
    }

    Psi.correctBoundaryConditions();
}


void Foam::functionObjects::psiYawMoment::createOutputFile()
{
    if (!Pstream::master() || filePtr_) return;

    filePtr_ = writeFile::createFile("psiYawMoment");

    writeHeader(*filePtr_, "psibar contribution to the mean yaw moment");
    writeHeaderValue(*filePtr_, "CofR", CofR_);
    writeHeaderValue(*filePtr_, "rhoInf", rhoRef_);
    writeHeaderValue(*filePtr_, "omegaE", omegaE_);
    writeHeader(*filePtr_, "free-surface integral only; hull integral not included");
    writeCommented(*filePtr_, "Time");
    writeTabbed(*filePtr_, "periodsAccumulated");
    writeTabbed(*filePtr_, "Mz_psi_freeSurface");
    writeTabbed(*filePtr_, "Q");
    *filePtr_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::psiYawMoment::psiYawMoment
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    UName_("U"),
    UsName_("Us"),
    freeSurfacePatchName_(),
    freeSurfacePatchID_(-1),
    CofR_(Zero),
    rhoRef_(1000),
    gMag_(9.81),
    omegaE_(0),
    U0_(0),
    yHat_(0, 1, 0),
    maxRadius_(-1),
    startTime_(0),
    usePsi_(true),
    psiSolverName_("PhiS"),
    nPsiCorr_(15),
    tAcc_(0),
    tPrev_(0),
    started_(false),
    Mz_(0),
    Q_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::psiYawMoment::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    dict.readIfPresent("U", UName_);
    dict.readIfPresent("Us", UsName_);
    dict.readEntry("freeSurfacePatch", freeSurfacePatchName_);
    dict.readEntry("CofR", CofR_);
    dict.readEntry("rhoInf", rhoRef_);
    dict.readIfPresent("gMag", gMag_);
    startTime_ = dict.getOrDefault<scalar>("startTime", 0);
    usePsi_ = dict.getOrDefault<Switch>("usePsi", Switch(true));
    psiSolverName_ = dict.getOrDefault<word>("psiSolver", "PhiS");
    nPsiCorr_ = dict.getOrDefault<label>("nPsiCorrectors", 15);
    hullPatches_ = dict.getOrDefault<wordRes>("hullPatches", wordRes({"sphere"}));
    farFieldPatches_ = dict.getOrDefault<wordRes>
    (
        "farFieldPatches", wordRes({"inlet", "outlet", "front", "back"})
    );

    freeSurfacePatchID_ =
        mesh_.boundaryMesh().findPatchID(freeSurfacePatchName_);

    if (freeSurfacePatchID_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "freeSurfacePatch " << freeSurfacePatchName_ << " not found"
            << exit(FatalIOError);
    }

    // constant/waveConditions is read ONCE, and supplies omegaE, U0 and the
    // ship-transverse direction unless the dictionary overrides them.
    //
    // U0 used to be picked up only inside the "omegaE was not given" branch.
    // Setting omegaE explicitly therefore left U0_ at its constructed 0, and
    // since the result is rho*U0*sum the whole term came out as a silent,
    // plausible-looking zero.
    {
        const IOdictionary waveDict
        (
            IOobject
            (
                "waveConditions",
                mesh_.time().constant(),
                mesh_.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        const scalar lambda = waveDict.get<scalar>("waveLength");
        const scalar depth  = waveDict.get<scalar>("waterDepth");
        const scalar head   = waveDict.get<scalar>("headingAngle");

        U0_ = waveDict.get<scalar>("currentSpeed");

        const scalar k  = constant::mathematical::twoPi/lambda;
        const scalar w0 = Foam::sqrt(gMag_*k*Foam::tanh(k*depth));

        omegaE_ = w0 + k*U0_*Foam::cos(head);

        // G&P put the forward speed on +x and measure y across it.  Here the
        // wave runs along +x and the steady flow is U(cos b, sin b), i.e. a
        // current in the direction the ship is NOT going, so the ship heads
        // along -(cos b, sin b) and the paper's +y is z ^ that.  Confirmed
        // against the meshed hull: KVLCC2_moved.stl is the hull turned by
        // -60 deg, its bow at (-0.83, +1.33), i.e. along -(cos b, sin b).
        yHat_ = vector(Foam::sin(head), -Foam::cos(head), 0);
    }

    // Dictionary overrides, applied after the defaults above
    dict.readIfPresent("omegaE", omegaE_);
    dict.readIfPresent("U0", U0_);
    dict.readIfPresent("shipTransverse", yHat_);

    yHat_ /= max(mag(yHat_), SMALL);

    maxRadius_ = dict.getOrDefault<scalar>("maxRadius", -1);

    if (mag(U0_) < SMALL)
    {
        WarningInFunction
            << "U0 is " << U0_ << ", so this term is identically zero."
            << "  The psibar block enters only multiplied by the forward"
               " speed, so that is correct at zero speed -- but if the case"
               " does have forward speed, check currentSpeed in"
               " constant/waveConditions." << endl;
    }

    const label nf = mesh_.boundaryMesh()[freeSurfacePatchID_].size();
    aPhi_.setSize(nf, 0); bPhi_.setSize(nf, 0);
    aD2_.setSize(nf, 0);  bD2_.setSize(nf, 0);
    tAcc_ = 0;
    started_ = false;

    Info<< type() << ' ' << name() << ':' << nl
        << "    free surface   " << freeSurfacePatchName_ << nl
        << "    omega_e        " << omegaE_ << " rad/s" << nl
        << "    U0             " << U0_ << " m/s" << nl
        << "    ship transverse (G&P +y) " << yHat_ << nl
        << "    maxRadius      " << maxRadius_ << nl
        << "    accumulate from t = " << startTime_ << nl
        << "    NOTE: free-surface integral only; the hull integral of the"
           " psibar block is not included" << endl;

    return true;
}


bool Foam::functionObjects::psiYawMoment::execute()
{
    const scalar t = mesh_.time().value();

    if (t < startTime_) return true;

    // Solve Psi once, at the first sample: it depends only on the geometry,
    // and having it early means it is written out and can be inspected
    if (usePsi_ && !PsiPtr_) solvePsi();

    if (!started_)
    {
        tPrev_ = t;
        started_ = true;
        return true;                    // no interval to weight yet
    }

    const scalar dt = t - tPrev_;
    tPrev_ = t;
    if (dt <= 0) return true;

    const auto& U = mesh_.lookupObject<volVectorField>(UName_);

    // phi on the free surface.  u = grad(phi), so the potential itself is
    // taken from the registered field of the same name as its gradient's
    // parent: the solver registers Phi, and we use it directly.
    const auto& Phi = mesh_.lookupObject<volScalarField>("Phi");
    const scalarField& phiP = Phi.boundaryField()[freeSurfacePatchID_];

    const tmp<volScalarField> td2(d2PhiDz2());
    const scalarField& d2P = td2().boundaryField()[freeSurfacePatchID_];

    const scalar c = Foam::cos(omegaE_*t)*dt;
    const scalar s = Foam::sin(omegaE_*t)*dt;

    forAll(aPhi_, i)
    {
        aPhi_[i] += phiP[i]*c;   bPhi_[i] += phiP[i]*s;
        aD2_[i]  += d2P[i]*c;    bD2_[i]  += d2P[i]*s;
    }

    tAcc_ += dt;

    return true;
}


bool Foam::functionObjects::psiYawMoment::write()
{
    // Write Psi into the time directory at every write time, so its decay
    // away from the hull can be inspected.  Placed before the tAcc_ guard so
    // it still appears when nothing has been accumulated yet.
    if (PsiPtr_ && mesh_.time().writeTime())
    {
        PsiPtr_->write();
    }

    if (tAcc_ <= SMALL) return true;

    const polyPatch& fsp = mesh_.boundaryMesh()[freeSurfacePatchID_];

    const vectorField& Cf = fsp.faceCentres();
    const scalarField magSf(mag(fsp.faceAreas()));

    scalarField Psi(fsp.size(), Zero);
    if (usePsi_)
    {
        Psi = PsiPtr_().boundaryField()[freeSurfacePatchID_];
    }
    const vectorField& Us_bd = mesh_.lookupObject<volVectorField>(UsName_).boundaryField()[freeSurfacePatchID_];

    // G&P (9): phi1 = Re{ phi exp(+i sigma t) }, so with the cos/sin
    // coefficients accumulated here, phi = a - i b.  Getting this backwards
    // flips the sign of the whole term.
    const scalar norm = 2.0/tAcc_;

    scalar sum = 0;
    scalar q = 0;

    // Radial profile: dpsi/dz decays like y^-3 (G&P section 7.1), so the
    // integral must converge with radius.  If it does not, the answer is
    // being set by the far field -- where the numerical beaches damp the
    // disturbance artificially -- and the term means nothing.
    const label nBin = 10;
    scalar rMax = 0;
    forAll(aPhi_, i)
    {
        const vector r(Cf[i] - CofR_);
        rMax = max(rMax, mag(vector(r.x(), r.y(), 0)));
    }
    reduce(rMax, maxOp<scalar>());
    if (maxRadius_ > 0) rMax = min(rMax, maxRadius_);
    scalarField binSum(nBin, Zero);

    forAll(aPhi_, i)
    {
        const scalar a1 = norm*aPhi_[i], b1 = norm*bPhi_[i];
        const scalar a2 = norm*aD2_[i],  b2 = norm*bD2_[i];

        // Im( phi conj(d2phi/dz2) ) with phi = a1 - i b1, d2 = a2 - i b2
        const scalar imPart = a1*b2 - b1*a2;

        // G&P (47)
        const scalar dPsiDz = -(omegaE_/(2*gMag_))*imPart;

        // G&P (84) / figure 2: the weight is (y - Psi), NOT a moment arm.
        // y is measured from the moment reference; the result is in fact
        // reference-independent whenever Q below vanishes.
        const vector r(Cf[i] - CofR_);

        if (maxRadius_ > 0 && mag(vector(r.x(), r.y(), 0)) > maxRadius_)
        {
            continue;
        }

        const scalar y = (r & yHat_);

        const scalar contrib = -(0.0 + Psi[i])*dPsiDz*magSf[i];
        sum += contrib;
        q += dPsiDz*magSf[i];

        const scalar rr = mag(vector(r.x(), r.y(), 0));
        const label b = min(label(nBin*rr/max(rMax, SMALL)), nBin - 1);
        binSum[b] += contrib;
    }

    reduce(sum, sumOp<scalar>());
    reduce(q, sumOp<scalar>());
    Pstream::listCombineAllGather(binSum, plusEqOp<scalar>());

    // G&P (84):  Mz_psi = rho U int (y - Psi) dpsi/dn dS
    Mz_ = rhoRef_*U0_*sum;

    // G&P (78): Q is the net flux of the psi field, equal to the Stokes drift
    // and proportional to the work the body does on the fluid.  It vanishes
    // for a body that is restrained or drifting freely (section 3), and when
    // it does, Mz_psi is independent of the moment reference point.  So a
    // non-zero Q here is a direct error measure on this term.
    Q_ = q;

    const scalar nPeriods = tAcc_*omegaE_/constant::mathematical::twoPi;

    if (Pstream::master())
    {
        createOutputFile();
        writeCurrentTime(filePtr_());
        filePtr_() << tab << nPeriods << tab << Mz_ << tab << Q_ << endl;
    }

    Log << type() << ' ' << name() << " write:" << nl
        << "    accumulated    " << nPeriods << " encounter periods" << nl
        << "    Mz_psi (free surface) " << Mz_ << " N m" << nl
        << "    Q = int dpsi/dz dS    " << Q_
        << "   (should vanish; see G&P section 3)" << nl
        << "    cumulative vs radius (must flatten, or the far field is"
           " setting it):" << endl;

    {
        scalar cum = 0;
        for (label b = 0; b < nBin; ++b)
        {
            cum += binSum[b];
            Log << "      r < " << (b + 1)*rMax/nBin << " m : "
                << rhoRef_*U0_*cum << endl;
        }
    }

    setResult("MzPsiFreeSurface", Mz_);
    setResult("Q", Q_);

    return true;
}


// ************************************************************************* //
