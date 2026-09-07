/*---------------------------------------------------------------------------*/

#include "nearFieldForm.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "uniformDimensionedFields.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(nearFieldForm, 0);
    addToRunTimeSelectionTable(functionObject, nearFieldForm, dictionary);
}
}

namespace
{

//- The segment two faces share, if they share exactly two points
bool sharedSegment
(
    const Foam::face& a,
    const Foam::face& b,
    const Foam::pointField& pts,
    Foam::point& p0,
    Foam::point& p1
)
{
    Foam::label s[2] = {-1, -1};
    Foam::label n = 0;

    forAll(a, i)
    {
        forAll(b, j)
        {
            if (a[i] == b[j])
            {
                if (n < 2) s[n] = a[i];
                ++n;
                break;
            }
        }
    }

    if (n != 2) return false;

    p0 = pts[s[0]];
    p1 = pts[s[1]];
    return true;
}

} // End anonymous namespace


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vector Foam::functionObjects::nearFieldForm::bodyVector
(
    const word& nm
) const
{
    // linBodyMotion publishes these.  Absent -- a restrained body, or a
    // postProcess run where the condition never fires -- the motion is zero,
    // which is the right answer rather than a fatal error.
    const auto* ptr = mesh_.findObject<uniformDimensionedVectorField>(nm);
    return ptr ? ptr->value() : vector::zero;
}


void Foam::functionObjects::nearFieldForm::buildSurfaceStencil() const
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    label n = 0;
    for (const label patchi : hullPatchIDs_) n += pbm[patchi].size();

    lsNbr_.setSize(n);
    lsInvA_.setSize(n, tensor::zero);
    lsPatchStart_.setSize(hullPatchIDs_.size(), 0);

    label off = 0;
    label nPoor = 0;

    forAll(hullPatchIDs_, k)
    {
        const label patchi = hullPatchIDs_[k];
        const polyPatch& pp = pbm[patchi];
        const labelListList& ff = pp.faceFaces();
        const vectorField& fc = pp.faceCentres();
        const vectorField nu(mesh_.boundary()[patchi].nf());

        lsPatchStart_[k] = off;

        forAll(pp, i)
        {
            lsNbr_[off + i] = ff[i];

            symmTensor A(Zero);
            for (const label j : ff[i])
            {
                const vector d(fc[j] - fc[i]);
                A += symm(d*d);
            }

            // the normal direction carries no surface data; regularise so the
            // matrix inverts, then project the result tangential afterwards
            const tensor Areg(tensor(A) + nu[i]*nu[i]);

            if (ff[i].size() < 2 || mag(det(Areg)) < SMALL)
            {
                ++nPoor;                       // leave zero -> fall back to U
            }
            else
            {
                lsInvA_[off + i] = inv(Areg);
            }
        }

        off += pp.size();
    }

    reduce(nPoor, sumOp<label>());
    reduce(n, sumOp<label>());

    if (Pstream::master())
    {
        Info<< "    surface-gradient stencil on " << n << " hull faces ("
            << nPoor << " with too few neighbours, falling back to U)" << endl;
    }

    lsValid_ = true;
}


void Foam::functionObjects::nearFieldForm::reset()
{
    quadForce_ = quadMoment_ = Zero;
    motForce_  = motMoment_  = Zero;
    wlForce_   = wlMoment_   = Zero;
}


void Foam::functionObjects::nearFieldForm::hullIntegral()
{
    const bool eq27_ = (formulation_ == "eq27");

    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const auto& p = mesh_.lookupObject<volScalarField>(pName_);
    const auto& Phi = mesh_.lookupObject<volScalarField>(PhiName_);

    if (surfaceTangential_ && !lsValid_) buildSurfaceStencil();

    // grad(p) is needed for the transfer term X1.grad(P1).  Its boundary
    // value has the normal row set from snGrad and the tangential rows
    // extrapolated from the cell, so the tangential part is first order --
    // the same limitation the m-terms have.
    const tmp<volVectorField> tGradP(fvc::grad(p));

    const vector xi(bodyVector("bodyDisp"));
    const vector th(bodyVector("bodyRot"));

    for (const label patchi : hullPatchIDs_)
    {
        const fvPatch& pp = mesh_.boundary()[patchi];

        const vectorField nu(pp.nf());
        const scalarField& area = pp.magSf();
        const vectorField& Cf = pp.Cf();

        vectorField u(U.boundaryField()[patchi]);

        if (surfaceTangential_)
        {
            // u_t = grad_s(phi) from the patch values; u.nu kept from U
            const label k = hullPatchIDs_.find(patchi);
            const label off = lsPatchStart_[k];

            const polyPatch& ppoly = mesh_.boundaryMesh()[patchi];
            const vectorField& fc = ppoly.faceCentres();
            const scalarField& phiB = Phi.boundaryField()[patchi];

            forAll(u, i)
            {
                if (lsInvA_[off + i] == tensor::zero) continue;   // fallback

                vector b(Zero);
                for (const label j : lsNbr_[off + i])
                {
                    b += (fc[j] - fc[i])*(phiB[j] - phiB[i]);
                }

                const vector G(lsInvA_[off + i] & b);
                const vector Gt(G - (G & nu[i])*nu[i]);

                u[i] = Gt + (u[i] & nu[i])*nu[i];
            }
        }
        const scalarField& pb = p.boundaryField()[patchi];
        const vectorField& gp = tGradP().boundaryField()[patchi];

        forAll(nu, i)
        {
            const vector r(Cf[i] - CofR_);
            const vector S(xi + (th ^ r));

            // Chen (18): P2 = -rho/2 |grad phi|^2, psibar dropped
            const vector fQ(-0.5*rhoRef_*magSqr(u[i])*nu[i]*area[i]);

            vector fM(Zero);
            vector mM(Zero);

            if (eq27_)
            {
                // Chen (25) first line: X1.grad(P1) and [P1 - g(X1.e3)]T1
                fM =
                (
                    rhoRef_*(S & gp[i])*nu[i]
                  + rhoRef_*(pb[i] - gMag_*S.z())*(th ^ nu[i])
                )*area[i];

                // and the lever arm itself moves: <X1 ^ P1 nu>
                mM = (r ^ fM) + (S ^ (rhoRef_*pb[i]*nu[i]*area[i]));
            }
            else
            {
                // Chen (29): only the NORMAL displacement enters, and the
                // rotation term is absorbed
                fM = rhoRef_*(S & nu[i])*gp[i]*area[i];
                mM = (r ^ fM);
            }

            quadForce_ += fQ;
            motForce_ += fM;

            quadMoment_ += (r ^ fQ);
            motMoment_ += mM;
        }
    }
}


void Foam::functionObjects::nearFieldForm::waterlineIntegral()
{
    const bool eq27_ = (formulation_ == "eq27");

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const polyPatch& fsPatch = pbm[freeSurfacePatchID_];

    const label fsStart = fsPatch.start();
    const label fsEnd = fsStart + fsPatch.size();

    const faceList& faces = mesh_.faces();
    const pointField& pts = mesh_.points();
    const cellList& cells = mesh_.cells();
    const labelList& owner = mesh_.faceOwner();

    const auto& zeta = mesh_.lookupObject<volScalarField>(zetaName_);
    const scalarField& zetaP = zeta.boundaryField()[freeSurfacePatchID_];

    const vector xi(bodyVector("bodyDisp"));
    const vector th(bodyVector("bodyRot"));

    for (const label patchi : hullPatchIDs_)
    {
        const polyPatch& hp = pbm[patchi];
        const vectorField nuAll(mesh_.boundary()[patchi].nf());

        forAll(hp, hf)
        {
            const label facei = hp.start() + hf;
            const label celli = owner[facei];

            for (const label facej : cells[celli])
            {
                if (facej < fsStart || facej >= fsEnd) continue;

                point p0(Zero), p1(Zero);
                if (!sharedSegment(faces[facej], faces[facei], pts, p0, p1))
                {
                    continue;
                }

                const scalar L = mag(p1 - p0);
                if (L < SMALL) continue;

                const point mid(0.5*(p0 + p1));
                const vector r(mid - CofR_);
                const vector S(xi + (th ^ r));

                // eq27 uses the RELATIVE elevation; eq29 has absorbed the
                // motion part into the hull integral and uses zeta alone
                const scalar zr =
                    eq27_
                  ? zetaP[facej - fsStart] - S.z()
                  : zetaP[facej - fsStart];

                // 1/cos(alpha) for the hull inclination, Chen (13)
                const vector& nu = nuAll[hf];
                const scalar cosA = Foam::sqrt(max(1 - sqr(nu.z()), SMALL));

                const vector f
                (
                    0.5*rhoRef_*gMag_*sqr(zr)*nu*(L/cosA)
                );

                wlForce_ += f;
                wlMoment_ += (r ^ f);
            }
        }
    }
}


void Foam::functionObjects::nearFieldForm::createFiles()
{
    if (!Pstream::master() || forceFilePtr_) return;

    forceFilePtr_ = createFile("force");
    momentFilePtr_ = createFile("moment");

    for (auto* os : {forceFilePtr_.get(), momentFilePtr_.get()})
    {
        writeHeader(*os, "Mean wave loads, nearfield formulation (Chen 2022)");
        writeHeaderValue(*os, "CofR", CofR_);
        writeHeaderValue(*os, "rhoInf", rhoRef_);
        writeCommented(*os, "Time");
        writeTabbed(*os, "total_x\ttotal_y\ttotal_z");
        writeTabbed(*os, "quadratic_x\tquadratic_y\tquadratic_z");
        writeTabbed(*os, "motion_x\tmotion_y\tmotion_z");
        writeTabbed(*os, "waterline_x\twaterline_y\twaterline_z");
        *os << endl;
    }
}


void Foam::functionObjects::nearFieldForm::writeFiles()
{
    if (!Pstream::master()) return;

    createFiles();

    writeCurrentTime(forceFilePtr_());
    forceFilePtr_() << tab << force() << tab << quadForce_
                    << tab << motForce_ << tab << wlForce_ << endl;

    writeCurrentTime(momentFilePtr_());
    momentFilePtr_() << tab << moment() << tab << quadMoment_
                     << tab << motMoment_ << tab << wlMoment_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::nearFieldForm::nearFieldForm
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    UName_("U"),
    pName_("p"),
    zetaName_("zeta"),
    PhiName_("Phi"),
    formulation_("eq29"),
    surfaceTangential_(false),
    lsValid_(false),
    freeSurfacePatchID_(-1),
    CofR_(Zero),
    rhoRef_(1000),
    gMag_(9.81)
{
    reset();
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::nearFieldForm::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    dict.readIfPresent("U", UName_);
    dict.readIfPresent("p", pName_);
    dict.readIfPresent("zeta", zetaName_);
    dict.readIfPresent("Phi", PhiName_);
    surfaceTangential_ =
        dict.getOrDefault<Switch>("surfaceTangential", Switch(false));
    dict.readEntry("CofR", CofR_);
    dict.readEntry("rhoInf", rhoRef_);
    dict.readIfPresent("gMag", gMag_);
    formulation_ = dict.getOrDefault<word>("formulation", "eq29");

    if (formulation_ != "eq27" && formulation_ != "eq29")
    {
        FatalIOErrorInFunction(dict)
            << "formulation must be eq27 or eq29, got " << formulation_
            << exit(FatalIOError);
    }

    hullPatchIDs_ =
        mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches")).sortedToc();

    if (hullPatchIDs_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "no hull patches matched" << exit(FatalIOError);
    }

    dict.readEntry("freeSurfacePatch", freeSurfacePatchName_);
    freeSurfacePatchID_ =
        mesh_.boundaryMesh().findPatchID(freeSurfacePatchName_);

    if (freeSurfacePatchID_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "freeSurfacePatch " << freeSurfacePatchName_ << " not found"
            << exit(FatalIOError);
    }

    Info<< type() << ' ' << name() << ':' << nl
        << "    hull patches   " << hullPatchIDs_.size() << nl
        << "    free surface   " << freeSurfacePatchName_ << nl
        << "    CofR           " << CofR_ << nl
        << "    formulation    Chen " << formulation_ << nl
        << "    tangential u   "
        << (surfaceTangential_ ? "surface gradient of Phi" : "boundary U")
        << nl
        << "    NOTE: psibar and the O(Fr^2) steady-flow blocks are dropped"
        << endl;

    return true;
}


bool Foam::functionObjects::nearFieldForm::execute()
{
    reset();

    hullIntegral();
    waterlineIntegral();

    reduce(quadForce_, sumOp<vector>());
    reduce(quadMoment_, sumOp<vector>());
    reduce(motForce_, sumOp<vector>());
    reduce(motMoment_, sumOp<vector>());
    reduce(wlForce_, sumOp<vector>());
    reduce(wlMoment_, sumOp<vector>());

    Log << type() << ' ' << name() << " write:" << nl
        << "    total     " << force() << nl
        << "    quadratic " << quadForce_ << nl
        << "    motion    " << motForce_ << nl
        << "    waterline " << wlForce_ << endl;

    setResult("force", force());
    setResult("moment", moment());

    return true;
}


bool Foam::functionObjects::nearFieldForm::write()
{
    if (writeToFile())
    {
        writeFiles();
    }

    return true;
}


Foam::vector Foam::functionObjects::nearFieldForm::force() const
{
    return quadForce_ + motForce_ + wlForce_;
}


Foam::vector Foam::functionObjects::nearFieldForm::moment() const
{
    return quadMoment_ + motMoment_ + wlMoment_;
}


// ************************************************************************* //
