/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "vector.H"
#include "linearizedRigidBodyFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "spatialVector.H"
#include "FixedList.H"
#include <algorithm> 
#include "IOdictionary.H"
#include "fvcGrad.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearizedRigidBodyFvPatchScalarField::
linearizedRigidBodyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
     fixedGradientFvPatchField<scalar>(p, iF),
    UName_("U"),
	lastUpdateTimeIndex_(0),
    a_new_(spatialVector::zero),
    Xb_new_(spatialVector::zero),
    Ub_new_(spatialVector::zero)
{}


Foam::linearizedRigidBodyFvPatchScalarField::
linearizedRigidBodyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<scalar>(p, iF, dict, IOobjectOption::NO_READ),  // ✅ correct
    UName_(dict.getOrDefault<word>("U", "U")),
    lastUpdateTimeIndex_(0),
    a_new_(spatialVector::zero),
    Xb_new_(spatialVector::zero),
    Ub_new_(spatialVector::zero)
{
    gradient() = Zero;  // safe initialization
}


Foam::linearizedRigidBodyFvPatchScalarField::
linearizedRigidBodyFvPatchScalarField
(
    const linearizedRigidBodyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
     fixedGradientFvPatchField<scalar>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    lastUpdateTimeIndex_(ptf.lastUpdateTimeIndex_),
    a_new_(ptf.a_new_),
    Xb_new_(ptf.Xb_new_),
    Ub_new_(ptf.Ub_new_)
{}


Foam::linearizedRigidBodyFvPatchScalarField::
linearizedRigidBodyFvPatchScalarField
(
    const linearizedRigidBodyFvPatchScalarField& wbppsf
)
:
     fixedGradientFvPatchField<scalar>(wbppsf),
    UName_(wbppsf.UName_),
    lastUpdateTimeIndex_(wbppsf.lastUpdateTimeIndex_),
    a_new_(wbppsf.a_new_),
    Xb_new_(wbppsf.Xb_new_),
    Ub_new_(wbppsf.Ub_new_)
{}


Foam::linearizedRigidBodyFvPatchScalarField::
linearizedRigidBodyFvPatchScalarField
(
    const linearizedRigidBodyFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
     fixedGradientFvPatchField<scalar>(wbppsf, iF),
    UName_(wbppsf.UName_),
    lastUpdateTimeIndex_(wbppsf.lastUpdateTimeIndex_),
    a_new_(wbppsf.a_new_),
    Xb_new_(wbppsf.Xb_new_),
    Ub_new_(wbppsf.Ub_new_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
namespace Foam
{
    using Vec6 = Foam::FixedList<Foam::scalar, 6>;
    using Mat6 = Foam::FixedList<Vec6, 6>;

    // init to zero
    static inline Vec6 zero6()
    {
        Vec6 a;
        for (int i=0;i<6;++i) a[i]=0;
        return a;
    }
    // add two Vec6
    static inline Vec6 add6(const Vec6& a, const Vec6& b)
    {
        Vec6 r;
        for (int i=0;i<6;++i) r[i]=a[i]+b[i];
        return r;
    }

    // subtract two Vec6
    static inline Vec6 sub6(const Vec6& a, const Vec6& b)
    {
        Vec6 r;
        for (int i=0;i<6;++i) r[i]=a[i]-b[i];
        return r;
    }

    // scale Vec6 by scalar
    static inline Vec6 scal6(const scalar s, const Vec6& a)
    {
        Vec6 r;
        for (int i=0;i<6;++i) r[i]=s*a[i];
        return r;
    }

    // matrix-vector multiplication
    static inline Vec6 matVec6(const Mat6& A, const Vec6& x)
    {
        Vec6 r = zero6();
        for (int i=0;i<6;++i)
            for (int j=0;j<6;++j)
                r[i] += A[i][j]*x[j];
        return r;
    }

    // init Mat6to zero
    static inline Mat6 zero66()
    {
        Mat6 A;
        for (int i=0;i<6;++i)
            for (int j=0;j<6;++j)
                A[i][j]=0;
        return A;
    }

    // add two Mat6
    static inline Mat6 add66(const Mat6& A, const Mat6& B)
    {
        Mat6 R = zero66();
        for (int i=0;i<6;++i)
            for (int j=0;j<6;++j)
                R[i][j] = A[i][j] + B[i][j];
        return R;
    }

    // scale Mat6 by scalar
    static inline Mat6 scal66(const scalar s, const Mat6& A)
    {
        Mat6 R = zero66();
        for (int i=0;i<6;++i)
            for (int j=0;j<6;++j)
                R[i][j] = s*A[i][j];
        return R;
    }

    // A small robust 6×6 solver (Gaussian elimination with pivoting)
    static Vec6 solve6x6(Mat6 A, Vec6 b)
    {
        // Gaussian elimination with partial pivoting
        for (int k=0; k<6; ++k)
        {
            // Pivot
            int piv = k;
            scalar amax = Foam::mag(A[k][k]);
            for (int i=k+1; i<6; ++i)
            {
                scalar v = Foam::mag(A[i][k]);
                if (v > amax) { amax = v; piv = i; }
            }

            if (amax < VSMALL)
            {
                FatalErrorInFunction
                    << "Singular 6x6 system at row " << k
                    << " (pivot ~ 0). Check M/C/K."
                    << Foam::abort(FatalError);
            }

            // Swap rows
            if (piv != k)
            {
                for (int j=k; j<6; ++j) std::swap(A[k][j], A[piv][j]);
                std::swap(b[k], b[piv]);
            }

            // Eliminate
            const scalar Akk = A[k][k];
            for (int i=k+1; i<6; ++i)
            {
                const scalar f = A[i][k]/Akk;
                A[i][k] = 0;
                for (int j=k+1; j<6; ++j) A[i][j] -= f*A[k][j];
                b[i] -= f*b[k];
            }
        }

        // Back substitution
        Vec6 x = zero6();
        for (int i=5; i>=0; --i)
        {
            scalar s = b[i];
            for (int j=i+1; j<6; ++j) s -= A[i][j]*x[j];
            x[i] = s/A[i][i];
        }
        return x;
    }

    static inline Vec6 q_from_spatial(const spatialVector& X)
    {
        // X.w() = (phi,theta,psi) (angular)
        // X.l() = (x,y,z)         (linear)
        const vector ang = X.w();
        const vector lin = X.l();

        Vec6 q;
        q[0]=lin.x(); q[1]=lin.y(); q[2]=lin.z();
        q[3]=ang.x(); q[4]=ang.y(); q[5]=ang.z();
        return q;
    }

    static inline spatialVector spatial_from_q(const Vec6& q)
    {
        const vector lin(q[0], q[1], q[2]);
        const vector ang(q[3], q[4], q[5]);

        // spatialVector(w,l) = (angular, linear)
        return spatialVector(ang, lin);
    }

    static inline Vec6 f_from_wrench(const spatialVector& W)
    {
        const vector M = W.w();
        const vector F = W.l();

        Vec6 f;
        f[0]=F.x(); f[1]=F.y(); f[2]=F.z();
        f[3]=M.x(); f[4]=M.y(); f[5]=M.z();
        return f;
    }
}

void Foam::linearizedRigidBodyFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const Time& runTime = db().time();
    const scalar dt = runTime.deltaTValue();

    const scalar aRelax_ = 0.3;
    const scalar aDamp_  = 1.0;

    // Newmark parameters
    const scalar gamma = 0.5;
    const scalar beta  = 0.25;

    // Read properties
    const IOdictionary rbDict
    (
        IOobject
        (
            "rigidBodyMotionProperties",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    mass_ = rbDict.lookupOrDefault<scalar>("Mass", 1.0);
    xG_   = rbDict.lookupOrDefault<scalar>("xG", 0.0);
    zG_   = rbDict.lookupOrDefault<scalar>("zG", 0.0);

    I4_   = rbDict.lookupOrDefault<scalar>("I4", 1.0);
    I5_   = rbDict.lookupOrDefault<scalar>("I5", 1.0);
    I6_   = rbDict.lookupOrDefault<scalar>("I6", 1.0);
    I46_  = rbDict.lookupOrDefault<scalar>("I46", 0.0);

    C33_  = rbDict.lookupOrDefault<scalar>("C33", 0.0);
    C44_  = rbDict.lookupOrDefault<scalar>("C44", 0.0);
    C55_  = rbDict.lookupOrDefault<scalar>("C55", 0.0);
    C35_  = rbDict.lookupOrDefault<scalar>("C35", 0.0);

    heading_ = rbDict.lookupOrDefault<scalar>("heading", 0.0);

    // Init on first step
    if (runTime.timeIndex() == 1)
    {
        a_old_     = spatialVector::zero;
        Xb_old_    = spatialVector::zero;
        Ub_old_    = spatialVector::zero;
        aPrevIter_ = spatialVector::zero;

        a_new_     = spatialVector::zero;
        Xb_new_    = spatialVector::zero;
        Ub_new_    = spatialVector::zero;

        lastUpdateTimeIndex_ = 0;
    }

    // Shift once per time step
    if (runTime.timeIndex() != lastUpdateTimeIndex_)
    {
        Xb_old_ = Xb_new_;
        Ub_old_ = Ub_new_;
        a_old_  = a_new_;

        writeMotion(Xb_old_);  // motion output
    }

    // --- Build 6x6 matrices (order: [x y z phi theta psi])
    Mat6 M = zero66();
    M[0][0] = mass_;
    M[0][4] = mass_*zG_;
    M[1][1] = mass_;
    M[1][3] = -mass_*zG_;
    M[2][2] = mass_;
    M[3][3] = I4_;
    M[3][1] = -mass_*zG_;
    M[4][0] = mass_*zG_;
    M[4][4] = I5_;
    M[5][5] = I6_;
    M[3][5] = -I46_;
    M[5][3] = -I46_;

    Mat6 K = zero66();
    K[2][2] = C33_;
    K[3][3] = C44_;
    K[4][4] = C55_;
    K[2][4] = C35_;
    K[4][2] = C35_;

    // - Soft Mooring
    K[0][0] = 60.0;
    K[1][1] = 60.0;
    K[5][5] = 180.0;

    Mat6 C = zero66();  // no damping yet
    // Add roll damping
    // C[0][0] = 50.0;
    // C[1][1] = 50.0;
    // C[2][2] = 50.0;
    // C[3][3] = 50.0;
    // C[4][4] = 50.0;
    // C[5][5] = 50.0;

    // --- External load at n+1: wrench W=(M,F) -> f=[F,M]
    const spatialVector W = computeForce();
    Vec6 f = f_from_wrench(W);

    // -- to account for heading, rotate forces/moments
    scalar cH = cos(heading_);
    scalar sH = sin(heading_);

    const scalar f0_old = f[0];
    const scalar f3_old = f[3];

    f[0] = f0_old*cH + f[1]*sH;
    f[1] = -f0_old*sH + f[1]*cH;

    f[3] = f3_old*cH + f[4]*sH;
    f[4] = -f3_old*sH + f[4]*cH;

    // --- Convert old states to arrays
    const Vec6 q_n   = q_from_spatial(Xb_old_);
    const Vec6 qd_n  = q_from_spatial(Ub_old_);
    Vec6 qdd_n       = q_from_spatial(a_old_);   // treat as same ordering

    // Predictor terms
    const Vec6 xStar = add6( add6(q_n, scal6(dt, qd_n)),
                             scal6(dt*dt*(0.5 - beta), qdd_n) );

    const Vec6 vStar = add6( qd_n, scal6(dt*(1.0 - gamma), qdd_n) );

    // Effective system: A*qdd_{n+1} = f - C*vStar - K*xStar
    Mat6 A = add66( add66(M, scal66(gamma*dt, C)),
                    scal66(beta*dt*dt, K) );

    Vec6 rhs = sub6( sub6(f, matVec6(C, vStar)),
                     matVec6(K, xStar) );

    // Solve for acceleration at n+1
    Vec6 qdd_np1 = solve6x6(A, rhs);

    // Optional relaxation/damping on acceleration (component-wise)
    Vec6 qdd_prev = q_from_spatial(aPrevIter_);
    for (int i=0;i<6;++i)
    {
        qdd_np1[i] = aDamp_*(aRelax_*qdd_np1[i] + (1.0 - aRelax_)*qdd_prev[i]);
    }

    // Update v and x
    const Vec6 qd_np1 = add6(vStar, scal6(gamma*dt, qdd_np1));
    const Vec6 q_np1  = add6(xStar, scal6(beta*dt*dt, qdd_np1));

    // Store back to spatialVectors
    a_new_     = spatial_from_q(qdd_np1);
    Ub_new_    = spatial_from_q(qd_np1);
    Xb_new_    = spatial_from_q(q_np1);
    aPrevIter_ = a_new_;

    // Debug
    Info<< "BC iteration " << runTime.timeIndex()
        << " | z=" << q_np1[2]
        << " | Fz=" << f[2]
        << " | My=" << f[4]
        << endl;

    // --- Apply Neumann BC using linear velocity only
    tmp<vectorField> tnf = patch().nf();
    const vectorField& n = tnf();
    const vectorField& Cf = patch().Cf();

    

    const vector xRef(xG_, 0.0, 0.0);

    // -- Ulin and omega with respect to global frame
    cH = cos(-heading_);
    sH = sin(-heading_);

    const vector Ulin(qd_np1[0] * cH + qd_np1[1] * sH, -qd_np1[0] * sH + qd_np1[1] * cH, qd_np1[2]);
    const vector omega(qd_np1[3] * cH + qd_np1[4] * sH, -qd_np1[3] * sH + qd_np1[4] * cH, qd_np1[5]);

    vectorField UbFace(Cf.size());
    forAll(Cf, i)
    {
        UbFace[i] = Ulin + (omega ^ (Cf[i] - xRef));
    }

    // mj term:
    // -- Load Ucur 
    const volVectorField& Ucur = db().lookupObject<volVectorField>("Ucur");

    // current patch index
    const label patchi = patch().index();

    // patch unit normals
    const vectorField np(patch().nf());

    // Volume gradient of Ucur, then take boundary values on this patch
    tmp<volTensorField> tGradU = fvc::grad(Ucur);
    const volTensorField& gradU = tGradU();

    const fvPatchTensorField& gradUpP = gradU.boundaryField()[patchi];
    tensorField gradUp(gradUpP);  // make owned copy

    // Displacement amplitudes (mean-surface): translation X and small rotation theta
    const vector X = Xb_new_.l();
    const vector th = Xb_new_.w();

    // Rotate to global frame
    vector Xg( X.x()*cH + X.y()*sH, -X.x()*sH + X.y()*cH, X.z() );
    vector thg( th.x()*cH + th.y()*sH, -th.x()*sH + th.y()*cH, th.z() );

    scalarField mj(patch().size(), 0.0);

    forAll(Cf, i)
    {
        const vector r = Cf[i] - xRef;

        // Total displacement at the face: S = X + theta x r
        const vector S = Xg + (thg ^ r);

        // Approximate grad(n·U) = (gradU)^T · n   (ignores curvature term U·∇n)
        const vector grad_nDotU = (gradUp[i].T() & np[i]);

        // mj = - S · grad(n·U)
        mj[i] = (S & grad_nDotU);
    }

    
	// Info << "mj sample: " << mj[0] << endl;
    scalarField rhsBC = -(n & UbFace) - mj;
    // scalarField rhsBC = -(n & UbFace);
    this->gradient() = rhsBC;
    // Info<< "Updated linearized rigid body BC on patch " << patch().name()
    //     << " Xb_new " << Xb_new_
    //     << " (timeIndex=" << runTime.timeIndex() << ")"
    //     << endl;

    // Pout << "print a few mj terms for debugging: " << mj[0] << " " << mj[mj.size()/2] << " " << mj[mj.size()-1] << endl;

    lastUpdateTimeIndex_ = runTime.timeIndex();

    fixedGradientFvPatchField<scalar>::updateCoeffs();
}







// * * * * * * * Helper Functions  * * * * * * * * * * * * * //



namespace Foam {  

    void Foam::linearizedRigidBodyFvPatchScalarField::createMotionFile()
    {
        
        if (!Pstream::master()) return;  // Only master writes

        if (motionFilePtr_.valid()) return;
        
        const Time& runTime = db().time();

        //fileName dir(runTime.path()/"postProcessing"/"rigidBodyMotion");
        fileName dir("./postProcessing/rigidBodyMotion");  
        mkDir(dir);

        fileName filePath(dir/"motionTurgut.dat");
        mkDir(filePath.path());

        Info<< "Creating motion output file: " << filePath << endl;

        motionFilePtr_.reset(new OFstream(filePath));
        writeMotionHeader(motionFilePtr_());
    }


    void Foam::linearizedRigidBodyFvPatchScalarField::writeMotionHeader(OFstream& os) const
    {
        os << "# Rigid Body Motion: " << patch().name() << endl;
        os << "# Time   eta1   eta2   eta3   eta4   eta5   eta6" << endl;
        
    }


    void Foam::linearizedRigidBodyFvPatchScalarField::writeMotion(const spatialVector& W)
    {
        if (!Pstream::master()) return;  // Only master writes

        
        if (!motionFilePtr_.valid())
            createMotionFile();

        const Time& runTime = db().time();

        const vector& lin = W.l();   // surge, sway, heave forces
        const vector& ang = W.w();   // roll, pitch, yaw moments

        (*motionFilePtr_)
            << runTime.timeOutputValue() << '\t'
            << lin.x() << '\t'
            << lin.y() << '\t'
            << lin.z() << '\t'
            << ang.x() << '\t'
            << ang.y() << '\t'
            << ang.z() << endl;
    }


    Foam::spatialVector Foam::linearizedRigidBodyFvPatchScalarField::computeForce()
    {
        const vectorField& Sfb = patch().Sf();
        const vectorField& Cf  = patch().Cf();

        const volScalarField& p = db().lookupObject<volScalarField>("p2");
        const scalarField& pPatch = p.boundaryField()[patch().index()];

        vector F = vector::zero;
        vector M = vector::zero;
        const scalar rho = 1000.0;  // water density
        const vector rRef(xG_, 0.0, 0.0);   // moment reference point

        forAll(Sfb, i)
        {
            const vector f = pPatch[i] * Sfb[i];
            F += f * rho;

            const vector r = Cf[i] - rRef;
            M += r ^ f * rho;
        }

        reduce(F, sumOp<vector>());
        reduce(M, sumOp<vector>());

        // spatialVector stores (angular, linear) = (moment, force)
        return spatialVector(M, F);
    }
}




void Foam::linearizedRigidBodyFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        linearizedRigidBodyFvPatchScalarField
    );
}

// ************************************************************************* //
