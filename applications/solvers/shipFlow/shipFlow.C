/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    potentialFoam  with Free Surface Developed by Turgut
Group
    grpBasicSolvers

Description
    Potential flow solver which solves for the velocity potential, to
    calculate the flux-field, from which the velocity field is obtained by
    reconstructing the flux.

    \heading Solver details
    The potential flow solution is typically employed to generate initial fields
    for full Navier-Stokes codes.  The flow is evolved using the equation:

    \f[
        \laplacian \Phi = \div(\vec{U})
    \f]

    Where:
    \vartable
        \Phi      | Velocity potential [m2/s]
        \vec{U}   | Velocity [m/s]
    \endvartable

    The corresponding pressure field could be calculated from the divergence
    of the Euler equation:

    \f[
        \laplacian p + \div(\div(\vec{U}\otimes\vec{U})) = 0
    \f]

    but this generates excessive pressure variation in regions of large
    velocity gradient normal to the flow direction.  A better option is to
    calculate the pressure field corresponding to velocity variation along the
    stream-lines:

    \f[
        \laplacian p + \div(\vec{F}\cdot\div(\vec{U}\otimes\vec{U})) = 0
    \f]
    where the flow direction tensor \f$\vec{F}\f$ is obtained from
    \f[
        \vec{F} = \hat{\vec{U}}\otimes\hat{\vec{U}}
    \f]

    \heading Required fields
    \plaintable
        U         | Velocity [m/s]
    \endplaintable

    \heading Optional fields
    \plaintable
        p         | Kinematic pressure [m2/s2]
        Phi       | Velocity potential [m2/s]
                  | Generated from p (if present) or U if not present
    \endplaintable

    \heading Options
    \plaintable
        -writep   | write the Euler pressure
        -writephi | Write the final volumetric flux
        -writePhi | Write the final velocity potential
        -initialiseUBCs | Update the velocity boundaries before solving for Phi
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
//#include "fvOptions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote // output when -help option used
    (
        "Potential flow solver which solves for the velocity potential"
    );

    argList::addOption // -pName option
    (
        "pName",
        "pName",
        "Name of the pressure field"
    );

   

    argList::addBoolOption
    (
        "writephi",
        "Write the final volumetric flux field"
    );

    argList::addBoolOption
    (
        "writePhi",
        "Write the final velocity potential field"
    );

  

    argList::addBoolOption
    (
        "withFunctionObjects",
        "Execute functionObjects"
    );

    #include "addRegionOption.H" // just another option like before
    #include "addCheckCaseOptions.H" // for dry-run options
    #include "setRootCaseLists.H"  // Checks if there is a controlDict. Also can use flags to print e.g. types of BCs etc. For user-friendliness.
    #include "createTime.H" // creates runTime
    #include "createMesh.H" // reads mesh
	
    
	
	// #include "readGravitationalAcceleration.H" // self-explanatory
	#include "createControl.H" // loads custom file by Turgut located in constant dir: turgutFlow. It reads in nNonOrtho
	pisoControl turgutFlow(mesh, "turgutFlow"); // create pisoControl object. Named: turgutFlow. turgutFlow solver named in fvSolutions
	pisoControl steadyFlow(mesh, "steadyFlow");

    #include "createFields.H" // read/init fields: U, Ucur, phi(flux), zeta, zetaDx, p, p2, p3, Phi, PhiCur, PhiDx, PhiDy, PhiCurDz2 (dphi/dz**2)
	// Also read in waveCurConditions dictionary here
	
	//#include "createFvOptions.H"
	
	
	
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "PHIWaveCurSph2 by Turgut (Phi wave only)with warm start flow around 3D sphere" << endl;
	
    // Since solver contains no time loop it would never execute function objects so do it ourselves
    runTime.functionObjects().start();
	
	// Load in parameters from waveCurConditions dictionary 
	scalar steepness = waveCurConditions.lookupOrDefault<scalar>("steepness", 0.01);
	scalar wavelength = waveCurConditions.lookupOrDefault<scalar>("wavelength", 2.0);
	scalar U0 = waveCurConditions.lookupOrDefault<scalar>("currentspeed", 0.0);
	scalar hdepth = waveCurConditions.lookupOrDefault<scalar>("waterdepth", 1.0);
    scalar xdamp = waveCurConditions.lookupOrDefault<scalar>("xdamp", 40.0);
    scalar currentspeed = waveCurConditions.lookupOrDefault<scalar>("currentspeed", 0.0);
    scalar head_ang = waveCurConditions.lookupOrDefault<scalar>("head_ang", 0.0);
    scalar rampperiod = waveCurConditions.lookupOrDefault<scalar>("rampperiod", 5.0);

	// derived parameters
	const scalar wavenumber(2.0 * Foam::constant::mathematical::pi / wavelength); 
	const scalar amp(0.5 * steepness * wavelength);
	const scalar w(Foam::sqrt(9.81 * wavenumber * Foam::tanh(wavenumber*hdepth)));
	const scalar celerity(w/wavenumber);
    const scalar T(2.0 * Foam::constant::mathematical::pi / w);
    const scalar ramp_time(rampperiod * T);
     
	

	Info << "Wavelength: " << wavelength << " steepness: " << steepness << endl;
	Info << "currentspeed: " << U0 << " water depth: " << hdepth << endl;
	Info << "amp: " << amp << endl;
	Info << "k (Wavenumber): " << wavenumber << endl;
	Info << "w (Angular Frequency): " << w << endl;
	Info << "C (Celerity) " << celerity << endl;
    Info << "Period T: " << T << endl;
	
	
	// ---------- Calculation of steady potential PhiCur ---------------------------------
	Info<< nl << "Calculating steady potential PhiCur" << endl;
	while (steadyFlow.correctNonOrthogonal())
    {
        fvScalarMatrix PhiCurEqn
        (
            fvm::laplacian(dimensionedScalar("1", dimless, 1), PhiCur)
         ==
            Zero
        );

        // PhiCurEqn.setReference(PhiRefCell, PhiRefValue);
        PhiCurEqn.solve();

        // if (steadyFlow.finalNonOrthogonalIter())
        // {
        //     phi -= PhiCurEqn.flux();
        // }
    }

    // Calculate d²PhiCur/dz²
    const volTensorField PhiCurD2 = fvc::grad(fvc::grad(PhiCur));
    //PhiCurDz2 = PhiCurD2.component(tensor::ZZ);

    // For the mj terms in the linearized BC
    gradUcur = fvc::grad(Ucur);
    gradUcur.write();

    // ---------- Set PhiCur instead using NKL everywhere in the domain (as a freestream undisturbed potential) ---------------------------------
    //scalar heading	=  -30.0;
    //scalar U0	=  0.3087;
    //scalar head_ang = heading*constant::mathematical::pi/180.0;

    //scalar Ux = U0*Foam::cos(head_ang);
    //scalar Uy = U0*Foam::sin(head_ang);

    //const volVectorField& C = mesh.C();
    //forAll(C, celli)
   // {
     //   const scalar xcor = C[celli].x();
       // const scalar ycor = C[celli].y();

       // PhiCur[celli] = -Ux*xcor - Uy*ycor;
        //Ucur[celli] = vector(Ux, Uy, 0.0);
   // }

   // forAll(PhiCur.boundaryFieldRef(), patchI)
   // {
     //   auto& pPhi = PhiCur.boundaryFieldRef()[patchI];

       // forAll(pPhi, faceI)
        //{
            // use face-centres for boundary faces
          //  const vector& cf = mesh.Cf().boundaryField()[patchI][faceI];
            //PhiCur.boundaryFieldRef()[patchI][faceI] = -Ux*cf.x() - Uy*cf.y();
            //Ucur.boundaryFieldRef()[patchI][faceI] = vector(Ux, Uy, 0.0);
        //}
    //}
    
	// Calculate Ucur
	Ucur=-fvc::grad(PhiCur);//U = fvc::reconstruct(phi);
	p3 = -0.5*(Ucur & Ucur);
	p3.write();
	PhiCur.write();
	Ucur.write();
    PhiCurDz2.write();

	// ---------- End of steady potential calculation ---------------------------------
    

    // -- Main time loop
    while (runTime.loop()) 
    {
        Info << "Time = " << runTime.timeName() << "  deltaT = " << runTime.deltaTValue() << endl; 
        
        #include "CourantNo.H"

        // =========================================================================
        // PHASE 1: ADVANCE INCIDENT POTENTIAL (BEFORE SOLVING AND WRITING!)
        // =========================================================================
        scalar t = runTime.value();
        
        // Safely get a mutable reference to the internal field array
        scalarField& PhiI_internal = PhiI.primitiveFieldRef();
        const vectorField& cellCenters = mesh.C();
        const scalar ramp_factor = 0.5 * (1.0 - Foam::cos(Foam::constant::mathematical::pi * min(1.0, t / ramp_time)));
        forAll(cellCenters, cellI)
        {
            const scalar xc = cellCenters[cellI].x();
            const scalar zc = cellCenters[cellI].z();
            
            const scalar coshTerm = Foam::cosh(wavenumber*(hdepth + zc)) / Foam::cosh(wavenumber*hdepth);
            const scalar sinhTerm = Foam::sinh(wavenumber*(hdepth + zc)) / Foam::cosh(wavenumber*hdepth);
            const scalar phase = wavenumber * xc - (w + wavenumber * currentspeed * Foam::cos(head_ang)) * t;

            PhiI_internal[cellI] = -amp * (9.81 / w)
                * coshTerm
                * Foam::sin(wavenumber * xc  - (w + wavenumber * currentspeed * Foam::cos(head_ang)) * t)
                * ramp_factor;
            

            UI[cellI] = vector(
                amp * (9.81 * wavenumber / w) * coshTerm * Foam::cos(phase) * ramp_factor,
                0.0,
                amp * (9.81 * wavenumber / w) * sinhTerm * Foam::sin(phase) * ramp_factor
            );
            
        }

        // Set values on boundaries!
        forAll(mesh.boundary(), patchI)
        {
            vectorField& UIPatch = UI.boundaryFieldRef()[patchI];
            scalarField& PhiIPatch = PhiI.boundaryFieldRef()[patchI];
            vectorField& zetaIPatch = zetaI.boundaryFieldRef()[patchI];
            
            const vectorField& CfPatch = mesh.boundary()[patchI].Cf();

            forAll(CfPatch, faceI)
            {
                const scalar xc = CfPatch[faceI].x();
                const scalar zc = CfPatch[faceI].z();
                const scalar coshTerm = Foam::cosh(wavenumber*(hdepth + zc)) / Foam::cosh(wavenumber*hdepth);
                const scalar sinhTerm = Foam::sinh(wavenumber*(hdepth + zc)) / Foam::cosh(wavenumber*hdepth);
                const scalar phase = wavenumber * xc - (w + wavenumber * currentspeed * Foam::cos(head_ang)) * t;

                PhiIPatch[faceI] = -amp * (9.81 / w) * coshTerm * Foam::sin(phase) * ramp_factor;
                
                UIPatch[faceI] = vector(
                    amp * (9.81 * wavenumber / w) * coshTerm * Foam::cos(phase) * ramp_factor,
                    0.0,
                    amp * (9.81 * wavenumber / w) * sinhTerm * Foam::sin(phase) * ramp_factor
                );

                // Only apply zeta to the upper patch (assuming patch 0 or check name)
                if (mesh.boundary()[patchI].name() == "upper") {
                    zetaIPatch[faceI] = vector(0.0, 0.0, amp * Foam::cos(phase) * ramp_factor);
                }
            }
        }

        forAll(cellCenters, cellI)
        {
            const scalar xc = cellCenters[cellI].x();
            const scalar zc = cellCenters[cellI].z();
            
            const scalar we = w + wavenumber * currentspeed * Foam::cos(head_ang);
            const scalar coshTerm = Foam::cosh(wavenumber*(hdepth + zc)) / Foam::cosh(wavenumber*hdepth);
            const scalar phase = wavenumber * xc - we * t;


            // Analytical time derivative of PhiI!
            ddtPhiI[cellI] = we/w*amp * 9.81 * coshTerm * Foam::cos(phase) * ramp_factor;
        }

        // Force it onto the boundaries just like you did for PhiI and UI
        const fvBoundaryMesh& boundary = mesh.boundary();
        forAll(boundary, patchI)
        {
            scalarField& ddtPhiIPatch = ddtPhiI.boundaryFieldRef()[patchI];
            
            // Get the face centers from the fvPatch, not the polyPatch!
            const vectorField& CfPatch = boundary[patchI].Cf();
            
            forAll(CfPatch, faceI)
            {
                const scalar xc = CfPatch[faceI].x();
                const scalar zc = CfPatch[faceI].z();
                const scalar we = w + wavenumber * currentspeed * Foam::cos(head_ang);
                const scalar coshTerm = Foam::cosh(wavenumber*(hdepth + zc)) / Foam::cosh(wavenumber*hdepth);
                const scalar phase = wavenumber * xc - we* t;
                
                ddtPhiIPatch[faceI] = we/w*amp * 9.81 * coshTerm * Foam::cos(phase) * ramp_factor;
            }
        }

        // CRITICAL: Synchronize parallel ghost cells and update other boundaries
        PhiI.correctBoundaryConditions();
        zetaI.correctBoundaryConditions();
        UI.correctBoundaryConditions();       
        ddtPhiI.correctBoundaryConditions(); 
        
        // =========================================================================
        // PHASE 2: SOLVER EQUATIONS (Your Scattered Potential)
        // =========================================================================
        while (turgutFlow.correctNonOrthogonal()) 
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), PhiD) 
            ==
                Zero 
            );

            PhiEqn.setReference(PhiRefCell, PhiRefValue); 
            PhiEqn.solve();

            // 1. Calculate UD first
            UD = -fvc::grad(PhiD);
            
            // 2. Sum Phi and U BEFORE taking gradients/fluxes!
            Phi = PhiI + PhiD;
            U = UI + UD;
            zeta = zetaI + zetaD;

            // Now the BC will see the updated pressure on the next iteration!
            p = ddtPhiI + fvc::ddt(PhiD) - (Ucur & U);
            p2 = ddtPhiI + fvc::ddt(PhiD) - (Ucur & U);
            

            // 3. Update total velocities and fluxes
            if (turgutFlow.finalNonOrthogonalIter())
            {
                // 3. NOW calculate the total flux perfectly!
                phi = -fvc::snGrad(Phi) * mesh.magSf();
            }
            
        }

        // U.correctBoundaryConditions();
        Phi.correctBoundaryConditions();
        zeta.correctBoundaryConditions();
        // Sum to total
        

        
        Info << "Iterative loop ended \n" << endl;
        Info << "Continuity error from phi = " << mag(fvc::div(phi))().weightedAverage(mesh.V()).value() << endl;
        
        // p = ddtPhiI + fvc::ddt(PhiD)-(Ucur & U); 
        // p2 = ddtPhiI + fvc::ddt(PhiD)-(Ucur & U);
        p3 = 0.5*(U & U);
        
        Info << "Continuity error from U = " << mag(fvc::div(U))().weightedAverage(mesh.V()).value() << endl;


        // =========================================================================
        // PHASE 3: WRITE TO DISK
        // =========================================================================
        if (args.found("writephi")) { phi.write(); }
        if (args.found("writePhi")) { Phi.write(); }

        if (runTime.writeTime())
        {
            U.write();
            zetaDx.write();
            PhiDx.write();
            PhiDy.write();
        }
        
        // PhiI and zetaI are AUTO_WRITE, so they will be correctly written to disk here!
        runTime.write(); 
        
        runTime.printExecutionTime(Info);
    }

	runTime.functionObjects().end();

    runTime.printExecutionTime(Info);


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
