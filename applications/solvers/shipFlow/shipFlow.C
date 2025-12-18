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
	
    
	
	#include "readGravitationalAcceleration.H" // self-explanatory
	#include "createControl.H" // loads custom file by Turgut located in constant dir: turgutFlow. It reads in nNonOrtho
	pisoControl turgutFlow(mesh, "turgutFlow"); // create pisoControl object. Named: turgutFlow. turgutFlow solver named in fvSolutions
	pisoControl steadyFlow(mesh, "steadyFlow");
	// TODO: add control for steady flow?
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
	// const scalar R0 = readScalar(waveCurConditions.lookup("R0")); // radius of sphere
	// const scalar xc = readScalar(waveCurConditions.lookup("xc")); // x coordinate of sphere center

	// derived parameters
	const scalar wavenumber(2.0 * Foam::constant::mathematical::pi / wavelength); 
	const scalar amp(0.5 * steepness * wavelength);
	const scalar w(U0 * wavenumber + Foam::sqrt(9.81 * wavenumber * Foam::tanh(wavenumber*hdepth)));
	const scalar celerity(w/wavenumber); 
	

	Info << "Wavelength: " << wavelength << " steepness: " << steepness << endl;
	Info << "currentspeed: " << U0 << " water depth: " << hdepth << endl;
	Info << "amp: " << amp << endl;
	Info << "k (Wavenumber): " << wavenumber << endl;
	Info << "w (Angular Frequency): " << w << endl;
	Info << "C (Celerity) " << celerity << endl;
	
	
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

	// Calculate Ucur
	Ucur=-fvc::grad(PhiCur);//U = fvc::reconstruct(phi);
	p3 = -0.5*(Ucur & Ucur);
	p3.write();
	PhiCur.write();
	Ucur.write();

	// ---------- End of steady potential calculation ---------------------------------

	while (runTime.run()) // main loop
    {
		
		

		++runTime; // why is this needed? -> because runtime.run() is used not .loop(). So only checks if it should finish.
		// no auto update of time
		
		Info << "currentspeed: " << U0 << " water depth: " << hdepth << endl;
		Info << " \n Current time: " << runTime.timeName() << endl;	
		
		// Non-orthogonal velocity potential corrector loop
		// while (turgutFlow.correctNonOrthogonal()) // Numer given in fvSolution for turgutFlow
		// {
			
			fvScalarMatrix PhiEqn
			(
				fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)  // div(Gamma*grad(Phi)) where Gamma("Name", dimension, value) = 1
			  ==
				Zero //dimensionedScalar("1", dimensionSet(0,-2,1,0,0,0,0), 1)*fvOptions(Phi) //fvc::div(phi)//
			);

			
			PhiEqn.setReference(PhiRefCell, PhiRefValue); // Set reference to fix the potential level at give sel in fvSolution
			PhiEqn.solve();

			
			
			//if (turgutFlow.finalNonOrthogonalIter())
			//{
				
			//	phi -= PhiEqn.flux();
		    //}
		// }
		
		Info << "Iterative loop ended \n" << endl;
		//MRF.makeAbsolute(phi);

		Info<< "Continuity error from phi = "
			<< mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
			<< endl;
		
		p = fvc::ddt(Phi)-(Ucur & U); // integrated U**2 term should be 0
		U=-fvc::grad(Phi);//U = fvc::reconstruct(phi);
		//U2=-fvc::snGrad(Phi);
		p2 = fvc::ddt(Phi)-(Ucur & U);
		p3 = 0.5*(U & U);
		
		phi = fvc::flux(U);
		
		Info<< "Continuity error from U = "
			<< mag(fvc::div(U))().weightedAverage(mesh.V()).value()
			<< endl;
		
		
		Info<< "Continuity error  = "
        << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
        << endl;
		
		
		
		// Optionally write the volumetric flux, phi
		if (args.found("writephi"))
		{
			phi.write();
		}

		// Optionally write velocity potential, Phi
		if (args.found("writePhi"))
		{
			Phi.write();
		}
		
		
		
		runTime.write(); // write time directory when needed
		
		runTime.printExecutionTime(Info);
		
    }

	runTime.functionObjects().end();

    runTime.printExecutionTime(Info);


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
