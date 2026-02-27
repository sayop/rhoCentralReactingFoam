// Read div scheme preference from system/fvSchemes
IOdictionary fvSchemesDict
(
    IOobject
    (
        "fvSchemes",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dictionary& divDict = fvSchemesDict.subDict("divSchemes");

word divScheme("div(phi,Yi)");
if (!divDict.found(divScheme))
{
    Info<< "No 'div(phi,Yi)' entry found in fvSchemes/divSchemes. "
        << "Using 'div(phi,default)'." << nl << endl;
    divScheme = "div(phi,default)";
}

forAll(Y, i)
{
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];
	fvScalarMatrix YiEqn
        (
	    fvm::ddt(rho, Yi)
          + fvm::div(phi, Yi, divScheme)
	  + thermophysicalTransport->divj(Yi)
        );

	// Optional: reduce relaxation for faster convergence
        //YiEqn.relax(); // or remove entirely with direct solver
        YiEqn.solve("Yi");
    }
}
// normalize once, after all species are updated
composition.normalise();
