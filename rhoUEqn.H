{
    // Hyperbolic (inviscid) update in conservative form for the momentum rhoU
    fvVectorMatrix rhoUEqn
    (
        fvm::ddt(rhoU) 
      + fvc::div(phiUp)
    );

    // Usually no relaxation here in a shock-capturing density-based update
    // rhoUEqn.relax();

    rhoUEqn.solve();

    // Reconstruct primitive velocity
    U.ref() = rhoU()/rho();
    U.correctBoundaryConditions();

    // Keep boundary consistency
    rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

    // Optional viscous correction applied to the primitive velocity U (parabolic terms), 
    // then it re-forms rhoU = rho*U
    if (!inviscid)
    {
	fvVectorMatrix UEqnVisc
        (
            fvm::ddt(rho, U) - fvc::ddt(rho, U)
          - fvm::laplacian(muEff, U)
          - fvc::div(tauMC)
        );

	UEqnVisc.relax();
	UEqnVisc.solve();

        rhoU = rho*U;
	rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
    }
}
