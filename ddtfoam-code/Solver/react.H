{
    omegaC_ign *= 0.0;

    const label tautableDim = 3;

    if(autoIgnitionTauModel != "off")
    {
	  volScalarField dc1 = scalar(1.0)-c;
	  //dc1.min(scalar(1.0)/ignRelax);
	  omegaC_ign = rho*dc1/runTime.deltaT() * Foam::pos(tau-scalar(0.99));
    }

    if(autoIgnitionTauModel == "standard") forAll(p,i)
    {
    	// induction time:
    	List<scalar> lookupCoordinates(tautableDim,0.0);
		lookupCoordinates[0]=fH[i];
		lookupCoordinates[1]=p[i];		// old values, no thermo.update()
		lookupCoordinates[2]=T[i];

    	const List<scalar> &LOG_t_ign_temp = tignTable.lookUpBounded(lookupCoordinates,false);
    	tIgn[i] = min(max(1e-14,pow(10,LOG_t_ign_temp[0])),1e9);
    }
    else if(autoIgnitionTauModel == "ann") forAll(p,i) {
        double coordinates[3];
		coordinates[0]=fH[i];
		coordinates[1]=p[i];		// old values, no thermo.update()
		coordinates[2]=T[i];
		tIgn[i] = min(max(1e-14, tign_model(coordinates)),1e9);
    }
    tIgn.correctBoundaryConditions(); // only for post-processing: update zeroGradient BC

    omegaTau = rho / tIgn;

}



