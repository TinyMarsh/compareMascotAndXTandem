# compareMascotAndXTandem

This program takes the RAW file name of the Thermo RAW file we want to compare as a single argument. The program looks in the location F:\\data\\CTAM\\selectedRuns\\RAW\\ for that RAW file.

4 CSV files are output with the following columns respectively, where each row is a common peptide;
	FWHMs.csv: Mascot FWHM isotope 1, Mascot FWHM isotope 2, Mascot FWHM isotope 3, X!Tandem FWHM isotope 1, X!Tandem FWHM isotope 2, X!Tandem FWHM isotope 3
	peakAreas.csv: Mascot peak area (peptide abundance), X!Tandem peak area (peptide abundance)
	ChiSquared.csv: Mascot chi-squared, X!Tandem chi-squared, Mascot degrees of freedom, X!Tandem degrees of freedom
	peptideInfo.csv: Mascot Pescal++ dbID, X!Tandem Pescal++ dbID, Mascot m/z, X!Tandem m/z, Mascot RT, X!Tandem RT, charge (z)
	bestPeakRts.csv: Mascot best peak RT, X!Tandem best peak RT
