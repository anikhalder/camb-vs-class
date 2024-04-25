# camb-vs-class
Scripts to compare quantities between CAMB and CLASS

First, please make sure to install the python packages of these two codes. See installation instructions in these links:

CAMB: https://camb.readthedocs.io/en/latest/

CLASS: https://github.com/lesgourg/class_public/tree/master 

## Power spectrum

**camb_vs_class_Pk.py** compares the power spectrum (linear or nonlinear) computed by the two codes (check description at the top of the file to see how to use it).

*Status*: CAMB and CLASS matter power spectra agree to within a sub-percent level for the linear power spectrum. For nonlinear power spectrum (esp. for cosmologies far away from LCDM) CAMB and CLASS start to differ (significantly) more than 1% on small scales (i.e. k > 1 [Mpc^-1]) 

Important note: CAMB by default has one massive neutrino species whereas CLASS sets that to zero. Therefore, need to explicitly turn off massive neutrinos in CAMB parameter settings for a fair comparison to CLASS.