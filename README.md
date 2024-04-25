# camb-vs-class
Scripts to compare quantities between CAMB and CLASS

First, please make sure to install the python packages of these two codes. See installation instructions in these links:

CAMB: https://camb.readthedocs.io/en/latest/

CLASS: https://github.com/lesgourg/class_public/tree/master 

## Power spectrum

**camb_vs_class_Pk.py** compares the power spectrum (linear or nonlinear)  computed by the two codes (check description at the top of the file to see how to use it).

*Status*: currently observing a large discrepancy (upto 5%) between CAMB and CLASS P(k) already at the linear level --> possible issues with the code which need to be fixed! 