# Script to compare the matter power spectrum (linear or nonlinear) between CAMB and CLASS at a desired cosmology and redshift
#
# Usage of this script:
#
# python camb_vs_class_Pk.py <cosmology_name> <lin_or_nl>
#
# where,
# <cosmology_name> = LCDM, test1, test2 or whatever the name of your desired cosmology as defined in the cosmo_dict (see below)
# <lin_or_nl> = lin or nl (equivalently, linear or nonlinear) --> for linear or nonlinear power spectrum
#
# Example call from the terminal:
#
# python camb_vs_class_Pk.py LCDM lin

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.interpolate import interp1d

import camb
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))

import classy
from classy import Class
print('Using CLASS %s installed at %s'%(classy.__version__,os.path.dirname(classy.__file__)))

########################################
# Read terminal arguments
########################################

cosmology_name = sys.argv[1]
lin_or_nl = sys.argv[2]

########################################
# Parameters
########################################

# Note: c_min in CLASS is called A_baryon in CAMB (one of the two HMCode 2015/2016 parameters)
# Order of the paramters for a given cosmology: [Omega_m, ln(10^10*As), w0, h, c_min, z]

cosmo_dict = {
    'LCDM' : [ 0.279,  3.044522437723423, -1.0,  0.7,  3.13, 0.5], # a standard LCDM like cosmology (with A_s = 2.1e-9)
    'test1' : [ 0.12836942,  3.22501109, -2.34029847,  0.57573668,  3.38315697, 2.04090028],
    'test2' : [ 0.13063248,  3.36214195, -2.71867168,  0.58195185,  3.70005089,  0.0433379 ],
    'test3' : [ 0.13051834,  4.16525975, -1.96813609,  0.56149989,  5.34963737,  0.25553606],
    'test4' : [ 0.13003067,  3.96725238, -1.91138158,  0.69027659,  1.66192719, 0.13685009],
}

cosmo_params = cosmo_dict[cosmology_name]

Omega_b = 0.046
Omega_m = cosmo_params[0]
Omega_cdm = Omega_m - Omega_b
A_s = np.exp(cosmo_params[1]) / 1e10
n_s = 0.965

w0 = cosmo_params[2]
h = cosmo_params[3]

def eta_0_val(c_min):
    return 1.03-0.11*c_min

c_min = cosmo_params[4]
eta_0 = eta_0_val(c_min)

z_test = cosmo_params[5]

omega_b = Omega_b*h**2
omega_cdm = Omega_cdm*h**2

k = np.logspace(-4, np.log10(1.5), num=100) #Mpc^-1

kmax = 50.0

########################################
#  CAMB
########################################

pars = camb.set_params(H0=100*h, ombh2=omega_b, omch2=omega_cdm, As=A_s, ns=n_s, halofit_version='mead', HMCode_A_baryon=c_min, HMCode_eta_baryon=eta_0 )
pars.set_dark_energy(w=w0, cs2=1.0, wa=0, dark_energy_model='fluid')
pars.set_matter_power(redshifts=[z_test], kmax=kmax+50) # manually setting kmax to a larger value for better interpolation

if (lin_or_nl == 'lin' or lin_or_nl == 'linear'):
    # linear power spectrum
    pars.NonLinear = camb.model.NonLinear_none

elif (lin_or_nl == 'nl' or lin_or_nl == 'nonlinear'):
    # nonlinear power spectrum
    pars.NonLinear = camb.model.NonLinear_both

results = camb.get_results(pars)
kh, z, pk_h = results.get_matter_power_spectrum(minkh=1e-5, maxkh=kmax+20.0, npoints=300) # manually setting minkh and maxkh to a wider range for better interpolation
camb_pk = interp1d(kh*h, pk_h[0]/h**3)
P_camb = np.array([camb_pk(ki) for ki in k])

########################################
#  CLASS
########################################

commonsettings_nl  = {
    'h':h,
    'omega_b':omega_b,
    'omega_cdm': omega_cdm,
    'A_s':A_s,
    'n_s':n_s,
    'Omega_Lambda':0.0,
    'fluid_equation_of_state':'CLP',
    'w0_fld':w0,
    'wa_fld':0.0,
    'output':'mPk',
    'P_k_max_1/Mpc':kmax,
    'z_max_pk':z_test + 0.5, # manually setting this to larger than the desired redshift for better interpolation
    'non linear':'hmcode',
    'eta_0':eta_0,
    'c_min':c_min,
    }
       
cosmo_class = Class()
cosmo_class.set(commonsettings_nl)
cosmo_class.compute()

if (lin_or_nl == 'lin'):
    # linear power spectrum
    P_class = np.array([cosmo_class.pk_lin(ki, z_test) for ki in k])
elif (lin_or_nl == 'nl'):
    # nonlinear power spectrum
    P_class = np.array([cosmo_class.pk(ki, z_test) for ki in k])

########################################
#  Plots
########################################

fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]}, sharex=True)

ax = axes[0]
ax.plot(k, P_camb, color='b', label='camb')
ax.plot(k, P_class, color='r', label='class')
ax.set_ylabel('P(k) [Mpc^3]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('Matter power at z=%s'%(z_test));
ax.legend()

ax = axes[1]
ax.plot(k, P_camb/P_class - 1, color='b', label='camb/class - 1')
ax.axhline(0.05, c='grey', ls='dashed')
ax.axhline(0.0, c='grey', ls='dotted')
ax.axhline(-0.05, c='grey', ls='dashed')
ax.set_xlabel('k [Mpc]')
ax.set_ylabel('frac. diff')
ax.set_xscale('log')
ax.set_ylim([-0.1,0.1])
ax.legend()

plt.savefig('./plots/'cosmology_name+'_'+lin_or_nl+'.png')