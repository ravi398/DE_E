import numpy as np
import matplotlib.pyplot as plt
from classy import Class
from scipy.interpolate import interp1d

# Define redshift grid
zeff = np.linspace(0.1, 3.0, 200)

params = {
        'output': 'tCl,pCl,lCl,mPk',
        'lensing': 'yes',
        'l_max_scalars': 2600,
        'n_s': 0.9678,
        'ln10^{10}A_s': 3.054,
        'tau_reio': 0.0586,
        'omega_b': 0.02240,
        #'theta_s_100':1.0413,
        'h': 0.70707,
        'P_k_max_h/Mpc': 100.0,
        'z_max_pk': 2.0,
        'omega_cdm': 0.12,
        'N_ur': 0.00441,
        'N_ncdm': 1,
        'deg_ncdm': 3,
        'Omega_Lambda':0,
        'Omega_scf':0,
        'fluid_equation_of_state':'EDE',
        'w0_fld':-0.9,
        'Omega_EDE':1e-3,
        'cs2_fld':1,
        'integral_type':'num',
        'background_verbose':4,


    }
cosmo = Class()
cosmo.set(params)
cosmo.compute()
print(cosmo.get_background().keys)
print(cosmo.h(),cosmo.sigma(8,0))
redshift=cosmo.get_background()['z']
rho_de=cosmo.get_background()["(.)rho_fld"]
rho_crit=cosmo.get_background()["(.)rho_crit"]
w_de=cosmo.get_background()["(.)w_fld"]

cosmo.struct_cleanup()
cosmo.empty()





plt.semilogx(redshift,w_de,label='equation of state $w_{ede}$')
plt.semilogx(redshift,rho_de/rho_crit,label='$\Omega_{ede}$')
#plt.semilogx(redshift,1e-3/rho_crit,label='omega')




plt.xlabel("redshift")
plt.ylabel(r"$w_{ede}$")
plt.title("Evolution of fluid density fraction")
plt.xlim(0.01,100000)
plt.legend()
plt.show()

plt.semilogx(redshift,rho_de/rho_crit,label='$\Omega_{ede}$')
plt.xlabel("redshift")
plt.ylabel(r"$\rho_{\rm fld} / \rho_{\rm crit}$")
plt.title("Evolution of fluid density fraction")
plt.xlim(0.01,100000)
plt.legend()
plt.show()

plt.savefig("rho_ratio_vs_a.png", dpi=150)



