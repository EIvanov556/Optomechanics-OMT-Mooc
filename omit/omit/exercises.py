import numpy as np
from cmath import *
from ipywidgets import interact, interactive, fixed, interact_manual
# https://ipywidgets.readthedocs.io/en/stable/examples/Using%20Interact.html
import matplotlib.pyplot as plt
import ipywidgets as widgets


class InteractiveFunction:
    params = dict()
    main_fun = None
    int_fun = None

    def run_raw(self):
        if self.main_fun is None:
            raise AttributeError('No function attributoed.')
        else:
            self.main_fun()

    def run(self):
        if self.int_fun is None:
            raise AttributeError('No function attributoed.')
        else:
            self.int_fun()


class OpticsOnly(InteractiveFunction):
    params = dict([('eta', 1)])

    def main_fun(self, eta=1):
        # eta is k0/kext
        # eta = params['eta']
        omegas = np.arange(-5, 5, 0.01)  # in units of kext
        susceptibility = 1 / (((1+eta)/2)+1j*omegas)
        r = np.abs(1-susceptibility)**2

        plt.plot(omegas, r)
        plt.ylim(0, 1.05)
        plt.xlabel('detuning from resonance')
        plt.ylabel('R')
        plt.show()

        if eta == 0:
            print('Extremely overcoupled: \nThe losses are negligible compared to the rate at which light can escape. '
                  '\nAll the light is then reflected back out of the cavity')
        elif eta < 1:
            print('Overcoupled: \nThe cavity losses are smaller than the external coupling rate.')
        elif eta == 1:
            print('Critical coupling:\nThe rate at which the light enters the cavity matches exactly the decay rate.'
                  '\nAt resonance, all photons will be dissipated')
        else:
            bad_cavity_text = """Undercoupled:\nThis is the so-called "bad cavity" regime. All excitations die out 
            faster than they can be extracted from the cavity.
            """
            print(bad_cavity_text)

    def int_fun(self):
        # NEEDS A PROPER GUI
        eta_slider = widgets.FloatSlider(min=0, max=10, step=0.1, continuous_update=False, value=1)
        interactive(self.main_fun, eta=eta_slider)


class OneToneExperiment(InteractiveFunction):
    params = dict([('eta', 1), ('g0', 1e2/1e5)])

    def main_fun(self, eta=1, g0=1e2/1e5, omega_m=8, t=300):
        """"
        ALL FREQUENCIES NORMALISED TO kext (arbitrarily taken to be 1e5)!!!
        eta is the effective coupling efficiency k0/kext
        g0 is the vaccum coupling rate
        omega_m is the mechanical frequency
        t is the temperature (in K)
        """

        kB = 1e-23
        hbar = 1e-34
        nth = (kB*t)/(hbar*omega_m*1e5)  # number of thermal phonons
        omega_m_max = 15
        omega_m_margin = 5

        omegas = np.arange(-omega_m_max-omega_m_margin, omega_m_max+omega_m_margin, 0.01)

        denominator_rsb = -1j*(omegas+omega_m) + (1+eta)/2
        denominator_bsb = -1j*(omegas-omega_m) + (1+eta)/2

        chi_ab = nth*(g0/2)**2*(1/denominator_rsb + 1/denominator_bsb)
        chi_aa = 1/(-1j*omegas + (1+eta)/2 + chi_ab)
        r_opt = np.abs(1 - chi_aa)**2

        r_rsb = nth*(g0/2)**2*np.abs(chi_aa/denominator_rsb)**2
        r_bsb = nth*(g0/2)**2*np.abs(chi_aa/denominator_bsb)**2

        plt.plot(omegas, r_opt)
        plt.plot(omegas, r_rsb, color='xkcd:red', linestyle='--')
        plt.plot(omegas, r_bsb, color='xkcd:blue', linestyle='--')
        plt.ylim(0, 1.01)
        plt.xlim(-omega_m_max-omega_m_margin, omega_m_max+omega_m_margin)
        plt.xlabel('detuning from resonance')
        plt.ylabel('R')

    def int_fun(self):
        # NEEDS A PROPER GUI
        eta_slider = widgets.FloatSlider(min=0, max=10, step=0.1, continuous_update=False, value=1)
        g0_slider = widgets.FloatSlider(min=0, max=2e2/1e5, step=50/1e5, continuous_update=False, value=1e2/1e5)
        omega_m_slider = widgets.FloatSlider(min=1, max=15, step=2, continuous_update=False, value=8)
        t_slider = widgets.FloatSlider(min=1, max=300, step=50, continuous_update=False, value=300)
        interactive(ote.main_fun, eta=eta_slider, g0=g0_slider, omega_m=omega_m_slider, t=t_slider)
