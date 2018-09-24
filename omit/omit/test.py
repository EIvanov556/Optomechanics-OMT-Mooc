import numpy as np
from cmath import *
from ipywidgets import interact, interactive, fixed, interact_manual
#https://ipywidgets.readthedocs.io/en/stable/examples/Using%20Interact.html
import ipywidgets as widgets


def optics_only(coupling = 1):
    # coupling is k0/kext
    omegas = np.arange(-5, 5, 0.01)  # in units of kext
    susceptibility = 1 / (((1+coupling)/2)+1j*omegas)
    R = np.abs(1-susceptibility)**2

    plot(omegas, R)
    ylim(0, 1.05)
    xlabel('detuning from resonance')
    ylabel('R')
    show()

    if coupling == 0:
        print('Extremely overcoupled: \nThe losses are negligible compared to the rate at which light can escape. \nAll the light is then reflected back out of the cavity')
    elif coupling < 1:
        print('Overcoupled: \nThe cavity losses are smaller than the external coupling rate.')
    elif coupling == 1:
        print('Critical coupling:\nThe rate at which the light enters the cavity matches exactly the decay rate.\nAt resonance, all photons will be dissipated')
    else:
        bad_cavity_text = """Undercoupled:\nThis is the so-called "bad cavity" regime. All excitations die out faster than they can be extracted from the cavity.
        """
        print(bad_cavity_text)
