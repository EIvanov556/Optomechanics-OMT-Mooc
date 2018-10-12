import numpy as np
from cmath import *
from ipywidgets import interact, interactive, fixed, interact_manual
# https://ipywidgets.readthedocs.io/en/stable/examples/Using%20Interact.html
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
#import ipywidgets as widgets


class InteractiveFunction:
    """
    Interactive GUI plotting Object.
    Is given a set of widgets that are to be shown on the interactive UI.
    To pass this onto the class, a dictionary called widgets must be created in the child, containing the information
    of the widgets to the plotted. This is organised in the following way:

    widgets = {
        'widget_name': [UIWidget, {*args}]}
    The UI widget is an object that must contain as attributes an init_fun, called at the initialisation, and a
    reset_fun, used when the reset button is pressed.
    The live updates of the widget are to be handled in the child class itself.
    """
    axcolor = 'lightgoldenrodyellow'
    window_size = [50, 100, 1000, 800]
    widgets = None
    reset_pos = [0,0,0,0]

    def __init__(self):
        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)
        self.fig_manager = plt.get_current_fig_manager()
        self.fig_manager.window.setGeometry(self.window_size[0], self.window_size[1],
                                            self.window_size[2], self.window_size[3])

        for key in self.widgets.keys():
            key_widget = self.widgets[key][0]
            key_widget.init_fun()

        self.resetax = plt.axes(self.reset_pos)
        self.button = Button(self.resetax, 'Reset', color=self.axcolor, hovercolor='0.975')
        self.button.on_clicked(self.reset_ui)

    def reset_ui(self,event):
        for key in self.widgets.keys():
            key_widget = self.widgets[key][0]
            reset_fun = key_widget.reset_fun
            if reset_fun is not None:
                reset_fun()


class MatplotlibUIWidget:
    """
    A simple Object class to contain all the attributed required by the Interactive Function __init__.
    """
    def __init__(self, init_fun=None, update_fun=None, reset_fun=None):
        self.init_fun = init_fun
        self.update_fun = update_fun
        self.reset_fun = reset_fun


class OpticsOnly(InteractiveFunction):
    """
    An InteractiveFunction for light sent into a cavity, with no other effects given. The transmission is then plotted
    as a function of the ratio between the loss coefficients (kext+k0)/kext
    Everything is given in units of kext.
    """
    reset_pos = [0.85, 0.33, 0.1, 0.04]

    def __init__(self):
        self.window_size = [50, 100, 1000, 800]
        self.omega_range = 5  # axis will go from - omega_range to +omega_range

        plot_xlength = 0.65
        plot_ylength = plot_xlength
        self.widgets = {
            'plot': [MatplotlibUIWidget(self.plot_fig, self.update_plot_fig, None),
                     {
                         'plot_coords': [0.05,1-plot_ylength-0.02,plot_xlength,plot_ylength],
                         'plot_range': [-self.omega_range, self.omega_range, 0, 1.1]
                     }],
            'axes': [MatplotlibUIWidget(self.init_axes, None, self.reset_axes),
                     {
                         'eta_pos': ['eta', [1, 0, 5, 0.1], [0.75, 1-plot_ylength-0.02+0.1-0.04, 0.2, 0.025]]
                     }],
            'explanatory_text': [MatplotlibUIWidget(self.show_text, self.update_text, None),
                                 {
                                     'text_pos': [0.01, 1 - 0.65 - 0.02 - 0.3, 1 - 0.02, 0.25]
                                 }]}
        self.params = {
            'eta': [1, 0, 5, 0.1]}
        super().__init__()

    def main_fun(self, *pars):
        """
        Calculated the reflection coefficient of a Fabry-Perot cavity around its resonance frequency.
        :param pars: Contains only one item: eta, the ratio k0+kext/kext
        :return: the frequency range and the reflection coefficient in this range
        """
        omegas = np.arange(-self.omega_range, self.omega_range, 0.01)  # in units of kext
        susceptibility = 1 / (((1 + pars[0]) / 2) + 1j * omegas)
        r = np.abs(1 - susceptibility) ** 2
        return omegas, r

    def plot_fig(self):
        default_vals = ()
        for par in self.params.keys():
            default_vals += (self.params[par][0],)
        x, y = self.main_fun(*default_vals)
        self.l, = plt.plot(x, y, lw=2, color='xkcd:black')
        plot_coords = self.widgets['plot'][1]['plot_coords']
        self.l.axes.set_position(plot_coords)
        plot_range = self.widgets['plot'][1]['plot_range']
        plt.axis(plot_range)

    def update_plot_fig(self):
        x, y = self.main_fun(*self.current_vals)
        self.l.set_ydata(y)
        self.fig.canvas.draw_idle()

    def init_axes(self):
        self.slider_axes = dict()
        self.sliders = dict()
        for i, par in enumerate(self.widgets['axes'][1].keys()):
            slider_positions = self.widgets['axes'][1][par][2]
            self.slider_axes[par] = plt.axes(slider_positions,
                                             facecolor=self.axcolor)

            slider_tag = self.widgets['axes'][1][par][0]
            slider_pars = self.widgets['axes'][1][par][1]
            self.sliders[par] = Slider(self.slider_axes[par], slider_tag,
                                       slider_pars[1],
                                       slider_pars[2],
                                       valinit=slider_pars[0],
                                       valstep=slider_pars[3])
            self.sliders[par].on_changed(self.update_ui)

    def reset_axes(self):
        for ax in self.widgets['axes'][1].keys():
            self.sliders[ax].reset()

    def show_text(self):
        default_vals = ()
        for par in self.params.keys():
            default_vals += (self.params[par][0],)
        text = self.get_explanatory_text(*default_vals)

        text_pos = self.widgets['explanatory_text'][1]['text_pos']
        self.textax = plt.axes(text_pos, facecolor=self.axcolor)  # an axis needs to be created for each interactive object

        text = self.get_explanatory_text(*default_vals)
        self.textbox = TextBox(self.textax, '', initial=text, color=self.axcolor, hovercolor=self.axcolor,
                               label_pad=0.01)

    def update_text(self):
        new_text = self.get_explanatory_text(*self.current_vals)
        self.textbox.set_val(new_text)

    @staticmethod
    def get_explanatory_text(*pars):
        """
        Gives a short explanatory text about the current selected situation
        :param pars:
        :return: text to be displayed
        """
        eta = pars[0]
        if eta == 0:
            text = 'Extremely overcoupled: \nThe losses are negligible compared to the rate at which light can escape.\n' \
                   'All the light is then reflected back out of the cavity'
        elif eta < 1:
            text = 'Overcoupled: \nThe cavity losses are smaller than the external coupling rate.'
        elif eta == 1:
            text = 'Critical coupling:\nThe rate at which the light enters the cavity matches exactly the decay rate.\n' \
                   'At resonance, all photons will be dissipated'
        else:
            text = 'Undercoupled:\nThis is the so-called "bad cavity" regime. All excitations die out ' \
                   'faster than they can be extracted from the cavity.'
        return text

    def update_ui(self, val):
        for widget in self.widgets.keys():
            update_fun = self.widgets[widget][0].update_fun
            if update_fun is not None:
                update_fun()

    @property
    def current_vals(self):
        current_vals = ()
        for par in self.widgets['axes'][1].keys():
            current_vals += (self.sliders[par].val,)
        return current_vals


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
