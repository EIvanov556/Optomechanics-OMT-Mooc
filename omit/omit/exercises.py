import numpy as np
from cmath import *
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox

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
    reset_pos = [0, 0, 0, 0]

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

    def reset_ui(self, event):
        """
        Note that the event input seems redundant but is absolutely necessary. Something within the functioning
        of the matplotlib widget's response to action requires it. [MORE THOROUGH EXPLANATION WELCOME IF KNOWN]
        """
        for key in self.widgets.keys():
            key_widget = self.widgets[key][0]
            reset_fun = key_widget.reset_fun
            if reset_fun is not None:
                reset_fun()

    def update_ui(self, val):
        """
        Will go through each of the widgets of self.widgets and, if present, do the update_function of that widget.
        """

        for widget in self.widgets.keys():
            update_fun = self.widgets[widget][0].update_fun
            if update_fun is not None:
                update_fun()

    def init_axes(self):
        """
        If the instance calls upon the creation of axes (that is, self.widgets contains an entry called 'axes',
        this function defines them and places them.
        """
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
                         'plot_coords': [0.05, 0.98-plot_ylength, plot_xlength, plot_ylength],
                         'plot_range': [-self.omega_range, self.omega_range, 0, 1.1]
                     }],
            'axes': [MatplotlibUIWidget(self.init_axes, None, self.reset_axes),
                     {
                         'eta_pos': ['$\eta$', [1, 0, 5, 0.01], [0.75, 1.04-plot_ylength, 0.2, 0.025]]
                     }],
            'explanatory_text': [MatplotlibUIWidget(self.show_text, self.update_text, None),
                                 {
                                     'text_pos': [0.01, 0.03, 0.98, 0.25]
                                 }]}
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
        x, y = self.main_fun(*self.default_vals)
        self.l, = plt.plot(x, y, lw=2, color='xkcd:black')
        plot_coords = self.widgets['plot'][1]['plot_coords']
        self.l.axes.set_position(plot_coords)
        plot_range = self.widgets['plot'][1]['plot_range']
        plt.axis(plot_range)

    def update_plot_fig(self):
        x, y = self.main_fun(*self.current_vals)
        self.l.set_ydata(y)
        self.fig.canvas.draw_idle()

    def show_text(self):
        default_vals = self.default_vals
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

    @property
    def current_vals(self):
        current_vals = ()
        for par in self.widgets['axes'][1].keys():
            current_vals += (self.sliders[par].val,)
        return current_vals

    @property
    def default_vals(self):
        default_vals = ()
        for par in self.widgets['axes'][1].keys():
            default_vals += (self.widgets['axes'][1][par][1][0],)
        return default_vals


class OneToneExperiment(InteractiveFunction):
    """
    An InteractiveFunction for light sent into a cavity, where it interacts with a mechanical degree of freedom, with an
    assumed small coupling. The transmission is then plotted.
    as a function of the ratio between the loss coefficients (kext+k0)/kext
    Everything is given in units of kext.
    """
    reset_pos = [0.85, 0.33, 0.1, 0.04]

    def __init__(self):
        self.window_size = [50, 100, 1000, 800]
        self.omega_m_max = 15
        self.omega_m_margin = 2
        self.omega_range = self.omega_m_max+self.omega_m_margin  # axis will go from - omega_range to +omega_range
        self.axes = (None, )  # initialises the axes

        plot_xlength = 0.65
        plot_ylength = plot_xlength
        self.widgets = {
            'plot': [MatplotlibUIWidget(self.plot_fig, self.update_plot_fig, None),
                     {
                         'plot_coords': [0.05, 0.98-plot_ylength, plot_xlength, plot_ylength],
                         'plot_range': [-self.omega_range, self.omega_range, 0, 1.1]
                     }],
            'axes': [MatplotlibUIWidget(self.init_axes, None, self.reset_axes),
                     {
                         'eta_pos': ['$\eta$', [1, 0, 5, 0.01], [0.75, 1.04-plot_ylength, 0.2, 0.025]],
                         'g0_pos': ['$g_0$', [5, 0, 100, 0.1], [0.75, 1.04-plot_ylength+0.1, 0.2, 0.025]],
                         'omega_m_pos': ['$\Omega_\mathsf{m}$', [5, 1, self.omega_m_max, 1], [0.75, 1.04-plot_ylength+0.2, 0.2, 0.025]],
                         't_pos': ['$T$', [300, 20, 500, 1], [0.75, 1.04-plot_ylength+0.3, 0.2, 0.025]]

                     }]
        }
        super().__init__()

    def main_fun(self, *pars):
        """
        Calculated the reflection coefficient of a Fabry-Perot cavity around its resonance frequency. pars takes (in
        that order), the following arguments
        eta is the effective coupling efficiency k0/kext
        g0 is the vaccum coupling rate
        omega_m is the mechanical frequency
        t is the temperature (in K)
        :param pars: a tuple containing (eta, g0, omega_m, t)
        :return: the frequency range, the optical reflection coefficient, and the r at the side_band freqs
        """
        eta, g0, omega_m, t = pars
        kb = 1e-23
        hbar = 1e-34
        g0 = g0/1e5  # normalised by kext
        nth = (kb*t)/(hbar*omega_m*1e5)  # number of thermal phonons
        omegas = np.arange(-self.omega_m_max-self.omega_m_margin, self.omega_m_max+self.omega_m_margin, 0.01)

        denominator_rsb = -1j*(omegas+omega_m) + (1+eta)/2
        denominator_bsb = -1j*(omegas-omega_m) + (1+eta)/2
        chi_ab = nth*(g0/2)**2*(1/denominator_rsb + 1/denominator_bsb)
        chi_aa = 1/(-1j*omegas + (1+eta)/2 + chi_ab)

        r_opt = np.abs(1 - chi_aa)**2
        r_rsb = nth*(g0/2)**2*np.abs(chi_aa/denominator_rsb)**2
        r_bsb = nth*(g0/2)**2*np.abs(chi_aa/denominator_bsb)**2

        return omegas, r_opt, r_rsb, r_bsb

    def plot_fig(self):
        x, y1, y2, y3 = self.main_fun(*self.default_vals)
        self.axes = plt.plot(x, y1, 'xkcd:black', x, y2, 'xkcd:red', x, y3, 'xkcd:blue', lw=2)
        plot_coords = self.widgets['plot'][1]['plot_coords']
        for ax in self.axes:
            ax.axes.set_position(plot_coords)
        plot_range = self.widgets['plot'][1]['plot_range']
        plt.axis(plot_range)

    def update_plot_fig(self):
        x, *ys = self.main_fun(*self.current_vals)
        for ax, y in zip(self.axes, ys):
            ax.set_ydata(y)
        self.fig.canvas.draw_idle()

    @property
    def current_vals(self):
        current_vals = ()
        for par in self.widgets['axes'][1].keys():
            current_vals += (self.sliders[par].val,)
        return current_vals

    @property
    def default_vals(self):
        default_vals = ()
        for par in self.widgets['axes'][1].keys():
            default_vals += (self.widgets['axes'][1][par][1][0],)
        return default_vals
