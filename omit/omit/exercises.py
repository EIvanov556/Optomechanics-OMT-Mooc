import numpy as np
from cmath import *
import matplotlib.pyplot as plt
import matplotlib.image as mimg
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import os

SCRIPT_DIR = os.path.dirname(__file__)
RESOURCES_DIR = 'exercises_resources'


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
    axcolor = 'xkcd:light pink'
    window_size = [50, 100, 1000, 800]
    widgets = None
    reset_pos = [0, 0, 0, 0]
    window_title = 'Some Title'

    def __init__(self):
        self.fig, self.ax = plt.subplots(num=self.window_title)
        plt.subplots_adjust(left=0.25, bottom=0.25)
        self.fig_manager = plt.get_current_fig_manager()
        self.fig_manager.window.setGeometry(self.window_size[0], 
                                            self.window_size[1],
                                            self.window_size[2], 
                                            self.window_size[3])

        for key in self.widgets.keys():
            key_widget = self.widgets[key][0]
            if key_widget.init_fun is not None:
                key_widget.init_fun()

        self.resetax = plt.axes(self.reset_pos)
        self.button = Button(self.resetax, 'Reset', 
                             color=self.axcolor, hovercolor='xkcd:orange')
        self.button.on_clicked(self.reset_ui)

    def reset_ui(self, event):
        """
        Note that the event input seems redundant but is absolutely necessary. 
        Something within the functioning
        of the matplotlib widget's response to action requires it. 
        [MORE THOROUGH EXPLANATION WELCOME IF KNOWN]
        """
        for key in self.widgets.keys():
            key_widget = self.widgets[key][0]
            reset_fun = key_widget.reset_fun
            if reset_fun is not None:
                reset_fun()

    def update_ui(self, val):
        """
        Will go through each of the widgets of self.widgets and, if present, 
        do the update_function of that widget.
        The value input is obligatory for the syntax but can be a dummy
        """
        for widget in self.widgets.keys():
            update_fun = self.widgets[widget][0].update_fun
            if update_fun is not None:
                update_fun()

    def init_axes(self):
        """
        If the instance calls upon the creation of axes (that is, self.widgets 
        contains an entry called 'axes', this function defines them and places 
        them.)
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

    def add_image(self, image_location):
        """
        Adds an image on the widget (it will be an axis-less plot).
        :param image_location: Give the file address of the image
        :param position: is a vector containing the position information 
        [xmin, xmax, ymin, ymax]
        """
        img = mimg.imread(image_location)
        self.img_ax = plt.subplot2grid((3,3), (0,2), colspan=1, rowspan=1)
        self.img_fig = self.img_ax.imshow(img)
        self.img_fig.axes.axis('off')


class MatplotlibUIWidget:
    """
    A simple Object class to contain all the attributed required by the 
    Interactive Function __init__.
    """
    def __init__(self, init_fun=None, update_fun=None, reset_fun=None):
        self.init_fun = init_fun
        self.update_fun = update_fun
        self.reset_fun = reset_fun


class OpticsOnly(InteractiveFunction):
    """
    An InteractiveFunction for light sent into a cavity, with no other effects 
    given. The transmission is then plotted
    as a function of the ratio between the loss coefficients (kext+k0)/kext
    Everything is given in units of kext.
    """
    reset_pos = [0.85, 0.33, 0.1, 0.04]
    window_title = 'One-Tone, Optics Only, Experiment'

    def __init__(self):
        self.window_size = [50, 100, 700, 700]
        self.omega_range = 5  # axis will go from - omega_range to +omega_range in MHz
    
        self.widgets = {
            'axes': [MatplotlibUIWidget(self.init_axes, 
                                    None, 
                                    self.reset_axes),
                 {
                     'eta_pos': ['$(\kappa_\mathsf{intr}+\kappa_\mathsf{ext})/\kappa_\mathsf{ext}$', 
                                 [2, 1, 5, 0.01], 
                                 [0.55, 0.34, 0.2, 0.025]],
                     'alpha_pos': ['$\kappa_\mathsf{out}/\kappa_\mathsf{in}$',
                                   [0, 0, 5, 0.01], 
                                   [0.55, 0.3, 0.2, 0.025]],
                     'kappa_intr': ['$\kappa_\mathsf{intr}$ (MHz)',
                                   [1, 0, 3, 0.01], 
                                   [0.2, 0.3, 0.1, 0.025]]
                 }],
            'plot': [MatplotlibUIWidget(self.init_fig, 
                                        self.update_fig, 
                                        None),
                     {                        
                         }],
            'cavity_image': [MatplotlibUIWidget(self.init_cavity_image, 
                                                self.update_cavity_image, 
                                                self.reset_cavity_image),
                             {
                                 'names': ('hanger unidirectional',
                                           'hanger bidirectional',
                                           'Fabry-Pérot'),
                                 'rax_pos': [0.675, 0.475, 0.3, 0.12],
    
                                 }],
            'explanatory_text': [MatplotlibUIWidget(self.init_text, 
                                                    self.update_text, 
                                                    None),
                                 {
                                     'text_pos': [0.01, 0.03, 0.98, 0.25]
                                     }],
            'plot_units_radio': [MatplotlibUIWidget(self.init_units_radio, 
                                                    None, 
                                                    self.reset_units_radio),
                                 {
                                         'rax_pos': [0.025, 0.7555, 0.12, 0.12],
                                         'choices': ('linear', 'dB')
                                         }],
        }
        super().__init__()

    def main_fun(self, *pars):
        """
        Calculated the reflection coefficient of a Fabry-Perot cavity around 
        its resonance frequency.
        :param pars: Contains : eta, the ratio (k_intr+k_ext)/k_ext; alpha, the ratio
        k_out/k_in; k_intr
        Frequencies in units of MHz
        :return: the frequency range and the reflection coefficient in this range
        """
        dw = -0.5 #asymmetry detuning factor
        eta = pars[0]
        k_intr = pars[2]
        omegas = np.arange(-self.omega_range, self.omega_range, 2*self.omega_range/100)  
        if self.cavity_choice == 'hanger unidirectional': 
            freq_factor = (eta-1)/k_intr
            susceptibility = (1 - 1j*dw*freq_factor) / ( eta/2 + 1j*omegas*freq_factor )
            s21 = np.abs(1 - susceptibility)
        elif self.cavity_choice == 'hanger bidirectional': 
            freq_factor = (eta-1)/k_intr
            susceptibility = (1 - 1j*dw*freq_factor) / ( eta/2 + 1j*omegas*freq_factor )
            s21 = np.abs(1 - susceptibility/2)
        elif self.cavity_choice == 'Fabry-Pérot':
            alpha = pars[1]
            freq_factor = (eta-1-alpha)/k_intr
            susceptibility = ( np.sqrt(alpha) - 1j*dw*freq_factor ) / ( eta/2 + 1j*omegas*freq_factor )
            s21 = np.abs(susceptibility)
        if self.plot_units=='dB':
            s21 = 20*np.log10(s21)
        return omegas, s21

    def init_fig(self):
        x, y = self.main_fun(*self.default_vals)
        self.plt_ax = plt.subplot2grid((3,3), (0,0), colspan=2, rowspan=2)
        self.l, = self.plt_ax.plot(x, y, lw=2, color='xkcd:black')
        self.set_plot_labels_and_range()

    def set_plot_labels_and_range(self):
        ylabel = '|S21|' if self.plot_units=='linear' else 'S21'
        self.plt_ax.set_ylabel('%s (%s)' %(ylabel, self.plot_units))
        self.plt_ax.set_xlabel('detuning (MHz)')
        self.plt_ax.set_xlim(-self.omega_range, +self.omega_range)
        if self.plot_units=='linear':
            self.plt_ax.set_ylim(0,1.1)
        else:
            self.plt_ax.set_ylim(-30, 0)
        
    def update_fig(self):
        x, y = self.main_fun(*self.current_vals)
        self.set_plot_labels_and_range()
        self.l.set_ydata(y)
        self.fig.canvas.draw_idle()

    def init_text(self):
        default_vals = self.default_vals
        text = self.get_explanatory_text(*default_vals)
        text_pos = self.widgets['explanatory_text'][1]['text_pos']
        # an axis needs to be created for each interactive object
        self.textax = plt.axes(text_pos, facecolor=self.axcolor)  
        text = self.get_explanatory_text(*default_vals)
        self.textbox = TextBox(self.textax, '', initial=text, 
                               color=self.axcolor, hovercolor=self.axcolor,
                               label_pad=0.01)
        
    @staticmethod
    def get_explanatory_text(*pars):
        """
        Gives a short explanatory text about the current selected situation
        :param pars:
        :return: text to be displayed
        """
        eta = pars[0]
        if eta == 0:
            text = 'Extremely overcoupled: \nThe losses are negligible' \
            ' compared to the rate at which light can escape.\n' \
            'All the light is then reflected back out of the cavity'
        elif eta < 1:
            text = 'Overcoupled: \nThe cavity losses are smaller than the' \
            'external coupling rate.'
        elif eta == 1:
            text = 'Critical coupling:\nThe rate at which the light enters' \
            ' the cavity matches exactly the decay rate.\n' \
            'At resonance, all photons will be dissipated'
        else:
            text = 'Undercoupled:\nThis is the so-called "bad cavity" regime.'\
            'All excitations die out ' \
            'faster \nthan they can be extracted from the cavity.'
        return text
    
    def update_text(self):
        new_text = self.get_explanatory_text(*self.current_vals)
        self.textbox.set_val(new_text)

    def init_cavity_image(self):
        self.update_cavity_image()
        
        imgrax_pos = self.widgets['cavity_image'][1]['rax_pos']
        imgrax = plt.axes(imgrax_pos, facecolor=self.axcolor)
        imgrax_choices = self.widgets['cavity_image'][1]['names']
        self.radio_image = RadioButtons(imgrax, imgrax_choices)
        self.radio_image.on_clicked(self.update_ui)
    
    def update_cavity_image(self):
        img_name = self.cavity_choice
        img_loc = os.path.join(SCRIPT_DIR, RESOURCES_DIR, img_name+'.png')
        self.add_image(img_loc)
    
    @property
    def cavity_choice(self):
        try:
            return self.radio_image.value_selected
        except AttributeError:
            return self.widgets['cavity_image'][1]['names'][0]
        
    def reset_cavity_image(self):
        self.radio_image.set_active(0)
        #self.update_ui(None)
    
    def init_units_radio(self):
        rax_pos = self.widgets['plot_units_radio'][1]['rax_pos']
        rax = plt.axes(rax_pos, facecolor=self.axcolor)
        rax_choices = self.widgets['plot_units_radio'][1]['choices']
        self.radio_units = RadioButtons(rax, rax_choices)
        self.radio_units.on_clicked(self.update_ui)

    @property
    def plot_units(self):
        try:
            return self.radio_units.value_selected
        #default choice if the radio button does not exist
        except AttributeError:
            return self.widgets['plot_units_radio'][1]['choices'][0]
        
    def reset_units_radio(self):
        self.radio_units.set_active(0)
        
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
    An InteractiveFunction for light sent into a cavity, where it interacts 
    with a mechanical degree of freedom, with an
    assumed small coupling. The transmission is then plotted.
    as a function of the ratio between the loss coefficients (kext+k0)/kext
    Everything is given in units of kext.
    """
    reset_pos = [0.85, 0.33, 0.1, 0.04]
    window_title = 'One-Tone Optomechanics Experiment'

    def __init__(self):
        self.window_size = [50, 100, 1000, 800]
        self.omega_m_max = 15
        self.omega_m_margin = 2
        self.omega_range = self.omega_m_max+self.omega_m_margin  # axis will go from - omega_range to +omega_range
        self.axes = (None, )  # initialises the axes

        plot_xlength = 0.65
        plot_ylength = plot_xlength
        self.widgets = {
            'plot': [MatplotlibUIWidget(self.plot_fig, self.update_plot_fig, 
                                        None),
                     {
                         'plot_coords': [0.05, 0.98-plot_ylength, 
                                         plot_xlength, plot_ylength],
                         'plot_range': [-self.omega_range, self.omega_range, 
                                        0, 1.1]
                     }],
            'axes': [MatplotlibUIWidget(self.init_axes, None, self.reset_axes),
                     {
                         'eta_pos': ['$\eta$', [1, 0, 5, 0.01], 
                                     [0.75, 1.04-plot_ylength, 0.2, 0.025]],
                         'g0_pos': ['$g_0$', [35, 0, 100, 0.1], 
                                    [0.75, 1.04-plot_ylength+0.1, 0.2, 0.025]],
                         'omega_m_pos':       ['$\Omega_\mathsf{m}$', 
                                               [5, 1, self.omega_m_max, 1], 
                                               [0.75, 1.04-plot_ylength+0.2, 0.2, 0.025]],
                         't_pos': ['$T$', [300, 20, 500, 1], 
                                   [0.75, 1.04-plot_ylength+0.3, 0.2, 0.025]]

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
        omegas = np.arange(-self.omega_m_max-self.omega_m_margin, 
                           self.omega_m_max+self.omega_m_margin, 0.01)

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
