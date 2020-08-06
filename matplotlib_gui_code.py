#! /usr/bin/env python
#
"""
matplotlib_gui_code.py is an interactive plotting tool using matplotlib.

This is intended to be a front-end for a number of MATPLOTLIB functions similar
in general usage to 'xmgrace" (although not intending to duplicate all the
xmgrace funcionality...).  From the interface one is able to read in data and
interactively change the plot properties (axes, symbols, lines, and so on).


** NOTE:  for this program to work properly it needs to be in a directory
listed in the $PATH environment variable, and also in the $PYTHON_PATH
environment variable.  These allow the code to be called from the command
line and from inside python, respectively.

** Note:  this code uses python 3, it will not work with python 2.


Various of the functions provided here are not available in the normal
.plot() interface of matplotlib.

The code uses the more general matplotlib functions and not the pyplot
interface for this purpose.

The code needs matplotlib and numpy to function.  Numpy loadtxt is used to
read in data.

The code uses TKinter for the widget generation.  It needs the "Tk" backend.
One can change the backend right at the start of the code if another such is
needed.

The code assumes that the "times new roman" font is installed along with the
default matplotlib fonts.  If this is not the case, one will get a fall-back
font instead (usually "DejaVe Sans").  One can install the Times New Roman
font if this is not already on the system.  Due to preferences of the author
the Times New Roman font is the default.  If one wishes to change this, search
for 'times new roman' and replace the string by, say, 'sans-serif'.  There are
commented out lines to make 'sans-serif' the default font, which can be
uncommented and used to replace the instances of 'times new roman' font as
the default, should that be needed.

************************************************************

Use from inside Python:

  While the code was written to read data from files as a stand-alone
interface, one can also call the routines from inside python and add sets
from the python prompt.  The commands needed are as follows, assuming that
data arrays "xvalues" and "yvalues" already exist within the Python session:

>>>> import matplotlib_gui_code
>>>> root, myplot = matplotlib_gui_code.startup()
>>>> myplot.add_set(xvalues, yvalues)

(repeat for multiple sets, if needed)

>>>> myplot.make_plot()

If one wishes to do all the steps manually the following somewhat longer
set of commands can be used:

(1) import the required packages tkinter and matplotlib_gui_code

>>> import tkinter as Tk
>>> import matplotlib_gui_code

The next step defines a window in tkinter to use for the plotting:

>>> root = Tk.Tk()

This produces a window with nothing in it.

(2) Generate x and y data as numpy arrays as needed

Assume that "xvalues" and "yvalues" hold the data values for a plot.  These
need to be numpy arrays and should not be lists.

One simple example would be

>>>> xvalues = numpy.arange(1, 101, 1)
>>>> yvalues = numpy.log10(xvalues * xvalues)

(3) Define the plotting object

>>> myplot=matplotlib_gui_code.PlotGUI(root)

This fills in the buttons/functions in the window.  The functionality is
then avaialble.

(4) Add a set or sets to the plot.  Any number of sets up a limit (currently
100) can be added.  See the add_set subroutine doc string for the possible
use of error bars in the set.  The value self.max_sets determines the maximum
number of sets that can be handled.

>>> myplot.add_set(xvalues, yvalues)

(5) Tell the code to plot or replot.  This call can be used as many times
as required.

>>> myplot.make_plot()

The add_set routine takes as input the x and y data values and makes a new
set.  Hence one can mix reading in values and generating values inside python
as needed by following the above set of steps.  There are additional
parameter values one can pass to this routine, for error values in the data
points.

When running in a script rather than at the python interpreter, one needs
to add

root.mainloop()

after making the call to make_plot with the sets, to keep the window in place.
Otherwise the window will disappear after the make_plot call if nothing else
is being done in the script.

"""
import math
from copy import deepcopy
import sys
import os
import bisect
import tkinter as Tk
import tkinter.ttk
import tkinter.filedialog
import tkinter.simpledialog
import tkinter.messagebox
from tkinter.colorchooser import askcolor
from tkinter.scrolledtext import ScrolledText
import numpy
from numpy.polynomial import polynomial, legendre, laguerre, chebyshev
from scipy.interpolate import UnivariateSpline
import matplotlib
import matplotlib.lines as mlines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Ellipse, FancyArrow
from matplotlib.ticker import MultipleLocator

# The following are "global" variables with line/marker information from
# matplotlib.  These are used, but not changed, in the code in more than one
# place hence I am using single variables here for the values.
matplotlib_symbol_list = [None, ".", ", ", "o", "v", "^", "<", ">", "1",
                          "2", "3", "4", "8", "s", "p", "P", "*", "h", "H",
                          "+", "x", "X", "D", "d", "|", "_", "histogram"]
matplotlib_symbol_name_list = ['None', 'point', 'pixel', 'circle',
                               'triangle down', 'triangle up',
                               'triangle left', 'triangle right', 'tri_down',
                               'tri_up', 'tri_left', 'tri_right', 'octagon',
                               'square', 'pentagon', 'plus (filled)', 'star',
                               'hexagon1', 'hexagon2', 'plus', 'x',
                               'x (filled)', 'diamond', 'thin_diamond',
                               'vline', 'hline', 'histogram']
matplotlib_line_list = ['-', '--', '-.', ':', None]
matplotlib_line_name_list = ['solid', 'dashed', 'dashdot', 'dotted', 'None']
matplotlib.use('TkAgg')
# define a background colour for windows
BGCOL = '#F8F8FF'


def startup():
    """
    Startup.py is a wrapper for starting the plot tool.

    This is a wrapper for making a plot from the Python command prompt.
    Assuming that numpy data arrays xvalues and yvalues exist one can use the
    commands

    >>> import matplotlib_gui_code
    >>> root, myplot = matplotlib_gui_code.startup()
    >>> myplot.add_set(xvalues, yvalues)
    >>> myplot.make_plot()

    to bring up the window and plot the data set from inside Python.  One can
    then add other data sets and use the functionality of the tool.

    Parameters
    ----------
        None.

    Returns
    -------
        newroot :   The Tkinter window class variable for the plot window.

        plotobject :  The matplotlib_gui_code plot GUI object variable.  This is the
                      variable one uses for making plots.

    """
    # Make a Tkinter window
    newroot = Tk.Tk()
    newroot.title("Plotting Tool")
    plotobject = PlotGUI(newroot)
    return newroot, plotobject


def save_values(xvalues, yvalues, labelstring=None):
    """
    Save plot values to an ascii output file.

    Parameters
    ----------
        xvalues :  A numpy array of float or integer values, the x plot values

        yvalues :  A numpy array of float or integer values, the y plot values

        labelstring :  An optional string to write to the top of the file.
                       If the value is None, nothing is written.

    Returns
    -------
        Nothing
    """
    outfilename = tkinter.filedialog.asksaveasfilename()
    outfile = open(outfilename, 'w')
    if labelstring is not None:
        print(labelstring, file=outfile)
    for loop in range(len(xvalues)):
        print('%20.6g %20.6g' % (xvalues[loop], yvalues[loop]), file=outfile)
    outfile.close()


def save_png_figure(fig):
    """
    Save the current plot as a PNG file.

    Parameters
    ----------
        fig :   A matplotlib Figure instance, the figure to be
                written out

    Returns
    -------
        Nothing

    """
    outfile = tkinter.filedialog.asksaveasfilename(filetypes=[('PNG',
                                                               '.png')])
    if isinstance(outfile, type('string')):
        s1 = outfile.split('.')
    if 'png' not in s1[-1]:
        outfile = outfile + '.png'
    fig.savefig(outfile, format="PNG")


def save_ps_figure(fig):
    """
    Save the current plot as a Postscript file.

    Parameters
    ----------
        fig :  A matplotlib Figure instance, the figure to be
               written out

    Returns
    -------
    Nothing

    """
    outfile = tkinter.filedialog.asksaveasfilename(filetypes=[('PS',
                                                               '.ps')])
    if isinstance(outfile, type('string')):
        s1 = outfile.split('.')
        if 'ps' not in s1[-1]:
            outfile = outfile + '.ps'
        fig.savefig(outfile, format="PS")


def hybrid_transform(datavalues):
    """
    Apply an IRAF-style hybrid log transformation.

    This routine transforms a set of input data values for the hybrid
    log plot:

    Numbers with absolute value > 10 have the sign times the base
    10 log of the absolute value i.e. 100 -> 2.0, -100 -> -2.0.

    Numbers with absolute value <= 10 are scaled down by a factor of
    10 to lie in the range from +1.0 to -1.0 i.e 8 -> 0.8 and -6 -> -0.6.

    Parameter
    ---------
        datavalues :   A numpy array of numbers (real or integer), or a
                       single float value

    Returns
    -------
        newdatavalues :   A numpy array of the same length as datavalues,
                          with the transformed numbers, or a single float
                          value with the transformed mumber

    """
    try:
        newdatavalues = datavalues.astype(numpy.float32)
        inds = numpy.where(numpy.abs(datavalues) < 10.)
        newdatavalues[inds] = datavalues[inds]/10.
        inds = numpy.where(datavalues >= 10.)
        newdatavalues[inds] = numpy.log10(newdatavalues[inds])
        inds = numpy.where(datavalues <= -10.)
        newdatavalues[inds] = -1.0*numpy.log10(numpy.abs(newdatavalues[inds]))
    except:
        if abs(datavalues) >= 10.:
            newdatavalues = math.log10(abs(datavalues))
            if datavalues < 0.:
                newdatavalues = -1.*newdatavalues
        else:
            newdatavalues = datavalues/10.
    return newdatavalues


def inverse_hybrid_transform(value):
    """
    Transform back from the IRAF-style hybrid log values.

    This takes the hybrid log value and transforms it back to the
    actual value.  That value is returned.  Unlike the hybrid_transform
    function, this works on single values not a numpy array.  That is because
    one is not going to have a data array in the hybrid transformation form.

    Parameters
    ----------
        value :  A real number to be transformed back from the hybrid
                 log scaling

    Returns
    -------
        newvalue :   The associated value that maps to the given hybrid
                     log value.

    """
    if value < 0.:
        workvalue = -1.0*value
        sign = -1.0
    else:
        workvalue = value
        sign = +1.0
    if workvalue < 1.0:
        newvalue = 10.*workvalue
    else:
        newvalue = 10.**workvalue
    newvalue = sign * newvalue
    return newvalue


def generate_labels(datavalues):
    """
    Create labels for the hybrid log case.

    This routine makes hybrid log labels for the range of values in a
    given data set: it assumes that the range is large enough that
    integer ticks between  10 and +10 are suitable.  There is no reason
    to use the hybridlog option if the range is smaller than this, or
    indeed if the range is less than a couple of orders of magnitude.

    Parameters
    ----------
        datavalues :     A one-dimensional numpy array of float data values.

    Returns
    -------
        rangeout :  A list of the positions for the tick marks (floats).

        ticksout :  A list of the tick labels (strings).

    Both the return values are used with set_xticks/set_xticklabels
    (main matplotlib) or xticks (pyplot) or the y axis equivalents.
    """
    ticklabels = []
    datamax = numpy.max(datavalues)
    datamin = numpy.min(datavalues)
    baserange = numpy.arange(-10, 12, 2)/10.0
    subvalues = numpy.arange(2, 11)
    subvalues = numpy.log10(subvalues)
    n1 = int(datamin)
    if datamin < 0.:
        n1 = n1 - 1
    n2 = int(datamax)
    if n1 < -1:
        for loop in range(-1, n1-1, -1):
            for n1 in range(len(subvalues)-1, -1, -1):
                baserange = numpy.insert(baserange, 0, loop-subvalues[n1])
    if n2 > 1:
        for loop in range(1, n2+1):
            for n1 in range(len(subvalues)):
                baserange = numpy.append(baserange, loop+subvalues[n1])
    ticklabels = []
    tickzero = None
    for loop in range(len(baserange)):
        if (tickzero is None) and (datamin < baserange[loop]) and \
           (datamax > baserange[loop]):
            tickzero = baserange[loop]
    for loop in range(len(baserange)):
        if baserange[loop] < datamin:
            baserange[loop] = tickzero
            ticklabels.append('')
        elif baserange[loop] > datamax:
            baserange[loop] = tickzero
            ticklabels.append('')
        else:
            if abs(baserange[loop]) < 1.1:
                ticklabels.append('%.0f' % (baserange[loop]*10.0))
            else:
                if (baserange[loop] > 1.1) and (math.floor(
                        baserange[loop]) == baserange[loop]):
                    ticklabels.append(r'$10^%d$' % (baserange[loop]))
                elif (baserange[loop] < -1.1) and (math.floor(
                        baserange[loop]) == baserange[loop]):
                    ticklabels.append(r'$-10^%d$' % (abs(baserange[loop])))
                else:
                    ticklabels.append('')
    rangeout = []
    ticksout = []
    for loop, value in enumerate(baserange):
        if not math.isnan(value):
            rangeout.append(value)
            ticksout.append(ticklabels[loop])
    return rangeout, ticksout


def parse_text(text):
    """
    This routine parses a list of text lines to extract the numerical
    values for a set, used with the text entry option.

    Parameters
    ----------
    text : a list of strings

    Returns
    -------
    xvalues : a list of float values, the x data values for a set

    dxvalues1 : a list of float values, the lower x error values for a set

    dxvalues2 : a list of float values, the upper x error values for a set

    yvalues : a list of float values, the y data values for a set

    dyvalues1 : a list of float values, the lower y error values for a set

    dyvalues2 : a list of float values, the upper y error values for a set

    errorflag : boolean value, flags whether the uncertaintes are defined

    """
    xvalues = []
    dxvalues1 = []
    dxvalues2 = []
    yvalues = []
    dyvalues1 = []
    dyvalues2 = []
    lines = text.split('\n')
    for line in lines:
        if '#' in line:
            pass
        else:
            values = line.split()
            numbers = []
            errorflag = False
            for loop in range(len(values)):
                try:
                    v1 = float(values[loop])
                    numbers.append(v1)
                except:
                    pass
            if len(numbers) == 2:
                xvalues.append(numbers[0])
                yvalues.append(numbers[1])
                dxvalues1.append(0.0)
                dxvalues2.append(0.0)
                dyvalues1.append(0.0)
                dyvalues2.append(0.0)
            elif len(numbers) == 4:
                xvalues.append(numbers[0])
                yvalues.append(numbers[2])
                dxvalues1.append(numbers[1])
                dxvalues2.append(numbers[1])
                dyvalues1.append(numbers[3])
                dyvalues2.append(numbers[3])
                errorflag = True
            elif len(numbers) == 6:
                xvalues.append(numbers[0])
                yvalues.append(numbers[3])
                dxvalues1.append(numbers[1])
                dxvalues2.append(numbers[2])
                dyvalues1.append(numbers[4])
                dyvalues2.append(numbers[5])
                errorflag = True
            elif len(numbers) > 2:
                xvalues.append(numbers[0])
                yvalues.append(numbers[1])
                dxvalues1.append(0.0)
                dxvalues2.append(0.0)
                dyvalues1.append(0.0)
                dyvalues2.append(0.0)
    return xvalues, dxvalues1, dxvalues2, yvalues, dyvalues1, \
        dyvalues2, errorflag


def slope_calculation(xdata, ydata, yerrors=None):
    """
    Calculate a standard least-squares linear fit to input data.

    This is for the case where the standard library function fails, as is
    found to be the case in some instances.

    The treatment is as in "Numerical Recipes in C" (second edition) section
    15.2.

    Parameters
    ----------

    xdata:  A one-dimensional numpy float or integer array of the x values

    ydata:  A one-dimensional numpy float or integer array of the y values,
            which must be the same length as the xdata array

    yerrors:  An optional one-dimensional numpy float array of the y value
              uncertainties, which must be the same length as the xdata array
              if defined.  If not defined, the values are set to a constant.
              The error values must be strictly positive, non-zero.

    Returns
    -------

    slope:           A float value, the best fit slope

    intercept:       A float value, the best fit intercept

    slopeerror:      A float value, the uncertainty in the best fit slope
                     (standard deviation estimate)

    intercepterror:  A float value, the uncertainty in the best fit intercept
                     (standard deviation estimate)

    covariance:      A float value, the covariance between the fit parameters

    correlation:     A float value, the correlation coefficient of the fit

    """
    if yerrors is None:
        yerrors = ydata*0. + 1.
    else:
        inds = numpy.where(yerrors > 0.)
        mean1 = numpy.mean(yerrors[inds])
        inds = numpy.where(yerrors <= 0.)
        yerrors[inds] = mean1
    if (xdata.shape != ydata.shape) or (len(xdata.shape) > 1) or \
       (xdata.shape != yerrors.shape):
        return None, None, None, None, None, None
    invvariance = 1./(yerrors*yerrors)
    sum1 = numpy.sum(invvariance)
    sum2 = numpy.sum(xdata*invvariance)
    sum3 = numpy.sum(ydata*invvariance)
    xmean = sum2/sum1
    t1 = (xdata - xmean) / yerrors
    sum4 = numpy.sum(t1*t1)
    slope = numpy.sum(t1*ydata/yerrors)/sum4
    intercept = (sum3 - (slope*sum2))/sum1
    interceptvariance = (1. + (sum2*sum2)/(sum1*sum4))/sum1
    slopevariance = 1./sum4
    covariance = -1.*sum2/(sum1*sum4)
    correlation = covariance/numpy.sqrt(slopevariance*interceptvariance)
    return slope, intercept, math.sqrt(slopevariance), \
        math.sqrt(interceptvariance), covariance, correlation


def list_fitpars(fit_type, fit_order, fitpars):
    """
    This code writes out the fit parameters to the fit_values.txt file.

    Parameters
    ----------
    fit_type : integer variable, flags the function used in the fit

    fit_order : integer variable, gives the order of the fit function

    fitpars : numpy float array, the fit parameters

    Returns
    -------
    None.

    """
    outfile = open('fit_values.txt', 'a')
    print('Order %d %s polynomial fit:' % (fit_order, fit_type), file=outfile)
    print('Parameter    Value', file=outfile)
    for loop in range(len(fitpars)):
        print('%3d %f' % (loop, fitpars[loop]), file=outfile)
    print(' ', file=outfile)
    outfile.close()


class PlotGUI(Tk.Frame):
    """
    PlotGUI is the object class for the plotting functions.

    The PlotGUI class is the object for the plotting window and all assocated
    variables.  Typical usage on the python command line would be

    root = Tk.Tk()
    plotobject = PlotGUI(root)

    The plotobject would then be used for all the plot functionality.

    Parameters
    ----------
       Tk.Frame :   A Tkinter root or Toplevel variable that is parent to the
                    GUI window.
    """

    def __init__(self, parent=None, **args):
        """
        Initialize variables for the class and start the widget.

        This routine sets a few variables for the interface and calls the
        main widget function.

        Parameters
        ----------
            parent      An optional parameter giving the parent TK root window
                        name.  If it is None the window is not created.  This
                        allows one to use the code non-interactively,
                        although right now that option has not been developed.

            **args      A possible list of additional arguments.  This is
                        currently not used.

        Returns
        -------
            No value is returned by this routine.

        """
        if sys.version_info[0] == 2:
            tkinter.messagebox.showinfo(
                "Error",
                "The code requires Python version 3.")
            return
        # The following are the default colours and symbols for 10 sets.  If
        # more than 10 sets are read in, the values repeat but the symbol sizes
        # are changed.
        self.colourset = ['black', 'blue', 'forestgreen', 'orange', 'red',
                          'cyan', 'lime', 'brown', 'violet', 'grey',
                          'select']
        self.altcolourset = ['none', 'black', 'blue', 'forestgreen', 'orange',
                             'red', 'cyan', 'lime', 'brown', 'violet',
                             'grey', 'select']
        self.markerset = ['o', 's', 'v', '^', '<', '>', 'p', '*', '+', 'x']
        # This variable marks the number of active data sets.
        self.nsets = 0
        # The code allows up to self.max_sets data sets.  The value is large
        # enough that normally one will not have all the slots taken.
        self.max_sets = 100
        # Ditto for labels
        self.max_labels = 100
        # Ditto for lines, ellipses, boxes
        self.max_ellipses = 100
        self.max_boxes = 100
        self.max_lines = 100
        self.max_vectors = 100
        self.equal_aspect = [False, ]
        # The following two variables hold the data values and some associated
        # information.
        self.xdata = []
        self.ydata = []
        self.original_range = []
        # The set_properties variable holds the set symbol and other such
        # parameters that affect how the values are plotted.
        self.set_properties = []
        # The plot_range value holds the current x and y axis range values.
        # These are the default values when no data is yet available.
        self.plot_range = [[0.0, 1.0, 0.0, 1.0], ]
        # The title string for the plot is kept in self.title.
        self.title = [' ', ]
        # For labels, define a variable to hold the position and the labels.
        self.plot_labels = []
        for loop in range(self.max_labels):
            self.plot_labels.append({'xposition': None, 'yposition': None,
                                     'labelstring': '', 'plot': 1,
                                     'colour': 'black', 'size': 12,
                                     'font': 'times new roman',
                                     # 'font': 'sans-serif',
                                     'fontweight': 'normal'})
        # Similar thing for other drawing objects
        self.plot_lines = []
        for loop in range(self.max_lines):
            self.plot_lines.append({'xstart': None, 'ystart': None,
                                    'xend': None, 'yend': None, 'plot': 1,
                                    'line_type': 'solid', 'line_colour':
                                    'black', 'line_thickness': 1.0})
        self.plot_vectors = []
        for loop in range(self.max_vectors):
            self.plot_vectors.append({'xstart': None, 'ystart': None,
                                      'xend': None, 'yend': None,
                                      'delx': None, 'dely': None,
                                      'plot': 1, 'line_type': 'solid',
                                      'line_colour': 'black',
                                      'line_thickness': 1.0, 'fill': True,
                                      'fill_colour': 'black'})
        self.plot_ellipses = []
        for loop in range(self.max_ellipses):
            self.plot_ellipses.append({'xposition': None, 'yposition': None,
                                       'major': None, 'minor': None,
                                       'rotation': 0.0, 'plot': 1,
                                       'line_type': 'solid', 'line_colour':
                                       'black', 'line_thickness': 1.0,
                                       'fill_colour': 'none'})
        self.plot_boxes = []
        for loop in range(self.max_boxes):
            self.plot_boxes.append({'xstart': None, 'ystart': None,
                                    'xend': None, 'yend': None,
                                    'rotation': 0.0, 'plot': 1,
                                    'line_type': 'solid', 'line_colour':
                                    'black', 'line_thickness': 1.0,
                                    'fill_colour': 'none'})
        # flags for setting the different objects, variables for the number
        # of objects of each type
        self.label_flag = False
        self.box_flag = False
        self.ellipse_flag = False
        self.line_flag = False
        self.vector_flag = False
        self.number_of_labels = 0
        self.number_of_lines = 0
        self.number_of_boxes = 0
        self.number_of_ellipses = 0
        self.number_of_vectors = 0
        self.positions = []
        self.datavalues = None
        self.legend_labels = [None, ]
        self.legend_handles = [None, ]
        self.legend_variable = [None, ]
        self.legend_frame = [None, ]
        self.legend_options = [None, ]
        self.legend_position = [None, ]
        self.legend_user_position = [None, ]
        # Load the set_properties variable with the parameters.  For the
        # xdata and ydata variables, load "None".
        for loop in range(self.max_sets):
            self.set_properties.append({
                'symbol': None, 'symbolsize': 4.0,
                'linestyle': 'None', 'linewidth': 1.0, 'colour': 'black',
                'label': '', 'xmin': 0.0, 'xmax': 1.0, 'ymin': 0.0,
                'ymax': 1.0, 'display': True, 'errors': False,
                'legend': True, 'plot': 1})
            self.xdata.append(None)
            self.ydata.append(None)
            self.original_range.append(True)
        # Axis parameter variables
        self.xparameters = [{
            'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
            'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
            'invert': 0, 'hide': 0, 'hideticks': 0,
            'hidelabels': 0, 'hybridlog': 0,
            'inverseticks': 0, 'ticklength': 6,
            'bothticks': 0, 'minorticks': 0, 'oppositeaxis': 0}, ]
        self.yparameters = [{
            'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
            'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
            'invert': 0, 'hide': 0, 'hideticks': 0,
            'hidelabels': 0, 'hybridlog': 0,
            'inverseticks': 0, 'ticklength': 6,
            'bothticks': 0, 'minorticks': 0, 'oppositeaxis': 0}, ]
        self.fontname = ['sans-serif', ]
#        self.fontname = ['times new roman', ]
        self.fontsize = ['12', ]
        self.fontweight = ['normal', ]
        # The following are window names, set to None before and after the
        # window is active to avoid duplicate windows.
        self.font_window = None
        self.label_font_window = None
        self.read_window = None
        self.data_set_window = None
        self.data_entry_window = None
        self.data_text = None
        self.data_set_transformation_window = None
        self.data_set_delete_window = None
        self.data_set_fitting_window = None
        self.data_set_sort_window = None
        self.box_window = None
        self.line_window = None
        self.ellipse_window = None
        self.vector_window = None
        self.plot_control_window = None
        self.tile_window = None
        self.setplot_window = None
        self.hideplot_window = None
        self.plot_margin = 0.
        self.plot_frame = [0., ]
        self.nxplots = 1
        self.nyplots = 1
        self.current_plot = 1
        self.number_of_plots = self.nxplots * self.nyplots
        self.plot_area_flag = True
        self.save_version = 1.0
        if parent is not None:
            # initialize the window and make the plot area.
            Tk.Frame.__init__(self, parent, args)
            self.root = parent
            self.make_widget()

    def make_widget(self):
        """
        Create the main GUI window.

        This routine makes the main GUI window.  It uses the root variable
        passed to the class when it is initialized (although I do not know
        how this works as the value is not passed to this routine...but it
        does work).

        Parameters
        ----------
            None

        Returns
        -------
            Nothing is returned from this routine.

        """
        menuframe = Tk.Frame(self.root)
        menuframe.pack(side=Tk.TOP, anchor=Tk.W)
        self.make_menus(menuframe)
        controlframe = Tk.Frame(self.root)
        controlframe.pack(side=Tk.LEFT, fill=Tk.Y, expand=1)
        self.make_controls(controlframe)
        sl = self.separator_line(self.root, 5, 750, 5, False, Tk.LEFT)
        self.plotframe = Tk.Frame(self.root)
        self.plotframe.pack(side=Tk.LEFT, fill=Tk.Y, expand=1)
        self.make_plot_area(self.plotframe)

    def make_menus(self, parent):
        """
        Create pull-down menus for plot functionality.

        Given a Tk Frame variable "parent" this routine makes a pull-down
        menu area within this frame.

        Parameters
        ----------
            parent     A Tk.Frame variable, that holds the menus

        Returns
        -------
            No values are returned by the routine.

        """
        menubutton1 = Tk.Menubutton(parent, text="Data")
        menubutton1.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
        menu1 = Tk.Menu(menubutton1)
        menubutton1['menu'] = menu1
        menu1.add_command(label='Read Data', command=self.read_data_set)
        menu1.add_command(label='Create Values by Formula',
                          command=self.create_data_set)
        menu1.add_command(label='Create Values in Widget',
                          command=self.create_data_set_by_editor)
        menu1.add_command(label='Write Data', command=self.write_data_sets)
        label1 = Tk.Label(parent, text="    ")
        label1.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
        menubutton2 = Tk.Menubutton(parent, text="Sets")
        menubutton2.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
        menu2 = Tk.Menu(menubutton2)
        menubutton2['menu'] = menu2
        menu2.add_command(label='Set Properties',
                          command=self.make_data_set_window)
        menu2.add_command(label='Fit Sets',
                          command=self.make_data_set_fitting_window)
        menu2.add_command(label='Set Statistics',
                          command=self.set_statistics)
        menu2.add_command(label='Transform Set',
                          command=self.make_data_set_transformation_window)
        menu2.add_command(label='Edit Set in Widget',
                          command=self.make_data_set_edit_window)
        menu2.add_command(label='Sort Set',
                          command=self.make_data_set_sort_window)
        menu2.add_command(label='Delete Set',
                          command=self.make_data_set_delete_window)
        label2 = Tk.Label(parent, text="    ")
        label2.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
        menubutton3 = Tk.Menubutton(parent, text="Plot")
        menubutton3.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
        menu3 = Tk.Menu(menubutton3)
        menubutton3['menu'] = menu3
        menu3.add_command(label='Plot Parameters',
                          command=self.make_plot_control_window)
        menu3.add_command(label='Clear Plot', command=self.clear_plot)
        menu3.add_command(label='Clear Current Plot',
                          command=self.clear_current_plot)
        menu3.add_command(label='Opposite Y Axis', command=self.set_opposite_y)
        menu3.add_command(label='Opposite X Axis', command=self.set_opposite_x)
        menu3.add_command(label='Tile Plots', command=self.tile_plots)
        menu3.add_command(label='Set Plot', command=self.set_plot_number)
        menu3.add_command(label='Hide/Show Plot', command=self.set_plot_hide)
        menu3.add_command(label='Toggle Equal Aspect',
                          command=self.toggle_equal_aspect)
        menu3.add_separator()
        menu3.add_command(label='Add a Label', command=self.set_label)
        menu3.add_command(label='Edit Labels', command=self.edit_labels)
        menu3.add_command(label='Clear All Labels', command=self.clear_labels)
        menu3.add_command(label='Set Font', command=self.set_font)
        label3 = Tk.Label(parent, text="    ")
        label3.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
        menubutton4 = Tk.Menubutton(parent, text="Plot Items")
        menubutton4.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
        menu4 = Tk.Menu(menubutton4)
        menubutton4['menu'] = menu4
        menu4.add_command(label='Add a line', command=self.add_line)
        menu4.add_command(label='Add an ellipse', command=self.add_ellipse)
        menu4.add_command(label='Add a box', command=self.add_box)
        menu4.add_command(label='Add a vector', command=self.add_vector)
        menu4.add_separator()
        menu4.add_command(label='Edit lines', command=self.edit_lines)
        menu4.add_command(label='Edit ellipses', command=self.edit_ellipses)
        menu4.add_command(label='Edit boxes', command=self.edit_boxes)
        menu4.add_command(label='Edit vectors', command=self.edit_vectors)
        menu4.add_separator()
        menu4.add_command(label='Remove last line', command=self.remove_line)
        menu4.add_command(label='Remove last ellipse',
                          command=self.remove_ellipse)
        menu4.add_command(label='Remove last box', command=self.remove_box)
        menu4.add_command(label='Remove last vector',
                          command=self.remove_vector)
        menu4.add_separator()
        menu4.add_command(label='Remove all lines',
                          command=self.remove_all_lines)
        menu4.add_command(label='Remove all ellipses',
                          command=self.remove_all_ellipses)
        menu4.add_command(label='Remove all boxes',
                          command=self.remove_all_boxes)
        menu4.add_command(label='Remove all vectors',
                          command=self.remove_all_vectors)
        label4 = Tk.Label(parent, text="    ")
        label4.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
        menubutton5 = Tk.Menubutton(parent, text="Save/Restore Plot")
        menubutton5.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
        menu5 = Tk.Menu(menubutton5)
        menubutton5['menu'] = menu5
        menu5.add_command(label='Save configuration', command=self.save_plot)
        menu5.add_command(label='Read configuration', command=self.load_plot)
        menu5.add_command(label='Save as PNG',
                          command=lambda: save_png_figure(self.figure))
        menu5.add_command(label='Save as postscript',
                          command=lambda: save_ps_figure(self.figure))

    def toggle_equal_aspect(self):
        """
        Toggle the equal aspect plot display flag.

        Returns
        -------
        None.

        """
        self.equal_aspect[self.current_plot-1] = not \
            self.equal_aspect[self.current_plot-1]
        self.make_plot()

    def set_opposite_y(self):
        """
        Arrange right side independent y axis on one of the plots.

        Returns
        -------
        None.

        """
        plot_number = self.current_plot
        self.subplot.append(self.figure.add_subplot(
            self.nxplots, self.nyplots, self.current_plot,
            sharex=self.subplot[self.current_plot-1], frameon=False))
#        self.subplot.append(
#            self.subplot[self.current_plot-1].secondary_yaxis(
#            "right"))
        self.share_axis.append(plot_number)
        self.xparameters.append({
            'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
            'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
            'invert': 0, 'hide': 0, 'hideticks': 0,
            'hidelabels': 0, 'hybridlog': 0,
            'inverseticks': 0, 'ticklength': 6,
            'bothticks': 0, 'minorticks': 0, 'oppositeaxis': 0})
        for key in self.xparameters[-1].keys():
            self.xparameters[-1][key] = self.xparameters[
                self.current_plot-1][key]
        self.yparameters.append({
            'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
            'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
            'invert': 0, 'hide': 0, 'hideticks': 0,
            'hidelabels': 0, 'hybridlog': 0,
            'inverseticks': 0, 'ticklength': 6,
            'bothticks': 0, 'minorticks': 0,
            'oppositeaxis': 1-self.yparameters[self.current_plot-1][
                'oppositeaxis']})
#        self.yparameters[-1]['oppositeaxis'] = \
#            1-self.yparameters[self.current_plot-1]['oppositeaxis']
        self.fontname.append(deepcopy(self.fontname[self.current_plot-1]))
        self.fontsize.append(deepcopy(self.fontsize[self.current_plot-1]))
        self.fontweight.append(deepcopy(self.fontweight[self.current_plot-1]))
        self.legend_variable.append(None)
        self.legend_frame.append(None)
        self.legend_options.append(None)
        self.legend_position.append(None)
        self.legend_user_position.append(None)
        self.plot_frame.append(deepcopy(self.plot_frame[self.current_plot-1]))
        self.plot_range.append(deepcopy(self.plot_range[self.current_plot-1]))
        self.bounding_box.append(deepcopy(self.bounding_box[
            self.current_plot-1]))
        self.title.append('')
        self.current_plot = len(self.subplot)
        self.number_of_plots = self.number_of_plots+1
        self.make_plot()

    def set_opposite_x(self):
        """
        Arrange top independent x axis on one of the plots.

        Returns
        -------
        None.

        """
        plot_number = self.current_plot
        self.subplot.append(self.figure.add_subplot(
            self.nxplots, self.nyplots, self.current_plot,
            sharey=self.subplot[self.current_plot-1], frameon=False))
#        self.subplot.append(
#            self.subplot[self.current_plot-1].secondary_xaxis(
#            "top"))
        self.share_axis.append(-1*plot_number)
        self.xparameters.append({
            'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
            'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
            'invert': 0, 'hide': 0, 'hideticks': 0,
            'hidelabels': 0, 'hybridlog': 0,
            'inverseticks': 0, 'ticklength': 6,
            'bothticks': 0, 'minorticks': 0,
            'oppositeaxis': 1-self.xparameters[self.current_plot-1][
                'oppositeaxis']})
        self.yparameters.append({
            'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
            'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
            'invert': 0, 'hide': 0, 'hideticks': 0,
            'hidelabels': 0, 'hybridlog': 0,
            'inverseticks': 0, 'ticklength': 6,
            'bothticks': 0, 'minorticks': 0,
            'oppositeaxis': 0})
        for key in self.yparameters[-1].keys():
            self.yparameters[-1][key] = self.yparameters[
                self.current_plot-1][key]
        self.xparameters[-1]['oppositeaxis'] = \
            1-self.xparameters[self.current_plot-1]['oppositeaxis']
        self.fontname.append(deepcopy(self.fontname[self.current_plot-1]))
        self.fontsize.append(deepcopy(self.fontsize[self.current_plot-1]))
        self.fontweight.append(deepcopy(self.fontweight[self.current_plot-1]))
        self.legend_variable.append(None)
        self.legend_frame.append(None)
        self.legend_options.append(None)
        self.legend_position.append(None)
        self.legend_user_position.append(None)
        self.plot_frame.append(deepcopy(self.plot_frame[self.current_plot-1]))
        self.plot_range.append(deepcopy(self.plot_range[self.current_plot-1]))
        self.bounding_box.append(deepcopy(self.bounding_box[
            self.current_plot-1]))
        self.title.append('')
        self.current_plot = len(self.subplot)
        self.number_of_plots = self.number_of_plots+1
        self.make_plot()

    def clear_current_plot(self):
        """
        Clear the current plot if the user OKs this action.

        This routine clears the sets, parameters and plot items for the
        currently selected plot.  All values return to their initial values.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if (self.nxplots == 1) and (self.nyplots == 1):
            self.clear_plot()
        else:
            response = tkinter.messagebox.askyesno(
                "Verify",
                "Do you want to abandon the current plot?")
            if not response:
                return
            if self.nsets > 0:
                for loop in range(self.max_sets):
                    if self.set_properties[loop]['plot'] == self.current_plot:
                        self.set_properties[loop] = {
                            'symbol': None, 'symbolsize': 4.0,
                            'linestyle': 'None', 'linewidth': 1.0,
                            'colour': 'black', 'label': '', 'xmin': 0.0,
                            'xmax': 1.0, 'ymin': 0.0, 'ymax': 1.0,
                            'display': True, 'errors': False, 'legend': True,
                            'plot': 1}
                        self.xdata[loop] = None
                        self.ydata[loop] = None
                        self.original_range[loop] = True
            self.plot_range[self.current_plot-1] = [0., 1., 0., 1.]
            self.original_range[self.current_plot-1] = True
            self.title[self.current_plot-1] = ' '
            self.xparameters[self.current_plot-1] = {
                'label': ' ', 'minimum': 0.0, 'maximum': 1.0, 'major': 0.1,
                'minor': 0.1, 'logarithmic': 0, 'invert': 0, 'hide': 0,
                'hideticks': 0, 'hidelabels': 0, 'hybridlog': 0,
                'inverseticks': 0, 'ticklength': 6, 'bothticks': 0,
                'minorticks': 0, 'oppositeaxis': 0}
            self.yparameters[self.current_plot-1] = {
                'label': ' ', 'minimum': 0.0, 'maximum': 1.0, 'major': 0.1,
                'minor': 0.1, 'logarithmic': 0, 'invert': 0, 'hide': 0,
                'hideticks': 0, 'hidelabels': 0, 'hybridlog': 0,
                'inverseticks': 0, 'ticklength': 6, 'bothticks': 0,
                'minorticks': 0, 'oppositeaxis': 0}
            try:
                self.legend_variable[self.current_plot-1].set(0)
            except:
                self.legend_variable[self.current_plot-1] = None
        self.legend_handles[self.current_plot-1] = None
        self.legend_labels[self.current_plot-1] = None
        self.legend_frame[self.current_plot-1] = None
        self.legend_position[self.current_plot-1] = None
        self.legend_user_position[self.current_plot-1] = None
        self.plot_frame[self.current_plot-1] = 0.0
        for loop in range(self.number_of_lines):
            if self.plot_lines[loop]['plot'] == self.current_plot:
                self.plot_lines[loop] = {
                    'xstart': None, 'ystart': None, 'xend': None,
                    'yend': None, 'plot': 1, 'line_type': 'solid',
                    'line_colour': 'black', 'line_thickness': 1.0}
        for loop in range(self.number_of_labels):
            if self.plot_labels[loop]['plot'] == self.current_plot:
                self.plot_labels = {
                    'xposition': None, 'yposition': None, 'labelstring': '',
                    'plot': 1, 'colour': 'black', 'size': 12,
                    # 'font': 'times new roman', 'fontweight': 'normal'}
                    'font': 'sans-serif', 'fontweight': 'normal'}
        for loop in range(self.number_of_vectors):
            if self.plot_vectors[loop]['plot'] == self.current_plot:
                self.plot_vectors[loop] = {
                    'xstart': None, 'ystart': None, 'xend': None,
                    'yend': None, 'delx': None, 'dely': None, 'plot': 1,
                    'line_type': 'solid', 'line_colour': 'black',
                    'line_thickness': 1.0, 'fill': True,
                    'fill_colour': 'black'}
        for loop in range(self.number_of_boxes):
            if self.plot_boxes[loop]['plot'] == self.current_plot:
                self.plot_boxes[loop] = {
                    'xstart': None, 'ystart': None, 'xend': None,
                    'yend': None, 'rotation': 0.0, 'plot': 1,
                    'line_type': 'solid', 'line_colour': 'black',
                    'line_thickness': 1.0, 'fill_colour': 'none'}
        for loop in range(self.number_of_ellipses):
            if self.plot_ellipses[loop]['plot'] == self.current_plot:
                self.plot_ellipses[loop] = {
                    'xposition': None, 'yposition': None, 'major': None,
                    'minor': None, 'rotation': 0.0, 'plot': 1,
                    'line_type': 'solid', 'line_colour': 'black',
                    'line_thickness': 1.0, 'fill_colour': 'none'}
        self.clear_sets()
        self.current_plot = 1
        self.make_plot()

    def clear_sets(self):
        """
        Clear inactive sets from the set variable.

        Reset the set values to remove ones that have been removed with a plot.
        These are marked by self.xdata[self.current_set-1]['values'] set to
        None.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing
        """
        nsets = 0
        properties = []
        xdata = []
        ydata = []
        original_range = []
        for loop in range(self.nsets):
            if self.xdata[loop]['values'] is not None:
                xdata.append(self.xdata[loop])
                ydata.append(self.ydata[loop])
                properties.append(self.set_properties[loop])
                original_range.append(self.original_range[loop])
                nsets = nsets + 1
        for loop in range(nsets, self.max_sets):
            xdata.append(None)
            ydata.append(None)
            original_range.append(True)
            properties.append({
                'symbol': None, 'symbolsize': 4.0,
                'linestyle': 'None', 'linewidth': 1.0, 'colour': 'black',
                'label': '', 'xmin': 0.0, 'xmax': 1.0, 'ymin': 0.0,
                'ymax': 1.0, 'display': True, 'errors': False,
                'legend': True, 'plot': 1})
        self.nsets = nsets
        self.xdata = xdata
        self.ydata = ydata
        self.set_properties = properties
        self.original_range = original_range

    def save_plot(self):
        """
        Save the plot state to an ascii file that can be loaded later.

        This routine writes the current sets and parameters to an ascii output
        file in a set format that can be read back in later.

        No values are passed to this routine or returned from this routine.
        """
        outfilename = tkinter.filedialog.asksaveasfilename()
        if isinstance(outfilename, type('string')):
            outfile = open(outfilename, 'w')
            print('# matplotlib_gui_code.py save file version 1.0', file=outfile)
            print('# number of sets: %d maximum: %d' %
                  (self.nsets, self.max_sets), file=outfile)
            for loop in range(self.nsets):
                print('# set number %d: %d points' %
                      (loop+1, len(self.xdata[loop]['values'])), file=outfile)
                for n1 in range(len(self.xdata[loop]['values'])):
                    print(self.xdata[loop]['values'][n1],
                          self.ydata[loop]['values'][n1],
                          self.xdata[loop]['lowerror'][n1],
                          self.ydata[loop]['lowerror'][n1],
                          self.xdata[loop]['higherror'][n1],
                          self.ydata[loop]['higherror'][n1], file=outfile)
                print('# set parameters', file=outfile)
                print(self.xdata[loop]['minimum'],
                      self.ydata[loop]['minimum'],
                      self.xdata[loop]['maximum'],
                      self.ydata[loop]['maximum'], file=outfile)
                print(self.xdata[loop]['errors'], self.ydata[loop]['errors'],
                      self.xdata[loop]['legend'], self.ydata[loop]['legend'],
                      file=outfile)
            print('# set properties', file=outfile)
            for loop in range(self.nsets):
                str1 = '%s\t%f\t%s\t%f\t%s\t%s\t%f\t%f\t%f\t%f' % (
                    self.set_properties[loop]['symbol'],
                    self.set_properties[loop]['symbolsize'],
                    self.set_properties[loop]['linestyle'],
                    self.set_properties[loop]['linewidth'],
                    self.set_properties[loop]['colour'],
                    self.set_properties[loop]['label'],
                    self.set_properties[loop]['xmin'],
                    self.set_properties[loop]['xmax'],
                    self.set_properties[loop]['ymin'],
                    self.set_properties[loop]['ymax'])
                str1 = str1 + '\t' \
                    + str(self.set_properties[loop]['display']) \
                    + '\t' + str(self.set_properties[loop]['errors']) \
                    + '\t' + str(self.set_properties[loop]['legend']) \
                    + '\t%d' % (self.set_properties[loop]['plot'])
                print(str1, file=outfile)
            print('# plot properties', file=outfile)
            print('# nxplots, nyplots, total: %d %d %d ' %
                  (self.nxplots, self.nyplots, self.number_of_plots),
                  file=outfile)
            str1 = '# hide plots: '
            for loop in range(self.number_of_plots):
                str1 = str1 + str(self.hide_subplot[loop]) + '\t'
            str1 = str1.rstrip('\t')
            print(str1, file=outfile)
            print('# margin: %f ' % (self.plot_margin), file=outfile)
            nplot = 0
            for n1 in range(self.nxplots):
                for n2 in range(self.nyplots):
                    print('# plot %d index values: %d %d' %
                          (nplot, n1, n2), file=outfile)
                    print('# title: %s ' % (self.title[nplot]), file=outfile)
                    print('# frame: %f ' %
                          (self.plot_frame[nplot]), file=outfile)
                    print('# font name: %s ' %
                          (self.fontname[nplot]), file=outfile)
                    print('# font size: %s ' %
                          (self.fontsize[nplot]), file=outfile)
                    print('# font weight: %s ' %
                          (self.fontweight[nplot]), file=outfile)
                    strformat = '# x parameters: \t%s\t%f\t%f\t%f\t%f\t%d' \
                                + '\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d'
                    print(strformat % (
                        self.xparameters[nplot]['label'],
                        self.xparameters[nplot]['minimum'],
                        self.xparameters[nplot]['maximum'],
                        self.xparameters[nplot]['major'],
                        self.xparameters[nplot]['minor'],
                        self.xparameters[nplot]['logarithmic'],
                        self.xparameters[nplot]['invert'],
                        self.xparameters[nplot]['hide'],
                        self.xparameters[nplot]['hideticks'],
                        self.xparameters[nplot]['hidelabels'],
                        self.xparameters[nplot]['hybridlog'],
                        self.xparameters[nplot]['inverseticks'],
                        self.xparameters[nplot]['ticklength'],
                        self.xparameters[nplot]['bothticks'],
                        self.xparameters[nplot]['minorticks'],
                        self.xparameters[nplot]['oppositeaxis']),
                          file=outfile)
                    strformat = '# y parameters: \t%s\t%f\t%f\t%f\t%f\t%d' \
                        + '\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d'
                    print(strformat % (
                        self.yparameters[nplot]['label'],
                        self.yparameters[nplot]['minimum'],
                        self.yparameters[nplot]['maximum'],
                        self.yparameters[nplot]['major'],
                        self.yparameters[nplot]['minor'],
                        self.yparameters[nplot]['logarithmic'],
                        self.yparameters[nplot]['invert'],
                        self.yparameters[nplot]['hide'],
                        self.yparameters[nplot]['hideticks'],
                        self.yparameters[nplot]['hidelabels'],
                        self.yparameters[nplot]['hybridlog'],
                        self.yparameters[nplot]['inverseticks'],
                        self.yparameters[nplot]['ticklength'],
                        self.yparameters[nplot]['bothticks'],
                        self.yparameters[nplot]['minorticks'],
                        self.yparameters[nplot]['oppositeaxis']),
                          file=outfile)
                    print('# plot range: %f %f %f %f %s' % (
                        self.plot_range[nplot][0],
                        self.plot_range[nplot][1],
                        self.plot_range[nplot][2],
                        self.plot_range[nplot][3],
                        self.original_range[nplot]), file=outfile)
                    nplot = nplot + 1
            if len(self.xparameters) > self.nxplots*self.nyplots:
                for n1 in range(nplot, len(self.xparameters)):
                    print('# plot %d index values: %d %d' %
                          (nplot, -1, -1), file=outfile)
                    print('# title: %s ' % (self.title[nplot]), file=outfile)
                    print('# frame: %f ' %
                          (self.plot_frame[nplot]), file=outfile)
                    print('# font name: %s ' %
                          (self.fontname[nplot]), file=outfile)
                    print('# font size: %s ' %
                          (self.fontsize[nplot]), file=outfile)
                    print('# font weight: %s ' %
                          (self.fontweight[nplot]), file=outfile)
                    strformat = '# x parameters: \t%s\t%f\t%f\t%f\t%f\t%d' \
                                + '\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d'
                    print(strformat % (
                        self.xparameters[nplot]['label'],
                        self.xparameters[nplot]['minimum'],
                        self.xparameters[nplot]['maximum'],
                        self.xparameters[nplot]['major'],
                        self.xparameters[nplot]['minor'],
                        self.xparameters[nplot]['logarithmic'],
                        self.xparameters[nplot]['invert'],
                        self.xparameters[nplot]['hide'],
                        self.xparameters[nplot]['hideticks'],
                        self.xparameters[nplot]['hidelabels'],
                        self.xparameters[nplot]['hybridlog'],
                        self.xparameters[nplot]['inverseticks'],
                        self.xparameters[nplot]['ticklength'],
                        self.xparameters[nplot]['bothticks'],
                        self.xparameters[nplot]['minorticks'],
                        self.xparameters[nplot]['oppositeaxis']),
                          file=outfile)
                    strformat = '# y parameters: \t%s\t%f\t%f\t%f\t%f\t%d' \
                        + '\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d'
                    print(strformat % (
                        self.yparameters[nplot]['label'],
                        self.yparameters[nplot]['minimum'],
                        self.yparameters[nplot]['maximum'],
                        self.yparameters[nplot]['major'],
                        self.yparameters[nplot]['minor'],
                        self.yparameters[nplot]['logarithmic'],
                        self.yparameters[nplot]['invert'],
                        self.yparameters[nplot]['hide'],
                        self.yparameters[nplot]['hideticks'],
                        self.yparameters[nplot]['hidelabels'],
                        self.yparameters[nplot]['hybridlog'],
                        self.yparameters[nplot]['inverseticks'],
                        self.yparameters[nplot]['ticklength'],
                        self.yparameters[nplot]['bothticks'],
                        self.yparameters[nplot]['minorticks'],
                        self.yparameters[nplot]['oppositeaxis']),
                          file=outfile)
                    print('# plot range: %f %f %f %f %s' % (
                        self.plot_range[nplot][0],
                        self.plot_range[nplot][1],
                        self.plot_range[nplot][2],
                        self.plot_range[nplot][3],
                        self.original_range[nplot]), file=outfile)
                    nplot = nplot + 1
            print('# number of labels: %d maximum: %d' %
                  (self.number_of_labels, self.max_labels), file=outfile)
            for loop in range(self.number_of_labels):
                print('%f\t%f\t%d\t%s\t%s\t%d\t%s\t%s' % (
                    self.plot_labels[loop]['xposition'],
                    self.plot_labels[loop]['yposition'],
                    self.plot_labels[loop]['plot'],
                    self.plot_labels[loop]['labelstring'],
                    self.plot_labels[loop]['colour'],
                    self.plot_labels[loop]['size'],
                    self.plot_labels[loop]['font'],
                    self.plot_labels[loop]['fontweight']), file=outfile)
            print('# number of lines: %d maximum: %d ' %
                  (self.number_of_lines, self.max_lines), file=outfile)
            for loop in range(self.number_of_lines):
                print('%f\t%f\t%f\t%f\t%d\t%s\t%s\t%f' % (
                    self.plot_lines[loop]['xstart'],
                    self.plot_lines[loop]['ystart'],
                    self.plot_lines[loop]['xend'],
                    self.plot_lines[loop]['yend'],
                    self.plot_lines[loop]['plot'],
                    self.plot_lines[loop]['line_type'],
                    self.plot_lines[loop]['line_colour'],
                    self.plot_lines[loop]['line_thickness']), file=outfile)
            print('# number of vectors: %d maximum: %d ' %
                  (self.number_of_vectors, self.max_vectors), file=outfile)
            for loop in range(self.number_of_vectors):
                print('%f\t%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%f\t%5.5s\t%s' % (
                    self.plot_vectors[loop]['xstart'],
                    self.plot_vectors[loop]['ystart'],
                    self.plot_vectors[loop]['xend'],
                    self.plot_vectors[loop]['yend'],
                    self.plot_vectors[loop]['delx'],
                    self.plot_vectors[loop]['dely'],
                    self.plot_vectors[loop]['plot'],
                    self.plot_vectors[loop]['line_type'],
                    self.plot_vectors[loop]['line_colour'],
                    self.plot_vectors[loop]['line_thickness'],
                    str(self.plot_vectors[loop]['fill']),
                    self.plot_vectors[loop]['fill_colour']), file=outfile)
            print('# number of ellipses: %d maximum: %d ' %
                  (self.number_of_ellipses, self.max_ellipses), file=outfile)
            for loop in range(self.number_of_ellipses):
                print('%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%f\t%s' % (
                    self.plot_ellipses[loop]['xposition'],
                    self.plot_ellipses[loop]['yposition'],
                    self.plot_ellipses[loop]['major'],
                    self.plot_ellipses[loop]['minor'],
                    self.plot_ellipses[loop]['rotation'],
                    self.plot_ellipses[loop]['plot'],
                    self.plot_ellipses[loop]['line_colour'],
                    self.plot_ellipses[loop]['line_type'],
                    self.plot_ellipses[loop]['line_thickness'],
                    self.plot_ellipses[loop]['fill_colour']), file=outfile)
            print('# number of boxes: %d maximum: %d ' %
                  (self.number_of_boxes, self.max_boxes), file=outfile)
            for loop in range(self.number_of_boxes):
                print('%f\t%f\t%f\t%f\t%f\t%d\t%s\t%s\t%f\t%s' % (
                    self.plot_boxes[loop]['xstart'],
                    self.plot_boxes[loop]['ystart'],
                    self.plot_boxes[loop]['xend'],
                    self.plot_boxes[loop]['yend'],
                    self.plot_boxes[loop]['rotation'],
                    self.plot_boxes[loop]['plot'],
                    self.plot_boxes[loop]['line_colour'],
                    self.plot_boxes[loop]['line_type'],
                    self.plot_boxes[loop]['line_thickness'],
                    self.plot_boxes[loop]['fill_colour']), file=outfile)
            print('# legend values:', file=outfile)
            nplot = 0
            for n1 in range(self.nxplots):
                for n2 in range(self.nyplots):
                    str1 = '# plot %d index values %d %d: \t' % (nplot, n1, n2)
                    if self.legend_variable[nplot] is None:
                        str1 = str1 + 'None\tNone\tNone\tNone\tNone\t'
                    else:
                        str1 = str1 + '%d\t%d\t%s\t%s\t' % (
                            self.legend_variable[nplot].get(),
                            self.legend_frame[nplot].get(),
                            self.legend_position[nplot],
                            self.legend_options[nplot].get())
                        try:
                            str1 = str1 + '%f\t%f' % (
                                self.legend_user_position[nplot][0],
                                self.legend_user_position[nplot][1])
                        except:
                            str1 = str1 + 'None\tNone'
                    print(str1, file=outfile)
                    nplot = nplot + 1
            if len(self.legend_variable) > self.nxplots*self.nyplots:
                for n1 in range(nplot, len(self.xparameters)):
                    str1 = '# plot %d index values %d %d: \t' % (nplot, -1, -1)
                    if self.legend_variable[nplot] is None:
                        str1 = str1 + 'None\tNone\tNone\tNone\tNone\t'
                    else:
                        str1 = str1 + '%d\t%d\t%s\t%s\t' % (
                            self.legend_variable[nplot].get(),
                            self.legend_frame[nplot].get(),
                            self.legend_position[nplot],
                            self.legend_options[nplot].get())
                        try:
                            str1 = str1 + '%f\t%f' % (
                                self.legend_user_position[nplot][0],
                                self.legend_user_position[nplot][1])
                        except:
                            str1 = str1 + 'None\tNone'
                    print(str1, file=outfile)
                    nplot = nplot + 1
            for n1 in range(len(self.share_axis)):
                if n1 == 0:
                    str1 = ' %d ' % (self.share_axis[n1])
                else:
                    str1 = str1 + '\t %d ' % (self.share_axis[n1])
            print('# share_axis', file=outfile)
            print(str1, file=outfile)
            str1 = ''
            for n1 in range(len(self.equal_aspect)):
                if self.equal_aspect[n1]:
                    str1 = str1 + ' True' + '\t'
                else:
                    str1 = str1 + 'False' + '\t'
            str1 = str1.rstrip('\t')
            print('# equal_aspect', file=outfile)
            print(str1, file=outfile)
            print('# end', file=outfile)
            outfile.close()
        else:
            tkinter.messagebox.showinfo(
                "Error",
                "File "+outfilename+" was not written properly.")

    def load_plot(self):
        """
        Read an ascii file of plot parameters to make a plot.

        This routine asks for the ascii file to load and calls the routine
        to parse the file and load in the data and parameters.

        No values are passed to this routine or returned from this routine.
        """
        savefilename = tkinter.filedialog.askopenfilename()
        if savefilename is None:
            return
        if not os.path.isfile(savefilename):
            tkinter.messagebox.showinfo(
                "Error",
                "The file %s was not found." % (savefilename))
            return
        savefile = open(savefilename, 'r')
        lines = savefile.readlines()
        savefile.close()
        flag = self.parse_save_file(lines, False)
        if not flag:
            tkinter.messagebox.showinfo(
                "Error",
                "The file %s is not formatted properly for a save file." %
                (savefilename))
            return
        response = tkinter.messagebox.askyesno(
            "Verify",
            "Do you want to abandon the current plot for the saved one?")
        if response:
            self.clear_plot(False)
            flag = self.parse_save_file(lines, True)
            save_plot_range = deepcopy(self.plot_range)
            save_xpars = deepcopy(self.xparameters)
            share_axis = deepcopy(self.share_axis)
            save_ypars = deepcopy(self.yparameters)
            nx1 = 1*self.nxplots
            ny1 = 1*self.nyplots
            nplots = 1*self.number_of_plots
            self.nxplots = 0
            self.nyplots = 0
            n1 = self.current_plot
            self.subplot = []
            self.hide_subplot = [False, ]
            self.make_plot_layout(nx1, ny1, 1)
            self.share_axis = deepcopy(share_axis)
            if self.number_of_plots < nplots:
                for loop in range(len(self.share_axis)):
                    if self.share_axis[loop] != 0:
                        n2 = abs(self.share_axis[loop])
                        self.current_plot = n2
                        if self.share_axis[loop] < 0:
                            self.subplot.append(self.figure.add_subplot(
                                self.nxplots, self.nyplots, n2,
                                sharey=self.subplot[n2-1],
                                frameon=False))
                            self.bounding_box.append(self.bounding_box[n2-1])
                            self.current_axis = len(self.subplot)
                        else:
                            self.subplot.append(self.figure.add_subplot(
                                self.nxplots, self.nyplots, n2,
                                sharex=self.subplot[n2-1],
                                frameon=False))
                            self.bounding_box.append(self.bounding_box[n2-1])
                            self.current_axis = len(self.subplot)
                self.number_of_plots = nplots
            self.xparameters = deepcopy(save_xpars)
            self.yparameters = deepcopy(save_ypars)
            for loop in range(len(self.subplot)):
                self.current_plot = loop+1
                self.make_plot()
            self.current_plot = n1
            self.plot_range = deepcopy(save_plot_range)
            self.make_plot()

    def parse_save_file(self, lines, flag):
        """
        Parse an ascii save file and optionally load the values.

        This program reads the lines from a matplotlib_gui_code.py save file.  It
        determines whether the file is structured properly.  If the
        flag value is True it also tries to set the parameters for the plot.

        Parameters
        ----------
            lines : A set of limes (assumed to be from the .readlines()
                    function) from a matplotlib_gui_code.py save file

            flag :  A boolean value, if True the code tries to assign the
                    values for the plot; if False it only tries to parse
                    the file

        Returns
        -------
            goodfile :  A boolean value for whether the lines are of the
                       expected structure for a matplotlib_gui_code.py save file

        """
        goodfile = False
        equal_aspect = None
        ndatasets = 0
        xdata = []
        ydata = []
        set_properties = []
        set_parameters = []
        plot_parameters = []
        for ind1 in range(len(lines)):
            line = lines[ind1]
            line = line.strip('\n')
            if '# matplotlib_gui_code.py save file' in line:
                values = line.split()
                try:
                    version = float(values[-1])
                except:
                    return goodfile
                if version != self.save_version:
                    return goodfile
            if '# number of sets:' in line:
                values = line.split()
                if len(values) != 7:
                    return goodfile
                try:
                    nsets = int(values[4])
                    maxsets = int(values[6])
                except:
                    return goodfile
                if (nsets < 1) or (nsets > maxsets):
                    return goodfile
            if '# set number' in line:
                ind2 = self.line_range(lines, ind1)
                npoints = ind2 - ind1 - 1
                xset = numpy.zeros((npoints), dtype=numpy.float32)
                yset = numpy.zeros((npoints), dtype=numpy.float32)
                xseterr1 = numpy.zeros((npoints), dtype=numpy.float32)
                yseterr1 = numpy.zeros((npoints), dtype=numpy.float32)
                xseterr2 = numpy.zeros((npoints), dtype=numpy.float32)
                yseterr2 = numpy.zeros((npoints), dtype=numpy.float32)
                for loop in range(ind1+1, ind2):
                    n1 = loop - ind1 - 1
                    line1 = lines[loop].strip('\n')
                    values = line1.split()
                    if len(values) == 6:
                        xset[n1] = values[0]
                        yset[n1] = values[1]
                        xseterr1[n1] = values[2]
                        yseterr1[n1] = values[3]
                        xseterr2[n1] = values[4]
                        yseterr2[n1] = values[5]
                    else:
                        return goodfile
                xdata.append(xset)
                xdata.append(xseterr1)
                xdata.append(xseterr2)
                ydata.append(yset)
                ydata.append(yseterr1)
                ydata.append(yseterr2)
                ndatasets = ndatasets + 1
            if '# set parameters' in line:
                p1 = {}
                line1 = lines[ind1+1].strip('\n')
                values = line1.split()
                if len(values) == 4:
                    p1['xminimum'] = float(values[0])
                    p1['yminimum'] = float(values[1])
                    p1['xmaximum'] = float(values[2])
                    p1['ymaximum'] = float(values[3])
                line2 = lines[ind1+2].strip('\n')
                values = line2.split()
                if len(values) == 4:
                    if values[0] == 'True':
                        p1['xerrors'] = True
                    elif values[0] == 'False':
                        p1['xerrors'] = False
                    else:
                        return goodfile
                    if values[1] == 'True':
                        p1['yerrors'] = True
                    elif values[1] == 'False':
                        p1['yerrors'] = False
                    else:
                        return goodfile
                    if values[2] == 'True':
                        p1['xlegend'] = True
                    elif values[2] == 'False':
                        p1['xlegend'] = False
                    else:
                        return goodfile
                    if values[3] == 'True':
                        p1['ylegend'] = True
                    elif values[3] == 'False':
                        p1['ylegend'] = False
                    else:
                        return goodfile
                set_parameters.append(p1)
            if '# set properties' in line:
                if ndatasets != nsets:
                    return goodfile
                ind2 = self.line_range(lines, ind1)
                npoints = ind2 - ind1 - 1
                if npoints == nsets:
                    for loop in range(ind1+1, ind2):
                        p2 = {}
                        line = lines[loop].strip('\n')
                        values = line.split('\t')
                        if len(values) != 14:
                            return goodfile
                        p2['symbol'] = values[0]
                        if p2['symbol'] == 'None':
                            p2['symbol'] = None
                        p2['symbolsize'] = float(values[1])
                        p2['linestyle'] = values[2]
                        if p2['linestyle'] == 'None':
                            p2['linestyle'] = None
                        p2['linewidth'] = float(values[3])
                        p2['colour'] = values[4]
                        p2['label'] = values[5]
                        p2['xmin'] = float(values[6])
                        p2['xmax'] = float(values[7])
                        p2['ymin'] = float(values[8])
                        p2['ymax'] = float(values[9])
                        if values[10] == 'True':
                            p2['display'] = True
                        elif values[10] == 'False':
                            p2['display'] = False
                        else:
                            return goodfile
                        if values[11] == 'True':
                            p2['errors'] = True
                        elif values[11] == 'False':
                            p2['errors'] = False
                        else:
                            return goodfile
                        if values[12] == 'True':
                            p2['legend'] = True
                        elif values[12] == 'False':
                            p2['legend'] = False
                        else:
                            return goodfile
                        p2['plot'] = int(values[13])
                        set_properties.append(p2)
                else:
                    return goodfile
            if ('# plot ' in line) and ('index values:' in line):
                pp = {}
            if '# nxplots, nyplots, total:' in line:
                values = line.split()
                nxplots = int(values[-3])
                nyplots = int(values[-2])
                nplots = int(values[-1])
                if (nxplots < 1) or (nyplots < 1):
                    return goodfile
                hide_subplot = []
                for loop in range(nplots):
                    hide_subplot.append(False)
            if '# hide plots' in line:
                line1 = line.replace('# hide plots: ', '')
                values = line.split('\t')
                if len(values) == nplots:
                    for loop in range(nplots):
                        if 'True' in values[loop]:
                            hide_subplot[loop] = True
                        else:
                            hide_subplot[loop] = False
            if '# margin' in line:
                values = line.split()
                if len(values) != 3:
                    return goodfile
                plot_margin = float(values[2])
                if plot_margin < 0.:
                    return goodfile
            if '# title:' in line:
                frag = line.replace('# title:', '')
                frag = frag.lstrip()
                frag = frag.rstrip()
                pp['title'] = frag
            if '# frame' in line:
                values = line.split()
                if len(values) != 3:
                    return goodfile
                pp['frame'] = float(values[2])
            if '# font name:' in line:
                frag = line.replace('# font name:', '')
                frag = frag.lstrip()
                frag = frag.rstrip()
                pp['fontname'] = frag
            if '# font size' in line:
                values = line.split()
                if len(values) != 4:
                    return goodfile
                pp['fontsize'] = int(values[3])
            if '# font weight:' in line:
                frag = line.replace('# font weight:', '')
                frag = frag.lstrip()
                frag = frag.rstrip()
                pp['fontweight'] = frag
            if '# x parameters:' in line:
                line = line.strip('\n')
                values = line.split('\t')
                if len(values) != 17:
                    return goodfile
                xp = {}
                xp['label'] = values[1]
                xp['minimum'] = values[2]
                xp['maximum'] = values[3]
                xp['major'] = float(values[4])
                xp['minor'] = float(values[5])
                xp['logarithmic'] = int(values[6])
                xp['invert'] = int(values[7])
                xp['hide'] = int(values[8])
                xp['hideticks'] = int(values[9])
                xp['hidelabels'] = int(values[10])
                xp['hybridlog'] = int(values[11])
                xp['inverseticks'] = int(values[12])
                xp['ticklength'] = int(values[13])
                xp['bothticks'] = int(values[14])
                xp['minorticks'] = float(values[15])
                xp['oppositeaxis'] = int(values[16])
                pp['xparameters'] = xp
            if '# y parameters:' in line:
                values = line.split('\t')
                if len(values) != 17:
                    return goodfile
                yp = {}
                yp['label'] = values[1]
                yp['minimum'] = values[2]
                yp['maximum'] = values[3]
                yp['major'] = float(values[4])
                yp['minor'] = float(values[5])
                yp['logarithmic'] = int(values[6])
                yp['invert'] = int(values[7])
                yp['hide'] = int(values[8])
                yp['hideticks'] = int(values[9])
                yp['hidelabels'] = int(values[10])
                yp['hybridlog'] = int(values[11])
                yp['inverseticks'] = int(values[12])
                yp['ticklength'] = int(values[13])
                yp['bothticks'] = int(values[14])
                yp['minorticks'] = float(values[15])
                yp['oppositeaxis'] = int(values[16])
                pp['yparameters'] = yp
            if '# plot range:' in line:
                values = line.split()
                if len(values) != 8:
                    return goodfile
                try:
                    xmin = float(values[3])
                    xmax = float(values[4])
                    ymin = float(values[5])
                    ymax = float(values[6])
                except:
                    return goodfile
                pp['plot_range'] = [0., 1., 0., 1.]
                pp['plot_range'][0] = 1.*xmin
                pp['plot_range'][1] = 1.*xmax
                pp['plot_range'][2] = 1.*ymin
                pp['plot_range'][3] = 1.*ymax
                if 'True' in values[7]:
                    pp['original_range'] = True
                elif 'False' in values[7]:
                    pp['original_range'] = False
                else:
                    return goodfile
                plot_parameters.append(pp)
            if '# number of labels:' in line:
                values = line.split()
                if len(values) != 7:
                    return goodfile
                labels = []
                if int(values[4]) > 0:
                    ind2 = self.line_range(lines, ind1)
                    npoints = ind2 - ind1 - 1
                    for loop in range(ind1+1, ind2):
                        lv = {}
                        line = lines[loop].strip('\n')
                        values = line.split('\t')
                        if len(values) != 8:
                            return goodfile
                        lv['xposition'] = float(values[0])
                        lv['yposition'] = float(values[1])
                        lv['plot'] = int(values[2])
                        lv['labelstring'] = values[3]
                        lv['colour'] = values[4]
                        lv['size'] = int(values[5])
                        lv['font'] = values[6]
                        lv['fontweight'] = values[7]
                        labels.append(lv)
            if '# number of lines:' in line:
                values = line.split()
                if len(values) != 7:
                    return goodfile
                plines = []
                if int(values[4]) > 0:
                    ind2 = self.line_range(lines, ind1)
                    npoints = ind2 - ind1 - 1
                    for loop in range(ind1+1, ind2):
                        lv = {}
                        line = lines[loop].strip('\n')
                        values = line.split('\t')
                        if len(values) != 8:
                            return goodfile
                        lv['xstart'] = float(values[0])
                        lv['ystart'] = float(values[1])
                        lv['xend'] = float(values[2])
                        lv['yend'] = float(values[3])
                        lv['plot'] = int(values[4])
                        lv['line_type'] = values[5]
                        lv['line_colour'] = values[6]
                        lv['line_width'] = float(values[7])
                        plines.append(lv)
            if '# number of vectors:' in line:
                values = line.split()
                if len(values) != 7:
                    return goodfile
                vectors = []
                if int(values[4]) > 0:
                    ind2 = self.line_range(lines, ind1)
                    npoints = ind2 - ind1 - 1
                    for loop in range(ind1+1, ind2):
                        v1 = {}
                        line = lines[loop].strip('\n')
                        values = line.split('\t')
                        if len(values) != 12:
                            return goodfile
                        v1['xstart'] = float(values[0])
                        v1['ystart'] = float(values[1])
                        v1['xend'] = float(values[2])
                        v1['yend'] = float(values[3])
                        v1['delx'] = float(values[4])
                        v1['dely'] = float(values[5])
                        v1['plot'] = int(values[6])
                        v1['line_type'] = values[7]
                        v1['line_colour'] = values[8]
                        v1['line_thickness'] = float(values[9])
                        if values[10] == ' True':
                            v1['fill'] = True
                        elif values[10] == 'False':
                            v1['fill'] = False
                        else:
                            return goodfile
                        v1['fill_colour'] = values[11]
                        vectors.append(v1)
            if '# number of ellipses:' in line:
                values = line.split()
                if len(values) != 7:
                    return goodfile
                ellipses = []
                if int(values[4]) > 0:
                    ind2 = self.line_range(lines, ind1)
                    npoints = ind2 - ind1 - 1
                    for loop in range(ind1+1, ind2):
                        e1 = {}
                        line = lines[loop].strip('\n')
                        values = line.split('\t')
                        if len(values) != 10:
                            return goodfile
                        e1['xstart'] = float(values[0])
                        e1['ystart'] = float(values[1])
                        e1['major'] = float(values[2])
                        e1['minor'] = float(values[3])
                        e1['rotation'] = float(values[4])
                        e1['plot'] = int(values[5])
                        e1['line_colour'] = values[6]
                        e1['line_type'] = values[7]
                        e1['line_thickness'] = float(values[8])
                        e1['fill_colour'] = values[9]
                        ellipses.append(e1)
            if '# number of boxes:' in line:
                values = line.split()
                if len(values) != 7:
                    return goodfile
                boxes = []
                if int(values[4]) > 0:
                    ind2 = self.line_range(lines, ind1)
                    npoints = ind2 - ind1 - 1
                    for loop in range(ind1+1, ind2):
                        b1 = {}
                        line = lines[loop].strip('\n')
                        values = line.split('\t')
                        if len(values) != 10:
                            return goodfile
                        b1['xstart'] = float(values[0])
                        b1['ystart'] = float(values[1])
                        b1['xend'] = float(values[2])
                        b1['yend'] = float(values[3])
                        b1['rotation'] = float(values[4])
                        b1['plot'] = int(values[5])
                        b1['line_colour'] = values[6]
                        b1['line_type'] = values[7]
                        b1['line_thickness'] = float(values[8])
                        b1['fill_colour'] = values[9]
                        boxes.append(b1)
            if '# legend values:' in line:
                legend = []
                for loop in range(nplots):
                    line1 = lines[ind1+1+loop].strip('\n')
                    values = line1.split('\t')
                    if len(values) != 7:
                        return goodfile
                    lg = {}
                    if values[1] == 'None':
                        lg['legend_variable_value'] = None
                        lg['legend_frame_value'] = None
                        lg['legend_position'] = None
                        lg['legend_option_value'] = None
                        lg['legend_user_position'] = None
                    else:
                        lg['legend_variable_value'] = int(values[1])
                        lg['legend_frame_value'] = int(values[2])
                        lg['legend_option_value'] = values[3]
                        lg['legend_position'] = values[4]
                        if values[5] == 'None':
                            lg['legend_user_position'] = None
                        else:
                            lg['legend_user_xposition'] = [float(values[5]),
                                                           float(values[6])]
                    legend.append(lg)
            if '# share_axis' in line:
                share_axis = []
                line1 = lines[ind1+1].strip('\n')
                values = line1.split('\t')
                for loop in range(len(values)):
                    share_axis.append(int(values[loop]))
            if '# equal_aspect' in line:
                equal_aspect = []
                line1 = lines[ind1+1].strip('\n')
                values = line1.split('\t')
                for loop in range(len(values)):
                    if 'True' in values[loop]:
                        equal_aspect.append(True)
                    else:
                        equal_aspect.append(False)
        goodfile = True
        if flag:
            self.maxsets = maxsets
            for loop in range(nsets):
                ind1 = 3 * loop
                self.add_set(xdata[ind1], ydata[ind1],
                             xlowerror=xdata[ind1+1],
                             xhigherror=xdata[ind1+2],
                             ylowerror=ydata[ind1+1],
                             yhigherror=ydata[ind1+2])
                self.xdata[loop]['minimum'] = set_parameters[loop]['xminimum']
                self.ydata[loop]['minimum'] = set_parameters[loop]['yminimum']
                self.xdata[loop]['maximum'] = set_parameters[loop]['xmaximum']
                self.ydata[loop]['maximum'] = set_parameters[loop]['ymaximum']
                self.xdata[loop]['errors'] = set_parameters[loop]['xerrors']
                self.ydata[loop]['errors'] = set_parameters[loop]['yerrors']
                self.xdata[loop]['legend'] = set_parameters[loop]['xlegend']
                self.ydata[loop]['legend'] = set_parameters[loop]['ylegend']
                self.set_properties[loop]['symbol'] = \
                    set_properties[loop]['symbol']
                self.set_properties[loop]['symbolsize'] = \
                    set_properties[loop]['symbolsize']
                self.set_properties[loop]['linestyle'] = \
                    set_properties[loop]['linestyle']
                self.set_properties[loop]['linewidth'] = \
                    set_properties[loop]['linewidth']
                self.set_properties[loop]['colour'] = \
                    set_properties[loop]['colour']
                self.set_properties[loop]['label'] = \
                    set_properties[loop]['label']
                self.set_properties[loop]['xmin'] = \
                    set_properties[loop]['xmin']
                self.set_properties[loop]['xmax'] = \
                    set_properties[loop]['xmax']
                self.set_properties[loop]['ymin'] = \
                    set_properties[loop]['ymin']
                self.set_properties[loop]['ymax'] = \
                    set_properties[loop]['ymax']
                self.set_properties[loop]['display'] = \
                    set_properties[loop]['display']
                self.set_properties[loop]['errors'] = \
                    set_properties[loop]['errors']
                self.set_properties[loop]['legend'] = \
                    set_properties[loop]['legend']
                self.set_properties[loop]['plot'] = \
                    set_properties[loop]['plot']
                self.plot_margin = plot_margin
            self.nsets = nsets
            self.title = []
            self.plot_frame = []
            self.plot_range = []
            self.fontname = []
            self.fontsize = []
            self.fontweight = []
            self.xparameters = []
            self.yparameters = []
            for loop in range(nplots):
                self.title.append(plot_parameters[loop]['title'])
                self.plot_frame.append(plot_parameters[loop]['frame'])
                self.fontname.append(plot_parameters[loop]['fontname'])
                self.fontsize.append(plot_parameters[loop]['fontsize'])
                self.fontweight.append(plot_parameters[loop]['fontweight'])
                self.xparameters.append({'label': ' ', 'minimum': 0.0,
                                         'maximum': 1.0, 'major': 0.1,
                                         'minor': 0.1, 'logarithmic': 0,
                                         'invert': 0, 'hide': 0,
                                         'hideticks': 0, 'hidelabels': 0,
                                         'hybridlog': 0, 'inverseticks': 0,
                                         'ticklength': 6, 'bothticks': 0,
                                         'minorticks': 0,
                                         'oppositeaxis': 0})
                self.yparameters.append({'label': ' ', 'minimum': 0.0,
                                         'maximum': 1.0, 'major': 0.1,
                                         'minor': 0.1, 'logarithmic': 0,
                                         'invert': 0, 'hide': 0,
                                         'hideticks': 0, 'hidelabels': 0,
                                         'hybridlog': 0, 'inverseticks': 0,
                                         'ticklength': 6, 'bothticks': 0,
                                         'minorticks': 0,
                                         'oppositeaxis': 0})
                self.plot_range.append(plot_parameters[loop]['plot_range'])
                self.original_range.append(True)
                self.xparameters[loop]['label'] = \
                    plot_parameters[loop]['xparameters']['label']
                self.xparameters[loop]['minimum'] = \
                    plot_parameters[loop]['xparameters']['minimum']
                self.xparameters[loop]['maximum'] = \
                    plot_parameters[loop]['xparameters']['maximum']
                self.xparameters[loop]['major'] = \
                    plot_parameters[loop]['xparameters']['major']
                self.xparameters[loop]['minor'] = \
                    plot_parameters[loop]['xparameters']['minor']
                self.xparameters[loop]['logarithmic'] = \
                    plot_parameters[loop]['xparameters']['logarithmic']
                self.xparameters[loop]['invert'] = \
                    plot_parameters[loop]['xparameters']['invert']
                self.xparameters[loop]['hide'] = \
                    plot_parameters[loop]['xparameters']['hide']
                self.xparameters[loop]['hideticks'] = \
                    plot_parameters[loop]['xparameters']['hideticks']
                self.xparameters[loop]['hidelabels'] = \
                    plot_parameters[loop]['xparameters']['hidelabels']
                self.xparameters[loop]['hybridlog'] = \
                    plot_parameters[loop]['xparameters']['hybridlog']
                self.xparameters[loop]['inverseticks'] = \
                    plot_parameters[loop]['xparameters']['inverseticks']
                self.xparameters[loop]['ticklength'] = \
                    plot_parameters[loop]['xparameters']['ticklength']
                self.xparameters[loop]['bothticks'] = \
                    plot_parameters[loop]['xparameters']['bothticks']
                self.xparameters[loop]['minorticks'] = \
                    plot_parameters[loop]['xparameters']['minorticks']
                self.xparameters[loop]['oppositeaxis'] = \
                    plot_parameters[loop]['xparameters']['oppositeaxis']
                self.yparameters[loop]['label'] = \
                    plot_parameters[loop]['yparameters']['label']
                self.yparameters[loop]['minimum'] = \
                    plot_parameters[loop]['yparameters']['minimum']
                self.yparameters[loop]['maximum'] = \
                    plot_parameters[loop]['yparameters']['maximum']
                self.yparameters[loop]['major'] = \
                    plot_parameters[loop]['yparameters']['major']
                self.yparameters[loop]['minor'] = \
                    plot_parameters[loop]['yparameters']['minor']
                self.yparameters[loop]['logarithmic'] = \
                    plot_parameters[loop]['yparameters']['logarithmic']
                self.yparameters[loop]['invert'] = \
                    plot_parameters[loop]['yparameters']['invert']
                self.yparameters[loop]['hide'] = \
                    plot_parameters[loop]['yparameters']['hide']
                self.yparameters[loop]['hideticks'] = \
                    plot_parameters[loop]['yparameters']['hideticks']
                self.yparameters[loop]['hidelabels'] = \
                    plot_parameters[loop]['yparameters']['hidelabels']
                self.yparameters[loop]['hybridlog'] = \
                    plot_parameters[loop]['yparameters']['hybridlog']
                self.yparameters[loop]['inverseticks'] = \
                    plot_parameters[loop]['yparameters']['inverseticks']
                self.yparameters[loop]['ticklength'] = \
                    plot_parameters[loop]['yparameters']['ticklength']
                self.yparameters[loop]['bothticks'] = \
                    plot_parameters[loop]['yparameters']['bothticks']
                self.yparameters[loop]['minorticks'] = \
                    plot_parameters[loop]['yparameters']['minorticks']
                self.yparameters[loop]['oppositeaxis'] = \
                    plot_parameters[loop]['yparameters']['oppositeaxis']
                self.original_range[loop] = \
                    plot_parameters[loop]['original_range']
            if len(labels) > 0:
                self.number_of_labels = len(labels)
                for loop in range(len(labels)):
                    self.plot_labels[loop]['xposition'] = \
                        labels[loop]['xposition']
                    self.plot_labels[loop]['yposition'] = \
                        labels[loop]['yposition']
                    self.plot_labels[loop]['plot'] = labels[loop]['plot']
                    self.plot_labels[loop]['labelstring'] = \
                        labels[loop]['labelstring']
                    self.plot_labels[loop]['colour'] = labels[loop]['colour']
                    self.plot_labels[loop]['size'] = labels[loop]['size']
                    self.plot_labels[loop]['font'] = labels[loop]['font']
                    self.plot_labels[loop]['fontweight'] = \
                        labels[loop]['fontweight']
            if len(plines) > 0:
                self.number_of_lines = len(plines)
                for loop in range(len(plines)):
                    self.plot_lines[loop]['xstart'] = plines[loop]['xstart']
                    self.plot_lines[loop]['ystart'] = plines[loop]['ystart']
                    self.plot_lines[loop]['xend'] = plines[loop]['xend']
                    self.plot_lines[loop]['yend'] = plines[loop]['yend']
                    self.plot_lines[loop]['plot'] = plines[loop]['plot']
                    self.plot_lines[loop]['line_type'] = \
                        plines[loop]['line_type']
                    self.plot_lines[loop]['line_colour'] = \
                        plines[loop]['line_colour']
                    self.plot_lines[loop]['line_width'] = \
                        plines[loop]['line_width']
            if len(vectors) > 0:
                self.number_of_vectors = len(vectors)
                for loop in range(len(vectors)):
                    self.plot_vectors[loop]['xstart'] = vectors[loop]['xstart']
                    self.plot_vectors[loop]['ystart'] = vectors[loop]['ystart']
                    self.plot_vectors[loop]['xend'] = vectors[loop]['xend']
                    self.plot_vectors[loop]['yend'] = vectors[loop]['yend']
                    self.plot_vectors[loop]['delx'] = vectors[loop]['delx']
                    self.plot_vectors[loop]['dely'] = vectors[loop]['dely']
                    self.plot_vectors[loop]['plot'] = vectors[loop]['plot']
                    self.plot_vectors[loop]['line_type'] = \
                        vectors[loop]['line_type']
                    self.plot_vectors[loop]['line_colour'] = \
                        vectors[loop]['line_colour']
                    self.plot_vectors[loop]['line_thickness'] = \
                        vectors[loop]['line_thickness']
            if len(ellipses) > 0:
                self.number_of_ellipses = len(ellipses)
                for loop in range(len(ellipses)):
                    self.plot_ellipses[loop]['xstart'] = \
                        ellipses[loop]['xstart']
                    self.plot_ellipses[loop]['ystart'] = \
                        ellipses[loop]['ystart']
                    self.plot_ellipses[loop]['major'] = ellipses[loop]['major']
                    self.plot_ellipses[loop]['minor'] = ellipses[loop]['minor']
                    self.plot_ellipses[loop]['rotation'] = \
                        ellipses[loop]['rotation']
                    self.plot_ellipses[loop]['plot'] = ellipses[loop]['plot']
                    self.plot_ellipses[loop]['line_type'] = \
                        ellipses[loop]['line_type']
                    self.plot_ellipses[loop]['line_colour'] = \
                        ellipses[loop]['line_colour']
                    self.plot_ellipses[loop]['line_thickness'] = \
                        ellipses[loop]['line_thickness']
                    self.plot_ellipses[loop]['fill_colour'] = \
                        ellipses[loop]['fill_colour']
            if len(boxes) > 0:
                self.number_of_boxes = len(boxes)
                for loop in range(len(boxes)):
                    self.plot_boxes[loop]['xstart'] = boxes[loop]['xstart']
                    self.plot_boxes[loop]['ystart'] = boxes[loop]['ystart']
                    self.plot_boxes[loop]['xend'] = boxes[loop]['xend']
                    self.plot_boxes[loop]['yend'] = boxes[loop]['yend']
                    self.plot_boxes[loop]['rotation'] = boxes[loop]['rotation']
                    self.plot_boxes[loop]['plot'] = boxes[loop]['plot']
                    self.plot_boxes[loop]['line_type'] = \
                        boxes[loop]['line_type']
                    self.plot_boxes[loop]['line_colour'] = \
                        boxes[loop]['line_colour']
                    self.plot_boxes[loop]['line_thickness'] = \
                        boxes[loop]['line_thickness']
                    self.plot_boxes[loop]['fill_colour'] = \
                        boxes[loop]['fill_colour']
#            nplots = self.nxplots * self.nyplots
            self.legend_position = []
            self.legend_user_position = []
            for loop in range(nplots):
                self.legend_position.append(legend[loop]['legend_position'])
                self.legend_user_position.append(legend[loop][
                    'legend_user_position'])
                if legend[loop]['legend_variable_value'] is None:
                    if loop < len(self.legend_variable):
                        try:
                            self.legend_variable[loop].set(0)
                        except:
                            self.legend_variable.append(None)
                    else:
                        self.legend_variable.append(None)
                else:
                    if loop < len(self.legend_variable):
                        try:
                            self.legend_variable[loop].set(
                                legend[loop]['legend_variable_value'])
                        except:
                            self.legend_variable[loop] = Tk.IntVar()
                            self.legend_variable[loop].set(
                                legend[loop]['legend_variable_value'])
                    else:
                        self.legend_variable.append(Tk.IntVar())
                        self.legend_variable[loop].set(
                            legend[loop]['legend_variable_value'])
                if legend[loop]['legend_frame_value'] is None:
                    if loop < len(self.legend_frame):
                        try:
                            self.legend_frame[loop].set(0)
                        except:
                            self.legend_frame.append(None)
                    else:
                        self.legend_frame.append(None)
                else:
                    if loop < len(self.legend_frame):
                        try:
                            self.legend_frame[loop].set(
                                legend[loop]['legend_frame_value'])
                        except:
                            self.legend_frame[loop] = Tk.IntVar()
                            self.legend_frame[loop].set(
                                legend[loop]['legend_frame_value'])
                    else:
                        self.legend_frame.append(Tk.IntVar())
                        self.legend_frame[loop].set(
                            legend[loop]['legend_frame_value'])
                if legend[loop]['legend_option_value'] is None:
                    if loop < len(self.legend_options):
                        try:
                            self.legend_options[loop].get()
                        except:
                            self.legend_options.append(None)
                    else:
                        self.legend_options.append(None)
                else:
                    if loop < len(self.legend_options):
                        try:
                            self.legend_options[loop].set(
                                legend[loop]['legend_option_value'])
                        except:
                            self.legend_options[loop] = Tk.StringVar()
                            self.legend_options[loop].set(
                                legend[loop]['legend_option_value'])
                    else:
                        self.legend_options.append(Tk.StringVar())
                        self.legend_options[loop].set(
                            legend[loop]['legend_option_value'])
            self.share_axis = deepcopy(share_axis)
            if equal_aspect is None:
                self.equal_aspect = []
                for loop in range(nplots):
                    self.equal_aspect.append(False)
            else:
                self.equal_aspect = deepcopy(equal_aspect)
            self.number_of_plots = 1*nplots
            self.hide_subplot = deepcopy(hide_subplot)
        return goodfile

    def line_range(self, lines, ind1):
        """
        Find a range of data lines within a line list.

        Given an input line list and a starting index, subsequent lines are
        examined to see where the next comment line is.  Comment lines are
        assumed to start with the # character.  Lines that are not comments
        are assumed to be data lines.  The index of the next comment line is
        returned, or the index that gives a range to the end of the line list
        where there is no such comment line after the index specified.

        Parameters
        ----------
            lines :  A list of input lines (assumed to be from the readlines()
                     function)

            ind1 :   A starting index in the list of lines

        Returns
        -------
            n1 : an integer value for the next comment line (assumed to
                 start with '#') in the list of input lines, or the index
                 for the next of the line list if no other comment line is
                 found
        """
        for n1 in range(ind1+1, len(lines)):
            if '#' in lines[n1][0:1]:
                return n1
        return len(lines)

    def create_data_set(self):
        """
        Open a window to define a function for making a data set.

        This routine makes the window that allows a limited capability to
        create data sets via a defined function.

        Parameters
        ----------
            None

        Returns
        -------
           Nothing

        """
        function_window = Tk.Toplevel()
        function_window.title('Create Data Set')
        holder = Tk.Frame(function_window)
        holder.pack(side=Tk.TOP)
        h1 = Tk.Frame(holder)
        h1.pack(side=Tk.TOP)
        label1 = Tk.Label(h1, text='    Start values at ')
        label1.pack(side=Tk.LEFT)
        self.start_value_field = Tk.Entry(h1, width=12)
        self.start_value_field.pack(side=Tk.LEFT)
        label2 = Tk.Label(h1, text='    Stop values at ')
        label2.pack(side=Tk.LEFT)
        self.stop_value_field = Tk.Entry(h1, width=12)
        self.stop_value_field.pack(side=Tk.LEFT)
        label3 = Tk.Label(h1, text=' Number of values/step ')
        label3.pack(side=Tk.LEFT)
        self.number_of_values_field = Tk.Entry(h1, width=6)
        self.number_of_values_field.pack(side=Tk.LEFT)
        h2 = Tk.Frame(holder)
        h2.pack()
        self.sequence_option = Tk.IntVar()
        lab1 = Tk.Label(
            h2,
            text='Spacing (logarithmic only works for postive range values): ')
        lab1.pack(side=Tk.LEFT)
        b1 = Tk.Radiobutton(
            h2, text='linear', variable=self.sequence_option,
            value=0)
        b1.pack(side=Tk.LEFT)
        b2 = Tk.Radiobutton(
            h2, text='logarithmic', variable=self.sequence_option,
            value=1)
        self.sequence_option.set(0)
        b2.pack(side=Tk.LEFT)
        h3 = Tk.Frame(holder)
        h3.pack(side=Tk.TOP)
        label1 = Tk.Label(h3, text=' x function: ')
        label1.grid(column=0, row=0)
        self.xfunction = Tk.Entry(h3, width=50)
        self.xfunction.grid(column=1, row=0)
        label2 = Tk.Label(h3, text=' y function: ')
        label2.grid(column=0, row=1)
        self.yfunction = Tk.Entry(h3, width=50)
        self.yfunction.grid(column=1, row=1)
        self.xfunction.insert(0, '$t')
        self.yfunction.insert(0, '$x*$x')
        label3 = Tk.Label(h3, text='Enter the function you want, '
                          + 'where $t represents the sequence of values'
                          + ' defined at top\nand either $x or $y refers '
                          + 'to the variables.  Note that while $x can be '
                          + 'used to define y \nand $y can be used to '
                          + 'define x one cannot use x or y in its own '
                          + 'function field.\n\nOne can use numpy and '
                          + 'math functions within the definition.',
                          justify=Tk.LEFT)
        label3.grid(column=0, row=2, columnspan=2)
        h4 = Tk.Frame(holder)
        h4.pack(side=Tk.TOP)
        button1 = Tk.Button(h4, text='Apply', command=self.parse_function)
        button1.pack(side=Tk.LEFT)
        label1 = Tk.Label(h4, text='   ')
        label1.pack(side=Tk.LEFT)
        button2 = Tk.Button(h4, text='Close', command=function_window.destroy)
        button2.pack(side=Tk.LEFT)

    def parse_function(self):
        """
        Read parameters and try to parse to make a data set.

        This routine reads the parameters to make a new set and attempts to
        evaluate them using a parser function.  Whilst one could use the
        eval() function this is considered as too dangerous.

        No values are passed to this routine or returned from it.
        """
        try:
            v1 = float(self.start_value_field.get())
            v2 = float(self.stop_value_field.get())
            if v1 > v2:
                temp = v1
                v1 = v2
                v2 = temp
            try:
                instring = self.number_of_values_field.get()
                if '.' in instring:
                    n1 = 0
                    step = float(instring)
                else:
                    n1 = int(instring)
                    step = 0.
            except:
                n1 = 0
                step = float(self.number_of_values_field.get())
            # if the value n1 is not an integer, assume it is a step size
            option = self.sequence_option.get()
            if (option == 0) or ((v1 <= 0.) or (v2 <= 0)):
                if n1 > 1:
                    step = (v2-v1)/(n1-1)
                else:
                    if step == 0.:
                        step = v2 - v1
                if (v1 == v2) or (step <= 0.):
                    tkinter.messagebox.showinfo(
                        'Error', 
                        'There was some error trying to generate the sets.  (1)')
                seq = numpy.arange(v1, v2+step, step)
                if (len(seq) > n1) and (n1 > 1):
                    seq = seq[0:n1]
            else:
                datarange = math.log10(v2/v1)
                if n1 > 1:
                    step = datarange/(n1-1)
                if n1 > 2:
                    seq = numpy.arange(0., datarange+step, step)
                else:
                    seq = numpy.asarray([0., datarange])
                seq = v1*numpy.power(10., seq)
                if (len(seq) > n1) and (n1 > 1):
                    seq = seq[0:n1]
            # seq is the starting sequence; read the function strings
            xstring = self.xfunction.get()
            ystring = self.yfunction.get()
            xstring = xstring.replace('$t', 'seq')
            ystring = ystring.replace('$t', 'seq')
            xstring = xstring.replace('$x', 'x')
            ystring = ystring.replace('$x', 'x')
            xstring = xstring.replace('$y', 'y')
            ystring = ystring.replace('$y', 'y')
            try:
                xvalues = self.my_eval(xstring, seq=seq, xvalues=None,
                                       yvalues=None)
                yvalues = self.my_eval(ystring, seq=seq, xvalues=xvalues,
                                       yvalues=None)
            except:
                yvalues = self.my_eval(ystring, seq=seq, xvalues=None,
                                       yvalues=None)
                xvalues = self.my_eval(xstring, seq=seq, xvalues=None,
                                       yvalues=yvalues)
            # deal with the case where one of x or y is entered as a constant
            try:
                n = len(xvalues)
            except:
                xvalues = numpy.asarray([xvalues, ])
            try:
                n = len(yvalues)
            except:
                yvalues = numpy.asarray([yvalues, ])
            if (len(xvalues) == 1) and (len(yvalues) > 1):
                xvalues = yvalues*0.+xvalues
            if (len(yvalues) == 1) and (len(xvalues) > 1):
                yvalues = xvalues*0.+yvalues
            if len(xvalues) != len(yvalues):
                tkinter.messagebox.showinfo(
                    'Error',
                    'There was some error trying to generate the sets. (2)')
                return
            if (xvalues is not None) and (yvalues is not None):
                self.add_set(xvalues, yvalues, current_plot=self.current_plot)
                self.make_plot()
            else:
                tkinter.messagebox.showinfo(
                    'Error',
                    'There was some error trying to generate the sets. (3)')
        except:
            tkinter.messagebox.showinfo(
                'Error',
                'There was some error trying to generate the sets. (4)')

    def my_eval(self, inputstring, seq, xvalues=None, yvalues=None):
        """
        Evaluate a string as an expression to make a data set.

        This routine attempts to evaluate a string as an expression.
        It uses the python "eval" function.  To guard against bad inputs,
        only numpy, math and builtin functions can be used in the
        transformation.

        Parameters
        ----------
            inputstring  a string that defines the new data set

            seq : a numpy vector of floating point or integer values,
                  nominally a sequence of values when the data creation
                  option is used, which could be another numpy array in
                  the transformation case

            xvalues :  optionally, the x data values in a set, a numpy
                       floating point vector

            yvalues :  optionally, the y data values in a set, a numpy
                       floating point vector

        Returns
        -------
            values :   a new numpy vector of floating point values calculated
                       from the input numpy arrays and the string defining the
                       function; or None if there is an issue

        Note: the three numpy arrays "seq", "xvalues", and "yvalues" need
        to be one dimensional and of the same lengths

        The python "eval" command is used here.  To avoid issues with this
        being used to run arbitrary commands, only the __builtin__, math,
        and numpy packages are available to the eval command upon execution.
        The assumption is that math and numpy have been imported in the main
        code (and that numpy is not abbreviated as "np" at import).

        """
        sh1 = seq.shape
        try:
            sh2 = xvalues.shape
        except:
            sh2 = seq.shape
        try:
            sh3 = yvalues.shape
        except:
            sh3 = seq.shape
        if (sh1 != sh2) or (sh2 != sh3) or (len(sh1) > 1):
            return None
        # check the input string for command elements that could cause issues
        if ('import' in inputstring) or ('os.' in inputstring) or \
           ('eval' in inputstring) or ('exec' in inputstring):
            return None
        str1 = inputstring.replace('np.', 'numpy.')
        try:
            # get the global environment, extract the three items allowed here
            global1 = globals()
            global2 = {}
            global2['__builtins__'] = global1['__builtins__']
            global2['math'] = global1['math']
            global2['numpy'] = global1['numpy']
            # define local variables, s, x, and y; only these will be
            # available in eval if they are actually defined in the call....
            local1 = {}
            s = numpy.copy(seq)
            local1['seq'] = s
            if xvalues is not None:
                x = numpy.copy(xvalues)
                local1['x'] = x
            if yvalues is not None:
                y = numpy.copy(yvalues)
                local1['y'] = y
            values = eval(str1, global2, local1)
            return values
        except:
            return None

    def create_data_set_by_editor(self):
        """
        Allow the user to make a set by entering numbers in a window.


        Parameters
        ----------
            None

        Returns
        -------
            No values are returned by the routine.

        """
        if self.data_entry_window is not None:
            return
        self.data_entry_window = Tk.Toplevel()
        self.data_entry_window.title('Enter Data Set Values')
        holder = Tk.Frame(self.data_entry_window)
        holder.pack(side=Tk.TOP)
        holder.config(bg='black')
        self.data_text = ScrolledText(holder, height=40, width=80,
                                      wrap=Tk.NONE, relief="solid")
        self.data_text.config(font=('courier', 16))
        self.data_text.pack(side=Tk.TOP, padx=10, pady=10)
        bframe = Tk.Frame(self.data_entry_window)
        bframe.pack()
        set_button = Tk.Button(
            bframe, text="Apply",
            command=self.apply_data_input)
        set_button.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            bframe, text="Close",
            command=lambda: self.close_data_window(
                self.data_entry_window))
        close_button.pack(side=Tk.LEFT)

    def apply_data_input(self):
        """
        Read the values from a data input text window and parse these to
        a data set for the plots.

        Returns
        -------
        None.

        """
        text = self.data_text.get("1.0", Tk.END)
        xvalues, dxvalues1, dxvalues2, yvalues, dyvalues1, dyvalues2,\
            errorflag = parse_text(text)
        try:
            xvalues = numpy.asarray(xvalues)
            yvalues = numpy.asarray(yvalues)
            dxvalues1 = numpy.asarray(dxvalues1)
            dxvalues2 = numpy.asarray(dxvalues2)
            dyvalues1 = numpy.asarray(dyvalues1)
            dyvalues2 = numpy.asarray(dyvalues2)
            if len(xvalues) < 1:
                tkinter.messagebox.showinfo(
                    'error',
                    'Unable to parse text from entry widget (2)')
                return
            self.add_set(xvalues, yvalues, xlowerror=dxvalues1,
                         xhigherror=dxvalues2, ylowerror=dyvalues1,
                         yhigherror=dyvalues2, xerrorflag=errorflag,
                         yerrorflag=errorflag,
                         labelstring='Set from editor',
                         current_plot=self.current_plot)
            self.make_plot()
        except:
            tkinter.messagebox.showinfo(
                'error',
                'Unable to parse text from entry widget')
            return

    def remove_line(self):
        """
        Remove the last line from the list of plot elements.

        This routine sets the line counter down by 1 so the last defined
        line is no longer plotted by the code.  The next new line
        subsequently defined will then overwrite the last current one.

        Parameters
        ----------
            None

        Returns
        -------
            No values are returned by the routine.  The self.number_of_lines
            variable is changed.

        """
        if self.number_of_lines > 0:
            self.number_of_lines = self.number_of_lines - 1
            self.make_plot()

    def remove_vector(self):
        """
        Remove the last vector from the list of plot elements.

        This routine sets the vector counter down by 1 so the last defined
        vector is no longer plotted by the code.  The next new vector
        subsequently defined will then overwrite the last current one.

        Parameters
        ----------
            None

        Returns
        -------
            No values are returned by the routine.  The self.number_of_vectors
            variable is changed.

        """
        if self.number_of_vectors > 0:
            self.number_of_vectors = self.number_of_vectors - 1
            self.make_plot()

    def remove_ellipse(self):
        """
        Remove the last ellipse from the list of plot elements.

        This routine sets the ellipse counter down by 1 so the last defined
        ellipse is no longer plotted by the code.  Any new ellipses
        subsquently defined will then overwrite the last current one.

        Parameters
        ----------
            None

        Returns
        -------
            No values are returned.  The self.number_of_ellipses variable
            is changed.

        """
        if self.number_of_ellipses > 0:
            self.number_of_ellipses = self.number_of_ellipses - 1
            self.make_plot()

    def remove_box(self):
        """
        Remove the last box from the list of plot elements.

        This routine sets the box counter down by 1 so the last box is no
        longer plotted by the code.  Any new box subsequently defined will
        then overwrite the last existing one.

        Parameters
        ----------
            None

        Returns
        -------
            No values are returned by the routine.  The self.number_of_boxes
            variable is changed.

        """
        if self.number_of_boxes > 0:
            self.number_of_boxes = self.number_of_boxes - 1
            self.make_plot()

    def remove_all_lines(self):
        """
        Remove all lines from the plot elements.

        This routine sets the line counter variable to zero so that any
        defined lines are not plotted, and subsequently defined lines will
        overwrite any that are currently in the self.plot_lines variable.

        Parameters
        ----------
            None

        Returns
        -------
            No values are returned by the routine.  The self.number_of_lines
            variable is changed.

        """
        if self.number_of_lines > 0:
            self.number_of_lines = 0
            self.positions = []
            self.make_plot()

    def remove_all_vectors(self):
        """
        Remove all vectors from the plot elements.

        This routine sets the vector counter variable to zero so that any
        defined vectors are not plotted, and subsequently defined vectors
        will overwrite any that are currently in the self.plot_vectors
        variable.

        Parameters
        ----------
            None

        Returns
        -------
            No values are returned by the routine.  The self.number_of_vectors
            variable is changed.

        """
        if self.number_of_vectors > 0:
            self.number_of_vectors = 0
            self.positions = []
            self.make_plot()

    def remove_all_ellipses(self):
        """
        Remove all ellipses from the plot elements.

        This routine sets the ellipse counter variable to zero so that any
        defined ellipses are not plotted, and subsequently defined ellipses
        will overwrite any that are currently in the self.plot_ellipses
        variable.

        Parameters
        ----------
            None

        Returns
        -------
            No values are returned.  The self.number_of_ellipses variable is
            changed.

        """
        if self.number_of_ellipses > 0:
            self.number_of_ellipses = 0
            self.positions = []
            self.make_plot()

    def remove_all_boxes(self):
        """
        Remove all boxes from the plot elements.

        This routine sets the box counter variable to zero so that any
        defined boxes are not plotted, and subsequently defined boxes
        will overwrite any that are currently in the self.plot_boxes variable.

        Parameters
        ----------
            None

        Returns
        -------
            No values are returned by the routine.  The self.number_of_boxes
            variable is changed.

        """
        if self.number_of_boxes > 0:
            self.number_of_boxes = 0
            self.positions = []
            self.make_plot()

    def add_line(self):
        """
        Set a flag for marking line points on the plot.

        This routine sets a flag so that mouse clicks are used for defining
        a line.

        Parameters
        ----------
            None

        Returns
        -------
           Nothing

        """
        self.line_flag = True

    def add_vector(self):
        """
        Set a flag for marking vector points on the plot.

        This routine sets a flag so that mouse clicks are used for defining
        a vector.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        self.vector_flag = True

    def add_ellipse(self):
        """
        Set a flag for marking ellipse points on the plot.

        This routine sets a flag so that mouse clicks are used for defining
        an ellipse.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        self.ellipse_flag = True

    def add_box(self):
        """
        Set a flag for marking box points on the plot.

        This routine sets a flag so that mouse clicks are used for defining
        a box.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        self.box_flag = True

    def set_font(self):
        """
        Bring up a window to set the plot font.

        This window allows one to set the font properties for the axis
        labels and titles.  The values are applied to all axis labels
        and the title in a given plot.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        font_window = Tk.Toplevel()
        font_window.title('Font Values')
        holder = Tk.Frame(font_window)
        holder.pack(side=Tk.TOP)
        fontnames = ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace',
                     'times new roman']
        fontsizes = ['8', '9', '10', '11', '12', '13', '14', '16', '18',
                     '20', '24', '30']
        fontweights = ['ultralight', 'light', 'normal', 'regular', 'book',
                       'medium', 'roman', 'semibold', 'demibold', 'demi',
                       'bold', 'heavy', 'extra bold', 'black']
        label = Tk.Label(holder, text='Font Name:')
        label.grid(column=0, row=0)
        label = Tk.Label(holder, text='Font Size:')
        label.grid(column=0, row=1)
        label = Tk.Label(holder, text='Font Weight:')
        label.grid(column=0, row=2)
        self.font_name_list = Tk.ttk.Combobox(holder, width=20)
        self.font_name_list.grid(column=1, row=0)
        self.font_name_list['values'] = fontnames
        self.font_name_list.set(self.fontname[self.current_plot-1])
        self.font_size_list = Tk.ttk.Combobox(holder, width=20)
        self.font_size_list.grid(column=1, row=1)
        self.font_size_list['values'] = fontsizes
        self.font_size_list.set(self.fontsize[self.current_plot-1])
        self.font_weight_list = Tk.ttk.Combobox(holder, width=20)
        self.font_weight_list.grid(column=1, row=2)
        self.font_weight_list['values'] = fontweights
        self.font_weight_list.set(self.fontweight[self.current_plot-1])
        bholder = Tk.Frame(font_window)
        bholder.pack(side=Tk.TOP)
        set_button = Tk.Button(bholder, text='Set',
                               command=self.set_font_values)
        set_button.pack()
        set_all_button = Tk.Button(bholder, text='Set for All',
                                   command=self.set_font_values_all)
        set_all_button.pack()
        close_button = Tk.Button(bholder, text='Close Window',
                                 command=font_window.destroy)
        close_button.pack()

    def set_font_values_all(self):
        """
        Read the font field and apply them to all plots.

        This routine reads the font fields and saves them to the internal
        variables.  It then calls the routine to re-plot.  This is done for
        all current plots (including hidden ones).

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        plotnumber = self.current_plot
        for loop in range(self.number_of_plots):
            self.current_plot = loop + 1
            self.set_font_values()
        self.current_plot = plotnumber

    def set_font_values(self):
        """
        Read the font fields and apply them to the plot.

        This routine reads the font fields and saves them to the internal
        variables.  It then calls the routine to re-plot.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        self.fontname[self.current_plot-1] = self.font_name_list.get()
        self.fontsize[self.current_plot-1] = self.font_size_list.get()
        self.fontweight[self.current_plot-1] = self.font_weight_list.get()
        self.make_plot()

    def edit_labels(self):
        """
        Bring up a window to edit the labels in the plot.

        This routine produces a text box in a window, within which one can
        edit the label values.  If no labels are defined the routine just
        exits with no action.

        Labels are presented one per line, with the parameter values
        separated by tab symbols.  One can edit the values within the text
        window and then these are applied when one clicks on the "Close
        Window" button.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.number_of_labels == 0:
            return
        str1 = 'Edit values below: fields are separated by tab ' \
            + 'characters\n x position    y position  plot   ' \
            + 'label\n--------------------------------------\n'
        for loop in range(self.number_of_labels):
            if self.xparameters[self.plot_labels[loop]['plot']-1]['hybridlog']:
                xpos1 = inverse_hybrid_transform(
                    self.plot_labels[loop]['xposition'])
            else:
                xpos1 = self.plot_labels[loop]['xposition']
            if self.yparameters[self.plot_labels[loop]['plot']-1]['hybridlog']:
                ypos1 = inverse_hybrid_transform(
                    self.plot_labels[loop]['yposition'])
            else:
                ypos1 = self.plot_labels[loop]['yposition']
            str1 = str1 + '%12.6g\t%12.6g\t%3d\t%s\t%s\t%d\t%s\t%s\n' % (
                xpos1,
                ypos1,
                self.plot_labels[loop]['plot'],
                self.plot_labels[loop]['labelstring'],
                self.plot_labels[loop]['colour'],
                self.plot_labels[loop]['size'],
                self.plot_labels[loop]['font'],
                self.plot_labels[loop]['fontweight'])
        label_window = Tk.Toplevel()
        label_window.title('Labels:')
        holder = Tk.Frame(label_window)
        holder.pack(side=Tk.TOP)
        label_message_text = ScrolledText(holder, height=40, width=90,
                                          wrap=Tk.NONE)
        label_message_text.config(font=('courier', 16, 'bold'))
        label_message_text.pack(side=Tk.TOP)
        label_message_text.insert(0.0, str1)
        bholder = Tk.Frame(label_window)
        bholder.pack(side=Tk.TOP)
        close_button = Tk.Button(
            bholder, text='Close Window',
            command=lambda: self.read_labels(label_message_text,
                                             label_window))
        close_button.pack()

    def read_labels(self, label_message_text, label_window):
        """
        Read and parse the label text field.

        This routine reads the label text field and makes the new set of
        labels and label positions.  It then applies these and closes the
        label window.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        The code does, however, change the self.plot_labels values as
        needed to match what is in the label text field.

        """
        labeltext = label_message_text.get(0.0, Tk.END)
        lines = labeltext.split('\n')
        newlabels = []
        nlabels = 0
        for line in lines:
            values = line.split('\t')
            if len(values) == 8:
                try:
                    x1 = float(values[0])
                    y1 = float(values[1])
                    nplot = int(values[2])
                    label = values[3]
                    colour = values[4]
                    size = int(values[5])
                    font = values[6]
                    fontweight = values[7]
                    label = label.strip('\n')
                    if self.xparameters[nplot-1]['hybridlog']:
                        x1 = hybrid_transform(x1)
                    if self.yparameters[nplot-1]['hybridlog']:
                        y1 = hybrid_transform(y1)
                    newlabels.append({'xposition': x1, 'yposition': y1,
                                      'labelstring': label, 'plot': nplot,
                                      'colour': colour, 'size': size,
                                      'font': font,
                                      'fontweight': fontweight})
                    nlabels = nlabels + 1
                except:
                    pass
        label_window.destroy()
        if nlabels > self.max_labels:
            self.max_labels = nlabels
        else:
            for loop in range(nlabels, self.max_labels):
                newlabels.append({'xposition': None, 'yposition': None,
                                  'labelstring': '', 'plot': 1,
                                  'colour': 'black', 'size': 12,
                                  'font': 'sans-serif',
                                  # 'font': 'times new roman',
                                  'fontweight': 'normal'})
            self.plot_labels = newlabels
            self.number_of_labels = nlabels
        self.make_plot()

    def set_label(self):
        """
        Set the flag so that a mouse click is used to position a label.

        This routine sets the label flag in response to the "Put Label"
        button.  When the flag is set, any key pressed leads to putting
        a label on the plot.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        self.label_flag = True

    def clear_labels(self):
        """
        Clear all the labels on the plot.

        This routine is called when the "Clear Labels" button is pressed.
        It asks whether the labels should be cleared, and if so all labels
        are removed.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        response = tkinter.messagebox.askyesno("Verify", "Delete all labels.")
        if response:
            self.label_flag = False
            self.plot_labels = []
            self.max_labels = 100
            for loop in range(self.max_labels):
                self.plot_labels.append({'xposition': None,
                                         'yposition': None,
                                         'labelstring': '', 'plot': 1,
                                         'colour': 'black', 'size': 12,
                                         'font': 'sans-serif',
                                         # 'font': 'times new roman',
                                         'fontweight': 'normal'})
            self.positions = []
            self.make_plot()

    def make_controls(self, parent):
        """
        Make the control area within the main window.

        This routine makes a control area within the main window, under
        frame "parent".  The overall root value is also passed here, but
        not currently used.  It may be neeed for orderly closing of the
        window depending on what is done in the main window, hence it is
        included here.

        Parameters
        ----------
            parent :   A Tk.Frame variable for the holder of the controla

        Returns
        -------
            No values are returned by this routine.

        """
        holder = Tk.Frame(parent)
        holder.pack(side=Tk.TOP)
        label1 = Tk.Label(holder, text=' ')
        label1.pack(side=Tk.TOP, fill=Tk.X)
        button1 = Tk.Button(holder, text='Plot', command=self.make_plot)
        button1.pack(side=Tk.TOP, fill=Tk.X)
        button2 = Tk.Button(holder, text='Auto scale',
                            command=self.autoscale_plot)
        button2.pack(side=Tk.TOP, fill=Tk.X)
        sl = self.separator_line(holder, 300, 25, 5, True)
        button3 = Tk.Button(holder, text='2-D Histogram',
                            command=self.make_hess_plot)
        button3.pack(side=Tk.TOP, fill=Tk.X)
        field = Tk.Frame(holder)
        field.pack(side=Tk.TOP)
        label1 = Tk.Label(field, text='Number of pixels: ')
        label1.pack(side=Tk.LEFT)
        self.npixelfield = Tk.Entry(field, width=5)
        self.npixelfield.pack(side=Tk.TOP)
        self.npixelfield.insert(0, '500')
        sl = self.separator_line(holder, 300, 25, 5, True)
        button4 = Tk.Button(holder, text='1-D Histogram',
                            command=self.make_histogram)
        button4.pack(side=Tk.TOP, fill=Tk.X)
        field = Tk.Frame(holder)
        field.pack(side=Tk.TOP)
        label1 = Tk.Label(field, text='Number of bins or bin size: ')
        label1.pack(side=Tk.LEFT)
        self.nbinfield = Tk.Entry(field, width=10)
        self.nbinfield.pack(side=Tk.TOP)
        self.nbinfield.insert(0, '500')
        self.histogramflag = Tk.IntVar()
        b1 = Tk.Frame(holder)
        self.put_yes_no(b1, self.histogramflag, ['x values', 'y values'], True)
        b1.pack(side=Tk.TOP)
        self.individualhistogramflag = Tk.IntVar()
        b1 = Tk.Frame(holder)
        self.put_yes_no(b1, self.individualhistogramflag,
                        ['all sets', 'individual sets'], True)
        b1.pack(side=Tk.TOP)
        sl = self.separator_line(holder, 300, 25, 5, True)
        self.matplotlib_rounding = Tk.IntVar()
        b1 = Tk.Frame(holder)
        label1 = Tk.Label(b1, text='Axis limits rounding algorithm: ')
        label1.pack(side=Tk.TOP)
        b1.pack(side=Tk.TOP)
        b1 = Tk.Frame(holder)
        self.put_yes_no(b1, self.matplotlib_rounding,
                        ['Matplotlib', 'Alternate'], True)
        b1.pack(side=Tk.TOP)
        sl = self.separator_line(holder, 300, 25, 5, True)
        button5 = Tk.Button(holder, text='Save as PNG',
                            command=lambda: save_png_figure(self.figure))
        button5.pack(side=Tk.TOP, fill=Tk.X)
        button6 = Tk.Button(holder, text='Save as PS',
                            command=lambda: save_ps_figure(self.figure))
        button6.pack(side=Tk.TOP, fill=Tk.X)
        button7 = Tk.Button(holder, text='Clear Current Plot',
                            command=self.clear_current_plot)
        button7.pack(side=Tk.TOP, fill=Tk.X)
        button8 = Tk.Button(holder, text='Tile Plots', command=self.tile_plots)
        button8.pack(side=Tk.TOP, fill=Tk.X)
        button9 = Tk.Button(holder, text='Set Plot',
                            command=self.set_plot_number)
        button9.pack(side=Tk.TOP, fill=Tk.X)
        button10 = Tk.Button(holder, text='Close Window',
                             command=self.root.destroy)
        button10.pack(side=Tk.TOP, fill=Tk.X)

    def make_histogram(self):
        """
        Create a histogram plot in a new window.

        This routine creates a new plot window within which a histogram
        plot is made for the current active main window plot.

        The new window has options for output of the histogram plot to a
        file.  The colours of the bars are the same as the colours of the
        points in the main plot.

        No parameters are passed to the routine, and no values are returned.
        """
        try:
            histogramwindow = Tk.Toplevel()
            histogramwindow.config(bg=BGCOL)
            optionflag = self.histogramflag.get()
            setoptionflag = self.individualhistogramflag.get()
            try:
                value = float(self.nbinfield.get())
                if value == 0:
                    tkinter.messagebox.showinfo(
                        "Error",
                        "Zero value for the number of bins/bin size."
                        + "  Check your inputs.")
                    return
                if value < 1.:
                    delx = value
                    nbins = 0
                else:
                    nbins = int(value+0.001)
                    delx = 0.
            except:
                tkinter.messagebox.showinfo(
                    "Error",
                    "Unable to read the number of bins/bin size.  "
                    + "Check your inputs.")
                return
            if optionflag == 1:
                xmin = self.plot_range[self.current_plot-1][0]
                xmax = self.plot_range[self.current_plot-1][1]
            else:
                xmin = self.plot_range[self.current_plot-1][2]
                xmax = self.plot_range[self.current_plot-1][3]
            xp = None
            for loop in range(self.nsets):
                if (self.set_properties[loop]['display']) and \
                   (self.set_properties[loop]['plot'] == self.current_plot):
                    mycolour = self.set_properties[loop]['colour']
                    if optionflag == 1:
                        values = numpy.copy(self.xdata[loop]['values'])
                    else:
                        values = numpy.copy(self.ydata[loop]['values'])
                    if xp is None:
                        xp = [values, ]
                        histcolours = [mycolour, ]
                    else:
                        if setoptionflag == 1:
                            oldvalues = numpy.copy(xp[0])
                            newvalues = numpy.append(oldvalues, values)
                            xp[0] = newvalues
                        else:
                            xp.append(values)
                            histcolours.append(mycolour)
            self.histogramLabelText = Tk.StringVar()
            self.histogramLabel = Tk.Label(
                histogramwindow,
                textvariable=self.histogramLabelText, anchor=Tk.N, width=70)
            self.histogramLabel.pack()
            self.histogramLabelText.set("Value:")
            self.p2 = Figure(figsize=(6, 6), dpi=100)
            sp1 = self.p2.add_subplot(1, 1, 1)
            c1 = FigureCanvasTkAgg(self.p2, master=histogramwindow)
            c1.mpl_connect("motion_notify_event", self.histogram_position)
            try:
                if delx > 0.:
                    npixels = int(abs((xmax-xmin)/delx))
                    delx = abs(delx)
                else:
                    npixels = nbins
                    delx = (xmax-xmin)/npixels
            except:
                nbins = 500
            histx = []
            histy = []
            for loop in range(len(xp)):
                histogramy, hxedges = numpy.histogram(
                    xp[loop], npixels, range=[xmin, xmax])
                histogramx = (hxedges[1:]+hxedges[0:-1])/2.
                sp1.bar(histogramx, histogramy,
                        color=histcolours[loop],
                        width=delx*0.9)
                histx.append(histogramx)
                histy.append(histogramy)
            if optionflag == 1:
                sp1.set_xlabel(self.xparameters[self.current_plot-1]['label'],
                               family=self.fontname[self.current_plot-1],
                               size=self.fontsize[self.current_plot-1],
                               weight=self.fontweight[self.current_plot-1])
            else:
                sp1.set_xlabel(self.yparameters[self.current_plot-1]['label'],
                               family=self.fontname[self.current_plot-1],
                               size=self.fontsize[self.current_plot-1],
                               weight=self.fontweight[self.current_plot-1])
            sp1.set_ylabel('Number of points per bin')
            invertxflag = self.xparameters[self.current_plot-1]['invert']
            invertyflag = self.yparameters[self.current_plot-1]['invert']
            if (invertxflag == 1) and (optionflag == 1):
                sp1.invert_xaxis()
            if (invertyflag == 1) and (optionflag == 0):
                sp1.invert_xaxis()
            c1.draw()
            c1.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES)
            h1 = Tk.Frame(histogramwindow)
            h1.pack(side=Tk.TOP)
            h1.config(bg=BGCOL)
            button = Tk.Button(h1, text="Save as PS",
                               command=lambda: save_ps_figure(self.p2))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(h1, text="Save as PNG",
                               command=lambda: save_png_figure(self.p2))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Write out values",
                command=lambda: self.writeHistogram(histx, histy))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Close", command=histogramwindow.destroy)
            button.pack()
            button.config(bg=BGCOL)
        except:
            pass

    def writeHistogram(self, xvalues, yvalues):
        """
        Write out the histgram values to a selected output file.

        Parameters
        ----------
            xvalues :  a list of vectors of x values for the histogram, each a
                       numpy float array

            yvalues :  a list of vector of y values for the histogram, each a
                       numpy int array

        Returns
        -------
            No values are returned by the routine

        """
        outfilename = tkinter.filedialog.asksaveasfilename()
        if isinstance(outfilename, type('string')):
            outfile = open(outfilename, 'w')
            for n1 in range(len(xvalues)):
                for loop in range(len(xvalues[n1])):
                    print('%13.6g %10d' % (xvalues[n1][loop],
                                           yvalues[n1][loop]),
                          file=outfile)
                print(' ', file=outfile)
            outfile.close()

    def set_plot_number(self):
        """
        Select the active plot, if there is more than one.

        This routine calls up a window so that the user can select the
        current active plot.  The routine is not useful unless more than
        one plot is defined.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.setplot_window is not None:
            return
        self.setplot_window = Tk.Toplevel()
        self.setplot_window.title('Set Plot Number')
        holder = Tk.Frame(self.setplot_window)
        holder.pack(side=Tk.TOP)
        label = Tk.Label(holder, text="Active plot:")
        label.grid(column=0, row=0)
        self.currentplot_select = Tk.Entry(holder, width=10)
        self.currentplot_select.grid(column=1, row=1)
        self.currentplot_select.insert(0, str(self.current_plot))
        buttonframe = Tk.Frame(self.setplot_window)
        buttonframe.pack(side=Tk.TOP)
        apply_button = Tk.Button(buttonframe, text="Select",
                                 command=self.set_current_plot)
        apply_button.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            buttonframe, text="Close",
            command=lambda: self.close_window(self.setplot_window))
        close_button.pack(side=Tk.LEFT)

    def set_plot_hide(self):
        """
        Hide one of the plots, if there is more than one.

        This routine calls up a window so that the user can select the
        plot to hide.  The routine is not useful unless more than
        one plot is defined.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.hideplot_window is not None:
            return
        self.hideplot_window = Tk.Toplevel()
        self.hideplot_window.title('Hide Plot Number')
        holder = Tk.Frame(self.hideplot_window)
        holder.pack(side=Tk.TOP)
        label = Tk.Label(holder, text="Plot to Toggle Hide/Show:")
        label.grid(column=0, row=0)
        self.hideplot_select = Tk.Entry(holder, width=10)
        self.hideplot_select.grid(column=1, row=1)
        self.hideplot_select.insert(0, str(self.current_plot))
        buttonframe = Tk.Frame(self.hideplot_window)
        buttonframe.pack(side=Tk.TOP)
        apply_button = Tk.Button(buttonframe, text="Select",
                                 command=self.hide_selected_plot)
        apply_button.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            buttonframe, text="Close",
            command=lambda: self.close_window(self.hideplot_window))
        close_button.pack(side=Tk.LEFT)

    def tile_plots(self):
        """
        Bring up a window to allow multiple plots in the main display.

        This subroutine creates a window for use in arranging multiple plots
        in a single display (i.e. canvas).

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.tile_window is not None:
            return
        self.tile_window = Tk.Toplevel()
        self.tile_window.title('Set Number of Plots')
        holder = Tk.Frame(self.tile_window)
        holder.pack(side=Tk.TOP)
        label = Tk.Label(holder, text="Number of plots in x direction:")
        label.grid(column=0, row=0)
        label = Tk.Label(holder, text="Number of plots in y direction:")
        label.grid(column=0, row=1)
        label = Tk.Label(holder, text="Active plot:")
        label.grid(column=0, row=2)
        self.xplots_set = Tk.Entry(holder, width=10)
        self.xplots_set.grid(column=1, row=0)
        self.yplots_set = Tk.Entry(holder, width=10)
        self.yplots_set.grid(column=1, row=1)
        self.xplots_set.insert(0, '1')
        self.yplots_set.insert(0, '1')
        self.currentplot_set = Tk.Entry(holder, width=10)
        self.currentplot_set.grid(column=1, row=2)
        self.currentplot_set.insert(0, '1')
        buttonframe = Tk.Frame(self.tile_window)
        buttonframe.pack(side=Tk.TOP)
        apply_button = Tk.Button(buttonframe, text="Apply Values",
                                 command=self.set_plot_layout)
        apply_button.pack(side=Tk.LEFT)
        label1 = Tk.Label(buttonframe, text="     ")
        label1.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            buttonframe, text="Close",
            command=lambda: self.close_window(self.tile_window))
        close_button.pack(side=Tk.LEFT)

    def set_plot_layout(self):
        """
        Read the plot layout parmeters and apply them.

        This sub-routine reads the plot layout parameters and attempts
        to apply them.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            nx1 = int(self.xplots_set.get())
            ny1 = int(self.yplots_set.get())
            n1 = int(self.currentplot_set.get())
            if (nx1 < 1) | (ny1 < 1) | (nx1 > 5) | (ny1 > 5) | (n1 < 1) \
               | (n1 > nx1*ny1):
                raise ValueError
            self.make_plot_layout(nx1, ny1, n1)
        except:
            tkinter.messagebox.showinfo(
                "Error",
                "There was some issue with the plot layout parameters.")

    def make_plot_layout(self, nx1, ny1, n1):
        """
        Layout multiple subplots in the main window.

        This routine lays out the plots in an xn1 by xy1 grid (values are 1
        or larger for both parameters) and sets the current plot to
        number n1.

        Parameters
        ----------
            nx1 :  an integer value of 1 or higher, the number of plots in the
                   x direction

            ny1 :  an integer value of 1 or higher, the number of plots in the
                   y direction

            n1 :   an integer value between 1 and nx1*yn1, the plot to set as
                   the current plot

        Returns
        -------
            No values are returned by this routine.

        If the xn1 and ny1 values match the current layout (1, 1, and 1
        upon startup) then nothing is done here.
        """
        if (nx1 != self.nxplots) | (ny1 != self.nyplots):
            self.nyplots = ny1
            self.nxplots = nx1
            self.number_of_plots = self.nxplots * self.nyplots
            self.current_plot = n1
            self.make_plot_area(self.plotframe)
            for loop in range(self.number_of_plots):
                try:
                    value = self.xparameters[loop]['xmin']
                except:
                    self.xparameters.append({'label': ' ',
                                             'minimum': 0.0,
                                             'maximum': 1.0, 'major': 0.1,
                                             'minor': 0.1, 'logarithmic': 0,
                                             'invert': 0, 'hide': 0,
                                             'hideticks': 0,
                                             'hidelabels': 0,
                                             'hybridlog': 0,
                                             'inverseticks': 0,
                                             'ticklength': 6,
                                             'bothticks': 0, 'minorticks': 0,
                                             'oppositeaxis': 0})
                    self.yparameters.append({'label': ' ',
                                             'minimum': 0.0,
                                             'maximum': 1.0, 'major': 0.1,
                                             'minor': 0.1, 'logarithmic': 0,
                                             'invert': 0, 'hide': 0,
                                             'hideticks': 0,
                                             'hidelabels': 0,
                                             'hybridlog': 0,
                                             'inverseticks': 0,
                                             'ticklength': 6,
                                             'bothticks': 0, 'minorticks': 0,
                                             'oppositeaxis': 0})
                    self.fontname.append('sans-serif')
                    self.fontsize.append('12')
                    self.fontweight.append('normal')
                    self.plot_range.append([0.0, 1.0, 0.0, 1.0])
                    self.original_range.append(True)
                    self.plot_frame.append(0.0)
                    self.title.append(' ')
            for loop in range(self.number_of_plots):
                try:
                    value = self.legend_position[loop]
                except:
                    self.legend_position.append(None)
                    self.legend_labels.append(None)
                    self.legend_handles.append(None)
                    self.legend_options.append(None)
                    self.legend_frame.append(None)
                    self.legend_variable.append(None)
                    self.legend_user_position.append(None)
            for loop in range(self.number_of_plots):
                try:
                    values = self.equal_aspect[loop]
                except:
                    self.equal_aspect.append(False)
        self.make_plot()

    def set_current_plot(self):
        """
        Set which plot is the current one, if there are multiple plots.

        The routine here brings up a window wherein the user can select
        the current plot number in the case where several plots are
        displayed.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            n1 = int(self.currentplot_select.get())
            if (n1 < 1) | (n1 > self.number_of_plots):
                raise ValueError
            if n1 != self.current_plot:
                self.current_plot = n1
        except:
            tkinter.messagebox.showinfo(
                "Error", "There was an error in the requested plot number.")

    def hide_selected_plot(self):
        """
        Set a plot to be hidden, if there are multiple plots.

        The routine here brings up a window wherein the user can select
        the current plot number in the case where several plots are
        displayed.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            n1 = int(self.hideplot_select.get())
            if (n1 < 1) | (n1 > self.number_of_plots):
                raise ValueError
            self.hide_subplot[n1-1] = not self.hide_subplot[n1-1]
            plotnumber = self.current_plot
            self.current_plot = n1
            self.make_plot()
            self.current_plot = plotnumber
        except:
            tkinter.messagebox.showinfo(
                "Error", "There was an error in the requested plot number.")

    def clear_plot(self, query=True):
        """
        Clear the plot area.

        This routine clears the data sets and resets the various parameters to
        the initial values.

        Parameters
        ----------
            query :  An optional Boolean value for whether the user is queried
                     before the plot is cleared, which defaults to True

        Returns
        -------
            Nothing

        """
        if query:
            response = tkinter.messagebox.askyesno(
                "Verify", "Do you want to abandon the plot?")
        else:
            response = True
        if not response:
            return
        self.xdata = []
        self.ydata = []
        self.set_properties = []
        for loop in range(self.max_sets):
            self.set_properties.append({
                'symbol': None,
                'symbolsize': 4.0, 'linestyle': 'None',
                'linewidth': 1.0, 'colour': 'black',
                'label': '', 'xmin': 0.0, 'xmax': 1.0, 'ymin': 0.0,
                'ymax': 1.0, 'display': True, 'errors': False,
                'legend': True, 'plot': 1})
            self.xdata.append(None)
            self.ydata.append(None)
        self.plot_range = [[0., 1., 0., 1.], ]
        self.original_range = []
        for loop in range(self.max_sets):
            self.original_range.append(True)
        self.nsets = 0
        self.title = [' ', ]
        self.xparameters = [{'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
                             'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
                             'invert': 0, 'hide': 0, 'hideticks': 0,
                             'hidelabels': 0, 'hybridlog': 0,
                             'inverseticks': 0, 'ticklength': 6,
                             'bothticks': 0, 'minorticks': 0,
                             'oppositeaxis': 0}, ]
        self.yparameters = [{'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
                             'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
                             'invert': 0, 'hide': 0, 'hideticks': 0,
                             'hidelabels': 0, 'hybridlog': 0,
                             'inverseticks': 0, 'ticklength': 6,
                             'bothticks': 0, 'minorticks': 0,
                             'oppositeaxis': 0}, ]
        self.position_stack = []
        self.label_flag = False
        self.line_flag = False
        self.ellipse_flag = False
        self.box_flag = False
        self.vector_flag = False
        self.number_of_lines = 0
        self.number_of_labels = 0
        self.number_of_ellipses = 0
        self.number_of_boxes = 0
        self.number_of_vectors = 0
        self.positions = []
        for loop in range(len(self.legend_variable)):
            try:
                self.legend_variable[loop].set(0)
            except:
                self.legend_variable[loop] = None
        self.legend_handles = [None, ]
        self.legend_labels = [None, ]
        self.legend_frame = [None, ]
        self.legend_position = [None, ]
        self.legend_user_position = [None, ]
        self.plot_margin = 0.0
        self.plot_frame = [0.0, ]
        self.nxplots = 1
        self.nyplots = 1
        self.number_of_plots = self.nxplots * self.nyplots
        self.current_plot = 1
        if len(self.subplot) > 1:
            for loop in range(1, len(self.subplot)):
                self.subplot[loop].clear()
            self.subplot = [self.subplot[0], ]
            self.hide_subplot = [self.hide_subplot[0], ]
        self.make_plot()

    def autoscale_plot(self):
        """
        Poll the data sets to autoscale the current plot.

        This routine uses the data sets to autoscale the plot.  It replaces
        the current plot limit values with the new ones calculated in this
        routine.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.nsets == 0:
            self.plot_range[self.current_plot-1][0] = 0.
            self.plot_range[self.current_plot-1][1] = 1.
            self.plot_range[self.current_plot-1][2] = 0.
            self.plot_range[self.current_plot-1][3] = 1.
        else:
            xmin = 0.
            xmax = 0.
            for loop in range(self.nsets):
                if (self.set_properties[loop]['display']) and \
                   (self.set_properties[loop]['plot'] == self.current_plot):
                    if xmin == xmax:
                        xmin = self.xdata[loop]['minimum']
                        xmax = self.xdata[loop]['maximum']
                        ymin = self.ydata[loop]['minimum']
                        ymax = self.ydata[loop]['maximum']
                    else:
                        xmin = min(xmin, self.xdata[loop]['minimum'])
                        xmax = max(xmax, self.xdata[loop]['maximum'])
                        ymin = min(ymin, self.ydata[loop]['minimum'])
                        ymax = max(ymax, self.ydata[loop]['maximum'])
        if self.matplotlib_rounding.get():
            self.original_range[self.current_plot-1] = True
        else:
            xmin1 = self.round_float(xmin, True)
            xmax1 = self.round_float(xmax, False)
            ymin1 = self.round_float(ymin, True)
            ymax1 = self.round_float(ymax, False)
            self.plot_range[self.current_plot-1][0] = xmin1
            self.plot_range[self.current_plot-1][1] = xmax1
            self.plot_range[self.current_plot-1][2] = ymin1
            self.plot_range[self.current_plot-1][3] = ymax1
        self.make_plot()

    def make_plot_area(self, parent):
        """
        Set up the main figure area for the plot.

        This routine makes the figure area and the sub-plot, and then
        sets up some event call-backs.

        Parameters
        ----------
            parent :   A Tk.Frame variable that holds the plot

        Returns
        -------
            No values are returned by this routine.

        """
        # The size of the plot area is determined here (values in inches).
        # The x/y ratio is 10 to 7 as in xmgrace.  One could also use the
        # golden ratio 1.618 to 1.  Take whatever values seem to be the
        # most aesthetically pleasing.  The DPI value should probably not
        # be set below 100.
        if self.plot_area_flag:
            self.figure = Figure(figsize=(7.142857, 5), dpi=150)
        else:
            self.figure.clf()
        # The code allows a grid of self.nxplots by self.nyplots panels
        self.bounding_box = []
        self.subplot = []
        self.hide_subplot = []
        self.share_axis = []
        plot_number = 1
        for xloop in range(self.nxplots):
            for yloop in range(self.nyplots):
                self.subplot.append(
                    self.figure.add_subplot(self.nyplots,
                                            self.nxplots, plot_number))
                self.share_axis.append(0)
                self.hide_subplot.append(False)
                bbox = self.subplot[-1].get_position()
                bound_values = bbox.bounds
                self.bounding_box.append(bound_values)
                plot_number = plot_number + 1
        if self.plot_area_flag:
            # The following sets up a label above the plot, for the plot
            # position.
            self.position_label_text = Tk.StringVar()
            self.position_label = Tk.Label(
                parent, textvariable=self.position_label_text)
            self.position_label.pack(side=Tk.TOP)
            self.position_label_text.set("Position:\n")
            self.canvas = FigureCanvasTkAgg(self.figure, master=parent)
            # Here are defined the events that the program responds to for
            # the plot area.
            self.canvas.mpl_connect("motion_notify_event", self.plot_position)
            self.canvas.mpl_connect("button_press_event", self.plot_marker_set)
            self.canvas.mpl_connect("button_release_event",
                                    self.plot_marker_release)
            self.canvas.mpl_connect("key_press_event", self.key_press_event)
            self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH,
                                             expand=Tk.YES)
            self.plot_area_flag = False
        self.canvas.draw()

    def key_press_event(self, event):
        """
        Change the plot focus if a key is pressed in a plot area.

        The "inaxes" parameter is compared to the available subplots and if
        they match then that subplot is set to be the current plot.

        Parameters
        ----------

        event:   A matplotlib KeyPress event

        Returns
        -------

        No values are returned.
        """
        for loop in range(self.number_of_plots):
            if event.inaxes == self.subplot[loop]:
                self.current_plot = loop + 1
        return

    def round_float(self, value, minimum_flag):
        """
        Round a floating point value to the nearest significant figure.

        If the self.matplotlib_rounding flag is False this performs rounding of
        floating point values to convenient values for such things as
        calculating plot ranges.  In some cases the matplotlib rounding is
        not very good so this is provided as a way to get better results.
        When one is using a logarithmic scaling the matplotlib ropunding is
        not always good, for example, and this routine can do better.  In
        other cases this routine does not work as well, such as when the
        total range of values on the plot is small.

        Parmeters
        ---------
            value :         A floating point number to be rounded off

            minimum_flag :  A boolean value value to determine whether the
                            value should be rounded down.  If True values
                            are rounded down as one wishes for the minimum
                            value in a plot, and if False values are rounded
                            up as one wishes for the maximum value in a plot.

        Returns
        -------
           rounded_value :    The rounded off floating point number calculated
                               from `value', or the input value in the case
                               where the self.matplotlib_routine flag is True.

        The code differs from use of the floor and ceil functions in that
        it tries to round to the nearest significant figure for a given
        value.  So passing in a value of 1.2e-08, for exmaple, with
        minimum_flag = False returns 2.0e-08, and if the flag is True it
        returns 1.0e-08.

        The code uses the math log10, floor, and ceil functions.

        """
        if self.matplotlib_rounding.get():
            return value
        if value == 0.:
            return value
        if value < 0.:
            sign = -1.
            value1 = -value
        else:
            sign = 1.
            value1 = value
        power = math.log10(value1)
        if power < 0.:
            exp = int(power-1.)
        else:
            exp = int(power)
        shift = 10.**exp
        x = value1/shift
        delx = 1.0
        if x < 1.7:
            x = x*10.
            shift = shift/10.
        elif x < 2.5:
            x = x*5.
            shift = shift/5.
        if (minimum_flag) and sign > 0.:
            x = math.floor(x)
        elif (minimum_flag) and sign < 0.:
            x = math.ceil(x)
        elif (not minimum_flag) and sign > 0.:
            x = math.ceil(x)
        elif (not minimum_flag) and sign < 0.:
            x = math.floor(x)
        rounded_value = x*shift*sign
        # If the rounded value is very close to the input value, offset
        # by one unit in x...not normally part of the routine, but needed
        # for matplotlib plots because of symbols close to the edges of
        # the plot.
        ratio = abs(value/rounded_value)
        if (ratio > 0.97) and (ratio < 1.03):
            if (minimum_flag) and sign > 0.:
                rounded_value = (x-delx)*shift*sign
            elif (minimum_flag) and sign < 0.:
                rounded_value = (x+delx)*shift*sign
            elif (not minimum_flag) and sign > 0.:
                rounded_value = (x+delx)*shift*sign
            elif (not minimum_flag) and sign < 0.:
                rounded_value = (x-delx)*shift*sign
        return rounded_value

    def plot_position(self, event):
        """
        Take a motion notify event and update the plot position field.

        This is an event handler routine that receives "motion notify events"
        within main the plot area.  It takes the position information and
        updates the position string at the top of the main plot area.

        Parameters
        ----------
            event :   A Tkinter event variable.  Only the position values
                      are used.

        Returns
        -------
            No values are returned by this routine.

        """
        if (event.xdata is None) or (event.ydata is None):
            return
        try:
            hlogxflag = self.xparameters[self.current_plot-1]['hybridlog']
            hlogyflag = self.yparameters[self.current_plot-1]['hybridlog']
            if hlogxflag == 0:
                xout = event.xdata
            else:
                xout = inverse_hybrid_transform(event.xdata)
            if hlogyflag == 0:
                yout = event.ydata
            else:
                yout = inverse_hybrid_transform(event.ydata)
            nset, npoint, xpoint, ypoint, distance = self.match_point(
                xout, yout)
            s1 = "Position: %.6g %.6g" % (xout, yout)
            if nset is not None:
                s1 = s1 + '\nNearest point: set %d point %d \n' % \
                     (nset+1, npoint) \
                     + 'data (%.6g, %.6g) distance %.6g' % (
                         xpoint, ypoint, distance)
            self.position_label_text.set(s1)
        except:
            pass

    def match_point(self, xposition, yposition):
        """
        Match a data point closest to the current position.

        Given an input plot position this routine finds the closest plot
        point and returns the point values.

        Parameters
        ----------
            xposition :  a floating point value, the x position in the
                         current plot

            yposition :  a floating point value, the y position in the
                         current plot

        Returns
        -------
            nset :    an integer value, the set number of the closest point

            ndata :   an integer value, the point number within the set

            xmin :    a float value, the x coordinate of the closest point

            ymin :    a float value, the y coordinate of the closest point

            dmin :    a float value, the distance of the closest point to the
                      input position

        The return values may all be None if the position is outside the
        plot area
        """
        if (xposition is None) or (yposition is None):
            return None, None, None, None, None
        if self.nsets == 0:
            return None, None, None, None, None
        try:
            dmin = -1.
            for loop in range(self.nsets):
                xdata = numpy.copy(self.xdata[loop]['values'])
                ydata = numpy.copy(self.ydata[loop]['values'])
                if (xdata is None) or (ydata is None):
                    break
                distances = numpy.sqrt(
                    (xdata - xposition)*(xdata - xposition)
                    + (ydata - yposition)*(ydata - yposition)
                )
                dargmin = numpy.argmin(distances)
                dminset = numpy.copy(distances[dargmin])
                try:
                    dminhere = dminset[0]
                except:
                    dminhere = dminset
                if (dmin < 0.) or (dminhere < dmin):
                    nset = loop
                    try:
                        ndata = dargmin[0]
                    except:
                        ndata = dargmin
                    dmin = dminhere
                    xmin = xdata[ndata]
                    ymin = ydata[ndata]
            if dmin >= 0.:
                return nset, ndata, xmin, ymin, dmin
        except:
            return None, None, None, None, None

    def plot_marker_set(self, event):
        """
        Respond to button press events.

        This is an event handler routine that responds to "button
        press events".  It normally notes the position for the drawing
        functions.  If the label_flag is set, this routine then gets
        the label information instead.

        Parameters
        ----------
            event :  A Tkinter event variable.  Only the position values
                     are used.

        Returns
        -------
            No values are returned by this routine.

        """
        xvalue = event.xdata
        yvalue = event.ydata
        if self.xparameters[self.current_plot-1]['hybridlog'] == 1:
            xvalue = inverse_hybrid_transform(xvalue)
        if self.yparameters[self.current_plot-1]['hybridlog'] == 1:
            yvalue = inverse_hybrid_transform(yvalue)
        position = [xvalue, yvalue]
        if self.label_flag:
            try:
                self.plot_labels[self.number_of_labels]['xposition'] = \
                    position[0]
                self.plot_labels[self.number_of_labels]['yposition'] = \
                    position[1]
                self.plot_labels[self.number_of_labels]['plot'] = \
                    self.current_plot
                self.plot_labels[self.number_of_labels]['labelstring'] = ''
                self.label_flag = False
                self.set_label_properties(self.number_of_labels)
                self.number_of_labels = self.number_of_labels + 1
                self.make_plot()
            except:
                tkinter.messagebox.showinfo(
                    "Label Error", "The label was not processed properly.")
            return
        if self.line_flag:
            self.positions.append(position)
            return
        if self.vector_flag:
            self.positions.append(position)
            return
        if self.box_flag:
            self.positions.append(position)
            return
        if self.ellipse_flag:
            self.positions.append(position)
            return

    def plot_marker_release(self, event):
        """
        Respond to button release events.

        This is an event handler routine that responds to "button
        release events".  It notes the position for line/box/vector/ellipse
        functions, and resets the flags for these functions.

        Parameters
        ----------
            event :  A Tkinter event variable.  Only the position values
                     are used.

        Returns
        -------
            No values are returned by this routine.

        """
        xvalue = event.xdata
        yvalue = event.ydata
        if self.xparameters[self.current_plot-1]['hybridlog'] == 1:
            xvalue = inverse_hybrid_transform(xvalue)
        if self.yparameters[self.current_plot-1]['hybridlog'] == 1:
            yvalue = inverse_hybrid_transform(yvalue)
        position = [xvalue, yvalue]
        if self.line_flag:
            self.positions.append(position)
            self.add_line_values()
        if self.vector_flag:
            self.positions.append(position)
            self.add_vector_values()
        if self.ellipse_flag:
            self.positions.append(position)
            self.add_ellipse_values()
        if self.box_flag:
            self.positions.append(position)
            self.add_box_values()
        self.make_plot()

    def add_box_values(self):
        """
        Create a box on the plot.

        This code is activated when the box definition option is selected.
        When a button press event and then a button release event are
        received then the positions are recorded in self.positions.  This
        routine reads these positions and presents a window with the box
        parameters for the user to change as they wish.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            self.positions[-2][0]
            self.positions[-2][1]
            self.positions[-1][0]
            self.positions[-1][1]
        except:
            tkinter.messagebox.showinfo(
                "Error",
                "The required start and stop positions are not available"
                + " to make a box.")
            return
        self.box_flag = False
        self.box_window = Tk.Toplevel()
        self.box_window.title('Set Box Properties')
        self.box_window.config(bg=BGCOL)
        frame1 = Tk.Frame(self.box_window)
        frame1.pack(side=Tk.TOP)
        label = Tk.Label(frame1, text='Corner 1 x')
        label.grid(column=0, row=0)
        label = Tk.Label(frame1, text='Corner 1 y')
        label.grid(column=0, row=1)
        label = Tk.Label(frame1, text='Corner 2 x')
        label.grid(column=0, row=2)
        label = Tk.Label(frame1, text='Corner 2 y')
        label.grid(column=0, row=3)
        label = Tk.Label(frame1, text='Orientation (degrees)')
        label.grid(column=0, row=4)
        label = Tk.Label(frame1, text='Line type')
        label.grid(column=0, row=5)
        label = Tk.Label(frame1, text='Line colour')
        label.grid(column=0, row=6)
        label = Tk.Label(frame1, text='Line thickness')
        label.grid(column=0, row=7)
        label = Tk.Label(frame1, text='Fill color')
        label.grid(column=0, row=8)
        # boxfields holds the box parameter entry/menu items
        # 0 to 3    positions
        # 4 orientation angle (degrees)
        # 5 line type (solid, dashed, etc)
        # 6 line colour
        # 7 line thickness
        # 8 interior colour (includes "none" for no colour, the default.
        self.boxfields = []
        self.boxfields.append(Tk.Entry(frame1, width=20))
        self.boxfields[-1].grid(column=1, row=0, sticky=Tk.W)
        self.boxfields.append(Tk.Entry(frame1, width=20))
        self.boxfields[-1].grid(column=1, row=1, sticky=Tk.W)
        self.boxfields.append(Tk.Entry(frame1, width=20))
        self.boxfields[-1].grid(column=1, row=2, sticky=Tk.W)
        self.boxfields.append(Tk.Entry(frame1, width=20))
        self.boxfields[-1].grid(column=1, row=3, sticky=Tk.W)
        self.boxfields.append(Tk.Entry(frame1, width=20))
        self.boxfields[-1].grid(column=1, row=4, sticky=Tk.W)
        self.boxfields[0].insert(0, str(self.positions[-2][0]))
        self.boxfields[1].insert(0, str(self.positions[-2][1]))
        self.boxfields[2].insert(0, str(self.positions[-1][0]))
        self.boxfields[3].insert(0, str(self.positions[-1][1]))
        self.boxfields[4].insert(0, '0.0')
        self.boxfields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.boxfields[-1].grid(column=1, row=5, sticky=Tk.W)
        self.boxfields[-1]['values'] = matplotlib_line_name_list
        self.boxfields[-1].current(0)
        self.boxfields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.boxfields[-1].grid(column=1, row=6, sticky=Tk.W)
        self.boxfields[-1]['values'] = self.colourset
        self.boxfields[-1].current(0)
        self.boxfields.append(Tk.Entry(frame1, width=15))
        self.boxfields[-1].grid(column=1, row=7, sticky=Tk.W)
        self.boxfields[-1].insert(0, '1.0')
        self.boxfields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.boxfields[-1].grid(column=1, row=8, sticky=Tk.W)
        self.boxfields[-1]['values'] = self.altcolourset
        self.boxfields[-1].current(0)
        frame2 = Tk.Frame(self.box_window)
        frame2.pack(side=Tk.TOP)
        apply_button = Tk.Button(frame2, text="Apply",
                                 command=self.apply_box_values)
        apply_button.pack(side=Tk.LEFT)
        label1 = Tk.Label(frame2, text="    ")
        label1.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            frame2, text="Close Window",
            command=lambda: self.close_window(self.box_window))
        close_button.pack(side=Tk.LEFT)

    def apply_line_values(self):
        """
        Create a line on the plot.

        This code reads the values in the line properties defintion window
        and applies them to the next available line.  The plot is then
        redone and the line properties window is removed.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            x1 = float(self.linefields[0].get())
            y1 = float(self.linefields[1].get())
            x2 = float(self.linefields[2].get())
            y2 = float(self.linefields[3].get())
            t1 = float(self.linefields[6].get())
            self.plot_lines[self.number_of_lines]['xstart'] = x1
            self.plot_lines[self.number_of_lines]['ystart'] = y1
            self.plot_lines[self.number_of_lines]['xend'] = x2
            self.plot_lines[self.number_of_lines]['yend'] = y2
            line_index = self.linefields[4].current()
            colour_index = self.linefields[5].current()
            self.plot_lines[self.number_of_lines]['line_thickness'] = t1
            if self.colourset[colour_index] == 'select':
                values = askcolor()
                self.plot_lines[self.number_of_lines]['line_colour'] = \
                    values[1]
            else:
                self.plot_lines[self.number_of_lines]['line_colour'] = \
                                self.colourset[colour_index]
            self.plot_lines[self.number_of_lines]['line_type'] = \
                matplotlib_line_name_list[line_index]
            self.plot_lines[self.number_of_lines]['plot'] = self.current_plot
            self.number_of_lines = self.number_of_lines + 1
            self.make_plot()
            self.close_data_window(self.line_window)
        except:
            return

    def apply_vector_values(self):
        """
        Create a vector for the plot.

        This code reads the values in the vector properties defintion window
        and applies them to the next available vector.  The plot is then
        redone and the vector properties window is removed.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            x1 = float(self.vectorfields[0].get())
            y1 = float(self.vectorfields[1].get())
            x2 = float(self.vectorfields[2].get())
            y2 = float(self.vectorfields[3].get())
            delx = float(self.vectorfields[4].get())
            dely = float(self.vectorfields[5].get())
            t1 = float(self.vectorfields[8].get())
            self.plot_vectors[self.number_of_vectors]['xstart'] = x1
            self.plot_vectors[self.number_of_vectors]['ystart'] = y1
            self.plot_vectors[self.number_of_vectors]['xend'] = x2
            self.plot_vectors[self.number_of_vectors]['yend'] = y2
            self.plot_vectors[self.number_of_vectors]['delx'] = delx
            self.plot_vectors[self.number_of_vectors]['dely'] = dely
            self.plot_vectors[self.number_of_vectors]['line_thickness'] = t1
            line_index = self.vectorfields[6].current()
            colour_index = self.vectorfields[7].current()
            if self.colourset[colour_index] == 'select':
                values = askcolor()
                self.plot_vectors[self.number_of_vectors]['line_colour'] = \
                    values[1]
            else:
                self.plot_vectors[self.number_of_vectors]['line_colour'] = \
                    self.colourset[colour_index]
            self.plot_vectors[self.number_of_vectors]['line_type'] = \
                matplotlib_line_name_list[line_index]
            self.plot_vectors[self.number_of_vectors]['plot'] = \
                self.current_plot
            flag = self.vectorfields[9].get()
            if flag == 1:
                self.plot_vectors[self.number_of_vectors]['fill'] = True
            else:
                self.plot_vectors[self.number_of_vectors]['fill'] = False
            colour_index = self.vectorfields[10].current()
            if self.colourset[colour_index] == 'select':
                values = askcolor()
                self.plot_vectors[self.number_of_vectors]['fill_colour'] = \
                    values[1]
            else:
                self.plot_vectors[self.number_of_vectors]['fill_colour'] = \
                    self.colourset[colour_index]
            self.number_of_vectors = self.number_of_vectors + 1
            self.make_plot()
            self.close_data_window(self.vector_window)
        except:
            return

    def apply_box_values(self):
        """
        Create a box for the plot.

        This code reads the values in the box properties defintion window and
        applies them to the next available box.  The plot is then redone and
        the box properties window is removed.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            x1 = float(self.boxfields[0].get())
            y1 = float(self.boxfields[1].get())
            x2 = float(self.boxfields[2].get())
            y2 = float(self.boxfields[3].get())
            angle = float(self.boxfields[4].get())
            angle = angle % 360.
            t1 = float(self.boxfields[7].get())
            self.plot_boxes[self.number_of_boxes]['xstart'] = x1
            self.plot_boxes[self.number_of_boxes]['ystart'] = y1
            self.plot_boxes[self.number_of_boxes]['xend'] = x2
            self.plot_boxes[self.number_of_boxes]['yend'] = y2
            self.plot_boxes[self.number_of_boxes]['rotation'] = angle
            line_index = self.boxfields[5].current()
            colour_index1 = self.boxfields[6].current()
            colour_index2 = self.boxfields[8].current()
            self.plot_boxes[self.number_of_boxes]['line_thickness'] = t1
            if self.colourset[colour_index1] == 'select':
                values = askcolor()
                self.plot_boxes[self.number_of_boxes]['line_colour'] = \
                    values[1]
            else:
                self.plot_boxes[self.number_of_boxes]['line_colour'] = \
                    self.colourset[colour_index1]
            self.plot_boxes[self.number_of_boxes]['line_type'] = \
                matplotlib_line_name_list[line_index]
            self.plot_boxes[self.number_of_boxes]['plot'] = self.current_plot
            self.plot_boxes[self.number_of_boxes]['fill_colour'] = \
                self.altcolourset[colour_index2]
            self.number_of_boxes = self.number_of_boxes + 1
            self.make_plot()
            self.close_window(self.box_window)
        except:
            return

    def add_ellipse_values(self):
        """
        Create an ellipse for the plot.

        This code is activated when the ellipse definition option is selected.
        When a button press event and then a button release event are received
        then the positions are recorded in self.positions.  This routine reads
        these positions and presents a window with the ellipse parameters for
        the user to change as they wish.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            self.positions[-2][0]
            self.positions[-2][1]
            self.positions[-1][0]
            self.positions[-1][1]
        except ValueError:
            tkinter.messagebox.showinfo(
                "Error", "The required start and"
                + " stop positions are not available to make an ellipse.")
            return
        self.ellipse_flag = False
        self.ellipse_window = Tk.Toplevel()
        self.ellipse_window.title('Set Ellipse Properties')
        self.ellipse_window.config(bg=BGCOL)
        frame1 = Tk.Frame(self.ellipse_window)
        frame1.pack(side=Tk.TOP)
        label = Tk.Label(frame1, text='Center x')
        label.grid(column=0, row=0)
        label = Tk.Label(frame1, text='Center y')
        label.grid(column=0, row=1)
        label = Tk.Label(frame1, text='Major Axis x')
        label.grid(column=0, row=2)
        label = Tk.Label(frame1, text='Minor Axis y')
        label.grid(column=0, row=3)
        label = Tk.Label(frame1, text='Orientation (degrees)')
        label.grid(column=0, row=4)
        label = Tk.Label(frame1, text='Line type')
        label.grid(column=0, row=5)
        label = Tk.Label(frame1, text='Line colour')
        label.grid(column=0, row=6)
        label = Tk.Label(frame1, text='Line thickness')
        label.grid(column=0, row=7)
        label = Tk.Label(frame1, text='Fill color')
        label.grid(column=0, row=8)
        # ellipsefields holds the ellipse parameter entry/menu items
        # 0 to 3    center and widths
        # 4 orientation angle (degrees)
        # 5 line type (solid, dashed, etc)
        # 6 line colour
        # 7 line thickness
        # 8 interior colour (includes "none" for no colour, the default.
        self.ellipsefields = []
        self.ellipsefields.append(Tk.Entry(frame1, width=20))
        self.ellipsefields[-1].grid(column=1, row=0, sticky=Tk.W)
        self.ellipsefields.append(Tk.Entry(frame1, width=20))
        self.ellipsefields[-1].grid(column=1, row=1, sticky=Tk.W)
        self.ellipsefields.append(Tk.Entry(frame1, width=20))
        self.ellipsefields[-1].grid(column=1, row=2, sticky=Tk.W)
        self.ellipsefields.append(Tk.Entry(frame1, width=20))
        self.ellipsefields[-1].grid(column=1, row=3, sticky=Tk.W)
        self.ellipsefields.append(Tk.Entry(frame1, width=20))
        self.ellipsefields[-1].grid(column=1, row=4, sticky=Tk.W)
        xcenter = (self.positions[-2][0]+self.positions[-1][0])/2.
        ycenter = (self.positions[-2][1]+self.positions[-1][1])/2.
        xwidth = abs(self.positions[-2][0]-self.positions[-1][0])
        ywidth = abs(self.positions[-2][1]-self.positions[-1][1])
        self.ellipsefields[0].insert(0, str(xcenter))
        self.ellipsefields[1].insert(0, str(ycenter))
        self.ellipsefields[2].insert(0, str(xwidth))
        self.ellipsefields[3].insert(0, str(ywidth))
        self.ellipsefields[4].insert(0, '0.0')
        self.ellipsefields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.ellipsefields[-1].grid(column=1, row=5, sticky=Tk.W)
        self.ellipsefields[-1]['values'] = matplotlib_line_name_list
        self.ellipsefields[-1].current(0)
        self.ellipsefields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.ellipsefields[-1].grid(column=1, row=6, sticky=Tk.W)
        self.ellipsefields[-1]['values'] = self.colourset
        self.ellipsefields[-1].current(0)
        self.ellipsefields.append(Tk.Entry(frame1, width=15))
        self.ellipsefields[-1].grid(column=1, row=7, sticky=Tk.W)
        self.ellipsefields[-1].insert(0, '1.0')
        self.ellipsefields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.ellipsefields[-1].grid(column=1, row=8, sticky=Tk.W)
        self.ellipsefields[-1]['values'] = self.altcolourset
        self.ellipsefields[-1].current(0)
        frame2 = Tk.Frame(self.ellipse_window)
        frame2.pack(side=Tk.TOP)
        apply_button = Tk.Button(frame2, text="Apply",
                                 command=self.apply_ellipse_values)
        apply_button.pack(side=Tk.LEFT)
        label1 = Tk.Label(frame2, text="    ")
        label1.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            frame2, text="Close Window",
            command=lambda: self.close_window(self.ellipse_window))
        close_button.pack(side=Tk.LEFT)

    def apply_ellipse_values(self):
        """
        Create an ellipse for the plot.

        This code reads the values in the ellipse properties defintion window
        and applies them to the next available ellipse.  The plot is then
        redone and the ellipse properties window is removed.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            x1 = float(self.ellipsefields[0].get())
            y1 = float(self.ellipsefields[1].get())
            x2 = float(self.ellipsefields[2].get())
            y2 = float(self.ellipsefields[3].get())
            angle = float(self.ellipsefields[4].get())
            angle = angle % 360.
            t1 = float(self.ellipsefields[7].get())
            self.plot_ellipses[self.number_of_ellipses]['xposition'] = x1
            self.plot_ellipses[self.number_of_ellipses]['yposition'] = y1
            self.plot_ellipses[self.number_of_ellipses]['major'] = x2
            self.plot_ellipses[self.number_of_ellipses]['minor'] = y2
            self.plot_ellipses[self.number_of_ellipses]['rotation'] = angle
            line_index = self.ellipsefields[5].current()
            colour_index1 = self.ellipsefields[6].current()
            colour_index2 = self.ellipsefields[8].current()
            self.plot_ellipses[self.number_of_ellipses]['line_thickness'] = t1
            if self.colourset[colour_index1] == 'select':
                values = askcolor()
                self.plot_ellipses[self.number_of_ellipses]['line_colour'] = \
                    values[1]
            else:
                self.plot_ellipses[self.number_of_ellipses]['line_colour'] = \
                    self.colourset[colour_index1]
            self.plot_ellipses[self.number_of_ellipses]['line_type'] = \
                matplotlib_line_name_list[line_index]
            self.plot_ellipses[self.number_of_ellipses]['plot'] = \
                self.current_plot
            self.plot_ellipses[self.number_of_ellipses]['fill_colour'] = \
                self.altcolourset[colour_index2]
            self.number_of_ellipses = self.number_of_ellipses + 1
            self.make_plot()
            self.close_window(self.ellipse_window)
        except:
            return

    def add_line_values(self):
        """
        Create a line for the plot.

        This code is activated when the line definition option is selected.
        When a button press event and then a button release event are
        received then the positions are recorded in self.positions.  This
        routine reads these positions and presents a window with the line
        parameters for the user to change as they wish.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            self.positions[-2][0]
            self.positions[-2][1]
            self.positions[-1][0]
            self.positions[-1][1]
        except ValueError:
            tkinter.messagebox.showinfo(
                "Error", "The required start and "
                + "stop positions are not available to make a line.")
            return
        self.line_flag = False
        self.line_window = Tk.Toplevel()
        self.line_window.title('Set Line Properties')
        self.line_window.config(bg=BGCOL)
        frame1 = Tk.Frame(self.line_window)
        frame1.pack(side=Tk.TOP)
        label = Tk.Label(frame1, text='Start x')
        label.grid(column=0, row=0)
        label = Tk.Label(frame1, text='Start y')
        label.grid(column=0, row=1)
        label = Tk.Label(frame1, text='End x')
        label.grid(column=0, row=2)
        label = Tk.Label(frame1, text='End y')
        label.grid(column=0, row=3)
        label = Tk.Label(frame1, text='Line type')
        label.grid(column=0, row=4)
        label = Tk.Label(frame1, text='Line colour')
        label.grid(column=0, row=5)
        label = Tk.Label(frame1, text='Line thickness')
        label.grid(column=0, row=6)
        # linefields holds the line parameter entry/menu items
        # 0 to 3    positions
        # 4 line type (solid, dashed, etc)
        # 5 line colour
        # 6 line thickness
        self.linefields = []
        self.linefields.append(Tk.Entry(frame1, width=20))
        self.linefields[-1].grid(column=1, row=0, sticky=Tk.W)
        self.linefields.append(Tk.Entry(frame1, width=20))
        self.linefields[-1].grid(column=1, row=1, sticky=Tk.W)
        self.linefields.append(Tk.Entry(frame1, width=20))
        self.linefields[-1].grid(column=1, row=2, sticky=Tk.W)
        self.linefields.append(Tk.Entry(frame1, width=20))
        self.linefields[-1].grid(column=1, row=3, sticky=Tk.W)
        self.linefields[0].insert(0, str(self.positions[-2][0]))
        self.linefields[1].insert(0, str(self.positions[-2][1]))
        self.linefields[2].insert(0, str(self.positions[-1][0]))
        self.linefields[3].insert(0, str(self.positions[-1][1]))
        self.linefields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.linefields[-1].grid(column=1, row=4, sticky=Tk.W)
        self.linefields[-1]['values'] = matplotlib_line_name_list[0:-1]
        self.linefields[-1].current(0)
        self.linefields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.linefields[-1].grid(column=1, row=5, sticky=Tk.W)
        self.linefields[-1]['values'] = self.colourset
        self.linefields[-1].current(0)
        self.linefields.append(Tk.Entry(frame1, width=15))
        self.linefields[-1].grid(column=1, row=6, sticky=Tk.W)
        self.linefields[-1].insert(0, '1.0')
        frame2 = Tk.Frame(self.line_window)
        frame2.pack(side=Tk.TOP)
        apply_button = Tk.Button(frame2, text="Apply",
                                 command=self.apply_line_values)
        apply_button.pack(side=Tk.LEFT)
        label1 = Tk.Label(frame2, text="    ")
        label1.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            frame2, text="Close",
            command=lambda: self.close_window(self.line_window))
        close_button.pack(side=Tk.LEFT)

    def add_vector_values(self):
        """
        Create a vector for the plot.

        This code is activated when the vector definition option is selected.
        When a button press event and then a button release event are received
        then the positions are recorded in self.positions.  This routine reads
        these positions and presents a window with the vector parameters for
        the user to change as they wish.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            self.positions[-2][0]
            self.positions[-2][1]
            self.positions[-1][0]
            self.positions[-1][1]
        except ValueError:
            tkinter.messagebox.showinfo(
                "Error", "The required start "
                + "and stop positions are not available to make a vector.")
            return
        self.vector_flag = False
        self.vector_window = Tk.Toplevel()
        self.vector_window.title('Set Vector Properties')
        self.vector_window.config(bg=BGCOL)
        frame1 = Tk.Frame(self.vector_window)
        frame1.pack(side=Tk.TOP)
        label = Tk.Label(frame1, text='Start x')
        label.grid(column=0, row=0)
        label = Tk.Label(frame1, text='Start y')
        label.grid(column=0, row=1)
        label = Tk.Label(frame1, text='End x')
        label.grid(column=0, row=2)
        label = Tk.Label(frame1, text='End y')
        label.grid(column=0, row=3)
        label = Tk.Label(frame1, text='Vector Head del x')
        label.grid(column=0, row=4)
        label = Tk.Label(frame1, text='Vector Head del y')
        label.grid(column=0, row=5)
        label = Tk.Label(frame1, text='Line type')
        label.grid(column=0, row=6)
        label = Tk.Label(frame1, text='Line colour')
        label.grid(column=0, row=7)
        label = Tk.Label(frame1, text='Line thickness')
        label.grid(column=0, row=8)
        label = Tk.Label(frame1, text='Fill head')
        label.grid(column=0, row=9)
        label = Tk.Label(frame1, text='Head colour')
        label.grid(column=0, row=10)
        # vectorfields holds the vector parameter entry/menu items
        # 0 to 3    positions
        # 4 and 5 vector head x and y sizes
        # 6 line type (solid, dashed, etc)
        # 7 line colour
        # 8 line thickness
        # 9 head fill (radio button)
        # 10 head fill colour
        self.vectorfields = []
        self.vectorfields.append(Tk.Entry(frame1, width=20))
        self.vectorfields[-1].grid(column=1, row=0, sticky=Tk.W)
        self.vectorfields.append(Tk.Entry(frame1, width=20))
        self.vectorfields[-1].grid(column=1, row=1, sticky=Tk.W)
        self.vectorfields.append(Tk.Entry(frame1, width=20))
        self.vectorfields[-1].grid(column=1, row=2, sticky=Tk.W)
        self.vectorfields.append(Tk.Entry(frame1, width=20))
        self.vectorfields[-1].grid(column=1, row=3, sticky=Tk.W)
        self.vectorfields.append(Tk.Entry(frame1, width=20))
        self.vectorfields[-1].grid(column=1, row=4, sticky=Tk.W)
        self.vectorfields.append(Tk.Entry(frame1, width=20))
        self.vectorfields[-1].grid(column=1, row=5, sticky=Tk.W)
        self.vectorfields[0].insert(0, str(self.positions[-2][0]))
        self.vectorfields[1].insert(0, str(self.positions[-2][1]))
        self.vectorfields[2].insert(0, str(self.positions[-1][0]))
        self.vectorfields[3].insert(0, str(self.positions[-1][1]))
        self.vectorfields[4].insert(0, str(0.1))
        self.vectorfields[5].insert(0, str(0.1))
        self.vectorfields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.vectorfields[-1].grid(column=1, row=6, sticky=Tk.W)
        self.vectorfields[-1]['values'] = matplotlib_line_name_list[0:-1]
        self.vectorfields[-1].current(0)
        self.vectorfields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.vectorfields[-1].grid(column=1, row=7, sticky=Tk.W)
        self.vectorfields[-1]['values'] = self.colourset
        self.vectorfields[-1].current(0)
        self.vectorfields.append(Tk.Entry(frame1, width=15))
        self.vectorfields[-1].grid(column=1, row=8, sticky=Tk.W)
        self.vectorfields[-1].insert(0, '1.0')
        self.vector_head_fill_flag = Tk.IntVar()
        self.vectorfields.append(self.vector_head_fill_flag)
        self.vectorfields.append(tkinter.ttk.Combobox(frame1, width=15))
        b1 = Tk.Frame(frame1, )
        b1.grid(column=1, row=9, sticky=Tk.W)
        self.put_yes_no(b1, self.vector_head_fill_flag, ['yes', 'no'], True)
        self.vectorfields[-1].grid(column=1, row=10, sticky=Tk.W)
        self.vectorfields[-1]['values'] = self.colourset
        self.vectorfields[-1].current(0)
        frame2 = Tk.Frame(self.vector_window)
        frame2.pack(side=Tk.TOP)
        apply_button = Tk.Button(frame2, text="Apply",
                                 command=self.apply_vector_values)
        apply_button.pack(side=Tk.LEFT)
        label1 = Tk.Label(frame2, text="    ")
        label1.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            frame2, text="Close",
            command=lambda: self.close_window(self.vector_window))
        close_button.pack(side=Tk.LEFT)

    def set_label_properties(self, ind1):
        """
        Make a window in which the label properties are assigned.

        Parameters
        ----------
            ind1 :  an integer value >= 0, the index of the label that is
                    having the properties set

        Returns
        -------
            Nothing

        The routine gets input from the window each time the "set properties"
        button is activated, and applies these to label ind1.
        """
        label_property_window = Tk.Toplevel()
        label_property_window.title('Label properties')
        holder = Tk.Frame(label_property_window)
        holder.pack(side=Tk.TOP)
        fontnames = ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace',
                     'times new roman']
        fontsizes = ['8', '9', '10', '11', '12', '13', '14', '16', '18', '20',
                     '24', '30']
        fontweights = ['ultralight', 'light', 'normal', 'regular', 'book',
                       'medium', 'roman', 'semibold', 'demibold', 'demi',
                       'bold', 'heavy', 'extra bold', 'black']
        label = Tk.Label(holder, text='Font Name:')
        label.grid(column=0, row=0)
        label = Tk.Label(holder, text='Font Size:')
        label.grid(column=0, row=1)
        label = Tk.Label(holder, text='Font Weight:')
        label.grid(column=0, row=2)
        label = Tk.Label(holder, text='Font Colour:')
        label.grid(column=0, row=3)
        label = Tk.Label(holder, text='Label Text:')
        label.grid(column=0, row=4)
        label = Tk.Label(holder, text='x position:')
        label.grid(column=0, row=5)
        label = Tk.Label(holder, text='y position:')
        label.grid(column=0, row=6)
        self.label_font_name_list = Tk.ttk.Combobox(holder, width=20)
        self.label_font_name_list.grid(column=1, row=0)
        self.label_font_name_list['values'] = fontnames
        self.label_font_name_list.set('sans-serif')
#        self.label_font_name_list.set('times new roman')
        for loop in range(len(fontnames)):
            if self.plot_labels[ind1]['font'] == fontnames[loop]:
                self.label_font_name_list.set(fontnames[loop])
        self.label_font_size_list = Tk.ttk.Combobox(holder, width=20)
        self.label_font_size_list.grid(column=1, row=1)
        self.label_font_size_list['values'] = fontsizes
        self.label_font_size_list.set('12')
        for loop in range(len(fontsizes)):
            if self.plot_labels[ind1]['size'] == int(fontsizes[loop]):
                self.label_font_size_list.set(fontsizes[loop])
        self.label_font_weight_list = Tk.ttk.Combobox(holder, width=20)
        self.label_font_weight_list.grid(column=1, row=2)
        self.label_font_weight_list['values'] = fontweights
        for loop in range(len(fontweights)):
            if self.plot_labels[ind1]['fontweight'] == fontweights[loop]:
                self.label_font_weight_list.set(fontweights[loop])
        self.label_font_colour_list = Tk.ttk.Combobox(holder, width=20)
        self.label_font_colour_list.grid(column=1, row=3)
        self.label_font_colour_list['values'] = self.colourset
        self.label_font_colour_list.set('black')
        for loop in range(len(self.colourset)):
            if self.plot_labels[ind1]['colour'] == self.colourset[loop]:
                self.label_font_colour_list.set(self.colourset[loop])
        self.label_text_entry = Tk.Entry(holder, width=20)
        self.label_text_entry.grid(column=1, row=4)
        self.label_text_entry.insert(0, self.plot_labels[ind1]['labelstring'])
        self.label_x_position = Tk.Entry(holder, width=20)
        self.label_x_position.grid(column=1, row=5)
        if self.xparameters[self.plot_labels[ind1]['plot']-1]['hybridlog']:
            xpos1 = inverse_hybrid_transform(
                self.plot_labels[ind1]['xposition'])
        else:
            xpos1 = self.plot_labels[ind1]['xposition']
        self.label_x_position.insert(0, str(xpos1))
        self.label_y_position = Tk.Entry(holder, width=20)
        self.label_y_position.grid(column=1, row=6)
        if self.yparameters[self.plot_labels[ind1]['plot']-1]['hybridlog']:
            xpos1 = inverse_hybrid_transform(
                self.plot_labels[ind1]['yposition'])
        else:
            ypos1 = self.plot_labels[ind1]['yposition']
        self.label_y_position.insert(0, str(ypos1))
        bholder = Tk.Frame(label_property_window)
        bholder.pack(side=Tk.TOP)
        set_button = Tk.Button(bholder, text='Set Properties',
                               command=lambda: self.set_label_values(ind1))
        set_button.pack()
        close_button = Tk.Button(bholder, text='Close Window',
                                 command=label_property_window.destroy)
        close_button.pack()

    def set_label_values(self, ind1):
        """
        Read parameters from the label input fields and apply them.

        Parameters
        ----------
            ind1 : an integer value >= 0, the index of the label that is
                   having the properties set

        Returns
        -------
            Nothing

        """
        try:
            font = self.label_font_name_list.get()
            fontweight = self.label_font_weight_list.get()
            fontsize = int(self.label_font_size_list.get())
            fontcolour = self.label_font_colour_list.get()
            labelstring = self.label_text_entry.get()
            xpos1 = float(self.label_x_position.get())
            if self.xparameters[self.plot_labels[ind1]['plot']-1]['hybridlog']:
                xpos1 = hybrid_transform(xpos1)
            ypos1 = float(self.label_y_position.get())
            if self.yparameters[self.plot_labels[ind1]['plot']-1]['hybridlog']:
                ypos1 = hybrid_transform(ypos1)
            self.plot_labels[ind1]['xposition'] = xpos1
            self.plot_labels[ind1]['yposition'] = ypos1
            self.plot_labels[ind1]['labelstring'] = labelstring
            self.plot_labels[ind1]['colour'] = fontcolour
            self.plot_labels[ind1]['size'] = fontsize
            self.plot_labels[ind1]['font'] = font
            self.plot_labels[ind1]['fontweight'] = fontweight
            self.make_plot()
        except:
            pass

    def make_plot(self):
        """
        Read the parameters and the data sets, and produce the plot.

        This is the main plotting routine where the graphs are produced.  It
        looks at the (many) parameters to determine how the plot should be
        made.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        myfont = {'family': self.fontname[self.current_plot-1],
                  'weight': self.fontweight[self.current_plot-1],
                  'size': self.fontsize[self.current_plot-1]}
        matplotlib.rc('font', **myfont)
        self.subplot[self.current_plot-1].clear()
        if self.hide_subplot[self.current_plot-1]:
            self.subplot[self.current_plot-1].axis('off')
        logxflag = self.xparameters[self.current_plot-1]['logarithmic']
        logyflag = self.yparameters[self.current_plot-1]['logarithmic']
        hlogxflag = self.xparameters[self.current_plot-1]['hybridlog']
        hlogyflag = self.yparameters[self.current_plot-1]['hybridlog']
        invertxflag = self.xparameters[self.current_plot-1]['invert']
        invertyflag = self.yparameters[self.current_plot-1]['invert']
        hidexflag = self.xparameters[self.current_plot-1]['hide']
        hideyflag = self.yparameters[self.current_plot-1]['hide']
        hidexticksflag = self.xparameters[self.current_plot-1]['hideticks']
        hideyticksflag = self.yparameters[self.current_plot-1]['hideticks']
        hidexlabelsflag = self.xparameters[self.current_plot-1]['hidelabels']
        hideylabelsflag = self.yparameters[self.current_plot-1]['hidelabels']
        inversexticksflag = \
            self.xparameters[self.current_plot-1]['inverseticks']
        inverseyticksflag = \
            self.yparameters[self.current_plot-1]['inverseticks']
        bothxticksflag = self.xparameters[self.current_plot-1]['bothticks']
        bothyticksflag = self.yparameters[self.current_plot-1]['bothticks']
        oppositexflag = self.xparameters[self.current_plot-1]['oppositeaxis']
        oppositeyflag = self.yparameters[self.current_plot-1]['oppositeaxis']
        try:
            xminorticks = \
                float(self.xparameters[self.current_plot-1]['minorticks'])
        except:
            xminorticks = 0.0
        try:
            yminorticks = \
                float(self.yparameters[self.current_plot-1]['minorticks'])
        except:
            yminorticks = 0.0
        for loop in range(self.nsets):
            if (self.set_properties[loop]['display']) \
               and (self.set_properties[loop]['plot'] == self.current_plot):
                if self.set_properties[loop]['errors']:
                    xerrors = numpy.zeros((2, len(self.xdata[loop]['values'])),
                                          dtype=numpy.float32)
                    xerrors[0, :] = numpy.copy(self.xdata[loop]['lowerror'])
                    xerrors[1, :] = numpy.copy(self.xdata[loop]['higherror'])
                    yerrors = numpy.zeros((2, len(self.xdata[loop]['values'])),
                                          dtype=numpy.float32)
                    yerrors[0, :] = numpy.copy(self.ydata[loop]['lowerror'])
                    yerrors[1, :] = numpy.copy(self.ydata[loop]['higherror'])
                try:
                    if self.set_properties[loop]['symbol'] == 'histogram':
                        barwidth = self.xdata[loop]['values'][1] - \
                            self.xdata[loop]['values'][0]
                        xoff = self.xdata[loop]['values'][1:] - \
                            self.xdata[loop]['values'][0:-1]
                        xoff = numpy.append(xoff, xoff[-1])
                        self.subplot[self.current_plot-1].bar(
                            self.xdata[loop]['values'] - xoff/2.,
                            self.ydata[loop]['values'],
                            width=barwidth,
                            color='white',
                            edgecolor=self.set_properties[loop]['colour'],
                            linewidth=1)
                    elif logyflag == 0 and logxflag == 0 and hlogxflag == 0 \
                            and hlogyflag == 0:
                        if self.set_properties[loop]['symbol'] is None:
                            self.subplot[self.current_plot-1].plot(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        elif self.set_properties[loop]['linestyle'] is None:
                            self.subplot[self.current_plot-1].plot(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle='none', markersize=
                                self.set_properties[loop]['symbolsize'])
                        else:
                            self.subplot[self.current_plot-1].plot(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                markersize=
                                self.set_properties[loop]['symbolsize'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        if self.set_properties[loop]['errors']:
                            self.subplot[self.current_plot-1].errorbar(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'], yerrors, xerrors,
                                fmt='None',
                                ecolor=self.set_properties[loop]['colour'],
                                elinewidth=
                                self.set_properties[loop]['linewidth'])
                        if xminorticks > 0:
                            self.subplot[self.current_plot-1].xaxis.\
                                set_minor_locator(MultipleLocator(xminorticks))
                            xtl = int(
                                self.xparameters[
                                    self.current_plot-1]['ticklength'] // 2)
                            self.subplot[self.current_plot-1].tick_params(
                                axis='x', which='minor', length=xtl)
                        if yminorticks > 0:
                            self.subplot[self.current_plot-1].yaxis.\
                                set_minor_locator(MultipleLocator(yminorticks))
                            ytl = int(
                                self.yparameters[
                                    self.current_plot-1]['ticklength'] // 2)
                            self.subplot[self.current_plot-1].tick_params(
                                axis='x', which='minor', length=ytl)
                    elif hlogyflag == 0 and hlogxflag == 1:
                        newxvalues = hybrid_transform(
                            self.xdata[loop]['values'])
                        if logyflag == 0:
                            if self.set_properties[loop]['symbol'] is None:
                                self.subplot[self.current_plot-1].plot(
                                    newxvalues, self.ydata[loop]['values'],
                                    color=self.set_properties[loop]['colour'],
                                    linestyle=
                                    self.set_properties[loop]['linestyle'],
                                    linewidth=
                                    self.set_properties[loop]['linewidth'])
                            elif self.set_properties[loop]['linestyle'] is \
                                    None:
                                self.subplot[self.current_plot-1].plot(
                                    newxvalues, self.ydata[loop]['values'],
                                    color=self.set_properties[loop]['colour'],
                                    marker=
                                    self.set_properties[loop]['symbol'],
                                    linestyle='none',
                                    markersize=
                                    self.set_properties[loop]['symbolsize'])
                            else:
                                self.subplot[self.current_plot-1].plot(
                                    newxvalues, self.ydata[loop]['values'],
                                    color=self.set_properties[loop]['colour'],
                                    marker=
                                    self.set_properties[loop]['symbol'],
                                    linestyle=
                                    self.set_properties[loop]['linestyle'],
                                    markersize=
                                    self.set_properties[loop]['symbolsize'],
                                    linewidth=
                                    self.set_properties[loop]['linewidth'])
                            if yminorticks > 0:
                                self.subplot[self.current_plot-1].yaxis.\
                                    set_minor_locator(
                                        MultipleLocator(yminorticks))
                                ytl = int(
                                    self.xparameters[self.current_plot-1]
                                    ['ticklength'] // 2)
                                self.subplot[self.current_plot-1].\
                                    tick_params(axis='y', which='minor',
                                                length=ytl)
                        else:
                            if self.set_properties[loop]['symbol'] is None:
                                self.subplot[self.current_plot-1].semilogy(
                                    newxvalues, self.ydata[loop]['values'],
                                    color=self.set_properties[loop]['colour'],
                                    linestyle=
                                    self.set_properties[loop]['linestyle'],
                                    linewidth=
                                    self.set_properties[loop]['linewidth'])
                            elif self.set_properties[loop]['linestyle'] is \
                                    None:
                                self.subplot[
                                    self.current_plot-1].semilogy(
                                        newxvalues, self.ydata[loop]['values'],
                                        color=
                                        self.set_properties[loop]['colour'],
                                        marker=
                                        self.set_properties[loop]['symbol'],
                                        linestyle='none', markersize=
                                        self.set_properties[loop]['symbolsize']
                                    )
                            else:
                                self.subplot[self.current_plot-1].semilogy(
                                    newxvalues, self.ydata[loop]['values'],
                                    color=self.set_properties[loop]['colour'],
                                    marker=
                                    self.set_properties[loop]['symbol'],
                                    linestyle=
                                    self.set_properties[loop]['linestyle'],
                                    markersize=
                                    self.set_properties[loop]['symbolsize'],
                                    linewidth=
                                    self.set_properties[loop]['linewidth'])
                    elif hlogyflag == 1 and hlogxflag == 0:
                        newyvalues = hybrid_transform(
                            self.ydata[loop]['values'])
                        if logxflag == 0:
                            if self.set_properties[loop]['symbol'] is None:
                                self.subplot[self.current_plot-1].plot(
                                    self.xdata[loop]['values'], newyvalues,
                                    color=self.set_properties[loop]['colour'],
                                    linestyle=
                                    self.set_properties[loop]['linestyle'],
                                    linewidth=
                                    self.set_properties[loop]['linewidth'])
                            elif self.set_properties[loop]['linestyle'] is \
                                    None:
                                self.subplot[self.current_plot-1].plot(
                                    self.xdata[loop]['values'], newyvalues,
                                    color=self.set_properties[loop]['colour'],
                                    marker=
                                    self.set_properties[loop]['symbol'],
                                    linestyle='none',
                                    markersize=
                                    self.set_properties[loop]['symbolsize'])
                            else:
                                self.subplot[self.current_plot-1].plot(
                                    self.xdata[loop]['values'], newyvalues,
                                    color=self.set_properties[loop]['colour'],
                                    marker=
                                    self.set_properties[loop]['symbol'],
                                    linestyle=
                                    self.set_properties[loop]['linestyle'],
                                    markersize=
                                    self.set_properties[loop]['symbolsize'],
                                    linewidth=
                                    self.set_properties[loop]['linewidth'])
                            if xminorticks > 0:
                                self.subplot[self.current_plot-1].xaxis.\
                                    set_minor_locator(MultipleLocator(
                                        xminorticks))
                                xtl = int(
                                    self.xparameters[
                                        self.current_plot-1]['ticklength']
                                    // 2)
                                self.subplot[self.current_plot-1].tick_params(
                                    axis='x', which='minor', length=xtl)
                        else:
                            if self.set_properties[loop]['symbol'] is None:
                                self.subplot[self.current_plot-1].semilogx(
                                    self.xdata[loop]['values'],
                                    newyvalues,
                                    color=self.set_properties[loop]['colour'],
                                    linestyle=
                                    self.set_properties[loop]['linestyle'],
                                    linewidth=
                                    self.set_properties[loop]['linewidth'])
                            elif self.set_properties[loop]['linestyle'] is \
                                    None:
                                self.subplot[self.current_plot-1].semilogx(
                                    self.xdata[loop]['values'],
                                    newyvalues,
                                    color=self.set_properties[loop]['colour'],
                                    marker=self.set_properties[loop]['symbol'],
                                    linestyle='none',
                                    markersize=
                                    self.set_properties[loop]['symbolsize'])
                            else:
                                self.subplot[self.current_plot-1].semilogx(
                                    self.xdata[loop]['values'],
                                    newyvalues,
                                    color=self.set_properties[loop]['colour'],
                                    marker=
                                    self.set_properties[loop]['symbol'],
                                    linestyle=
                                    self.set_properties[loop]['linestyle'],
                                    markersize=
                                    self.set_properties[loop]['symbolsize'],
                                    linewidth=
                                    self.set_properties[loop]['linewidth'])
                    elif hlogyflag == 1 and hlogxflag == 1:
                        newxvalues = hybrid_transform(
                            self.xdata[loop]['values'])
                        newyvalues = hybrid_transform(
                            self.ydata[loop]['values'])
                        if self.set_properties[loop]['symbol'] is None:
                            self.subplot[self.current_plot-1].plot(
                                newxvalues, newyvalues,
                                color=self.set_properties[loop]['colour'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        elif self.set_properties[loop]['linestyle'] is None:
                            self.subplot[self.current_plot-1].plot(
                                newxvalues, newyvalues,
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle='none',
                                markersize=
                                self.set_properties[loop]['symbolsize'])
                        else:
                            self.subplot[self.current_plot-1].plot(
                                newxvalues, newyvalues,
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                markersize=
                                self.set_properties[loop]['symbolsize'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                    elif logyflag == 0 and logxflag == 1 and hlogxflag == 0:
                        if self.set_properties[loop]['symbol'] is None:
                            self.subplot[self.current_plot-1].semilogx(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        elif self.set_properties[loop]['linestyle'] is None:
                            self.subplot[self.current_plot-1].semilogx(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle='None',
                                markersize=
                                self.set_properties[loop]['symbolsize'])
                        else:
                            self.subplot[self.current_plot-1].semilogx(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                markersize=
                                self.set_properties[loop]['symbolsize'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        if self.set_properties[loop]['errors']:
                            self.subplot[self.current_plot-1].errorbar(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'], xerrors, yerrors,
                                fmt=None, ecolor=
                                self.set_properties[loop]['colour'],
                                elinewidth=
                                self.set_properties[loop]['linewidth'])
                        if yminorticks > 0:
                            self.subplot[self.current_plot-1].\
                                yaxis.set_minor_locator(
                                    MultipleLocator(yminorticks))
                            ytl = int(
                                self.xparameters[self.current_plot-1][
                                    'ticklength'] // 2)
                            self.subplot[self.current_plot-1].tick_params(
                                axis='y', which='minor', length=ytl)
                    elif logxflag == 0 and logyflag == 1 and hlogyflag == 0:
                        if self.set_properties[loop]['symbol'] is None:
                            self.subplot[self.current_plot-1].semilogy(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        elif self.set_properties[loop]['linestyle'] is None:
                            self.subplot[self.current_plot-1].semilogy(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle='None',
                                markersize=
                                self.set_properties[loop]['symbolsize'])
                        else:
                            self.subplot[self.current_plot-1].semilogy(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                markersize=
                                self.set_properties[loop]['symbolsize'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        if self.set_properties[loop]['errors']:
                            self.subplot[self.current_plot-1].errorbar(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'], xerrors, yerrors,
                                fmt=None,
                                ecolor=self.set_properties[loop]['colour'],
                                elinewidth=
                                self.set_properties[loop]['linewidth'])
                        if xminorticks > 0:
                            self.subplot[self.current_plot-1].xaxis.\
                                set_minor_locator(MultipleLocator(xminorticks))
                            xtl = int(
                                self.xparameters[
                                    self.current_plot-1]['ticklength'] // 2)
                            self.subplot[self.current_plot-1].tick_params(
                                axis='x', which='minor', length=xtl)
                    elif logxflag == 1 and logyflag == 1:
                        if self.set_properties[loop]['symbol'] is None:
                            self.subplot[self.current_plot-1].loglog(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=
                                self.set_properties[loop]['colour'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        elif self.set_properties[loop]['linestyle'] is None:
                            self.subplot[self.current_plot-1].loglog(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle='None',
                                markersize=
                                self.set_properties[loop]['symbolsize'])
                        else:
                            self.subplot[self.current_plot-1].loglog(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                markersize=
                                self.set_properties[loop]['symbolsize'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        if self.set_properties[loop]['errors']:
                            self.subplot[self.current_plot-1].errorbar(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                xerrors, yerrors, fmt=None,
                                ecolor=self.set_properties[loop]['colour'],
                                elinewidth=
                                self.set_properties[loop]['linewidth'])
                    else:
                        # One should not get here, but if the flags do not
                        # correspond to the expected options use a linear
                        # plot as the default
                        if self.set_properties[loop]['symbol'] is None:
                            self.subplot[self.current_plot-1].plot(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        elif self.set_properties[loop]['linestyle'] is None:
                            self.subplot[self.current_plot-1].plot(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle='none',
                                markersize=
                                self.set_properties[loop]['symbolsize'])
                        else:
                            self.subplot[self.current_plot-1].plot(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'],
                                color=self.set_properties[loop]['colour'],
                                marker=self.set_properties[loop]['symbol'],
                                linestyle=
                                self.set_properties[loop]['linestyle'],
                                markersize=
                                self.set_properties[loop]['symbolsize'],
                                linewidth=
                                self.set_properties[loop]['linewidth'])
                        if self.set_properties[loop]['errors']:
                            self.subplot[self.current_plot-1].errorbar(
                                self.xdata[loop]['values'],
                                self.ydata[loop]['values'], yerrors, xerrors,
                                fmt='None',
                                ecolor=self.set_properties[loop]['colour'],
                                elinewidth=
                                self.set_properties[loop]['linewidth'])
                        if xminorticks > 0:
                            self.subplot[self.current_plot-1].xaxis.\
                                set_minor_locator(MultipleLocator(xminorticks))
                            xtl = int(self.xparameters[
                                self.current_plot-1]['ticklength'] // 2)
                            self.subplot[self.current_plot-1].tick_params(
                                axis='x', which='minor', length=xtl)
                        if yminorticks > 0:
                            self.subplot[
                                self.current_plot-1].yaxis.set_minor_locator(
                                    MultipleLocator(yminorticks))
                            ytl = int(self.yparameters[
                                self.current_plot-1]['ticklength'] // 2)
                            self.subplot[self.current_plot-1].tick_params(
                                axis='y', which='minor', length=ytl)
                except:
                    pass
        if bothyticksflag == 1:
            self.subplot[self.current_plot-1].tick_params(
                left=True, right=True, which='both')
        else:
            self.subplot[self.current_plot-1].tick_params(
                left=True, right=False, which='both')
        if bothxticksflag == 1:
            self.subplot[self.current_plot-1].tick_params(
                bottom=True, top=True, which='both')
        else:
            self.subplot[self.current_plot-1].tick_params(
                bottom=True, top=False, which='both')
        for n1 in range(self.number_of_lines):
            if self.plot_lines[n1]['plot'] == self.current_plot:
                xlvalues = numpy.asarray([self.plot_lines[n1]['xstart'],
                                          self.plot_lines[n1]['xend']])
                if hlogxflag == 1:
                    xlvalues[0] = hybrid_transform(xlvalues[0])
                    xlvalues[1] = hybrid_transform(xlvalues[1])
                ylvalues = numpy.asarray([self.plot_lines[n1]['ystart'],
                                          self.plot_lines[n1]['yend']])
                if hlogyflag == 1:
                    ylvalues[0] = hybrid_transform(ylvalues[0])
                    ylvalues[1] = hybrid_transform(ylvalues[1])
                self.subplot[self.current_plot-1].plot(
                    xlvalues, ylvalues,
                    color=self.plot_lines[n1]['line_colour'],
                    linestyle=self.plot_lines[n1]['line_type'],
                    linewidth=self.plot_lines[n1]['line_thickness'])
        patches = []
        for n1 in range(self.number_of_boxes):
            if self.plot_boxes[n1]['plot'] == self.current_plot:
                xcorner = self.plot_boxes[n1]['xstart']
                ycorner = self.plot_boxes[n1]['ystart']
                delx = self.plot_boxes[n1]['xend'] \
                    - self.plot_boxes[n1]['xstart']
                dely = self.plot_boxes[n1]['yend'] \
                    - self.plot_boxes[n1]['ystart']
                if delx < 0.:
                    xcorner = xcorner + delx
                    delx = -delx
                if dely < 0.:
                    ycorner = ycorner + dely
                    dely = -dely
                if hlogxflag == 1:
                    xc1 = xcorner
                    xc2 = xcorner+delx
                    xcorner = hybrid_transform(xc1)
                    delx = hybrid_transform(xc2) - xcorner
                if hlogyflag == 1:
                    yc1 = ycorner
                    yc2 = ycorner+dely
                    ycorner = hybrid_transform(yc1)
                    dely = hybrid_transform(yc2) - ycorner
                rect = Rectangle(
                    (xcorner, ycorner), delx, dely,
                    angle=self.plot_boxes[n1]['rotation'],
                    edgecolor=self.plot_boxes[n1]['line_colour'],
                    linestyle=self.plot_boxes[n1]['line_type'],
                    linewidth=self.plot_boxes[n1]['line_thickness'],
                    facecolor=self.plot_boxes[n1]['fill_colour'])
                self.subplot[self.current_plot-1].add_artist(rect)
                patches.append(rect)
        if len(patches) > 0:
            pc = PatchCollection(patches)
            self.subplot[self.current_plot-1].add_collection(pc)
        patches = []
        for n1 in range(self.number_of_ellipses):
            if self.plot_ellipses[n1]['plot'] == self.current_plot:
                xcenter = self.plot_ellipses[n1]['xposition']
                ycenter = self.plot_ellipses[n1]['yposition']
                delx = self.plot_ellipses[n1]['major']
                dely = self.plot_ellipses[n1]['minor']
                if hlogxflag == 1:
                    xc1 = xcorner
                    xc2 = xcorner+delx
                    xcorner = hybrid_transform(xc1)
                    delx = hybrid_transform(xc2) - xcorner
                if hlogyflag == 1:
                    yc1 = ycorner
                    yc2 = ycorner+dely
                    ycorner = hybrid_transform(yc1)
                    dely = hybrid_transform(yc2) - ycorner
                ellip = Ellipse(
                    (xcenter, ycenter), delx, dely,
                    angle=self.plot_ellipses[n1]['rotation'],
                    edgecolor=self.plot_ellipses[n1]['line_colour'],
                    linestyle=self.plot_ellipses[n1]['line_type'],
                    linewidth=self.plot_ellipses[n1]['line_thickness'],
                    facecolor=self.plot_ellipses[n1]['fill_colour'])
                self.subplot[self.current_plot-1].add_artist(ellip)
                patches.append(ellip)
        if len(patches) > 0:
            pc = PatchCollection(patches)
            self.subplot[self.current_plot-1].add_collection(pc)
        patches = []
        for n1 in range(self.number_of_vectors):
            if self.plot_vectors[n1]['plot'] == self.current_plot:
                x1 = self.plot_vectors[n1]['xstart']
                y1 = self.plot_vectors[n1]['ystart']
                xlength = self.plot_vectors[n1]['xend'] \
                    - self.plot_vectors[n1]['xstart']
                ylength = self.plot_vectors[n1]['yend'] \
                    - self.plot_vectors[n1]['ystart']
                delx = self.plot_vectors[n1]['delx']
                dely = self.plot_vectors[n1]['dely']
                if hlogxflag == 1:
                    xs1 = hybrid_transform(x1)
                    xs2 = hybrid_transform(self.plot_vectors[n1]['xend'])
                    xlength = xs2 - xs1
                    delx = xs2 - hybrid_transform(
                        self.plot_vectors[n1]['xend']-delx)
                if hlogyflag == 1:
                    ys1 = hybrid_transform(y1)
                    ys2 = hybrid_transform(self.plot_vectors[n1]['yend'])
                    ylength = ys2 - ys1
                    delx = ys2 - hybrid_transform(
                        self.plot_vectors[n1]['yend']-dely)
                arrow = FancyArrow(
                    x1, y1, xlength, ylength, head_width=delx,
                    head_length=dely, fill=self.plot_vectors[n1]['fill'],
                    edgecolor=self.plot_vectors[n1]['line_colour'],
                    linestyle=self.plot_vectors[n1]['line_type'],
                    linewidth=self.plot_vectors[n1]['line_thickness'],
                    facecolor=self.plot_vectors[n1]['fill_colour'])
                self.subplot[self.current_plot-1].add_artist(arrow)
                patches.append(arrow)
        if len(patches) > 0:
            pc = PatchCollection(patches)
            self.subplot[self.current_plot-1].add_collection(pc)
        if self.plot_frame[self.current_plot-1] > 0.:
            self.subplot[self.current_plot-1].spines['bottom'].set_linewidth(
                self.plot_frame[self.current_plot-1])
            self.subplot[self.current_plot-1].spines['top'].set_linewidth(
                self.plot_frame[self.current_plot-1])
            self.subplot[self.current_plot-1].spines['left'].set_linewidth(
                self.plot_frame[self.current_plot-1])
            self.subplot[self.current_plot-1].spines['right'].set_linewidth(
                self.plot_frame[self.current_plot-1])
        else:
            self.subplot[self.current_plot-1].spines['bottom'].set_linewidth(
                0.5)
            self.subplot[self.current_plot-1].spines['top'].set_linewidth(
                0.5)
            self.subplot[self.current_plot-1].spines['left'].set_linewidth(
                0.5)
            self.subplot[self.current_plot-1].spines['right'].set_linewidth(
                0.5)
        xisinverted = self.subplot[self.current_plot-1].xaxis_inverted()
        yisinverted = self.subplot[self.current_plot-1].yaxis_inverted()
        if (invertxflag == 1 and (not xisinverted)) or \
                (invertxflag == 0 and xisinverted):
            self.subplot[self.current_plot-1].invert_xaxis()
        if (invertyflag == 1 and (not yisinverted)) or \
                (invertyflag == 0 and yisinverted):
            self.subplot[self.current_plot-1].invert_yaxis()
        if self.matplotlib_rounding.get():
            xmin1, xmax1 = self.subplot[self.current_plot-1].get_xbound()
            ymin1, ymax1 = self.subplot[self.current_plot-1].get_ybound()
            if self.original_range[self.current_plot-1]:
                self.plot_range[self.current_plot-1][0] = xmin1
                self.plot_range[self.current_plot-1][1] = xmax1
                self.plot_range[self.current_plot-1][2] = ymin1
                self.plot_range[self.current_plot-1][3] = ymax1
                self.original_range[self.current_plot-1] = False
        if hlogxflag == 0:
            self.subplot[self.current_plot-1].set_xbound(
                self.plot_range[self.current_plot-1][0],
                self.plot_range[self.current_plot-1][1])
        else:
            xrange = numpy.asarray([self.plot_range[self.current_plot-1][0],
                                    self.plot_range[self.current_plot-1][1]],
                                   dtype=numpy.float32)
            xrange1 = hybrid_transform(xrange)
            self.subplot[self.current_plot-1].set_xbound(xrange1[0],
                                                         xrange1[1])
            newxvalues = hybrid_transform(self.xdata[loop]['values'])
            tickmarks, ticklabels = generate_labels(xrange1)
            self.subplot[self.current_plot-1].set_xticks(tickmarks)
            self.subplot[self.current_plot-1].set_xticklabels(ticklabels)
        if hlogyflag == 0:
            self.subplot[self.current_plot-1].set_ybound(
                self.plot_range[self.current_plot-1][2],
                self.plot_range[self.current_plot-1][3])
        else:
            yrange = numpy.asarray([self.plot_range[self.current_plot-1][2],
                                    self.plot_range[self.current_plot-1][3]],
                                   dtype=numpy.float32)
            yrange1 = hybrid_transform(yrange)
            self.subplot[self.current_plot-1].set_ybound(yrange1[0],
                                                         yrange1[1])
            newyvalues = hybrid_transform(self.ydata[loop]['values'])
            tickmarks, ticklabels = generate_labels(yrange1)
            self.subplot[self.current_plot-1].set_yticks(tickmarks)
            self.subplot[self.current_plot-1].set_yticklabels(ticklabels)
        self.subplot[self.current_plot-1].set_xlabel(
            self.xparameters[self.current_plot-1]['label'],
            family=self.fontname[self.current_plot-1],
            size=self.fontsize[self.current_plot-1],
            weight=self.fontweight[self.current_plot-1])
        self.subplot[self.current_plot-1].set_ylabel(
            self.yparameters[self.current_plot-1]['label'],
            family=self.fontname[self.current_plot-1],
            size=self.fontsize[self.current_plot-1],
            weight=self.fontweight[self.current_plot-1])
        try:
            self.subplot[self.current_plot-1].set_title(
                self.title[self.current_plot-1],
                family=self.fontname[self.current_plot-1],
                size=self.fontsize[self.current_plot-1],
                weight=self.fontweight[self.current_plot-1])
        except:
            pass
        # Adjust the margins if the plot_margin value is set.
        left = self.bounding_box[self.current_plot-1][0] + self.plot_margin
        right = self.bounding_box[self.current_plot-1][0] \
            + self.bounding_box[self.current_plot-1][2] - self.plot_margin
        bottom = self.bounding_box[self.current_plot-1][1] + self.plot_margin
        top = self.bounding_box[self.current_plot-1][1] \
            + self.bounding_box[self.current_plot-1][3] - self.plot_margin
        self.subplot[self.current_plot-1].set_position(
            [left, bottom, right-left, top-bottom], which='both')
        if self.number_of_labels > 0:
            for loop in range(self.number_of_labels):
                if self.current_plot == self.plot_labels[loop]['plot']:
                    xpos1 = self.plot_labels[loop]['xposition']
                    if self.xparameters[self.current_plot-1]['hybridlog'] == 1:
                        xpos1 = hybrid_transform(xpos1)
                    ypos1 = self.plot_labels[loop]['yposition']
                    if self.xparameters[self.current_plot-1]['hybridlog'] == 1:
                        ypos1 = hybrid_transform(ypos1)
                    self.subplot[self.current_plot-1].text(
                        xpos1, ypos1,
                        self.plot_labels[loop]['labelstring'],
                        {'color': self.plot_labels[loop]['colour'],
                         'fontsize': self.plot_labels[loop]['size'],
                         'fontfamily': self.plot_labels[loop]['font'],
                         'fontweight': self.plot_labels[loop]['fontweight']})
        try:
            legend_flag = self.legend_variable[self.current_plot-1].get()
        except:
            legend_flag = 0
        if legend_flag:
            self.generate_legend(None)
            legend_option = self.legend_options[self.current_plot-1].get()
            legend_position = None
            if legend_option == 'user':
                try:
                    str1 = self.legend_position_field.get()
                    values = str1.split()
                    if len(values) == 2:
                        xlpos = float(values[0])
                        ylpos = float(values[1])
                        legend_position = [xlpos, ylpos]
                        self.legend_user_position[self.current_plot-1] = \
                            legend_position
                except:
                    legend_option = 'best'
            else:
                legend_option = self.legend_position[self.current_plot-1]
                if legend_option is None:
                    legend_position = 'best'
            try:
                legend_frame = self.legend_frame[self.current_plot-1].get()
            except:
                legend_frame = 0
            if legend_position is None:
                self.subplot[self.current_plot-1].legend(
                    handles=self.legend_handles[self.current_plot-1],
                    labels=self.legend_labels[self.current_plot-1],
                    loc=legend_option, frameon=legend_frame)
            else:
                self.subplot[self.current_plot-1].legend(
                    handles=self.legend_handles[self.current_plot-1],
                    labels=self.legend_labels[self.current_plot-1],
                    loc=legend_position, frameon=legend_frame)
        if oppositexflag == 1:
            self.subplot[self.current_plot-1].get_xaxis().set_ticks_position(
                "top")
            self.subplot[self.current_plot-1].get_xaxis().set_label_position(
                "top")
        if oppositeyflag == 1:
            self.subplot[self.current_plot-1].get_yaxis().set_ticks_position(
                "right")
            self.subplot[self.current_plot-1].get_yaxis().set_label_position(
                "right")
        if inversexticksflag == 1:
            self.subplot[self.current_plot-1].tick_params(
                axis='x', direction='in', length=
                self.xparameters[self.current_plot-1]['ticklength'])
            self.subplot[self.current_plot-1].tick_params(
                axis='x', direction='in', which='minor')
        else:
            self.subplot[self.current_plot-1].tick_params(
                axis='x', direction='out', length=
                self.xparameters[self.current_plot-1]['ticklength'])
            self.subplot[self.current_plot-1].tick_params(
                axis='x', direction='out', which='minor')
        if inverseyticksflag == 1:
            self.subplot[self.current_plot-1].tick_params(
                axis='y', direction='in', length=
                self.yparameters[self.current_plot-1]['ticklength'])
            self.subplot[self.current_plot-1].tick_params(axis='y',
                                                          direction='in',
                                                          which='minor')
        else:
            self.subplot[self.current_plot-1].tick_params(
                axis='y', direction='out', length=
                self.yparameters[self.current_plot-1]['ticklength'])
            self.subplot[self.current_plot-1].tick_params(axis='y',
                                                          direction='out',
                                                          which='minor')
        if hidexticksflag == 1:
            self.subplot[self.current_plot-1].get_xaxis().set_ticks([])
        if hideyticksflag == 1:
            self.subplot[self.current_plot-1].get_yaxis().set_ticks([])
        if hidexlabelsflag == 1:
            self.subplot[self.current_plot-1].get_xaxis().set_ticklabels([])
        if hideylabelsflag == 1:
            self.subplot[self.current_plot-1].get_yaxis().set_ticklabels([])
        if hidexflag == 1:
            self.subplot[self.current_plot-1].get_xaxis().set_visible(False)
        else:
            self.subplot[self.current_plot-1].get_xaxis().set_visible(True)
        if hideyflag == 1:
            self.subplot[self.current_plot-1].get_yaxis().set_visible(False)
        else:
            self.subplot[self.current_plot-1].get_yaxis().set_visible(True)
        if self.equal_aspect[self.current_plot-1]:
            self.figure.gca().set_aspect('equal', adjustable='box')
        else:
            self.figure.gca().set_aspect('auto')
        self.canvas.draw()
        # The follow is duplicate information, but may be used to set tick
        # intervals....
        if not self.matplotlib_rounding.get():
            self.xparameters[self.current_plot-1]['minimum'] = \
                self.plot_range[self.current_plot-1][0]
            self.xparameters[self.current_plot-1]['maximum'] = \
                self.plot_range[self.current_plot-1][1]
            self.yparameters[self.current_plot-1]['minimum'] = \
                self.plot_range[self.current_plot-1][2]
            self.yparameters[self.current_plot-1]['maximum'] = \
                self.plot_range[self.current_plot-1][3]
        else:
            xr1, xr2 = self.subplot[self.current_plot-1].get_xbound()
            yr1, yr2 = self.subplot[self.current_plot-1].get_ybound()
            self.xparameters[self.current_plot-1]['minimum'] = xr1
            self.xparameters[self.current_plot-1]['maximum'] = xr2
            self.yparameters[self.current_plot-1]['minimum'] = yr1
            self.yparameters[self.current_plot-1]['maximum'] = yr2

    def make_hess_plot(self):
        """
        Make a two-dimensional histogram plot.

        This routine creates a new plot window within which a Hess plot (i.e.
        a two-dimensional histogram) is made for the current active plot.

        The new window has options for control of the two-dimensional
        histogram.

        No parameters are passed to the routine, and no values are returned.
        """
        try:
            hesswindow = Tk.Toplevel()
            hesswindow.config(bg=BGCOL)
            xmin = self.plot_range[self.current_plot-1][0]
            xmax = self.plot_range[self.current_plot-1][1]
            ymin = self.plot_range[self.current_plot-1][2]
            ymax = self.plot_range[self.current_plot-1][3]
            xp = None
            yp = None
            for loop in range(self.nsets):
                if (self.set_properties[loop]['display']) and \
                   (self.set_properties[loop]['plot'] == self.current_plot):
                    if xp is None:
                        xp = numpy.copy(self.xdata[loop]['values'])
                        yp = numpy.copy(self.ydata[loop]['values'])
                    else:
                        xp = numpy.append(xp, self.xdata[loop]['values'])
                        yp = numpy.append(yp, self.ydata[loop]['values'])
            self.hessLabelText = Tk.StringVar()
            self.hessLabel = Tk.Label(
                hesswindow,
                textvariable=self.hessLabelText, anchor=Tk.N, width=70)
            self.hessLabel.pack()
            self.hessLabelText.set("Value:")
            self.p1 = Figure(figsize=(6, 6), dpi=100)
            sp1 = self.p1.add_subplot(1, 1, 1)
            c1 = FigureCanvasTkAgg(self.p1, master=hesswindow)
            c1.mpl_connect("motion_notify_event", self.hess_position)
            try:
                npixels = int(self.npixelfield.get())
                if (npixels < 50) or (npixels > 5000):
                    print('Error: too many or too few pixels requested '
                          + '(%d).  Using the default value of 500.' %
                          (npixels))
                    npixels = 500
                    self.npixelfield.delete(0, Tk.END)
                    self.npixelfield.insert(0, '500')
            except:
                npixels = 500
            self.hist2d, self.xedges, self.yedges, image = sp1.hist2d(
                xp, yp, npixels, range=[[xmin, xmax], [ymin, ymax]],
                norm=LogNorm())
            sp1.set_xlabel(self.xparameters[self.current_plot-1]['label'],
                           family=self.fontname[self.current_plot-1],
                           size=self.fontsize[self.current_plot-1],
                           weight=self.fontweight[self.current_plot-1])
            sp1.set_ylabel(self.yparameters[self.current_plot-1]['label'],
                           family=self.fontname[self.current_plot-1],
                           size=self.fontsize[self.current_plot-1],
                           weight=self.fontweight[self.current_plot-1])
            invertxflag = self.xparameters[self.current_plot-1]['invert']
            invertyflag = self.yparameters[self.current_plot-1]['invert']
            hidexflag = self.xparameters[self.current_plot-1]['hide']
            hideyflag = self.yparameters[self.current_plot-1]['hide']
            hidexticksflag = self.xparameters[self.current_plot-1]['hideticks']
            hideyticksflag = self.yparameters[self.current_plot-1]['hideticks']
            hidexlabelsflag = self.xparameters[
                self.current_plot-1]['hidelabels']
            hideylabelsflag = self.yparameters[
                self.current_plot-1]['hidelabels']
            inversexticksflag = self.xparameters[
                self.current_plot-1]['inverseticks']
            inverseyticksflag = self.yparameters[
                self.current_plot-1]['inverseticks']
            bothxticksflag = self.xparameters[self.current_plot-1]['bothticks']
            bothyticksflag = self.yparameters[self.current_plot-1]['bothticks']
            try:
                xminorticks = float(self.xparameters[
                    self.current_plot-1]['minorticks'])
            except:
                xminorticks = 0.0
            try:
                yminorticks = float(self.yparameters[
                    self.current_plot-1]['minorticks'])
            except:
                yminorticks = 0.0
            if invertxflag == 1:
                sp1.invert_xaxis()
            if invertyflag == 1:
                sp1.invert_yaxis()
            if bothyticksflag == 1:
                sp1.tick_params(left=True, right=True, which='both')
            else:
                sp1.tick_params(left=True, right=False, which='both')
            if bothxticksflag == 1:
                sp1.tick_params(bottom=True, top=True, which='both')
            else:
                sp1.tick_params(bottom=True, top=False, which='both')
            if xminorticks > 0:
                sp1.xaxis.set_minor_locator(MultipleLocator(xminorticks))
                xtl = int(self.xparameters[self.current_plot-1]['ticklength']
                          // 2)
                sp1.tick_params(axis='x', which='minor', length=xtl)
            if yminorticks > 0:
                sp1.yaxis.set_minor_locator(MultipleLocator(yminorticks))
                ytl = int(self.yparameters[self.current_plot-1]['ticklength']
                          // 2)
                sp1.tick_params(axis='y', which='minor', length=ytl)
            if self.plot_frame[self.current_plot-1] > 0.:
                sp1.spines['bottom'].set_linewidth(
                    self.plot_frame[self.current_plot-1])
                sp1.spines['top'].set_linewidth(
                    self.plot_frame[self.current_plot-1])
                sp1.spines['left'].set_linewidth(
                    self.plot_frame[self.current_plot-1])
                sp1.spines['right'].set_linewidth(
                    self.plot_frame[self.current_plot-1])
            else:
                sp1.spines['bottom'].set_linewidth(0.5)
                sp1.spines['top'].set_linewidth(0.5)
                sp1.spines['left'].set_linewidth(0.5)
                sp1.spines['right'].set_linewidth(0.5)
            if inversexticksflag == 1:
                sp1.tick_params(
                    axis='x', direction='in',
                    length=
                    self.xparameters[self.current_plot-1]['ticklength'])
                sp1.tick_params(axis='x', direction='in', which='minor')
            else:
                sp1.tick_params(
                    axis='x', direction='out',
                    length=self.xparameters[self.current_plot-1]['ticklength'])
                sp1.tick_params(axis='x', direction='out', which='minor')
            if inverseyticksflag == 1:
                sp1.tick_params(
                    axis='y', direction='in',
                    length=self.yparameters[self.current_plot-1]['ticklength'])
                sp1.tick_params(axis='y', direction='in', which='minor')
            else:
                sp1.tick_params(
                    axis='y', direction='out',
                    length=self.yparameters[self.current_plot-1]['ticklength'])
                sp1.tick_params(axis='y', direction='out', which='minor')
            if hidexticksflag == 1:
                sp1.get_xaxis().set_ticks([])
            if hideyticksflag == 1:
                sp1.get_yaxis().set_ticks([])
            if hidexlabelsflag == 1:
                sp1.get_xaxis().set_ticklabels([])
            if hideylabelsflag == 1:
                sp1.get_yaxis().set_ticklabels([])
            if hidexflag == 0:
                sp1.get_xaxis().set_visible(True)
            else:
                sp1.get_xaxis().set_visible(False)
            if hideyflag == 0:
                sp1.get_yaxis().set_visible(True)
            else:
                sp1.get_yaxis().set_visible(False)
            c1.draw()
            c1.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES)
            h1 = Tk.Frame(hesswindow)
            h1.pack(side=Tk.TOP)
            h1.config(bg=BGCOL)
            button = Tk.Button(h1, text="Save as PS",
                               command=lambda: save_ps_figure(self.p1))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(h1, text="Save as PNG",
                               command=lambda: save_png_figure(self.p1))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(h1, text="Save as FITS", command=self.makeFits)
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(h1, text="Close", command=hesswindow.destroy)
            button.pack()
            button.config(bg=BGCOL)
        except:
            pass

    def makeFits(self):
        """
        Write out the two-dimensional histogram image as a FITS file.

        This routine makes a FITS output file of the two-dimensional histogram
        values in the plot, for use by other codes.

        No values are passed to this routine or returned from it.

        """
        filename = tkinter.filedialog.asksaveasfilename(
            filetypes=[('FITS', '*.fits')])
        st1 = filename.split('.')
        if 'fits' not in st1[-1]:
            filename = filename + '.fits'
        newimage = numpy.transpose(self.hist2d)
        hdu = fits.PrimaryHDU(newimage)
        hdulist = fits.HDUList(hdu)
        primary = hdulist[0].header
        primary['CRPIX1'] = (1.0, 'Axis 1 reference pixel')
        primary['CRPIX2'] = (1.0, 'Axis 2 reference pixel')
        x1 = (self.xedges[1]+self.xedges[0])/2.
        delx = (self.xedges[1]-self.xedges[0])
        y1 = (self.yedges[1]+self.yedges[0])/2.
        dely = (self.yedges[1]-self.yedges[0])
        primary['CRVAL1'] = (x1, 'mean value at reference pixel')
        primary['CRVAL2'] = (y1, 'sigma value at reference pixel')
        primary['CDELT1'] = (delx, 'change in mean value per pixel')
        primary['CDELT2'] = (dely, 'change in sigma value per pixel')
        primary['CTYPE1'] = (' ', 'axis 1 type')
        primary['CTYPE2'] = (' ', 'axis 2 type')
        primary['CUNIT1'] = (' ', 'axis 1 unit')
        primary['CUNIT2'] = (' ', 'axis 2 unit')
        hdulist.writeto(filename, overwrite=True)

    def histogram_position(self, event):
        """
        Apply the cursor position to the histogram plot.

        When a normal histogram plot exists, this routine takes the mouse
        position events and updates the position values at the top of
        the window.

        Parameters
        ----------
            event :   a standard Tkinter event variable.

        Returns
        -------
            No values are returned by this routine.

        """
        try:
            s1 = 'Position: [%f, %f]' % (event.xdata, event.ydata)
            self.histogramLabelText.set(s1)
        except:
            pass

    def hess_position(self, event):
        """
        Apply the cursor position to the field in the 2-D histogram plot.

        When a two-dimensional histogram plot exists, this
        routine takes the mouse position events and updates
        the position values at the top of the window.

        Parameters
        ----------
            event :   a standard Tkinter event variable.

        Returns
        -------
            No values are returned by this routine.

        """
        try:
            ix = bisect.bisect(self.xedges, event.xdata)
            iy = bisect.bisect(self.yedges, event.ydata)
            s1 = 'Position: [%f, %f] Value: %d' % (event.xdata, event.ydata,
                                                   self.hist2d[ix, iy])
            self.hessLabelText.set(s1)
        except:
            pass

    def write_data_sets(self):
        """
        Write the data values to an ascii output file.

        This routine writes the current set values (x, y) out to an ascii
        output file.  If no sets are defined, the routine simply returns.

        Parameters
        ----------
            event :   a standard Tkinter event variable.

        Returns
        -------
            No values are returned by this routine.

        The output is a simple (x, y) ascii dump with a blank line between
        data sets.
        """
        if self.nsets == 0:
            return
        outfilename = tkinter.filedialog.asksaveasfilename()
        outfile = open(outfilename, "w")
        for loop in range(self.nsets):
            xboth = False
            yboth = False
            for n in range(len(self.xdata[loop]['values'])):
                if n == 0:
                    if self.xdata[loop]['errors']:
                        for n in range(len(self.xdata[loop]['lowerror'])):
                            if self.xdata[loop]['lowerror'][n] != \
                               self.xdata[loop]['higherror'][n]:
                                xboth = True
                    if self.ydata[loop]['errors']:
                        for n in range(len(self.ydata[loop]['lowerror'])):
                            if self.ydata[loop]['lowerror'][n] != \
                               self.ydata[loop]['higherror'][n]:
                                yboth = True
                    headerstr = '# Set %d: %s' % (
                        loop+1, self.set_properties[loop]['label'])
                    print(headerstr, file=outfile)
                    headerstr = '# X Value  |'
                    if self.xdata[loop]['errors']:
                        if xboth:
                            headerstr = headerstr \
                                + 'X Error minus   | X Error plus   |'
                        else:
                            headerstr = headerstr + 'X Error   |}'
                    headerstr = headerstr + ' Y Value  |'
                    if self.ydata[loop]['errors']:
                        if yboth:
                            headerstr = headerstr \
                                + 'y Error minus   | X Error plus   |'
                        else:
                            headerstr = headerstr + 'Y Error   |'
                    print(headerstr, file=outfile)
                str1 = ''
                str1 = str1 + self.format(self.xdata[loop]['values'][n])
                if self.xdata[loop]['errors']:
                    str1 = str1 + self.format(self.xdata[loop]['lowerror'][n])
                    if xboth:
                        str1 = str1 + self.format(
                            self.xdata[loop]['higherror'][n])
                str1 = str1 + self.format(self.ydata[loop]['values'][n])
                if self.ydata[loop]['errors']:
                    str1 = str1 + self.format(self.ydata[loop]['lowerror'][n])
                    if yboth:
                        str1 = str1 + self.format(
                            self.ydata[loop]['higherror'][n])
                print(str1, file=outfile)
            print(' ', file=outfile)
        outfile.close()

    def format(self, value):
        """
        Apply a format to a real value for writing out the data values.

        This routine is used to format an input real value either in
        exponential format or in floating point format depending on
        the magnitude of the input value.
        This works better for constant width columns than the Python g format.

        Parameters
        ----------
            value :   a real number value

        Returns
        -------
            outstr :  a format string segment

        """
        if (abs(value) > 1.e+07) or (abs(value) < 1.e-06):
            outstr = '%20.12e ' % (value)
            if value == 0.:
                outstr = '%20.12f ' % (value)
        else:
            outstr = '%20.12f ' % (value)
        return outstr

    def set_statistics(self):
        """
        Calculate set statistics.

        This routine prints some statistics about the different sets to
        a pop-up window.  If no sets are defined, then the routine just
        returns.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.nsets == 0:
            return
        outstr = 'Set  Number of points    range (minimum, maximum)' \
                 + '       mean    standard deviation\n'
        for loop in range(self.nsets):
            npoints = len(self.xdata[loop]['values'])
            xdmin = numpy.min(self.xdata[loop]['values'])
            xdmax = numpy.max(self.xdata[loop]['values'])
            xdmean = numpy.mean(self.xdata[loop]['values'])
            xdsigma = numpy.std(self.xdata[loop]['values'])
            ydmin = numpy.min(self.ydata[loop]['values'])
            ydmax = numpy.max(self.ydata[loop]['values'])
            ydmean = numpy.mean(self.ydata[loop]['values'])
            ydsigma = numpy.std(self.ydata[loop]['values'])
            outstr = outstr + '%3d %10d  x:   %13.6g %13.6g %13.6g %13.6g\n'\
                % (loop+1, npoints, xdmin, xdmax, xdmean, xdsigma) \
                + '                y:   %13.6g %13.6g %13.6g %13.6g\n' \
                % (ydmin, ydmax, ydmean, ydsigma)
        stat_window = Tk.Toplevel()
        stat_window.title('Set Statistics')
        holder = Tk.Frame(stat_window)
        holder.pack(side=Tk.TOP)
        stat_message_text = ScrolledText(holder, height=40, width=90,
                                         wrap=Tk.NONE)
        stat_message_text.config(font=('courier', 16, 'bold'))
        stat_message_text.pack(side=Tk.TOP)
        stat_message_text.insert(0.0, outstr)
        bholder = Tk.Frame(stat_window)
        bholder.pack(side=Tk.TOP)
        close_button = Tk.Button(bholder, text='Close Window',
                                 command=stat_window.destroy)
        close_button.pack()

    def read_data_set(self):
        """
        Make a window for reading in data sets.

        This routine produces a window to read in ascii data sets from a file.
        The window stays until one clicks on the "Close" button.  The window
        is regenerated each time the parent button is clicked.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.read_window is not None:
            return
        self.labelstring = ' '
        self.read_window = Tk.Toplevel()
        self.read_window.title('Read in Data')
        holder = Tk.Frame(self.read_window)
        holder.pack(side=Tk.TOP)
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Data File Filter')
        label.pack(side=Tk.LEFT)
        self.file_template_field = Tk.Entry(field1, width=20)
        self.file_template_field.pack(side=Tk.LEFT)
        self.file_template_field.insert(0, '*')
        field2 = Tk.Frame(holder)
        field2.pack(side=Tk.TOP)
        label = Tk.Label(field2, text='Set Type')
        label.pack(side=Tk.LEFT)
        self.set_options = Tk.StringVar()
        self.set_options.set('XY')
        self.set_option_list = ['XY', 'XYdY', 'XdXY', 'XdXYdY',
                                'XYdYdY', 'XdXdXY', 'XdXdXYdYdY']
        menu1 = Tk.OptionMenu(field2, self.set_options, *self.set_option_list,
                              command=self.error_fields)
        menu1.config(width=10)
        menu1.pack(side=Tk.LEFT)
        field3 = Tk.Frame(holder)
        field3.pack(side=Tk.TOP)
        label = Tk.Label(field3, text='Autoscale')
        label.pack(side=Tk.LEFT)
        self.autoscale_options = Tk.StringVar()
        self.autoscale_options.set('XY')
        self.autoscale_option_list = ['XY', 'X', 'Y', 'None']
        menu1 = Tk.OptionMenu(field3, self.autoscale_options,
                              *self.autoscale_option_list)
        menu1.config(width=10)
        menu1.pack(side=Tk.LEFT)
        holder = Tk.Frame(self.read_window)
        holder.pack(side=Tk.TOP)
        self.read_window_label = Tk.Label(holder, text=self.labelstring)
        self.read_window_label.pack(side=Tk.TOP)
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='X data from column (0 for index): ')
        label.pack(side=Tk.LEFT)
        self.xdata_field = Tk.Entry(field1, width=10)
        self.xdata_field.pack(side=Tk.LEFT)
        self.xdata_field.insert(0, '1')
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Y data from column (0 for index): ')
        label.pack(side=Tk.LEFT)
        self.ydata_field = Tk.Entry(field1, width=10)
        self.ydata_field.pack(side=Tk.LEFT)
        self.ydata_field.insert(0, '2')
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='X uncertainties from column(s): ')
        label.pack(side=Tk.LEFT)
        self.dxdata_field = Tk.Entry(field1, width=10)
        self.dxdata_field.pack(side=Tk.LEFT)
        self.dxdata_field.insert(0, '3')
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Y uncertainties from column(s): ')
        label.pack(side=Tk.LEFT)
        self.dydata_field = Tk.Entry(field1, width=10)
        self.dydata_field.pack(side=Tk.LEFT)
        self.dydata_field.insert(0, '4')
        self.error_fields(None)
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='All other columns as Y: ')
        label.pack(side=Tk.LEFT)
        self.ally = Tk.IntVar()
        Tk.Radiobutton(field1, text='Yes', variable=self.ally, value=1).pack()
        Tk.Radiobutton(field1, text='No', variable=self.ally, value=0).pack()
        self.ally.set(0)
        field1 = Tk.Frame(holder)
        label = Tk.Label(field1, text=' ')
        label.pack(side=Tk.TOP)
        field1.pack(side=Tk.TOP)
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Sort data points: ')
        label.pack(side=Tk.LEFT)
        self.sortdata = Tk.IntVar()
        Tk.Radiobutton(field1, text='Sort on x',
                       variable=self.sortdata, value=2).pack()
        Tk.Radiobutton(field1, text='Sort on y',
                       variable=self.sortdata, value=1).pack()
        Tk.Radiobutton(field1, text='Do not sort',
                       variable=self.sortdata, value=0).pack()
        self.sortdata.set(0)
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        self.labelsstring = ' '
        field4 = Tk.Frame(holder)
        field4.pack(side=Tk.TOP)
        select_button = Tk.Button(field4, text="Select File",
                                  command=self.select_file_window)
        select_button.pack(side=Tk.LEFT)
        label1 = Tk.Label(field4, text="    ")
        label1.pack(side=Tk.LEFT)
        select_button = Tk.Button(
            field4, text="Get Values",
            command=lambda: self.get_set(self.datavalues, self.labelstring))
        select_button.pack(side=Tk.LEFT)
        label2 = Tk.Label(field4, text="    ")
        label2.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            field4, text="Cancel/Close",
            command=lambda: self.close_window(self.read_window))
        close_button.pack(side=Tk.LEFT)

    def error_fields(self, event):
        """
        Set the data error bar values from the window parameters.

        This routine processes events from the set type options menu field.
        When that menu is changed, this routine is called.  The event is
        actually ignored and the menu value used to set the error bar
        properties.

        Parameters
        ----------
            event :   a standard Tkinter event variable (not used)

        Returns
        -------
            No values are returned by the routine

        """
        set_option_string = self.set_options.get()
        if set_option_string is None:
            set_option_string = 'None'
        option1 = self.set_option_list.index(set_option_string)
        self.dxdata_field.config(state=Tk.NORMAL)
        self.dydata_field.config(state=Tk.NORMAL)
        if (option1 < 2) | (option1 == 4):
            self.dxdata_field.config(state=Tk.DISABLED)
        if (option1 == 0) | (option1 == 5):
            self.dydata_field.config(state=Tk.DISABLED)

    def close_window(self, windowname):
        """
        Close a window and set the associated variable to None.

        This routine removes a window from the systen and replaces the
        old window variable with None to mark that the window is not active.

        Parameters
        ----------
            windowname :  A Tk.Toplevel variable for an existing window
                          that is to be closed.

        Returns
        -------
            No value is returned by this routine.

        """
        windowname.destroy()
        # One needs to set the window name to None to signal that a new
        # window of this type can be opened.  For some reason one cannot
        # just set windowname to None here.  If one does, it has no affect
        # on the object variable that is passed into the routine.
        if self.read_window == windowname:
            self.read_window = None
        if self.data_set_window == windowname:
            self.data_set_window = None
        if self.data_set_transformation_window == windowname:
            self.data_set_transformation_window = None
        if self.plot_control_window == windowname:
            self.plot_control_window = None
        if self.plot_control_window == windowname:
            self.plot_control_window = None
        if self.box_window == windowname:
            self.box_window = None
        if self.vector_window == windowname:
            self.vector_window = None
        if self.ellipse_window == windowname:
            self.ellipse_window = None
        if self.line_window == windowname:
            self.line_window = None
        if self.setplot_window == windowname:
            self.setplot_window = None
        if self.tile_window == windowname:
            self.tile_window = None
        if self.hideplot_window == windowname:
            self.hideplot_window = None
        windowname = None

    def close_data_window(self, windowname):
        """
        Remove a data window and set the associated variables to None.

        This routine removes a window from the system and replaces the
        old window variable with None to mark that the window is not active.
        This version also sets the self.datavaulues variable to None, so it
        specific to the data set read window.

        Parameters
        ----------
            windowname :  A Tk.Toplevel variable for an existing window
                          that is to be closed.

        Returns
        -------
            No value is returned by this routine.

        """
        self.datavalues = None
        windowname.destroy()
        # One needs to set the window name to None to signal that a new
        # window of this type can be opened.  For some reason one cannot
        # just set windowname to None here.  If one does, it has no affect
        # on the object variable that is passed into the routine.
        if self.read_window == windowname:
            self.read_window = None
        if self.data_set_window == windowname:
            self.data_set_window = None
        if self.data_set_transformation_window == windowname:
            self.data_set_trasnformation_window = None
        if self.data_set_fitting_window == windowname:
            self.data_set_fitting_window = None
        if self.data_set_delete_window == windowname:
            self.data_set_delete_window = None
        if self.data_entry_window == windowname:
            self.data_entry_window = None
        windowname = None

    def select_file_window(self):
        """
        Select a file and read data columns from it for input.

        This window queries the user for an input file name and tries to
        read the data therein.  If this is successful, some of the fields
        in the "Read Data" window are updated.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            pattern = self.file_template_field.get()
            values = pattern.split(', ')
            filetypes = []
            for loop in range(len(values)):
                label1 = values[loop].replace('*.', '')
                label2 = values[loop].replace('*', '')
                label1 = label1.upper()
                template1 = (label1, label2)
                filetypes.append(template1)
        except:
            tkinter.messagebox.showinfo(
                'Error', 'There was '
                + 'some error reading the file pattern.  Using the default.')
            pattern = '*'
        if (pattern == '*') | (pattern == ''):
            self.filename = tkinter.filedialog.askopenfilename()
        else:
            try:
                self.filename = tkinter.filedialog.askopenfilename(
                    filetypes=filetypes)
            except:
                self.filename = tkinter.filedialog.askopenfilename()
        try:
            self.shortname = None
            inds, ncols = self.check_columns(self.filename)
            if inds is None:
                tkinter.messagebox.showinfo(
                    'Error', 'There was '
                    + 'some error reading file %s' % (self.filename))
            else:
                self.datavalues = numpy.loadtxt(
                    self.filename,
                    comments=['#', '\\', '|', 'x_or_RA'], usecols=inds)
                values = self.filename.split('/')
                self.shortname = values[-1]
                self.datashape = self.datavalues.shape
                if len(self.datashape) > 1:
                    newlabel = '%d values ' % (self.datashape[0])
                    newlabel = newlabel + '(%s) ' % (self.shortname)
                    self.labelstring = newlabel
                    self.ncolumns = self.datashape[1]
                else:
                    newlabel = '1 value '
                    newlabel = newlabel + '(%s) ' % (self.shortname)
                    self.labelstring = '(%s) ' % (self.shortname)
                    self.ncolumns = self.datashape[0]
                labeltext = '%d data columns, %d values ' % (
                    self.datashape[1], self.datashape[0])
                self.read_window_label['text'] = labeltext
                self.data_indexes = inds
                self.ndatacolumns = ncols
        except:
            try:
                if self.shortname is None:
                    tkinter.messagebox.showinfo(
                        'Error',
                        'There was some error reading file %s' % (
                            self.filename))
                else:
                    tkinter.messagebox.showinfo(
                        'Error',
                        'There was some error reading file %s' % (
                            self.shortname))
            except:
                tkinter.messagebox.showinfo(
                    'Error',
                    'There was some error trying to read a file')

    def check_columns(self, filename):
        """
        Check the columns in an input file for numerical values.

        This routine takes as input a file name.  It reads in the lines
        from the file and determines the number of columns plus the index
        values of the numerical columns.  It is assumed that the first
        column without '#', '|' or '\' is the template for all the lines, so
        only the first "data" line is split up and checked for
        numerical values in the columns.  The indices are then
        used with numpy.loadtxt.

        Parameters
        ----------
            filename :   the name of the file with regular columns to check

        Returns
        -------
            inds :       a list of the numerical columns in the file

            ncols :      the total number of columns in the file

        If an error occurs, inds and ncols are returned as None.

        """
        try:
            infile = open(filename, 'r')
            lines = infile.readlines()
            infile.close()
            ncols = 0
            inds = []
            for line in lines:
                if ncols == 0:
                    line = line.strip('\n')
                    if ('#' in line[0:1]) | ('\\' in line[0:1]) | \
                       ('|' in line[0:1]) | ('x_or_RA' in line):
                        pass
                    else:
                        values = line.split('#')
                        subline = values[0]
                        values = subline.split('|')
                        subline = values[0]
                        values = subline.split('\\')
                        subline = values[0]
                        values = subline.split()
                        ncols = len(values)
                        for loop in range(ncols):
                            try:
                                x = float(values[loop])
                                inds.append(loop)
                            except:
                                pass
            return inds, ncols
        except:
            return None, None

    def get_set(self, datavalues, labelstring):
        """
        Extract data values to a new data set.

        This routine tries to make a data set (x, y) pair according to the
        specified values from the read sets window.

        Parameters
        ----------
           datavalues :  This is a numpy array from the "loadtxt" function.

           labelstring : This is a string variable, used to set the labels
                          for the set(s)

        Returns
        -------
           No value is returned.

        """
        if datavalues is None:
            tkinter.messagebox.showinfo(
                'Error',
                'You must specify a file before reading data.')
            return
        if len(datavalues.shape) == 1:
            datavalues = numpy.expand_dims(datavalues, 0)
        try:
            set_option_string = self.set_options.get()
            if set_option_string is None:
                set_option_string = 'None'
            option1 = self.set_option_list.index(set_option_string)
            autoscale_option_string = self.autoscale_options.get()
            option2 = self.autoscale_option_list.index(
                autoscale_option_string)
            option3 = self.ally.get()
            option4 = self.sortdata.get()
            shape1 = datavalues.shape
            xerrorinds = [2, 3, 5, 6]
            yerrorinds = [1, 3, 4, 6]
            if option3 == 0:
                try:
                    inds = numpy.zeros((6), dtype=numpy.int8) - 1
                    str1 = self.xdata_field.get()
                    n1 = int(str1) - 1
                    inds[0] = n1
                    str1 = self.ydata_field.get()
                    n1 = int(str1) - 1
                    inds[1] = n1
                    if option1 in xerrorinds:
                        str1 = self.dxdata_field.get()
                        values = str1.split(', ')
                        n1 = int(values[0]) - 1
                        inds[2] = n1
                        if len(values) > 1:
                            n1 = int(values[1])-1
                            inds[3] = n1
                        else:
                            inds[3] = -1
                    if option1 in yerrorinds:
                        str1 = self.dydata_field.get()
                        values = str1.split(', ')
                        n1 = int(values[0]) - 1
                        inds[4] = n1
                        if len(values) > 1:
                            n1 = int(values[1])-1
                            inds[5] = n1
                        else:
                            inds[5] = -1
                except:
                    raise ValueError
                newinds = inds*0 - 1
                for loop in range(len(inds)):
                    for n in range(len(self.data_indexes)):
                        if inds[loop] == self.data_indexes[n]:
                            newinds[loop] = n
                if max(newinds) > self.ndatacolumns:
                    tkinter.messagebox.showinfo(
                        'Error',
                        'The requested columns do not match the shape'
                        + ' of the input data values.')
                    return
                if newinds[0] >= 0:
                    xvalues = numpy.squeeze(datavalues[:, newinds[0]])
                else:
                    xvalues = numpy.arange(1, shape1[0]+1)
                xerrorflag = False
                yerrorflag = False
                if option1 in xerrorinds:
                    xlowerror = numpy.squeeze(datavalues[:, newinds[2]])
                    if newinds[3] > 0:
                        xhigherror = numpy.squeeze(
                            datavalues[:, newinds[3]])
                    else:
                        xhigherror = numpy.squeeze(
                            datavalues[:, newinds[2]])
                    xerrorflag = True
                else:
                    xlowerror = xvalues * 0.
                    xhigherror = xvalues * 0.
                if newinds[1] >= 0:
                    yvalues = numpy.squeeze(datavalues[:, newinds[1]])
                else:
                    yvalues = numpy.arange(1, shape1[0]+1)
                if option1 in yerrorinds:
                    ylowerror = numpy.squeeze(datavalues[:, newinds[4]])
                    if newinds[5] > 0:
                        yhigherror = numpy.squeeze(
                            datavalues[:, newinds[5]])
                    else:
                        yhigherror = numpy.squeeze(
                            datavalues[:, newinds[4]])
                    yerrorflag = True
                else:
                    ylowerror = yvalues * 0.
                    yhigherror = yvalues * 0.
                flags = numpy.logical_not(numpy.isnan(xvalues))
                xvalues = xvalues[flags]
                xlowerror = xlowerror[flags]
                xhigherror = xhigherror[flags]
                yvalues = yvalues[flags]
                ylowerror = ylowerror[flags]
                yhigherror = yhigherror[flags]
                flags = numpy.logical_not(numpy.isnan(yvalues))
                xvalues = xvalues[flags]
                xlowerror = xlowerror[flags]
                xhigherror = xhigherror[flags]
                yvalues = yvalues[flags]
                ylowerror = ylowerror[flags]
                yhigherror = yhigherror[flags]
                if option4 > 0:
                    if option4 == 2:
                        inds = numpy.argsort(xvalues)
                    else:
                        inds = numpy.argsort(yvalues)
                    if newinds[0] >= 0:
                        xvalues = xvalues[inds]
                        xlowerror = xlowerror[inds]
                        xhigherror = xhigherror[inds]
                    else:
                        xlowerror = xvalues * 0.
                        xhigherror = xvalues * 0.
                    if newinds[1] >= 0:
                        yvalues = yvalues[inds]
                        ylowerror = ylowerror[inds]
                        xhigherror = yhigherror[inds]
                    else:
                        ylowerror = yvalues * 0.
                        yhigherror = yvalues * 0.
                newlabelstring = 'columns %d/%d %s' % (
                    newinds[0]+1, newinds[1]+1, labelstring)
                self.add_set(xvalues, yvalues, xlowerror, xhigherror,
                             ylowerror, yhigherror, xerrorflag, yerrorflag,
                             option2, newlabelstring, self.current_plot)
            else:
                str1 = self.xdata_field.get()
                n1 = int(str1) - 1
                for loop in range(len(self.data_indexes)):
                    if n1 == self.data_indexes[loop]:
                        xind = loop
                xvalues = numpy.squeeze(datavalues[:, xind])
                xlowerror = xvalues * 0.
                xhigherror = xvalues * 0.
                xerrorflag = False
                yerrorflag = False
                for loop in range(len(self.data_indexes)):
                    if loop == n1:
                        pass
                    else:
                        yvalues = numpy.squeeze(datavalues[:, loop])
                        ylowerror = yvalues * 0.
                        yhigherror = yvalues * 0.
                        newlabelstring = 'columns %d and %d from file %s' % (
                            n1+1, loop+1, labelstring)
                        self.add_set(xvalues, yvalues, xlowerror, xhigherror,
                                     ylowerror, yhigherror, xerrorflag,
                                     yerrorflag, option2, labelstring,
                                     self.current_plot)
            self.filename = None
            self.make_plot()
            return
        except:
            self.datavalues = None
            self.filename = None
            tkinter.messagebox.showinfo(
                'Error',
                'There was some error in finding the data values.')
            return

    def add_set(self, xvalues, yvalues, xlowerror=None, xhigherror=None,
                ylowerror=None, yhigherror=None, xerrorflag=False,
                yerrorflag=False, rangeoption=0, labelstring=" ",
                current_plot=1):
        """
        Add a data set to the plot object.

        Given some input arrays this routine loads the values into the
        internal data structure for plotting.

        At minimum one needs--

        Parameters
        ----------
            xvalues :  a numpy (floating point or integer) array of the x
                      values for the set

            yvalues :  a numpy (floating point or integer) array of the y
                       values for the set

            The remaining parameters are optional:

            xlowerror :   a numpy array of the error values on the low side
                          for x, defaults to 0.0

            xhigherror :  a numpy array of the error values on the high side
                          for x, defaults to 0.0

            ylowerror :   a numpy array of the error values on the low side
                          for y, defaults to 0.0

            yhigherror :  a numpy array of the error values on the high side
                          for y, defaults to 0.0

            xerrorflag :  a boolean value for whether the x error bars are
                          present, defaults to False

            yerrorflag :  a boolean value for whether the y error bars are
                          present, defaults to False

            rangeoption :  an integer flag for applying the set value to
                           present range

                           0: apply x and y range (default)
                           1: apply x range only
                           2: apply y range only

            labelstring :  a label to apply to the set, defaults to " "

            current_plot :  the plot number where the set is used/displayed
                            (default 1)

        Returns
        -------
            No values are returned by this routine.

        The minimumn call to the routine is to pass in a x data array and
        a y data array.  The error arrays and the flags are set if no
        values are assigned.

        Note: the code does not check whether the data arrays are of
        consistent lengths.

        """
        if xlowerror is None:
            xlowerror = xvalues * 0.
        if ylowerror is None:
            ylowerror = yvalues * 0.
        if xhigherror is None:
            xhigherror = xvalues * 0.
        if yhigherror is None:
            yhigherror = yvalues * 0.
        xmin, xmax, ymin, ymax = self.get_range(
            xvalues, yvalues, xlowerror, xhigherror, ylowerror, yhigherror,
            rangeoption)
        self.xdata[self.nsets] = {'values': xvalues, 'lowerror': xlowerror,
                                  'higherror': xhigherror, 'minimum': xmin,
                                  'maximum': xmax, 'errors': xerrorflag,
                                  'legend': True}
        self.ydata[self.nsets] = {'values': yvalues, 'lowerror': ylowerror,
                                  'higherror': yhigherror, 'minimum': ymin,
                                  'maximum': ymax, 'errors': yerrorflag,
                                  'legend': True}
        self.original_range[self.nsets] = True
        if (current_plot > self.number_of_plots) | (current_plot < 1):
            self.set_properties[self.nsets]['plot'] = 1
        else:
            self.set_properties[self.nsets]['plot'] = current_plot
        m = self.nsets % 10
        n = int(math.floor(self.nsets / 10))
        self.set_properties[self.nsets]['symbol'] = self.markerset[m]
        self.set_properties[self.nsets]['colour'] = self.colourset[m]
        self.set_properties[self.nsets]['symbolsize'] = 4.0 + 0.3*n
        self.set_properties[self.nsets]['label'] = labelstring
        if xerrorflag | yerrorflag:
            self.set_properties[self.nsets]['errors'] = True
        else:
            self.set_properties[self.nsets]['errors'] = False
        self.nsets = self.nsets + 1

    def get_range(self, xvalues, yvalues, xlowerror, xhigherror, ylowerror,
                  yhigherror, rangeoption):
        """
        Find the data range and round off the values.

        Given a set of values, this routine finds the data range and
        rounds the values off somewhat to make nice plot limits.

        Parameters
        ----------
            xvalues :  a numpy array of x data values (float or integer
                       values)

            yvalues :  a numpy array of y data values (float or integer
                       values)

            xlowerror : a numpy array of x error values (going negative)

            xhigherror :  a numpy array of x error values (going posative)

            ylowerror :  a numpy array of y error values (going negative)

            yhigherror : a numpy array of y error values (going posative)

            rangeoption : an integer flag for the values to update (0 for
                          both x and y, 1 for x only, 2 for y only)

        Returns
        -------
            xmin :    the minimum x value, a floating point number

            xmax :    the maximum x value, a floating point number

            ymin :    the minimum y value, a floating point number

            ymax :    the maximum y value, a floating point number

        Note that all four values are returned even if the update of the
        range parameters only affects x or y.

        """
        xmin = numpy.min(xvalues - xlowerror)
        xmax = numpy.max(xvalues + xhigherror)
        ymin = numpy.min(yvalues - ylowerror)
        ymax = numpy.max(yvalues + yhigherror)
        xmin1 = self.round_float(xmin, True)
        xmax1 = self.round_float(xmax, False)
        ymin1 = self.round_float(ymin, True)
        ymax1 = self.round_float(ymax, False)
        if self.nsets == 0:
            if rangeoption == 0:
                self.plot_range[self.current_plot-1][0] = xmin1
                self.plot_range[self.current_plot-1][1] = xmax1
                self.plot_range[self.current_plot-1][2] = ymin1
                self.plot_range[self.current_plot-1][3] = ymax1
            elif rangeoption == 1:
                self.plot_range[self.current_plot-1][0] = xmin1
                self.plot_range[self.current_plot-1][1] = xmax1
            elif rangeoption == 2:
                self.plot_range[self.current_plot-1][2] = ymin1
                self.plot_range[self.current_plot-1][3] = ymax1
            else:
                pass
        else:
            if rangeoption == 0:
                self.plot_range[self.current_plot-1][0] = min(
                    self.plot_range[self.current_plot-1][0], xmin1)
                self.plot_range[self.current_plot-1][1] = max(
                    self.plot_range[self.current_plot-1][1], xmax1)
                self.plot_range[self.current_plot-1][2] = min(
                    self.plot_range[self.current_plot-1][2], ymin1)
                self.plot_range[self.current_plot-1][3] = max(
                    self.plot_range[self.current_plot-1][3], ymax1)
            elif rangeoption == 1:
                self.plot_range[self.current_plot-1][0] = min(
                    self.plot_range[self.current_plot-1][0], xmin1)
                self.plot_range[self.current_plot-1][1] = max(
                    self.plot_range[self.current_plot-1][1], xmax1)
            elif rangeoption == 2:
                self.plot_range[self.current_plot-1][2] = min(
                    self.plot_range[self.current_plot-1][2], ymin1)
                self.plot_range[self.current_plot-1][3] = max(
                    self.plot_range[self.current_plot-1][3], ymax1)
        return xmin, xmax, ymin, ymax

    def make_data_set_edit_window(self):
        """
        Create a new text window in which to edit a data set.

        Returns
        -------
        None.

        """
        if self.nsets == 0:
            return
        if self.nsets > 1:
            str1 = "Which set do you want to edit (1 to %d)?" % self.nsets
            nset = tkinter.simpledialog.askinteger(
                "Input", str1, parent=self.root)
            if (nset is None) or (nset < 1) or (nset > self.nsets):
                return
        else:
            nset = self.nsets
        if self.data_entry_window is not None:
            return
        self.data_entry_window = Tk.Toplevel()
        self.data_entry_window.title('Edit Data Set Values')
        holder = Tk.Frame(self.data_entry_window)
        holder.pack(side=Tk.TOP)
        holder.config(bg='black')
        self.data_text = ScrolledText(holder, height=40, width=80,
                                      wrap=Tk.NONE, relief="solid")
        self.data_text.config(font=('courier', 16))
        self.data_text.pack(side=Tk.TOP, padx=10, pady=10)
        str1 = ''
        for loop in range(len(self.xdata[nset-1]['values'])):
            str1 = str1 + '%g %g %g %g %g %g\n' % (
                self.xdata[nset-1]['values'][loop],
                self.xdata[nset-1]['lowerror'][loop],
                self.xdata[nset-1]['higherror'][loop],
                self.ydata[nset-1]['values'][loop],
                self.ydata[nset-1]['lowerror'][loop],
                self.ydata[nset-1]['higherror'][loop]
            )
        self.data_text.insert(1.0, str1)
        bframe = Tk.Frame(self.data_entry_window)
        bframe.pack()
        set_button = Tk.Button(
            bframe, text="Apply",
            command=lambda: self.apply_data_edits(nset, self.data_text))
        set_button.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            bframe, text="Close",
            command=lambda: self.close_data_window(
                self.data_entry_window))
        close_button.pack(side=Tk.LEFT)

    def apply_data_edits(self, nset, data_text):
        text = data_text.get("1.0", Tk.END)
        xvalues, dxvalues1, dxvalues2, yvalues, dyvalues1, dyvalues2,\
            errorflag = parse_text(text)
        try:
            xvalues = numpy.asarray(xvalues)
            yvalues = numpy.asarray(yvalues)
            dxvalues1 = numpy.asarray(dxvalues1)
            dxvalues2 = numpy.asarray(dxvalues2)
            dyvalues1 = numpy.asarray(dyvalues1)
            dyvalues2 = numpy.asarray(dyvalues2)
            if len(xvalues) < 1:
                tkinter.messagebox.showinfo(
                    'error',
                    'Unable to parse text from edit widget (2)')
                return
            self.xdata[nset-1]['values'] = xvalues
            self.xdata[nset-1]['lowerror'] = dxvalues1
            self.xdata[nset-1]['higherror'] = dxvalues2
            self.ydata[nset-1]['values'] = yvalues
            self.ydata[nset-1]['lowerror'] = dyvalues1
            self.ydata[nset-1]['higherror'] = dyvalues2
            self.make_plot()
        except:
            tkinter.messagebox.showinfo(
                'error',
                'Unable to parse text from edit widget')
            return

    def make_data_set_sort_window(self):
        """
        Create the data set sorting window.

        This routine brings up a window within which one can sort a data set.
        The sorted data values replace the current data in the set.

        No values are passed to this routine or returned from the routine.
        """
        if self.data_set_sort_window is not None:
            return
        if self.nsets == 0:
            return
        self.data_set_sort_window = Tk.Toplevel()
        self.data_set_sort_window.title('Data Sets Sort Window')
        holder = Tk.Frame(self.data_set_sort_window)
        holder.pack(side=Tk.TOP)
        self.set_sort_list_area = tkinter.ttk.Combobox(holder, width=50,
                                                       height=10)
        self.set_sort_list_area.pack(side=Tk.TOP)
        setlist = []
        for loop in range(self.nsets):
            label = (
                'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                (loop+1, len(self.xdata[loop]['values']),
                 self.xdata[loop]['minimum'], self.xdata[loop]['maximum']))
            label = (
                label + 'y range %13.6g to %13.6g' %
                (self.ydata[loop]['minimum'], self.ydata[loop]['maximum'])
            )
            setlist.append(label)
        self.set_sort_list_area['values'] = setlist
        self.set_sort_list_area.current(0)
        frame1 = Tk.Frame(holder)
        frame1.pack()
        self.sortflag1 = Tk.IntVar()
        self.sortflag2 = Tk.IntVar()
        lb = Tk.Label(frame1, text="Sort on: ")
        lb.pack(side=Tk.LEFT)
        b1 = Tk.Frame(frame1)
        b1.pack()
        self.put_yes_no(b1, self.sortflag1, ['x values', 'y values'], True)
        lb = Tk.Label(frame1, text="Sort order: ")
        lb.pack(side=Tk.LEFT)
        b1 = Tk.Frame(frame1)
        b1.pack()
        self.put_yes_no(b1, self.sortflag2, ['Ascending', 'Decending'], True)
        frame2 = Tk.Frame(holder)
        frame2.pack(side=Tk.TOP)
        sort_button = Tk.Button(frame2, text="Sort", command=self.sort_set)
        sort_button.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            frame2, text="Close",
            command=lambda: self.close_data_window(self.data_set_sort_window))
        close_button.pack(side=Tk.LEFT)

    def sort_set(self):
        """
        Sort the set values (x, y) according to the options.

        This routine sorts a data set either on the x values or on the y
        values.  Sorting can be done in ascending or decending order.

        No values are passed to this routine or returned from this routine.
        """
        try:
            set_number = self.set_sort_list_area.current()
            # The following should not happen....but check just to be sure.
            if (set_number < 0) or (set_number >= self.nsets):
                tkinter.messagebox.showinfo('Error', 'The set was not sorted.')
                return
            option1 = self.sortflag1.get()
            option2 = self.sortflag2.get()
            flag = self.do_sort_set(set_number, option1, option2)
            if not flag:
                tkinter.messagebox.showinfo('Error', 'The set was not sorted.')
        except:
            tkinter.messagebox.showinfo('Error', 'The set was not sorted.')

    def do_sort_set(self, set_number, xyoption=1, direction_option=0):
        """
        Do the actual sorting of the set values.

        This routine carries out the sorting of values in a set.  It is
        provided sepaerately so it can be used from the Python command
        line if needed.

        Parameters
        ----------
            set_number : An integer value 0 or higher for the set to be sorted

            xyoption :   An optional integer value, 1 for sorting in x and
                         0 for sorting in y (defaults to soeting in x)

            direction_option :   An optional integer value, 1 for sorting in
                                 ascending order and 0 for sorting in
                                 decending order (defaults to ascending order)

        Values of 1 are the default for both options.  Actually any xyoption
        value but 1 will cause sorting in y and any direction_option value
        but zero will cause sorting in ascending order.

        Returns
        -------
            flag :    A boolean value, True if the sorting was successful,
                      False otherwise

        """
        try:
            xvalues = numpy.copy(self.xdata[set_number]['values'])
            xlowerror = numpy.copy(self.xdata[set_number]['lowerror'])
            xhigherror = numpy.copy(self.xdata[set_number]['higherror'])
            yvalues = numpy.copy(self.ydata[set_number]['values'])
            ylowerror = numpy.copy(self.ydata[set_number]['lowerror'])
            yhigherror = numpy.copy(self.ydata[set_number]['higherror'])
            if xyoption == 1:
                inds1 = numpy.argsort(xvalues)
            else:
                inds1 = numpy.argsort(yvalues)
            xvalues = xvalues[inds1]
            xlowerror = xlowerror[inds1]
            xhigherror = xhigherror[inds1]
            yvalues = yvalues[inds1]
            ylowerror = ylowerror[inds1]
            yhigherror = yhigherror[inds1]
            if direction_option == 0:
                xvalues = numpy.flip(xvalues)
                xlowerror = numpy.flip(xlowerror)
                xhigherror = numpy.flip(xhigherror)
                yvalues = numpy.flip(yvalues)
                ylowerror = numpy.flip(ylowerror)
                yhigherror = numpy.flip(yhigherror)
            self.xdata[set_number]['values'] = xvalues
            self.xdata[set_number]['lowerror'] = xlowerror
            self.xdata[set_number]['higherror'] = xhigherror
            self.ydata[set_number]['values'] = yvalues
            self.ydata[set_number]['lowerror'] = ylowerror
            self.ydata[set_number]['higherror'] = yhigherror
            self.make_plot()
            return True
        except:
            return False

    def make_data_set_fitting_window(self):
        """
        Create the window for data set fitting.

        This routine creates a window for the data set fitting functions.
        One can set the type of fit (i.e. polynomial, legendre, laguerre)
        and the order of the fit.  When a fit is made, the new fit is added
        as a set in the list.  One can clear the fits with the "Cancel Fits"
        button.  Otherwise once one closes the window any fits that have
        been made are retained in the list of data sets.

        If there are no sets to fit, the routine returns without doing
        anything.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.data_set_fitting_window is not None:
            return
        if self.nsets == 0:
            return
        self.data_set_fitting_window = Tk.Toplevel()
        self.data_set_fitting_window.title('Data Sets Fitting Window')
        holder = Tk.Frame(self.data_set_fitting_window)
        holder.pack(side=Tk.TOP)
        self.fit_text = ScrolledText(holder, height=6, width=50, bd=1,
                                     relief=Tk.RIDGE, wrap=Tk.NONE)
        self.fit_text.config(font=('courier', 12, 'bold'))
        self.fit_text.pack()
        self.set_fitting_list_area = tkinter.ttk.Combobox(holder, width=50,
                                                          height=10)
        self.set_fitting_list_area.pack(side=Tk.TOP)
        setlist = []
        for loop in range(self.nsets):
            label = (
                'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                (loop+1, len(self.xdata[loop]['values']),
                 self.xdata[loop]['minimum'], self.xdata[loop]['maximum'])
                + 'y range %13.6g to %13.6g' %
                (self.ydata[loop]['minimum'], self.ydata[loop]['maximum']))
            setlist.append(label)
        self.set_fitting_list_area['values'] = setlist
        self.set_fitting_list_area.current(0)
        nsets_start = self.nsets
        frame1 = Tk.Frame(holder)
        frame1.pack()
        labels = ['Function', 'Order/Smoothing']
        for loop in range(len(labels)):
            label = Tk.Label(frame1, text=labels[loop])
            label.grid(column=0, row=loop, sticky=Tk.W)
        self.set_fitting_fields = []
        self.set_fitting_fields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.set_fitting_fields[-1].grid(column=1, row=0, sticky=Tk.W)
        self.set_fitting_fields[-1]['values'] = [
            'Polynomial', 'Legendre',
            'Chebyshev', 'Laguerre', 'Spline', 'Cubic Spline',
            'Internal Linear Fit']
        self.set_fitting_fields[-1].current(1)
        self.set_fitting_fields.append(Tk.Entry(frame1, width=15, text='4'))
        self.set_fitting_fields[-1].grid(column=1, row=1, sticky=Tk.W)
        self.set_fitting_fields[-1].delete(0, Tk.END)
        self.set_fitting_fields[-1].insert(0, '4')
        frame2 = Tk.Frame(holder)
        frame2.pack()
        self.fit_option = Tk.IntVar()
        label = Tk.Label(frame2, text='Fit Type: ')
        label.pack(side=Tk.LEFT)
        Tk.Radiobutton(
            frame2, text='y = f(x)',
            variable=self.fit_option, value=0).pack(side=Tk.LEFT)
        Tk.Radiobutton(
            frame2, text='x = f(y)',
            variable=self.fit_option, value=1).pack(side=Tk.LEFT)
        self.fit_option.set(0)
        frame2 = Tk.Frame(holder)
        frame2.pack()
        cancel_button = Tk.Button(
            frame2, text="Cancel Fittings",
            command=lambda: self.cancel_fitting_fields(nsets_start))
        cancel_button.pack(side=Tk.LEFT)
        set_button = Tk.Button(
            frame2, text="Fit Set",
            command=self.apply_fitting_fields)
        set_button.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            frame2, text="Close",
            command=lambda: self.close_data_window(
                self.data_set_fitting_window))
        close_button.pack(side=Tk.LEFT)

    def cancel_fitting_fields(self, nsets):
        """
        Clear the set fitting results.

        This routine will "clear" the fitting data sets.  This is done by
        reverting the self.nsets value to what it was before the fitting
        was started.  Hence all the fitting data sets are still present but
        are not plotted, and if more data sets are read in the fitting sets
        will be overwritten.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        self.nsets = nsets
        self.make_plot()

    def apply_fitting_fields(self):
        """
        Apply the set fitting parameters from the window.

        This routine reads the fitting parameters from the window and then
        applies the fitting.  Each time it is called a new data set should
        be generated, unless the fit order is too large for the number of
        points.  The fit order must be at most 1 less than the number
        of points, otherwise the routine just shows an error message
        pop-up and returns.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        fit_type = self.set_fitting_fields[0].get()
        if 'Cubic Spline' in fit_type:
            fit_order = float(self.set_fitting_fields[1].get())
        else:
            try:
                fit_order = int(self.set_fitting_fields[1].get())
            except:
                str1 = 'Error: bad fit order (%s).  Settng to 4.' % (
                    self.set_fitting_fields[1].get())
                self.fit_text.insert(Tk.END, str1)
                self.fit_text.see(Tk.END)
                fit_order = 4
        set_number = self.set_fitting_list_area.current()
        fit_flag = self.fit_option.get()
        if fit_flag == 0:
            xvalues = numpy.copy(self.xdata[set_number]['values'])
            yvalues = numpy.copy(self.ydata[set_number]['values'])
        else:
            xvalues = numpy.copy(self.ydata[set_number]['values'])
            yvalues = numpy.copy(self.xdata[set_number]['values'])
        inds = numpy.argsort(xvalues)
        xvalues = xvalues[inds]
        yvalues = yvalues[inds]
        npoints = len(xvalues)
        if ('Cubic' not in fit_type) and ('Internal' not in fit_type):
            if npoints + 1 <= fit_order:
                tkinter.messagebox.showinfo(
                    'Error',
                    'The number of points is too few for the fit order.'
                    + '  Please check your inputs.')
                return
        xmin = numpy.min(xvalues)
        xmax = numpy.max(xvalues)
        delx = xmax - xmin
        xstep = ((xmax - xmin) + 2 * delx)/1001.
        xout = numpy.arange(xmin-delx, xmax+delx, xstep)
        if 'Internal' in fit_type:
            if npoints < 2:
                tkinter.messagebox.showinfo(
                    'Error',
                    'The number of points is too few for a linear fit.'
                    + '  Please check your inputs.')
                return
            if fit_flag == 0:
                yerrors = (self.ydata[set_number]['lowerror'] +
                           self.ydata[set_number]['higherror'])/2.
            else:
                yerrors = (self.xdata[set_number]['lowerror'] +
                           self.xdata[set_number]['higherror'])/2.
            if (numpy.min(yerrors) == 0.) and (numpy.max(yerrors) == 0.):
                yerrors = yerrors + 1.
            slope, intercept, slope_error, intercept_error, covariance, \
                correlation = slope_calculation(xvalues, yvalues, yerrors)
            if slope is None:
                tkinter.messagebox.showinfo(
                    'Error',
                    'Error in the standard slope fit.\n')
                return
            yfit = intercept + xvalues*slope
            yout = intercept + xout*slope
            errorterm1 = xout*0. + intercept_error
            errorterm2 = xout*slope_error
            youterror = numpy.sqrt(errorterm1*errorterm1 +
                                   errorterm2*errorterm2)
            labelstring = 'Standard linear fit'
            rms = numpy.sqrt(numpy.mean((yvalues-yfit)*(yvalues-yfit)))
            str1 = 'Regression calculation results:\n'
            str1 = str1 + 'Slope: %f +/- %f\n' % (slope, slope_error)
            str1 = str1 + 'Intercept: %f +/- %f\n' % (
                intercept, intercept_error)
            str1 = str1 + 'Covariance: %f\n' % (covariance)
            str1 = str1 + 'Correlation: %f\n' % (correlation)
            str1 = str1 + 'RMS deviation: %f\n' % (rms)
            tkinter.messagebox.showinfo('Information', str1)
            outfile = open('fit_values.txt', 'a')
            print(str1, file=outfile)
            print(' ', file=outfile)
            outfile.close()
        if fit_type == 'Polynomial':
            fitpars = polynomial.polyfit(xvalues, yvalues, fit_order)
            yout = polynomial.polyval(xout, fitpars)
            yfit = polynomial.polyval(xvalues, fitpars)
            labelstring = 'Order %d polynomial fit' % (fit_order)
            list_fitpars(fit_type, fit_order, fitpars)
        if fit_type == 'Legendre':
            fitpars = legendre.legfit(xvalues, yvalues, fit_order)
            yout = legendre.legval(xout, fitpars)
            yfit = legendre.legval(xvalues, fitpars)
            labelstring = 'Order %d Legendre polynomial fit' % (fit_order)
            list_fitpars(fit_type, fit_order, fitpars)
        if fit_type == 'Laguerre':
            fitpars = laguerre.lagfit(xvalues, yvalues, fit_order)
            yout = laguerre.lagval(xout, fitpars)
            yfit = laguerre.lagval(xvalues, fitpars)
            labelstring = 'Order %d Laguerre polynomial fit' % (fit_order)
            list_fitpars(fit_type, fit_order, fitpars)
        if fit_type == 'Chebyshev':
            fitpars = chebyshev.chebfit(xvalues, yvalues, fit_order)
            yout = chebyshev.chebval(xout, fitpars)
            yfit = chebyshev.chebval(xvalues, fitpars)
            labelstring = 'Order %d Chebyshev polynomial fit' % (fit_order)
            list_fitpars(fit_type, fit_order, fitpars)
        if fit_type == 'Spline':
            fitpars = UnivariateSpline(
                xvalues, yvalues, k=3, s=None, bbox=[xmin - delx,
                                                     xmax + delx])
            labelstring = 'Default spline fit'
            yout = fitpars(xout)
            yfit = fitpars(xvalues)
        if fit_type == 'Cubic Spline':
            if fit_order < 0.:
                str1 = 'Error: smoothing value %f (< 0) is not allowed.'\
                       + '  Settng to 0.0' % (fit_order)
                self.fit_text.insert(Tk.END, str1)
                self.fit_text.see(Tk.END)
                fit_order = 0.0
            fitpars = UnivariateSpline(
                xvalues, yvalues, k=3, bbox=[xmin - delx,
                                             xmax + delx], s=fit_order)
            yout = fitpars(xout)
            yfit = fitpars(xvalues)
            labelstring = 'Cubic spline fit, smoothing = %f' % (fit_order)
        rms = self.fit_statistics(yvalues, yfit)
        if 'Internal' in fit_type:
            if fit_flag == 0:
                xlowerror = youterror * 0.
                xhigherror = youterror * 0.
                ylowerror = youterror
                yhigherror = youterror
            else:
                xlowerror = youterror
                xhigherror = youterror
                ylowerror = youterror * 0.
                yhigherror = youterror * 0.
        else:
            xlowerror = xvalues * 0.
            xhigherror = xvalues * 0.
            ylowerror = yvalues * 0.
            yhigherror = yvalues * 0.
        xmin = numpy.min(xout)
        xmax = numpy.max(xout)
        ymin = numpy.min(yfit)
        ymax = numpy.max(yfit)
        if rms is not None:
            str1 = 'Fit: RMS = %g for %d points\n' % (rms, len(yfit))
            self.fit_text.insert(Tk.END, str1)
            self.fit_text.see(Tk.END)
        if fit_flag == 0:
            self.xdata[self.nsets] = {'values': xout, 'lowerror': xlowerror,
                                      'higherror': xhigherror,
                                      'minimum': xmin, 'maximum': xmax,
                                      'errors': False, 'legend': True}
            self.ydata[self.nsets] = {'values': yout, 'lowerror': ylowerror,
                                      'higherror': yhigherror,
                                      'minimum': ymin, 'maximum': ymax,
                                      'errors': False, 'legend': True}
        else:
            self.xdata[self.nsets] = {'values': yout, 'lowerror': ylowerror,
                                      'higherror': yhigherror,
                                      'minimum': ymin, 'maximum': ymax,
                                      'errors': False, 'legend': True}
            self.ydata[self.nsets] = {'values': xout, 'lowerror': xlowerror,
                                      'higherror': xhigherror,
                                      'minimum': xmin, 'maximum': xmax,
                                      'errors': False, 'legend': True}
        m = self.nsets % 10
        n = int(math.floor(self.nsets / 10))
        self.set_properties[self.nsets]['symbol'] = None
        self.set_properties[self.nsets]['linestyle'] = '-'
        self.set_properties[self.nsets]['colour'] = self.colourset[m]
        self.set_properties[self.nsets]['symbolsize'] = 4.0 + 0.3*n
        self.set_properties[self.nsets]['label'] = labelstring
        self.nsets = self.nsets + 1
        self.make_plot()

    def fit_statistics(self, yvalues, yfit):
        """
        Calculate the data set fitting statistics.

        This routine calculates the root mean square deviation of the fit
        values from the input values.

        Parameters
        ----------
            yvalues :   a numpy array of data values

            yfit :      a numpy array of fit values

        Both the arrays need to have the same dimensions.

        Returns
        -------
            rms:       the root mean square deviation between the two arrays,
                       or None if the arrays are of different shapes.

        """
        sh1 = yvalues.shape
        sh2 = yfit.shape
        if sh1 != sh2:
            return None
        deviations = yvalues - yfit
        npoints = yvalues.size
        rms = math.sqrt(numpy.sum(deviations*deviations)/npoints)
        return rms

    def make_data_set_transformation_window(self):
        """
        Create the data set transformation window.

        This routine makes a window within which a data set transformation
        can be defined.  The user can either transform an existing set or
        use the transformation to make a new set.

        No values are passed to this routine or returned from this routine.

        """
        if self.data_set_transformation_window is not None:
            return
        self.data_set_transformation_window = Tk.Toplevel()
        self.data_set_transformation_window.title(
            'Data Set Transformation Window')
        holder = Tk.Frame(self.data_set_transformation_window)
        holder.pack(side=Tk.TOP)
        self.set_list_area1 = tkinter.ttk.Combobox(holder, width=50, height=10)
        self.set_list_area1.pack(side=Tk.TOP)
        if self.nsets == 0:
            setlist = [' ']
        else:
            setlist = []
            for loop in range(self.nsets):
                label = (
                    'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                    (loop+1, len(self.xdata[loop]['values']),
                     self.xdata[loop]['minimum'],
                     self.xdata[loop]['maximum'])
                    + ' y range %13.6g to %13.6g' %
                    (self.ydata[loop]['minimum'],
                     self.ydata[loop]['maximum']))
                setlist.append(label)
        self.set_list_area1['values'] = setlist
        self.set_list_area1.current(0)
#        self.set_list_area1.bind('<<ComboboxSelected>>', self.set_fields1)
        frame1 = Tk.Frame(holder)
        frame1.pack()
        label1 = Tk.Label(frame1, text='Set Transformation:')
        label1.grid(row=0, column=0, columnspan=2)
        label2 = Tk.Label(frame1, text='x transformation')
        label2.grid(row=1, column=0)
        self.set_x_transformation_entry_field = Tk.Entry(frame1, width=30)
        self.set_x_transformation_entry_field.grid(row=1, column=1)
        label2 = Tk.Label(frame1, text='y transformation')
        label2.grid(row=2, column=0)
        self.set_y_transformation_entry_field = Tk.Entry(frame1, width=30)
        self.set_y_transformation_entry_field.grid(row=2, column=1)
        frame1 = Tk.Frame(holder)
        frame1.pack()
        label1 = Tk.Label(frame1, text='Create new set?')
        label1.grid(row=0, column=0)
        self.new_set = Tk.IntVar()
        b1 = Tk.Frame(frame1)
        self.put_yes_no(b1, self.new_set, ['Yes', 'No'], True)
        b1.grid(column=1, row=0, sticky=Tk.W)
        label1 = Tk.Label(holder, text="Enter the function you wish to apply.'\
                          + '  Use $x for the x values of the \n'\
                          + 'current set and $y for the y values of the '\
                          + 'current set selected above.  \nOne can also '\
                          + 'use $i for an index starting at 1.")
        label1.pack(side=Tk.TOP)
        frame2 = Tk.Frame(holder)
        frame2.pack()
        apply_button = Tk.Button(frame2, text="Apply",
                                 command=self.apply_transformation)
        apply_button.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            frame2, text="Close",
            command=lambda: self.close_data_window(
                self.data_set_transformation_window))
        close_button.pack(side=Tk.LEFT)

    def apply_transformation(self):
        """
        Apply the defined set transformation.

        This routine reads the strings defining a data transformation
        and attempts to apply them.  It uses the x and y values of the
        current set as inputs.  One can use either math or numpy functions
        in the expression that is defined, but nothing else aside from the
        builtin functions.

        No values are passed to this routine or returned from this routine.
        """
        try:
            ind1 = self.set_list_area1.current()
            xvalues = numpy.copy(self.xdata[ind1]['values'])
            yvalues = numpy.copy(self.ydata[ind1]['values'])
            seq = numpy.arange(len(xvalues)) + 1.
            xstr = self.set_x_transformation_entry_field.get()
            xstr = xstr.replace('$x', 'x')
            xstr = xstr.replace('$i', 'seq')
            ystr = self.set_y_transformation_entry_field.get()
            ystr = ystr.replace('$y', 'y')
            xstr = ystr.replace('$i', 'seq')
            x1 = self.my_eval(xstr, seq, xvalues, yvalues)
            y1 = self.my_eval(ystr, seq, xvalues, yvalues)
            if (x1 is None) or (y1 is None):
                return
            newopt = self.new_set.get()
            if newopt == 0:
                self.xdata[ind1]['values'] = x1
                self.ydata[ind1]['values'] = y1
                self.xdata[ind1]['lowerror'] = x1 * 0.
                self.xdata[ind1]['higherror'] = x1 * 0.
                self.ydata[ind1]['lowerror'] = y1 * 0.
                self.ydata[ind1]['higherror'] = y1 * 0.
            else:
                self.add_set(x1, y1, current_plot=self.current_plot)
            self.make_plot()
        except:
            return

    def make_data_set_delete_window(self):
        """
        Create the window for deleting data sets.

        This routine makes a window within which a data set can be
        removed from the data list.

        No values are passed to this routine or returned from this routine.

        """
        if self.data_set_delete_window is not None:
            return
        self.data_set_delete_window = Tk.Toplevel()
        self.data_set_delete_window.title('Data Set Delete Window')
        holder = Tk.Frame(self.data_set_delete_window)
        holder.pack(side=Tk.TOP)
        self.set_list_area2 = tkinter.ttk.Combobox(holder, width=50, height=10)
        self.set_list_area2.pack(side=Tk.TOP)
        if self.nsets == 0:
            setlist = [' ']
        else:
            setlist = []
            for loop in range(self.nsets):
                label = (
                    'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                    (loop+1, len(self.xdata[loop]['values']),
                     self.xdata[loop]['minimum'],
                     self.xdata[loop]['maximum'])
                    + ' y range %13.6g to %13.6g' %
                    (self.ydata[loop]['minimum'],
                     self.ydata[loop]['maximum']))
                setlist.append(label)
        self.set_list_area2['values'] = setlist
        self.set_list_area2.current(0)
#        self.set_list_area2.bind('<<ComboboxSelected>>', self.set_fields1)
        label1 = Tk.Label(
            holder,
            text="Select a set above and remove with the 'delete' button.\n"
            + " Note that this action cannot be undone so be careful.")
        label1.pack(side=Tk.TOP)
        frame2 = Tk.Frame(holder)
        frame2.pack()
        apply_button = Tk.Button(frame2, text="Delete Set",
                                 command=self.delete_set)
        apply_button.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            frame2, text="Close",
            command=lambda: self.close_data_window(
                self.data_set_delete_window))
        close_button.pack(side=Tk.LEFT)

    def delete_set(self):
        """
        Carry out clearing a data set from the variables.

        This is a work routine that takes the input for deleting a set, checks
        with the user, and if requested removes a set from the list.

        No values are passed to this routine or returned from it.
        """
        ind1 = self.set_list_area2.current()
        if (ind1 >= 0) and (ind1 < self.nsets):
            response = tkinter.messagebox.askyesno(
                'Verify',
                'Do you really want to remove set %d?' % (ind1+1))
            if response:
                if self.nsets > 0:
                    self.nsets = self.nsets - 1
                    for loop in range(ind1, self.nsets):
                        self.set_properties[loop] = deepcopy(
                            self.set_properties[loop+1])
                        self.xdata[loop] = deepcopy(self.xdata[loop+1])
                        self.ydata[loop] = deepcopy(self.ydata[loop+1])
                    self.xdata[self.nsets+1] = None
                    self.ydata[self.nsets+1] = None
                    self.set_properties[self.nsets+1]['symbol'] = None
                    self.set_properties[self.nsets+1]['symbolsize'] = 4.0
                    self.set_properties[self.nsets+1]['linestyle'] = 'None'
                    self.set_properties[self.nsets+1]['linewidth'] = 1.0
                    self.set_properties[self.nsets+1]['colour'] = 'black'
                    self.set_properties[self.nsets+1]['label'] = ''
                    self.set_properties[self.nsets+1]['xmin'] = 0.0
                    self.set_properties[self.nsets+1]['xmax'] = 0.0
                    self.set_properties[self.nsets+1]['ymin'] = 0.0
                    self.set_properties[self.nsets+1]['ymax'] = 0.0
                    self.set_properties[self.nsets+1]['display'] = True
                    self.set_properties[self.nsets+1]['errors'] = False
                    self.set_properties[self.nsets+1]['legend'] = True
                    self.set_properties[self.nsets+1]['plot'] = 1
                    self.make_plot()
                    # revise the set list area
                    setlist = []
                    for loop in range(self.nsets):
                        label = (
                            'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                            (loop+1, len(self.xdata[loop]['values']),
                             self.xdata[loop]['minimum'],
                             self.xdata[loop]['maximum'])
                            + ' y range %13.6g to %13.6g' %
                            (self.ydata[loop]['minimum'],
                             self.ydata[loop]['maximum']))
                        setlist.append(label)
                    self.set_list_area2['values'] = setlist

    def make_data_set_window(self):
        """
        Create the data set properties window (symbol, colour, etc).

        The makes a window within which one can set the display parameters for
        the data sets (symbol, colour, line properties, legend labels).  When
        one exits the window is destroyed, so it gets regenerated each time
        the parent button is clicked.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.data_set_window is not None:
            return
        self.data_set_window = Tk.Toplevel()
        self.data_set_window.title('Data Sets Window')
        holder = Tk.Frame(self.data_set_window)
        holder.pack(side=Tk.TOP)
        self.set_list_area = tkinter.ttk.Combobox(holder, width=50, height=10)
        self.set_list_area.pack(side=Tk.TOP)
        if self.nsets == 0:
            setlist = [' ']
        else:
            setlist = []
            for loop in range(self.nsets):
                label = (
                    'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                    (loop+1, len(self.xdata[loop]['values']),
                     self.xdata[loop]['minimum'], self.xdata[loop]['maximum'])
                    + ' y range %13.6g to %13.6g' %
                    (self.ydata[loop]['minimum'], self.ydata[loop]['maximum'])
                )
                setlist.append(label)
        self.set_list_area['values'] = setlist
        self.set_list_area.current(0)
        self.set_list_area.bind('<<ComboboxSelected>>', self.set_fields)
        frame1 = Tk.Frame(holder)
        frame1.pack()
        self.set_window_label = Tk.Label(frame1, text='Set Properties:')
        self.set_window_label.grid(row=0, column=0, columnspan=2)
        labels = ['Symbol', 'Symbol Size', 'Line', 'Line Width', 'Colour',
                  'Show', 'Error Bars', 'Label', 'Legend', 'Plot']
        for loop in range(len(labels)):
            label = Tk.Label(frame1, text=labels[loop])
            label.grid(column=0, row=loop+1, sticky=Tk.W)
        self.set_entry_fields = []
        # item 0: menu for the symbol
        self.set_entry_fields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.set_entry_fields[-1].grid(column=1, row=1, sticky=Tk.W)
        self.set_entry_fields[-1]['values'] = matplotlib_symbol_name_list
        self.set_entry_fields[-1].current(0)
        # item 1: entry field for the symbol size value
        self.set_entry_fields.append(Tk.Entry(frame1, width=15))
        self.set_entry_fields[-1].grid(column=1, row=2, sticky=Tk.W)
        self.set_entry_fields[-1].insert(0, '1.0')
        # item 2: menu for the line style
        self.set_entry_fields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.set_entry_fields[-1].grid(column=1, row=3, sticky=Tk.W)
        self.set_entry_fields[-1]['values'] = matplotlib_line_name_list
        self.set_entry_fields[-1].current(1)
        # item 3: entry field for the line width value
        self.set_entry_fields.append(Tk.Entry(frame1, width=15))
        self.set_entry_fields[-1].grid(column=1, row=4, sticky=Tk.W)
        # item 4: menu for the colour value
        self.set_entry_fields.append(tkinter.ttk.Combobox(frame1, width=15))
        self.set_entry_fields[-1].grid(column=1, row=5, sticky=Tk.W)
        self.set_entry_fields[-1]['values'] = self.colourset
        self.set_entry_fields[-1].current(0)
        # item 5: variable for the show/hide set option
        self.set_entry_fields.append(Tk.IntVar())
        b1 = Tk.Frame(frame1)
        self.put_yes_no(b1, self.set_entry_fields[-1], ['Yes', 'No'], True)
        b1.grid(column=1, row=6, sticky=Tk.W)
        # item 6: variable for the error bar display option
        self.set_entry_fields.append(Tk.IntVar())
        b1 = Tk.Frame(frame1)
        self.put_yes_no(b1, self.set_entry_fields[-1], ['Yes', 'No'], False)
        b1.grid(column=1, row=7, sticky=Tk.W)
        # item 7: entry field for the set label (for the legend, if any)
        self.set_entry_fields.append(Tk.Entry(frame1, width=35))
        self.set_entry_fields[-1].grid(column=1, row=8, sticky=Tk.W)
        # item 8: flag for whether to include the set in the legend
        self.set_entry_fields.append(Tk.IntVar())
        b1 = Tk.Frame(frame1)
        self.put_yes_no(b1, self.set_entry_fields[-1], ['Yes', 'No'], True)
        b1.grid(column=1, row=9, sticky=Tk.W)
        # item 9: Entry field for which plot the set belongs to
        self.set_entry_fields.append(Tk.Entry(frame1, width=5))
        self.set_entry_fields[-1].grid(column=1, row=10, sticky=Tk.W)
        self.set_entry_fields[-1].insert(0, str(self.current_plot))
        frame2 = Tk.Frame(holder)
        frame2.pack()
        apply_button = Tk.Button(
            frame2, text="Apply",
            command=self.apply_data_set_fields)
        apply_button.pack(side=Tk.LEFT)
        label1 = Tk.Label(frame2, text="    ")
        label1.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            frame2, text="Close",
            command=lambda: self.close_data_window(self.data_set_window))
        close_button.pack(side=Tk.LEFT)
        if self.nsets > 0:
            self.set_property_fields(0)

    def put_yes_no(self, root, var, labels, flag):
        """
        Create a Tkineter yes/no radio button.

        This is a utility routine to make a yes/no radio button in Tkinter.
        The required variables are passed to the code and it produces the
        button pair.

        Parameters
        ----------
            root :  The Tk frame variable to hold the buttons

            var :   The Tk IntVar that is used to communicate with the
                    buttons

            labels : A two element list with strings for the two states of
                     the buttons, first "yes" then "no".

            flag :  A boolean value that causes the code to set the yes field
                    (the first of the two) if it is True.  If the value is
                    False the no field (the second of the two) is set instead.

        Returns
        -------
            No values are returned by this routine.

        """
        yesfield = Tk.Radiobutton(root, text=labels[0], variable=var, value=1)
        yesfield.grid(row=0, column=0, sticky=Tk.W)
        nofield = Tk.Radiobutton(root, text=labels[1], variable=var, value=0)
        nofield.grid(row=0, column=1, sticky=Tk.W)
        if flag:
            yesfield.select()
        else:
            nofield.select()

    def set_fields(self, event):
        """
        Update the set selected in the set parameters window.

        This is a callback routine.  It determines the value of the set that
        is currently selected in the Combobox, and if this is an actual set
        it calls the routine to change the set parameter values in the set
        properties window.

        A check is needed because one might get an event even before any sets
        have been read in.  The code then calls the routine to change the
        parameter values.

        Parameters
        ----------
            event:  The Tkinter ComboboxSelected event (required, but not
                    used)

        Returns
        -------
            No values are returned by this routine.

        """
        index = self.set_list_area.current()
        if index < 0:
            tkinter.messagebox.showinfo('Error', 'You need to select a set.')
            return
        if index >= self.nsets:
            tkinter.messagebox.showinfo('Error', 'You need to select a set.')
            return
        self.set_property_fields(index)

    def set_property_fields(self, index):
        """
        Update the set properties window fields for the selected set.

        This routine attempts to change the set properties in the
        self.data_set_window based on which set is currently selected.
        It replaces the current values with those from a particular set
        using the self.set_properties values.

        Parameters
        ----------
        index :  The index of the currently selected set in the Combobox.
                 The value can be -1 if no line is selected, in which case
                 the routine does not take any action.

        Returns
        -------
            No values are returned by this routine.

        """
        if (index < 0) | (index >= self.nsets):
            return
        try:
            ind1 = matplotlib_symbol_list.index(
                self.set_properties[index]['symbol'])
            self.set_entry_fields[0].current(ind1)
            self.set_entry_fields[1].delete(0, Tk.END)
            self.set_entry_fields[1].insert(
                0, str(self.set_properties[index]['symbolsize']))
            try:
                ind1 = matplotlib_line_list.index(
                    self.set_properties[index]['linestyle'])
                self.set_entry_fields[2].current(ind1)
            except:
                self.set_entry_fields[2].current(4)
            self.set_entry_fields[3].delete(0, Tk.END)
            self.set_entry_fields[3].insert(
                0, str(self.set_properties[index]['linewidth']))
            try:
                self.set_entry_fields[4].set(
                    self.set_properties[index]['colour'])
            except:
                self.set_entry_fields[4].current(10)
            if self.set_properties[index]['display']:
                self.set_entry_fields[5].set(1)
            else:
                self.set_entry_fields[5].set(0)
            if self.set_properties[index]['errors']:
                self.set_entry_fields[6].set(1)
            else:
                self.set_entry_fields[6].set(0)
            self.set_entry_fields[7].delete(0, Tk.END)
            self.set_entry_fields[7].insert(
                0, self.set_properties[index]['label'])
            if self.set_properties[index]['legend']:
                self.set_entry_fields[8].set(1)
            else:
                self.set_entry_fields[8].set(0)
            self.set_entry_fields[9].delete(0, Tk.END)
            self.set_entry_fields[9].insert(
                0, str(self.set_properties[index]['plot']))
            self.make_plot()
        except:
            pass

    def apply_data_set_fields(self):
        """
        Apply the values in the data set window.

        Thie routine reads the values in the self.data_set_window and applies
        these to the selected data set on the list.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            index = self.set_list_area.current()
            symbol_index = self.set_entry_fields[0].current()
            line_index = self.set_entry_fields[2].current()
            colour_index = self.set_entry_fields[4].current()
            label = self.set_entry_fields[7].get()
            symbol_size = float(self.set_entry_fields[1].get())
            line_width = float(self.set_entry_fields[3].get())
            show_flag = self.set_entry_fields[5].get()
            error_bar_flag = self.set_entry_fields[6].get()
            legend_flag = self.set_entry_fields[8].get()
            target_plot = int(self.set_entry_fields[9].get())
            if (target_plot < 1) or (target_plot > self.number_of_plots):
                target_plot = self.set_properties[index]['plot']
            self.set_properties[index]['symbol'] = \
                matplotlib_symbol_list[symbol_index]
            self.set_properties[index]['linestyle'] = \
                matplotlib_line_list[line_index]
            if self.colourset[colour_index] == 'select':
                values = askcolor()
                self.set_properties[index]['colour'] = values[1]
            else:
                self.set_properties[index]['colour'] = \
                    self.colourset[colour_index]
            self.set_properties[index]['label'] = label
            self.set_properties[index]['linewidth'] = line_width
            self.set_properties[index]['symbolsize'] = symbol_size
            if show_flag == 1:
                self.set_properties[index]['display'] = True
            else:
                self.set_properties[index]['display'] = False
            if error_bar_flag == 1:
                self.set_properties[index]['errors'] = True
            else:
                self.set_properties[index]['errors'] = False
            if legend_flag == 1:
                self.set_properties[index]['legend'] = True
            else:
                self.set_properties[index]['legend'] = False
            self.set_properties[index]['plot'] = target_plot
            self.make_plot()
        except:
            tkinter.messagebox.showinfo(
                'Error',
                'There was some error applying the values.'
                + '  You may need to retry the settings.')

    def make_plot_control_window(self):
        """
        Create the plot control window (tick marks, axis labels, etc).

        This routine produces the plot control window wherein one sets the plot
        properties such as the axis labels and the title.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.plot_control_window is not None:
            return
        self.plot_control_window = Tk.Toplevel()
        self.plot_control_window.title('Plot Parameters')
        outframe = Tk.Frame(self.plot_control_window)
        outframe.pack(side=Tk.TOP)
        holder = Tk.Frame(outframe)
        holder.pack(side=Tk.LEFT)
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Plot Title:')
        label.pack(side=Tk.LEFT)
        self.title_field = Tk.Entry(field1, width=20)
        self.title_field.pack(side=Tk.LEFT)
        self.title_field.insert(0, self.title[self.current_plot-1])
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='X Axis Label:')
        label.pack(side=Tk.LEFT)
        self.xlabel_field = Tk.Entry(field1, width=20)
        self.xlabel_field.pack(side=Tk.LEFT)
        self.xlabel_field.insert(
            0, self.xparameters[self.current_plot-1]['label'])
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Y Axis Label:')
        label.pack(side=Tk.LEFT)
        self.ylabel_field = Tk.Entry(field1, width=20)
        self.ylabel_field.pack(side=Tk.LEFT)
        self.ylabel_field.insert(
            0, self.yparameters[self.current_plot-1]['label'])
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='X Axis Minimum:')
        label.pack(side=Tk.LEFT)
        self.xmin_field = Tk.Entry(field1, width=20)
        self.xmin_field.pack(side=Tk.LEFT)
        if self.xparameters[self.current_plot-1]['hybridlog'] == 1:
            xlim = inverse_hybrid_transform(
                self.xparameters[self.current_plot-1]['minimum'])
            self.xmin_field.insert(0, xlim)
        else:
            self.xmin_field.insert(
                0, self.xparameters[self.current_plot-1]['minimum'])
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='X Axis Maximum:')
        label.pack(side=Tk.LEFT)
        self.xmax_field = Tk.Entry(field1, width=20)
        self.xmax_field.pack(side=Tk.LEFT)
        if self.xparameters[self.current_plot-1]['hybridlog'] == 1:
            xlim = inverse_hybrid_transform(
                self.xparameters[self.current_plot-1]['maximum'])
            self.xmax_field.insert(0, xlim)
        else:
            self.xmax_field.insert(
                0, self.xparameters[self.current_plot-1]['maximum'])
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Y Axis Minimum:')
        label.pack(side=Tk.LEFT)
        self.ymin_field = Tk.Entry(field1, width=20)
        self.ymin_field.pack(side=Tk.LEFT)
        if self.yparameters[self.current_plot-1]['hybridlog'] == 1:
            ylim = inverse_hybrid_transform(
                self.yparameters[self.current_plot-1]['minimum'])
            self.ymin_field.insert(0, ylim)
        else:
            self.ymin_field.insert(
                0, self.yparameters[self.current_plot-1]['minimum'])
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Y Axis Maximum:')
        label.pack(side=Tk.LEFT)
        self.ymax_field = Tk.Entry(field1, width=20)
        self.ymax_field.pack(side=Tk.LEFT)
        if self.yparameters[self.current_plot-1]['hybridlog'] == 1:
            ylim = inverse_hybrid_transform(
                self.yparameters[self.current_plot-1]['maximum'])
            self.ymax_field.insert(0, ylim)
        else:
            self.ymax_field.insert(
                0, self.yparameters[self.current_plot-1]['maximum'])
        holder1 = Tk.Frame(outframe)
        holder1.pack(side=Tk.LEFT)
        # Note it is not possible to put in the line height
        # automatically as the frames it is intended to match have
        # not been packed yet .
        sl = self.separator_line(holder1, 5, 450, 5, False)
        holder2 = Tk.Frame(outframe)
        holder2.pack(side=Tk.LEFT)
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="logarithmic x axis")
        label.pack(side=Tk.LEFT)
        self.logx_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.xparameters[self.current_plot-1]['logarithmic'] == 1
        self.put_yes_no(b1, self.logx_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="hybrid log x axis")
        label.pack(side=Tk.LEFT)
        self.hlogx_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.xparameters[self.current_plot-1]['hybridlog'] == 1
        self.put_yes_no(b1, self.hlogx_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="logarithmic y axis")
        label.pack(side=Tk.LEFT)
        self.logy_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.yparameters[self.current_plot-1]['logarithmic'] == 1
        self.put_yes_no(b1, self.logy_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="hybrid log y axis")
        label.pack(side=Tk.LEFT)
        self.hlogy_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.yparameters[self.current_plot-1]['hybridlog'] == 1
        self.put_yes_no(b1, self.hlogy_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="invert x axis")
        label.pack(side=Tk.LEFT)
        self.invertx_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.xparameters[self.current_plot-1]['invert'] == 1
        self.put_yes_no(b1, self.invertx_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="invert y axis")
        label.pack(side=Tk.LEFT)
        self.inverty_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.yparameters[self.current_plot-1]['invert'] == 1
        self.put_yes_no(b1, self.inverty_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="hide x axis")
        label.pack(side=Tk.LEFT)
        self.hidex_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.xparameters[self.current_plot-1]['hide'] == 1
        self.put_yes_no(b1, self.hidex_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="hide y axis")
        label.pack(side=Tk.LEFT)
        self.hidey_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.yparameters[self.current_plot-1]['hide'] == 1
        self.put_yes_no(b1, self.hidey_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="hide x ticks")
        label.pack(side=Tk.LEFT)
        self.hidexticks_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.xparameters[self.current_plot-1]['hideticks'] == 1
        self.put_yes_no(b1, self.hidexticks_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="hide y ticks")
        label.pack(side=Tk.LEFT)
        self.hideyticks_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.yparameters[self.current_plot-1]['hideticks'] == 1
        self.put_yes_no(b1, self.hideyticks_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="hide x labels")
        label.pack(side=Tk.LEFT)
        self.hidexlabels_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.xparameters[self.current_plot-1]['hidelabels'] == 1
        self.put_yes_no(b1, self.hidexlabels_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="hide y labels")
        label.pack(side=Tk.LEFT)
        self.hideylabels_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.yparameters[self.current_plot-1]['hidelabels'] == 1
        self.put_yes_no(b1, self.hideylabels_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="x ticks both sides")
        label.pack(side=Tk.LEFT)
        self.bothxticks_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.xparameters[self.current_plot-1]['bothticks'] == 1
        self.put_yes_no(b1, self.bothxticks_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="y ticks both sides")
        label.pack(side=Tk.LEFT)
        self.bothyticks_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.yparameters[self.current_plot-1]['bothticks'] == 1
        self.put_yes_no(b1, self.bothyticks_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="opposite x axis")
        label.pack(side=Tk.LEFT)
        self.oppositex_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.xparameters[self.current_plot-1]['oppositeaxis'] == 1
        self.put_yes_no(b1, self.oppositex_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="opposite y axis")
        label.pack(side=Tk.LEFT)
        self.oppositey_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.yparameters[self.current_plot-1]['oppositeaxis'] == 1
        self.put_yes_no(b1, self.oppositey_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        b1 = Tk.Frame(field1)
        label = Tk.Label(field1, text="invert x ticks")
        label.pack(side=Tk.LEFT)
        self.inversexticks_variable = Tk.IntVar()
        flag = self.xparameters[self.current_plot-1]['inverseticks'] == 1
        self.put_yes_no(b1, self.inversexticks_variable, ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder2)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text="invert y ticks")
        label.pack(side=Tk.LEFT)
        self.inverseyticks_variable = Tk.IntVar()
        b1 = Tk.Frame(field1)
        flag = self.yparameters[self.current_plot-1]['inverseticks'] == 1
        self.put_yes_no(b1, self.inverseyticks_variable, ['Yes', 'No'], flag)
        b1.pack()
        #
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Tick length:')
        label.pack(side=Tk.LEFT)
        self.ticklengthfield = Tk.Entry(field1, width=3)
        self.ticklengthfield.pack()
        try:
            self.ticklengthfield.insert(
                0, str(self.xparameters[self.current_plot-1]['ticklength']))
        except:
            self.ticklengthfield.insert(0, '6')
            self.xparameters[self.current_plot-1]['ticklength'] = 6
            self.yparameters[self.current_plot-1]['ticklength'] = 6
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='X Minor tick interval:')
        label.pack(side=Tk.LEFT)
        self.xminortickfield = Tk.Entry(field1, width=10)
        self.xminortickfield.pack()
        self.xminortickfield.insert(
            0, str(self.xparameters[self.current_plot-1]['minorticks']))
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Y Minor tick interval:')
        label.pack(side=Tk.LEFT)
        self.yminortickfield = Tk.Entry(field1, width=10)
        self.yminortickfield.pack()
        self.yminortickfield.insert(
            0, str(self.yparameters[self.current_plot-1]['minorticks']))
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Plot Legend:')
        label.pack(side=Tk.LEFT)
        if self.legend_variable[self.current_plot-1] is None:
            self.legend_variable[self.current_plot-1] = Tk.IntVar()
            flag = False
        else:
            flag = self.legend_variable[self.current_plot-1].get()
        b1 = Tk.Frame(field1)
        self.put_yes_no(b1, self.legend_variable[self.current_plot-1],
                        ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Legend type:')
        label.pack(side=Tk.LEFT)
        if self.legend_options[self.current_plot-1] is None:
            self.legend_options[self.current_plot-1] = Tk.StringVar()
            self.legend_options[self.current_plot-1].set('best')
            self.legend_position[self.current_plot-1] = 'best'
        else:
            self.legend_options[self.current_plot-1].set(
                self.legend_position[self.current_plot-1])
        self.legend_option_list = ['user', 'best', 'upper right',
                                   'upper left', 'lower left',
                                   'lower right', 'right',
                                   'center left', 'center right',
                                   'lower center', 'upper center',
                                   'center']
        menu1 = Tk.OptionMenu(field1, self.legend_options[self.current_plot-1],
                              *self.legend_option_list,
                              command=self.generate_legend)
        menu1.config(width=10)
        menu1.pack(side=Tk.LEFT)
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='User legend position:')
        label.pack(side=Tk.LEFT)
        self.legend_position_field = Tk.Entry(field1, width=20)
        self.legend_position_field.pack(side=Tk.LEFT)
        if self.legend_user_position[self.current_plot-1] is None:
            self.legend_position_field.insert(0, '0.1 0.1')
        else:
            str1 = '%.3f %.3f' % (
                self.legend_user_position[self.current_plot-1][0],
                self.legend_user_position[self.current_plot-1][1])
            self.legend_position_field.insert(0, str1)
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Legend frame:')
        label.pack(side=Tk.LEFT)
        if self.legend_frame[self.current_plot-1] is None:
            self.legend_frame[self.current_plot-1] = Tk.IntVar()
            self.legend_frame[self.current_plot-1].set(0)
            flag = False
        else:
            flag = self.legend_frame[self.current_plot-1].get()
        b1 = Tk.Frame(field1)
        self.put_yes_no(b1, self.legend_frame[self.current_plot-1],
                        ['Yes', 'No'], flag)
        b1.pack()
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Plot margin:')
        label.pack(side=Tk.LEFT)
        self.plot_margin_field = Tk.Entry(field1, width=10)
        self.plot_margin_field.pack(side=Tk.LEFT)
        if self.plot_margin is None:
            self.plot_margin_field.insert(0, '0.0')
            self.plot_margin = 0.0
        else:
            self.plot_margin_field.insert(0, '%f' % (self.plot_margin))
        field1 = Tk.Frame(holder)
        field1.pack(side=Tk.TOP)
        label = Tk.Label(field1, text='Frame Width:')
        label.pack(side=Tk.LEFT)
        self.plot_frame_field = Tk.Entry(field1, width=10)
        self.plot_frame_field.pack(side=Tk.LEFT)
        if (self.plot_frame[self.current_plot-1] is None) or \
           (self.plot_frame[self.current_plot-1] <= 0.):
            self.plot_frame_field.insert(0, '0.0')
            self.plot_frame[self.current_plot-1] = 0.0
        else:
            self.plot_frame_field.insert(0, '%f' % (
                self.plot_frame[self.current_plot-1]))
        holder = Tk.Frame(self.plot_control_window)
        holder.pack(side=Tk.TOP)
        # find the width for the separator line and apply it...
        outframe.update()
        sl = self.separator_line(holder, outframe.winfo_width(), 5, 5, True)
        field1 = Tk.Frame(self.plot_control_window)
        field1.pack(side=Tk.TOP)
        apply_button = Tk.Button(field1, text="Apply",
                                 command=self.apply_plot_parameters)
        apply_button.pack(side=Tk.LEFT)
        label1 = Tk.Label(field1, text="    ")
        label1.pack(side=Tk.LEFT)
        apply_all_button = Tk.Button(field1, text="Apply to All",
                                     command=self.apply_plot_parameters_all)
        apply_all_button.pack(side=Tk.LEFT)
        label1 = Tk.Label(field1, text="    ")
        label1.pack(side=Tk.LEFT)
        close_button = Tk.Button(
            field1, text="Close",
            command=lambda: self.close_window(self.plot_control_window))
        close_button.pack(side=Tk.LEFT)

    def separator_line(self, parent, w1, h1, pad, flag, packvalue=Tk.TOP):
        """
        Create the Tkinter canvas object making a separator line.

        This is a utility routine to make a separator line within a Tkinter
        frame.

        Parameters
        ----------
            parent :  A Tkinter Frame variable, that will contain the line

            w1 :      An integer value for the line canvas width (pixels)

            h1 :      An integer value for the line canvas height (pixels)

            pad :     An integer value for the line padding (pixels)

            flag :    A Boolean value, the line is horizontal if the
                      value is True, vertical otherwise

            packvalue :  An optional Tkinter pack direction (one of Tk.LEFT,
                         Tk.RIGHT, Tk.TOP, or Tk.BOTTOM) with default value
                         of Tk.TOP

        Returns
        -------
            linecanvas:  the Tkinter Canvas variable for the line, in case
                         the user wishes to modify it later

        For a vertical line normally the height will be much larger than the
        width, as in

        sl = self.separator_line(frame, 10, 300, 5, False)

        while for a horizontal line normally the width will be much larger
        than the height as in

        sl = self.separator_line(frame, 500, 5, 5, True)

        """
        lincanvas = Tk.Canvas(parent, height=h1, width=w1)
        try:
            lincanvas.pack(side=packvalue, fill=Tk.BOTH, expand=Tk.YES)
        except:
            pass
        if flag:
            lincanvas.create_line(pad, h1/2, w1-pad, h1/2)
        else:
            lincanvas.create_line(w1/2, pad, w1/2, h1-pad)
        return lincanvas

    def generate_legend(self, event):
        """
        Generate the legend values from the sets of the plot.

        Parameters
        ----------
            event:   A Tkinter event or None.  The value is not used.

        Returns
        -------
           No value is returned by the code.

        The routine sets the self.legend_handles and self.legend_labels
        variables according to which sets have the 'legend' flag set to True.

        """
        if self.nsets == 0:
            return
        self.legend_handles[self.current_plot-1] = []
        self.legend_labels[self.current_plot-1] = []
        for loop in range(self.nsets):
            if (self.set_properties[loop]['legend']) and \
               (self.set_properties[loop]['plot'] == self.current_plot):
                if self.set_properties[loop]['symbol'] is None:
                    proxy = mlines.Line2D(
                        [], [],
                        color=self.set_properties[loop]['colour'],
                        linestyle=self.set_properties[loop]['linestyle'],
                        linewidth=self.set_properties[loop]['linewidth'])
                else:
                    proxy = mlines.Line2D(
                        [], [],
                        color=self.set_properties[loop]['colour'],
                        marker=self.set_properties[loop]['symbol'],
                        markersize=self.set_properties[loop]['symbolsize'],
                        linestyle='')
                self.legend_handles[self.current_plot-1].append(proxy)
                self.legend_labels[self.current_plot-1].append(
                    self.set_properties[loop]['label'])

    def apply_plot_parameters_all(self):
        """
        Apply the plot parameters to all active plots.

        Parameters
        ----------

        Nothing

        Returns
        -------

        None
        """
        plotnumber = self.current_plot
        for loop in range(self.number_of_plots):
            self.current_plot = loop+1
            if self.legend_options[self.current_plot-1] is None:
                self.legend_options[self.current_plot-1] = Tk.StringVar()
                self.legend_options[self.current_plot-1].set('best')
                self.legend_position[self.current_plot-1] = 'best'
            self.apply_plot_parameters()
        self.current_plot = plotnumber

    def apply_plot_parameters(self):
        """
        Apply the values from the plot parameters window.

        This routine is called when the values in the plot control window
        are to be applied to the plot.  It reads the values and places them
        in the internal variables, and then replots the data.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        try:
            self.title[self.current_plot-1] = self.title_field.get()
            instring = self.xlabel_field.get()
            # To do: test the string is OK as a label here....
            self.xparameters[self.current_plot-1]['label'] = instring
            instring = self.ylabel_field.get()
            self.yparameters[self.current_plot-1]['label'] = instring
            self.xparameters[self.current_plot-1]['logarithmic'] = \
                self.logx_variable.get()
            self.yparameters[self.current_plot-1]['logarithmic'] = \
                self.logy_variable.get()
            self.xparameters[self.current_plot-1]['hybridlog'] = \
                self.hlogx_variable.get()
            self.yparameters[self.current_plot-1]['hybridlog'] = \
                self.hlogy_variable.get()
            self.xparameters[self.current_plot-1]['invert'] = \
                self.invertx_variable.get()
            self.yparameters[self.current_plot-1]['invert'] = \
                self.inverty_variable.get()
            self.xparameters[self.current_plot-1]['hide'] = \
                self.hidex_variable.get()
            self.yparameters[self.current_plot-1]['hide'] = \
                self.hidey_variable.get()
            self.xparameters[self.current_plot-1]['hideticks'] = \
                self.hidexticks_variable.get()
            self.yparameters[self.current_plot-1]['hideticks'] = \
                self.hideyticks_variable.get()
            self.xparameters[self.current_plot-1]['hidelabels'] = \
                self.hidexlabels_variable.get()
            self.yparameters[self.current_plot-1]['hidelabels'] = \
                self.hideylabels_variable.get()
            self.xparameters[self.current_plot-1]['inverseticks'] = \
                self.inversexticks_variable.get()
            self.yparameters[self.current_plot-1]['inverseticks'] = \
                self.inverseyticks_variable.get()
            self.xparameters[self.current_plot-1]['bothticks'] = \
                self.bothxticks_variable.get()
            self.yparameters[self.current_plot-1]['bothticks'] = \
                self.bothyticks_variable.get()
            self.xparameters[self.current_plot-1]['oppositeaxis'] = \
                self.oppositex_variable.get()
            self.yparameters[self.current_plot-1]['oppositeaxis'] = \
                self.oppositey_variable.get()
            try:
                self.xparameters[self.current_plot-1]['minorticks'] = \
                     float(self.xminortickfield.get())
            except:
                self.xparameters[self.current_plot-1]['minorticks'] = 0.0
            if self.xparameters[self.current_plot-1]['minorticks'] < 0.0:
                self.xparameters[self.current_plot-1]['minorticks'] = 0.0
            try:
                self.yparameters[self.current_plot-1]['minorticks'] = \
                    float(self.yminortickfield.get())
            except:
                self.yparameters[self.current_plot-1]['minorticks'] = 0.0
            if self.yparameters[self.current_plot-1]['minorticks'] < 0.0:
                self.yparameters[self.current_plot-1]['minorticks'] = 0.0
            try:
                ticklength = int(self.ticklengthfield.get())
                if ticklength < 1:
                    ticklength = 1
            except:
                ticklength = 6
            self.xparameters[self.current_plot-1]['ticklength'] = ticklength
            self.yparameters[self.current_plot-1]['ticklength'] = ticklength
            try:
                frame = float(self.plot_frame_field.get())
                if (frame <= 0.) or (frame > 5.):
                    self.plot_frame[self.current_plot-1] = 0.
                else:
                    self.plot_frame[self.current_plot-1] = frame
            except:
                self.plot_frame[self.current_plot-1] = 0.
            try:
                margin = float(self.plot_margin_field.get())
                self.plot_margin = margin
                if margin > 0.3:
                    tkinter.messagebox.showinfo(
                        'Warning',
                        'Margins are limited to the range -0.1 to +0.3.')
                    margin = 0.3
                    self.plot_margin_field.delete(0, Tk.END)
                    self.plot_margin_field.insert(0, '0.3')
                if margin < -0.1:
                    tkinter.messagebox.showinfo(
                        'Warning',
                        'Margins are limited to the range -0.1 to +0.3.')
                    margin = -0.1
                    self.plot_margin_field.delete(0, Tk.END)
                    self.plot_margin_field.insert, (0, '-0.1')
                self.plot_margin = margin
            except:
                tkinter.messagebox.showinfo(
                    'Warning',
                    'Margin value could not be read, setting to zero.')
                self.plot_margin_field.delete(0, Tk.END)
                self.plot_margin_field.insert, (0, '0.0')
            try:
                xmin = float(self.xmin_field.get())
                xmax = float(self.xmax_field.get())
                ymin = float(self.ymin_field.get())
                ymax = float(self.ymax_field.get())
                if (xmin >= xmax) | (ymin >= ymax):
                    raise ValueError
                self.plot_range[self.current_plot-1][0] = xmin
                self.plot_range[self.current_plot-1][1] = xmax
                self.plot_range[self.current_plot-1][2] = ymin
                self.plot_range[self.current_plot-1][3] = ymax
            except:
                tkinter.messagebox.showinfo(
                    'Error',
                    'There was some error in the axis range values.'
                    + '  Please check your inputs.')
            self.legend_position[self.current_plot-1] = \
                self.legend_options[self.current_plot-1].get()
            for loop in range(len(self.share_axis)):
                if self.current_plot == abs(self.share_axis[loop]):
                    if self.share_axis[loop] < 0:
                        for key in self.yparameters[-1].keys():
                            self.yparameters[loop][key] = self.yparameters[
                                self.current_plot-1][key]
                    else:
                        for key in self.xparameters[-1].keys():
                            self.xparameters[loop][key] = self.xparameters[
                                self.current_plot-1][key]
            self.make_plot()
        except:
            tkinter.messagebox.showinfo(
                'Error',
                'There was some error in the input values.  '
                + 'Please check your inputs.')

    def edit_lines(self):
        """
        Make window for editing the line values.

        This routine produces a text box in a window, within which one can
        edit the line values.  If no lines are defined the routine just
        exits with no action.

        Line values are presented one per text line, with the start x
        position, the start y position, the end x position, the end y
        position, the plot number, the line type, the line colour, and the
        line thickness separated by tab symbols.  One can edit the values
        within the text window and then these are applied when one clicks
        on the "Close Window" button.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.number_of_lines == 0:
            return
        str1 = 'Edit values below: fields are separated by tab characters.\n'\
               + '     start x       start y       end x         end y   ' \
               + 'plot        colour   line type  thickness\n------------' \
               + '-------------------------------------------------------' \
               + '-----------------------------\n'
        for loop in range(self.number_of_lines):
            str1 = str1 + '%12.6g\t%12.6g\t%12.6g\t%12.6g\t%6d\t' % (
                self.plot_lines[loop]['xstart'],
                self.plot_lines[loop]['ystart'],
                self.plot_lines[loop]['xend'],
                self.plot_lines[loop]['yend'],
                self.plot_lines[loop]['plot'])
            str1 = str1 + '%15s\t%8s\t%7.3f\n' % (
                self.plot_lines[loop]['line_colour'],
                self.plot_lines[loop]['line_type'],
                self.plot_lines[loop]['line_thickness'])
        line_window = Tk.Toplevel()
        line_window.title('Lines:')
        holder = Tk.Frame(line_window)
        holder.pack(side=Tk.TOP)
        line_message_text = ScrolledText(holder, height=40, width=100,
                                         wrap=Tk.NONE)
        line_message_text.config(font=('courier', 16, 'bold'))
        line_message_text.pack(side=Tk.TOP)
        line_message_text.insert(0.0, str1)
        bholder = Tk.Frame(line_window)
        bholder.pack(side=Tk.TOP)
        close_button = Tk.Button(
            bholder, text='Close Window',
            command=lambda: self.read_lines(line_message_text, line_window))
        close_button.pack()

    def read_lines(self, line_message_text, line_window):
        """
        Read and apply the lines text field.

        This routine reads the line text field and makes the new set of lines
        and line positions.  It then applies these and closes the line window.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        The code does, however, change the self.plot_lines values as needed
        to match what is in the line text field.

        """
        linetext = line_message_text.get(0.0, Tk.END)
        lines = linetext.split('\n')
        newlines = []
        nlines = 0
        for line in lines:
            values = line.split('\t')
            if len(values) == 8:
                try:
                    x1 = float(values[0])
                    y1 = float(values[1])
                    x2 = float(values[2])
                    y2 = float(values[3])
                    nplot = int(values[4])
                    colour = values[5].strip(' ')
                    linetype = values[6].strip(' ')
                    thickness = float(values[7])
                    newlines.append({'xstart': x1, 'ystart': y1, 'xend': x2,
                                     'yend': y2, 'plot': nplot,
                                     'line_type': linetype,
                                     'line_colour': colour,
                                     'line_thickness': thickness})
                    nlines = nlines + 1
                except:
                    pass
        line_window.destroy()
        if nlines > self.max_lines:
            self.max_lines = nlines
        else:
            for loop in range(nlines, self.max_lines):
                newlines.append({'xstart': None, 'ystart': None,
                                 'xend': None, 'yend': None,
                                 'plot': 1, 'line_type': 'solid',
                                 'line_colour': 'black',
                                 'line_thickness': 1.0})
        self.plot_lines = newlines
        self.number_of_lines = nlines
        self.make_plot()

    def edit_boxes(self):
        """
        Make a window for editing the box values.

        This routine produces a text box in a window, within which one can
        edit the box values.  If no boxes are defined the routine just exits
        with no action.

        Box values are presented one per line, with the start x position, the
        start y position, the end x position, the end y position, the
        orientation, the plot number, the line type, the line colour, the
        line thickness, and the fill colour separated by tab symbols.
        One can edit the values within the text window and then these are
        applied when one clicks on the "Close Window" button.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.number_of_boxes == 0:
            return
        str1 = 'Edit values below: fields are separated by tab characters.\n' \
               + '     start x       start y       end x      end y      ' \
               + 'orient.   plot        colour   type  thickness  fill \n' \
               + '-------------------------------------------------------' \
               + '--------------------------------------------------------\n'
        for loop in range(self.number_of_boxes):
            str1 = str1 + '%12.6g\t%12.6g\t%12.6g\t%12.6g\t%8.3f\t%6d' % (
                self.plot_boxes[loop]['xstart'],
                self.plot_boxes[loop]['ystart'],
                self.plot_boxes[loop]['xend'],
                self.plot_boxes[loop]['yend'],
                self.plot_boxes[loop]['rotation'],
                self.plot_boxes[loop]['plot'])
            str1 = str1 + '\t%15s\t%8s\t%7.3f\t%15s\n' % (
                self.plot_boxes[loop]['line_colour'],
                self.plot_boxes[loop]['line_type'],
                self.plot_boxes[loop]['line_thickness'],
                self.plot_boxes[loop]['fill_colour'])
        box_window = Tk.Toplevel()
        box_window.title('Boxes:')
        holder = Tk.Frame(box_window)
        holder.pack(side=Tk.TOP)
        box_message_text = ScrolledText(holder, height=40, width=115,
                                        wrap=Tk.NONE)
        box_message_text.config(font=('courier', 16, 'bold'))
        box_message_text.pack(side=Tk.TOP)
        box_message_text.insert(0.0, str1)
        bholder = Tk.Frame(box_window)
        bholder.pack(side=Tk.TOP)
        close_button = Tk.Button(
            bholder, text='Close Window',
            command=lambda: self.read_boxes(box_message_text, box_window))
        close_button.pack()

    def read_boxes(self, box_message_text, box_window):
        """
        Read and apply the box text field.

        This routine reads the box text field and makes the new set of boxes
        and box positions.  It then applies these and closes the box window.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        The code does, however, change the self.plot_boxes values as needed
        to match what is in the box text field.

        """
        boxtext = box_message_text.get(0.0, Tk.END)
        lines = boxtext.split('\n')
        newboxes = []
        nboxes = 0
        for line in lines:
            values = line.split('\t')
            if len(values) == 10:
                try:
                    x1 = float(values[0])
                    y1 = float(values[1])
                    x2 = float(values[2])
                    y2 = float(values[3])
                    theta = float(values[4])
                    nplot = int(values[5])
                    colour = values[6].strip(' ')
                    linetype = values[7].strip(' ')
                    thickness = float(values[8])
                    fillcolour = values[9].strip(' ')
                    newboxes.append({'xstart': x1, 'ystart': y1,
                                     'xend': x2, 'yend': y2,
                                     'rotation': theta, 'plot': nplot,
                                     'line_type': linetype,
                                     'line_colour': colour,
                                     'line_thickness': thickness,
                                     'fill_colour': fillcolour})
                    nboxes = nboxes + 1
                except:
                    pass
        box_window.destroy()
        if nboxes > self.max_boxes:
            self.max_boxes = nboxes
        else:
            for loop in range(nboxes, self.max_boxes):
                newboxes.append({'xstart': None, 'ystart': None,
                                 'xend': None, 'yend': None,
                                 'rotation': 0.0, 'plot': 1,
                                 'line_type': 'solid',
                                 'line_colour': 'black',
                                 'line_thickness': 1.0,
                                 'fill_colour': 'none'})
        self.plot_boxes = newboxes
        self.number_of_boxes = nboxes
        self.make_plot()

    def edit_ellipses(self):
        """
        Make a window for editing the ellipse values.

        This routine produces a text box in a window, within which one can
        edit the ellipse values.  If no ellipses are defined the routine
        just exits with no action.

        Ellipse values are presented one per line, with the center x position,
        the center y position, the major axis length (x), the minor axis
        length (y), the orientation, the plot number, the line type, the
        line colour, the line thickness, and the fill colour separated by
        tab symbols.  One can edit the values within the text window and
        then these are applied when one clicks on the "Close Window" button.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        """
        if self.number_of_ellipses == 0:
            return
        str1 = 'Edit values below: fields are separated by tab characters.\n' \
               + '     center x       center y       major      minor      ' \
               + 'orient.   plot        colour   type  thickness  fill \n--' \
               + '---------------------------------------------------------' \
               + '----------------------------------------------------\n'
        for loop in range(self.number_of_ellipses):
            str1 = str1 + '%12.6g\t%12.6g\t%12.6g\t%12.6g\t%8.3f\t%6d' % (
                self.plot_ellipses[loop]['xposition'],
                self.plot_ellipses[loop]['yposition'],
                self.plot_ellipses[loop]['major'],
                self.plot_ellipses[loop]['minor'],
                self.plot_ellipses[loop]['rotation'],
                self.plot_ellipses[loop]['plot'])
            str1 = str1 + '\t%15s\t%8s\t%7.3f\t%15s\n' % (
                self.plot_ellipses[loop]['line_colour'],
                self.plot_ellipses[loop]['line_type'],
                self.plot_ellipses[loop]['line_thickness'],
                self.plot_ellipses[loop]['fill_colour'])
        ellipse_window = Tk.Toplevel()
        ellipse_window.title('Ellipses:')
        holder = Tk.Frame(ellipse_window)
        holder.pack(side=Tk.TOP)
        ellipse_message_text = ScrolledText(holder, height=40, width=115,
                                            wrap=Tk.NONE)
        ellipse_message_text.config(font=('courier', 16, 'bold'))
        ellipse_message_text.pack(side=Tk.TOP)
        ellipse_message_text.insert(0.0, str1)
        bholder = Tk.Frame(ellipse_window)
        bholder.pack(side=Tk.TOP)
        close_button = Tk.Button(
            bholder, text='Close Window',
            command=lambda: self.read_ellipses(ellipse_message_text,
                                               ellipse_window))
        close_button.pack()

    def read_ellipses(self, ellipse_message_text, ellipse_window):
        """
        Read and apply the ellipse text field.

        This routine reads the ellipse text field and makes the new set of
        ellipses and ellipse positions.  It then applies these and closes
        the ellipse window.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        The code does, however, change the self.plot_ellipses values as
        needed to match what is in the ellipse text field.

        """
        ellipsetext = ellipse_message_text.get(0.0, Tk.END)
        lines = ellipsetext.split('\n')
        newellipses = []
        nellipses = 0
        for line in lines:
            values = line.split('\t')
            if len(values) == 10:
                try:
                    x1 = float(values[0])
                    y1 = float(values[1])
                    x2 = float(values[2])
                    y2 = float(values[3])
                    theta = float(values[4])
                    nplot = int(values[5])
                    colour = values[6].strip(' ')
                    linetype = values[7].strip(' ')
                    thickness = float(values[8])
                    fillcolour = values[9].strip(' ')
                    newellipses.append({'xposition': x1, 'yposition': y1,
                                        'major': x2, 'minor': y2,
                                        'rotation': theta, 'plot': nplot,
                                        'line_type': linetype,
                                        'line_colour': colour,
                                        'line_thickness': thickness,
                                        'fill_colour': fillcolour})
                    nellipses = nellipses + 1
                except:
                    pass
        ellipse_window.destroy()
        if nellipses > self.max_ellipses:
            self.max_ellipses = nellipses
        else:
            for loop in range(nellipses, self.max_ellipses):
                newellipses.append({'xposition': None, 'yposition': None,
                                    'major': None, 'minor': None,
                                    'rotation': 0.0, 'plot': 1,
                                    'line_type': 'solid',
                                    'line_colour': 'black',
                                    'line_thickness': 1.0,
                                    'fill_colour': 'none'})
        self.plot_ellipses = newellipses
        self.number_of_ellipses = nellipses
        self.make_plot()

    def edit_vectors(self):
        """
        Make a window for editing the vector values.

        This routine produces a text box in a window, within which one
        can edit the vector values.  If no vectors are defined the routine
        just exits with no action.

        Vectors are presented one per line, with the start x position, the
        start y position, the end x position, the end y position, the head
        width, the head length, the plot number, the line type, the line
        colour, and the line thickness separated by tab symbols.  One can
        edit the values within the text window and then these are applied
        when one clicks on the "Close Window" button.

        Parameters:     None

        Return values:  None

        """
        if self.number_of_vectors == 0:
            return
        str1 = 'Edit values below: fields are separated by tab characters.\n' \
               + '     start x       start y       end x         end y   ' \
               + 'head width head length   plot     line colour   line type' \
               + '  thickness head fill head colour\n----------------------' \
               + '---------------------------------------------------------' \
               + '---------------------------------------------------------' \
               + '------\n'
        for loop in range(self.number_of_vectors):
            if self.plot_vectors[loop]['fill']:
                flag1 = 'True'
            else:
                flag1 = 'False'
            str1 = str1 + '%12.6g\t%12.6g\t%12.6g\t%12.6g\t%12.6g' % (
                self.plot_vectors[loop]['xstart'],
                self.plot_vectors[loop]['ystart'],
                self.plot_vectors[loop]['xend'],
                self.plot_vectors[loop]['yend'],
                self.plot_vectors[loop]['delx'])
            str1 = str1 + '\t%12.6g\t%6d\t%15s\t%8s\t%10.3f' % (
                self.plot_vectors[loop]['dely'],
                self.plot_vectors[loop]['plot'],
                self.plot_vectors[loop]['line_colour'],
                self.plot_vectors[loop]['line_type'],
                self.plot_vectors[loop]['line_thickness'])
            str1 = str1 + '\t     %s\t    %s\n' % (
                flag1,
                self.plot_vectors[loop]['fill_colour'])
        vector_window = Tk.Toplevel()
        vector_window.title('Vectors:')
        holder = Tk.Frame(vector_window)
        holder.pack(side=Tk.TOP)
        vector_message_text = ScrolledText(holder, height=40, width=145,
                                           wrap=Tk.NONE)
        vector_message_text.config(font=('courier', 16, 'bold'))
        vector_message_text.pack(side=Tk.TOP)
        vector_message_text.insert(0.0, str1)
        bholder = Tk.Frame(vector_window)
        bholder.pack(side=Tk.TOP)
        close_button = Tk.Button(
            bholder, text='Close Window',
            command=lambda: self.read_vectors(vector_message_text,
                                              vector_window))
        close_button.pack()

    def read_vectors(self, vector_message_text, vector_window):
        """
        Read and apply the vectors text field.

        This routine reads the vector text field and makes the new set of
        vectors and vector positions.  It then applies these and closes
        the vector window.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing

        The code does, however, change the self.plot_vectors values as
        needed to match what is in the vector text field.

        """
        vectortext = vector_message_text.get(0.0, Tk.END)
        vectors = vectortext.split('\n')
        newvectors = []
        nvectors = 0
        for vector in vectors:
            values = vector.split('\t')
            if len(values) == 12:
                try:
                    x1 = float(values[0])
                    y1 = float(values[1])
                    x2 = float(values[2])
                    y2 = float(values[3])
                    delx = float(values[4])
                    dely = float(values[5])
                    nplot = int(values[6])
                    colour = values[7].strip(' ')
                    linetype = values[8].strip(' ')
                    thickness = float(values[9])
                    if 'true' in values[10].lower():
                        flag = True
                    else:
                        flag = False
                    hcolour = values[11].strip(' ')
                    newvectors.append({'xstart': x1, 'ystart': y1,
                                       'xend': x2, 'yend': y2,
                                       'delx': delx, 'dely': dely,
                                       'plot': nplot, 'line_type': linetype,
                                       'line_colour': colour,
                                       'line_thickness': thickness,
                                       'fill': flag, 'fill_colour': hcolour})
                    nvectors = nvectors + 1
                except:
                    pass
        vector_window.destroy()
        if nvectors > self.max_vectors:
            self.max_vectors = nvectors
        else:
            for loop in range(nvectors, self.max_vectors):
                newvectors.append({'xstart': None, 'ystart': None,
                                   'xend': None, 'yend': None,
                                   'plot': 1, 'line_type': 'solid',
                                   'line_colour': 'black',
                                   'line_thickness': 1.0, 'fill': True,
                                   'fill_colour': 'black'})
        self.plot_vectors = newvectors
        self.number_of_vectors = nvectors
        self.make_plot()


if __name__ == "__main__":
    #    This is the main program.  Assuming that the program is called from
    #    the command line one defines the window and starts the plot widget.
    #
    #    One can run this code from within Python as well.  See the main help
    #    text for an example.
    for arg in sys.argv:
        if '--fast' in arg:
            import matplotlib.style as mplstyle
            mplstyle.use('fast')
    root, gui_window = startup()
    root.mainloop()
