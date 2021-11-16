#! /usr/bin/env python
#
"""
Matplotlib_user_interface.py is an interactive plotting tool using matplotlib.

This is intended to be a front-end for a number of MATPLOTLIB functions similar
in general usage to "xmgrace" (although not intending to duplicate all the
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

The code also includes some simple FITS image display functionality to have
this in the same file for reference, but that is independent of the main
plotting functionality.

The code assumes that the "times new roman" font is installed along with the
default matplotlib fonts.  If this is not the case, one will get a fall-back
font instead (usually "DejaVe Sans").  One can install the Times New Roman
font if this is not already on the system.  Due to preferences of the author
the Times New Roman font is the default.  If one wishes to change this, search
for 'times new roman' and replace the string by the default that is wanted, 
say 'sans-serif' for example.  There are commented out lines to make 
'sans-serif' the default font, which can be uncommented and used to replace 
the instances of 'times new roman' font as the default, should that be needed.

************************************************************

Use from inside Python:

  While the code was written to read data from files as a stand-alone
interface, one can also call the routines from inside python and add sets
from the python prompt.  The commands needed are as follows, assuming that
data arrays "xvalues" and "yvalues" already exist within the Python session:

>>>> import matplotlib_user_interface as mui
>>>> root, myplot = mui.startup()
>>>> myplot.add_set(xvalues, yvalues)

(repeat for multiple sets, if needed)

>>>> myplot.show_plot()

If one wishes to do all the steps manually the following somewhat longer
set of commands can be used:

(1) import the required packages tkinter and matplotlib_user_interface

>>> import tkinter as Tk
>>> import matplotlib_user_interface as mui

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

>>> myplot = mui.PlotGUI(root)

This fills in the buttons/functions in the window.  The functionality is
then available.

(4) Add a set or sets to the plot.  Any number of sets up a limit (currently
100) can be added.  See the add_set subroutine doc string for the possible
use of error bars in the set.  The value self.max_sets determines the maximum
number of sets that can be handled.

>>> myplot.add_set(xvalues, yvalues)

(5) Tell the code to plot or replot.  This call can be used as many times
as required.

>>> myplot.show_plot()

The add_set routine takes as input the x and y data values and makes a new
set.  Hence one can mix reading in values and generating values inside python
as needed by following the above set of steps.  There are additional
parameter values one can pass to this routine, for error values in the data
points.

When running in a script rather than at the python interpreter, one needs
to add

root.mainloop()

after making the call to show_plot with the sets, to keep the window in place.
Otherwise the window will disappear after the show_plot call if nothing else
is being done in the script.

************************************************************

Using the image display functionality

  In a similar way one can use the image display functionality from the python
command line.  There is a utility routine "showimage" for this purpose.

>>> from astropy.io import fits
>>>> imagefilename = 'test.fits'
>>> image = fits.getdata(imagefilename)
>>> import matplotlib_user_interface as mui
>>> root, myplot = mui.showimage(image, imagefilename, dpi=300)

The file name is optional.  The number of dots per inch is also optional.
The default number of dots per inch is 100.

If one wishes to make the calls explicitly they are (parameter self.indpi is
the number of dots per inch):

>>>> root = Tk.Tk()
>>>> plotobject = mui.ImageGUI(root)
>>>> plotobject.image = image
>>>> plotobject.imagefilename = 'some_name.fits'
>>>> plotobject.indpi = 300
>>>> plotobject.make_image_window()

As above, if done in a script one may need to add a call

>>>> root.mainloop()

to keep the window active.

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
import astropy.io.fits as fits
import general_utilities
import object_utilities
import label_utilities
import make_plot
import edit_objects
import save_and_restore_plot
import fits_image_display
import histogram_utilities
import data_set_utilities
import data_set_operations
import plot_flag_utilities
import window_creation
import plot_controls
import non_linear_fitting

# The following are "global" variables with line/marker information from
# matplotlib.  These are used, but not changed, in the code in more than one
# place hence I am using single variables here for the values.
matplotlib.use('TkAgg')
# define a background colour for windows
BGCOL = '#F8F8FF'


def startup():
    """
    Startup.py is a wrapper for starting the plot tool.

    This is a wrapper for making a plot from the Python command prompt.
    Assuming that numpy data arrays xvalues and yvalues exist one can use the
    commands

    >>> import matplotlib_user_interface as mui
    >>> root, myplot = mui.startup()
    >>> myplot.add_set(xvalues, yvalues)
    >>> myplot.show_plot()

    to bring up the window and plot the data set from inside Python.  One can
    then add other data sets and use the functionality of the tool.

    Parameters
    ----------
        None.

    Returns
    -------
        newroot :   The Tkinter window class variable for the plot window.

        plotobject :  The matplotlib_user_interface plot GUI object variable.  
                      This is the variable one uses for making plots.

    """
    # Make a Tkinter window
    newroot = Tk.Tk()
    newroot.title("Plotting Tool")
    plotobject = PlotGUI(newroot)
    return newroot, plotobject


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
        self.mpfit_values = []
        # The title string for the plot is kept in self.title.
        object_utilities.initialize_plot_objects(self)
        self.title = [' ', ]
        self.positions = []
        self.datavalues = None
        self.legend_labels = [None, ]
        self.legend_handles = [None, ]
        self.legend_variable = [None, ]
        self.legend_frame = [None, ]
        self.legend_options = [None, ]
        self.legend_position = [None, ]
        self.legend_user_position = [None, ]
        self.grid_colour = [None, ]
        self.grid_colour_variable = None
        self.grid_linetype = [None, ]
        self.grid_linetype_variable = None
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
            'bothticks': 0, 'minorticks': 0, 'oppositeaxis': 0,
            'majorgridlines': 0, 'minorgridlines': 0}, ]
        self.yparameters = [{
            'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
            'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
            'invert': 0, 'hide': 0, 'hideticks': 0,
            'hidelabels': 0, 'hybridlog': 0,
            'inverseticks': 0, 'ticklength': 6,
            'bothticks': 0, 'minorticks': 0, 'oppositeaxis': 0,
            'majorgridlines': 0, 'minorgridlines': 0}, ]
#        self.fontname = ['sans-serif', ]
        self.fontname = ['times new roman', ]
        self.fontsize = ['12', ]
        self.fontweight = ['normal', ]
        # The following are window names, set to None before and after the
        # window is active to avoid duplicate windows.
#        self.font_window = None
#        self.label_font_window = None
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
        self.non_linear_set_fitting_window = None
        self.plot_margin = 0.
        self.plot_frame = [0., ]
        self.nxplots = 1
        self.nyplots = 1
        self.current_plot = 1
        self.number_of_plots = self.nxplots * self.nyplots
        self.plot_area_flag = True
        self.zoom_flag = False
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
        window_creation.make_menus(self, menuframe)
        controlframe = Tk.Frame(self.root)
        controlframe.pack(side=Tk.LEFT, fill=Tk.Y, expand=1)
        window_creation.make_controls(self, controlframe)
        sl = general_utilities.separator_line(
            self.root, 5, 750, 5, False, Tk.LEFT)
        self.plotframe = Tk.Frame(self.root)
        self.plotframe.pack(side=Tk.LEFT, fill=Tk.Y, expand=1)
        window_creation.make_plot_area(self, self.plotframe)

    def show_plot(self):
        """
        Wrapper routine for displaying the plot.

        Parameters
        ----------

        None

        Returns
        -------

        None
        """
        make_plot.make_plot(self)
        
    def read_image(self):
        """
        Read in a FITS 2-D image for display.

        This routine allows the user to read in an image.  The image is then
        displayed in a new window.

        No values are passed to this routine or returned from this routine.
        """
        self.imagefilename = tkinter.filedialog.askopenfilename()
        try:
            self.image = fits.getdata(self.imagefilename)
        except:
            try:
                self.image = fits.getdata(self.imagefilename, ext=1)
            except:
                tkinter.messagebox.showinfo(
                    "Error", "File " + self.imagefilename
                    + " could not be read in.")
                self.imagefilename = None
                return
        sh1 = self.image.shape
        if len(sh1) < 2:
            tkinter.messagebox.showinfo(
                "Error", "The image from file " + self.imagefilename
                + " is not of dimension 2 or higher.")
            self.imagefilename = None
            return
        if len(sh1) == 3:
            self.image = numpy.squeeze(self.image[0, :, :])
        if len(sh1) == 4:
            self.image = numpy.squeeze(self.image[0, 0, :, :])
        if len(sh1) == 5:
            self.image = numpy.squeeze(self.image[0, 0, 0, :, :])
        if len(sh1) == 6:
            self.image = numpy.squeeze(self.image[0, 0, 0, 0, :, :])
        image_object = fits_image_display.ImageGUI()
        image_object.image = self.image
        image_object.imagefilename = self.imagefilename
        image_object.make_image_window()

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
            xmin1 = general_utilities.round_float(self, xmin, True)
            xmax1 = general_utilities.round_float(self, xmax, False)
            ymin1 = general_utilities.round_float(self, ymin, True)
            ymax1 = general_utilities.round_float(self, ymax, False)
            self.plot_range[self.current_plot-1][0] = xmin1
            self.plot_range[self.current_plot-1][1] = xmax1
            self.plot_range[self.current_plot-1][2] = ymin1
            self.plot_range[self.current_plot-1][3] = ymax1
        make_plot.make_plot(self)

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
        xpos = event.xdata
        ypos = event.ydata
        for loop in range(self.number_of_plots):
            if event.inaxes == self.subplot[loop]:
                self.current_plot = loop + 1
        return

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
                xout = general_utilities.inverse_hybrid_transform(event.xdata)
            if hlogyflag == 0:
                yout = event.ydata
            else:
                yout = general_utilities.inverse_hybrid_transform(event.ydata)
            nset, npoint, xpoint, ypoint, distance = self.match_point(
                xout, yout)
            s1 = "Position: %.6g %.6g" % (xout, yout)
            if nset is not None:
                s1 = s1 + '\nNearest point: set %d point %d \n' % \
                     (nset+1, npoint) \
                     + 'data (%.6g, %.6g) distance %.6g' % (
                         xpoint, ypoint, distance)
            self.position_label_text.set(s1)
        except ValueError:
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
                except IndexError:
                    dminhere = dminset
                if (dmin < 0.) or (dminhere < dmin):
                    nset = loop
                    try:
                        ndata = dargmin[0]
                    except IndexError:
                        ndata = dargmin
                    dmin = dminhere
                    xmin = xdata[ndata]
                    ymin = ydata[ndata]
            if dmin >= 0.:
                return nset, ndata, xmin, ymin, dmin
        except Exception:
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
            xvalue = general_utilities.inverse_hybrid_transform(xvalue)
        if self.yparameters[self.current_plot-1]['hybridlog'] == 1:
            yvalue = general_utilities.inverse_hybrid_transform(yvalue)
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
                label_utilities.set_label_properties(
                    self, self.number_of_labels)
                self.number_of_labels = self.number_of_labels + 1
                make_plot.make_plot(self)
            except Exception:
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
        self.zoom_flag = True
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
            xvalue = general_utilities.inverse_hybrid_transform(xvalue)
        if self.yparameters[self.current_plot-1]['hybridlog'] == 1:
            yvalue = general_utilities.inverse_hybrid_transform(yvalue)
        position = [xvalue, yvalue]
        if self.line_flag:
            self.positions.append(position)
            object_utilities.add_line_values(self)
        if self.vector_flag:
            self.positions.append(position)
            object_utilities.add_vector_values(self)
        if self.ellipse_flag:
            self.positions.append(position)
            object_utilities.add_ellipse_values(self)
        if self.box_flag:
            self.positions.append(position)
            object_utilities.add_box_values(self)
        if self.zoom_flag:
            self.positions.append(position)
            self.zoom_flag = False
            if xvalue > self.positions[-2][0]:
                xmin = self.positions[-2][0]
                xmax = xvalue
            else:
                xmin = xvalue
                xmax = self.positions[-2][0]
            if yvalue > self.positions[-2][1]:
                ymin = self.positions[-2][1]
                ymax = yvalue
            else:
                ymin = yvalue
                ymax = self.positions[-2][1]
            xratio = abs(xmax-xmin)/(self.plot_range[self.current_plot-1][1]-
                                     self.plot_range[self.current_plot-1][0])
            yratio = abs(ymax-ymin)/(self.plot_range[self.current_plot-1][3]-
                                  self.plot_range[self.current_plot-1][2])
            if (xmin != xmax) and (ymin != ymax) and \
               (xratio > 0.01) and (yratio > 0.01):
#                print(self.plot_range[self.current_plot-1])
#                print(xmin, xmax, ymin, ymax)
                self.original_range[self.current_plot-1] = False
                self.plot_range[self.current_plot-1] = [xmin, xmax, ymin, ymax]
        make_plot.make_plot(self)

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
            s1 = 'Position: [%g, %g]' % (event.xdata, event.ydata)
            self.histogramLabelText.set(s1)
        except ValueError:
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
            s1 = 'Position: [%g, %g] Value: %d' % (event.xdata, event.ydata,
                                                   self.hist2d[ix, iy])
            self.hessLabelText.set(s1)
        except (ValueError, TypeError):
            pass

    def close_window(self, windowname, label, clear_data_input=False):
        """
        Close a window and set the associated variable to None.

        This routine removes a window from the systen and replaces the
        old window variable with None to mark that the window is not active.

        Parameters
        ----------
           windowname :  A Tk.Toplevel variable for an existing window
                         that is to be closed.

           label :  A string identifying which window is being closed

           clear_data_input :   A boolean flag for whether the datavalues 
                                variable needs to be cleared.

        Returns
        -------
            No value is returned by this routine.

        """
        windowname.destroy()
        windowname.update()
        windowname = None
        if clear_data_input:
            self.datavalues = None
        # It is not clear why this cannot be put in a loop, but when it is
        # attempted the variables are not set to None.
        if label == 'read_window':
            self.read_window = None
        if label == 'data_set_window':
            self.data_set_window = None
        if label == 'data_set_transformation_window':
            self.data_set_transformation_window = None
        if label == 'data_set_delete_window':
            self.data_set_delete_window = None
        if label == 'data_entry_window':
            self.data_entry_window = None
        if label == 'data_set_sort_window':
            self.data_set_sort_window = None
        if label == 'data_set_fitting_window':
            self.data_set_fitting_window = None
        if label == 'plot_control_window':
            self.plot_control_window = None
        if label == 'box_window':
            self.box_window = None
        if label == 'vector_window':
            self.vector_window = None
        if label == 'ellipse_window':
            self.ellipse_window = None
        if label == 'line_window':
            self.line_window = None
        if label == 'setplot_window':
            self.setplot_window = None
        if label == 'hideplot_window':
            self.hideplot_window = None
        if label == 'tile_window':
            self.tile_window = None
        if label == 'non_linear_set_fitting_window':
            self.non_linear_set_fitting_window = None

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
            None

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
        self.grid_colour.append(None)
        self.grid_linetype.append(None)
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
        xmin1 = general_utilities.round_float(self, xmin, True)
        xmax1 = general_utilities.round_float(self, xmax, False)
        ymin1 = general_utilities.round_float(self, ymin, True)
        ymax1 = general_utilities.round_float(self, ymax, False)
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
        matplotlib_symbol_list = [None, ".", ", ", "o", "v", "^", "<", ">",
                                  "1", "2", "3", "4", "8", "s", "p", "P",
                                  "*", "h", "H", "+", "x", "X", "D", "d",
                                  "|", "_", "histogram"]
        matplotlib_symbol_name_list = ['None', 'point', 'pixel', 'circle',
                               'triangle down', 'triangle up',
                               'triangle left', 'triangle right', 'tri_down',
                               'tri_up', 'tri_left', 'tri_right', 'octagon',
                               'square', 'pentagon', 'plus (filled)', 'star',
                               'hexagon1', 'hexagon2', 'plus', 'x',
                               'x (filled)', 'diamond', 'thin_diamond',
                               'vline', 'hline', 'histogram']
        matplotlib_line_list = ['-', '--', '-.', ':', None]
        matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                     'dotted', 'None']
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
            except ValueError:
                self.set_entry_fields[2].current(4)
            self.set_entry_fields[3].delete(0, Tk.END)
            self.set_entry_fields[3].insert(
                0, str(self.set_properties[index]['linewidth']))
            try:
                self.set_entry_fields[4].set(
                    self.set_properties[index]['colour'])
            except Exception:
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
            make_plot.make_plot(self)
        except Exception:
            pass

    def set_grid_colour(self, event):
        """
        Set the grid line colour for the current plot.

        Parameters
        ----------

        event   A tkinter event variable (not used...)

        Returns
        -------

        None
        """
        colour = self.grid_colour_variable.get()
        self.grid_colour[self.current_plot-1] = colour
        
    def set_grid_linetype(self, event):
        """
        Set the grid line colour for the current plot.

        Parameters
        ----------

        event   A tkinter event variable (not used...)

        Returns
        -------

        None
        """
        linetype = self.grid_linetype_variable.get()
        self.grid_linetype[self.current_plot-1] = linetype
        
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
     # Uncomment the following command to set the foreground/background colours
     #root.tk_setPalette(background='#f8f8ff', foreground='black',
     #    activeBackground='black', activeForeground='#f8f8ff')
    root.mainloop()
