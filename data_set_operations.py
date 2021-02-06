import math
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
from scipy.interpolate import UnivariateSpline, make_lsq_spline
import matplotlib
import matplotlib.lines as mlines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Ellipse, FancyArrow
from matplotlib.ticker import MultipleLocator
import general_utilities
import make_plot
import data_set_utilities

def make_data_set_edit_window(plotgui):
    """
    Create a new text window in which to edit a data set.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if plotgui.nsets == 0:
        return
    if plotgui.nsets > 1:
        str1 = "Which set do you want to edit (1 to %d)?" % plotgui.nsets
        nset = tkinter.simpledialog.askinteger(
            "Input", str1, parent=plotgui.root)
        if (nset is None) or (nset < 1) or (nset > plotgui.nsets):
            return
    else:
        nset = plotgui.nsets
    if plotgui.data_entry_window is not None:
        return
    plotgui.data_entry_window = Tk.Toplevel()
    plotgui.data_entry_window.title('Edit Data Set Values')
    holder = Tk.Frame(plotgui.data_entry_window)
    holder.pack(side=Tk.TOP)
    holder.config(bg='black')
    plotgui.data_text = ScrolledText(holder, height=40, width=80,
                                     wrap=Tk.NONE, relief="solid")
    plotgui.data_text.config(font=('courier', 16))
    plotgui.data_text.pack(side=Tk.TOP, padx=10, pady=10)
    str1 = ''
    for loop in range(len(plotgui.xdata[nset-1]['values'])):
        str1 = str1 + '%g %g %g %g %g %g\n' % (
            plotgui.xdata[nset-1]['values'][loop],
            plotgui.xdata[nset-1]['lowerror'][loop],
            plotgui.xdata[nset-1]['higherror'][loop],
            plotgui.ydata[nset-1]['values'][loop],
            plotgui.ydata[nset-1]['lowerror'][loop],
            plotgui.ydata[nset-1]['higherror'][loop]
        )
    plotgui.data_text.insert(1.0, str1)
    bframe = Tk.Frame(plotgui.data_entry_window)
    bframe.pack()
    set_button = Tk.Button(
        bframe, text="Apply",
        command=lambda: apply_data_edits(plotgui, nset, plotgui.data_text))
    set_button.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        bframe, text="Close",
        command=lambda: plotgui.close_data_window(
            plotgui.data_entry_window))
    close_button.pack(side=Tk.LEFT)

def apply_data_edits(plotgui, nset, data_text):
    """
    Read the data edit window and apply the values.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

        nset:      an integer variable for the set number that is being edited

        data_text:  a tkinter text field variable which is read for the values

    Returns
    -------

       None

    """
    text = data_text.get("1.0", Tk.END)
    xvalues, dxvalues1, dxvalues2, yvalues, dyvalues1, dyvalues2,\
        errorflag = general_utilities.parse_data_input_text(text)
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
        plotgui.xdata[nset-1]['values'] = xvalues
        plotgui.xdata[nset-1]['lowerror'] = dxvalues1
        plotgui.xdata[nset-1]['higherror'] = dxvalues2
        plotgui.ydata[nset-1]['values'] = yvalues
        plotgui.ydata[nset-1]['lowerror'] = dyvalues1
        plotgui.ydata[nset-1]['higherror'] = dyvalues2
        make_plot.make_plot(plotgui)
    except ValueError:
        tkinter.messagebox.showinfo(
            'error',
            'Unable to parse text from edit widget')
        return

def make_data_set_sort_window(plotgui):
    """
    Create the data set sorting window.

    This routine brings up a window within which one can sort a data set.
    The sorted data values replace the current data in the set.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if plotgui.data_set_sort_window is not None:
        return
    if plotgui.nsets == 0:
        return
    plotgui.data_set_sort_window = Tk.Toplevel()
    plotgui.data_set_sort_window.title('Data Sets Sort Window')
    holder = Tk.Frame(plotgui.data_set_sort_window)
    holder.pack(side=Tk.TOP)
    plotgui.set_sort_list_area = tkinter.ttk.Combobox(holder, width=50,
                                                      height=10)
    plotgui.set_sort_list_area.pack(side=Tk.TOP)
    setlist = []
    for loop in range(plotgui.nsets):
        label = (
            'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
            (loop+1, len(plotgui.xdata[loop]['values']),
             plotgui.xdata[loop]['minimum'], plotgui.xdata[loop]['maximum']))
        label = (
            label + 'y range %13.6g to %13.6g' %
            (plotgui.ydata[loop]['minimum'], plotgui.ydata[loop]['maximum'])
        )
        setlist.append(label)
    plotgui.set_sort_list_area['values'] = setlist
    plotgui.set_sort_list_area.current(0)
    frame1 = Tk.Frame(holder)
    frame1.pack()
    plotgui.sortflag1 = Tk.IntVar()
    plotgui.sortflag2 = Tk.IntVar()
    lb = Tk.Label(frame1, text="Sort on: ")
    lb.pack(side=Tk.LEFT)
    b1 = Tk.Frame(frame1)
    b1.pack()
    general_utilities.put_yes_no(b1, plotgui.sortflag1, ['x values', 'y values'], True)
    lb = Tk.Label(frame1, text="Sort order: ")
    lb.pack(side=Tk.LEFT)
    b1 = Tk.Frame(frame1)
    b1.pack()
    general_utilities.put_yes_no(b1, plotgui.sortflag2, ['Ascending', 'Decending'], True)
    frame2 = Tk.Frame(holder)
    frame2.pack(side=Tk.TOP)
    sort_button = Tk.Button(frame2, text="Sort",
                            command=lambda: sort_set(plotgui))
    sort_button.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        frame2, text="Close",
        command=lambda: plotgui.close_data_window(plotgui.data_set_sort_window))
    close_button.pack(side=Tk.LEFT)

def sort_set(plotgui):
    """
    Sort the set values (x, y) according to the options.

    This routine sorts a data set either on the x values or on the y
    values.  Sorting can be done in ascending or decending order.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    try:
        set_number = plotgui.set_sort_list_area.current()
        # The following should not happen....but check just to be sure.
        if (set_number < 0) or (set_number >= plotgui.nsets):
            tkinter.messagebox.showinfo('Error', 'The set was not sorted.')
            return
        option1 = plotgui.sortflag1.get()
        option2 = plotgui.sortflag2.get()
        flag = do_sort_set(plotgui, set_number, option1, option2)
        if not flag:
            tkinter.messagebox.showinfo('Error', 'The set was not sorted.')
    except ValueError:
        tkinter.messagebox.showinfo('Error', 'The set was not sorted.')

def do_sort_set(plotgui, set_number, xyoption=1, direction_option=0):
    """
    Do the actual sorting of the set values.

    This routine carries out the sorting of values in a set.  It is
    provided sepaerately so it can be used from the Python command
    line if needed.

    Parameters
    ----------
        plotgui:   by assumption a matplotlib_user_interface object

        set_number : An integer value 0 or higher for the set to be sorted

        xyoption :   An optional integer value, 1 for sorting in x and
                     0 for sorting in y (defaults to sorting in x)

        direction_option :   An optional integer value, 1 for sorting in
                             ascending order and 0 for sorting in
                             descending order (defaults to ascending order)

    Values of 1 are the default for both options.  Actually any xyoption
    value but 1 will cause sorting in y and any direction_option value
    but zero will cause sorting in ascending order.

    Returns
    -------
        flag :    A boolean value, True if the sorting was successful,
                  False otherwise

    """
    try:
        xvalues = numpy.copy(plotgui.xdata[set_number]['values'])
        xlowerror = numpy.copy(plotgui.xdata[set_number]['lowerror'])
        xhigherror = numpy.copy(plotgui.xdata[set_number]['higherror'])
        yvalues = numpy.copy(plotgui.ydata[set_number]['values'])
        ylowerror = numpy.copy(plotgui.ydata[set_number]['lowerror'])
        yhigherror = numpy.copy(plotgui.ydata[set_number]['higherror'])
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
        plotgui.xdata[set_number]['values'] = xvalues
        plotgui.xdata[set_number]['lowerror'] = xlowerror
        plotgui.xdata[set_number]['higherror'] = xhigherror
        plotgui.ydata[set_number]['values'] = yvalues
        plotgui.ydata[set_number]['lowerror'] = ylowerror
        plotgui.ydata[set_number]['higherror'] = yhigherror
        make_plot.make_plot(plotgui)
        return True
    except Exception:
        return False

def make_data_set_fitting_window(plotgui):
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

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if plotgui.data_set_fitting_window is not None:
        return
    if plotgui.nsets == 0:
        return
    plotgui.data_set_fitting_window = Tk.Toplevel()
    plotgui.data_set_fitting_window.title('Data Sets Fitting Window')
    holder = Tk.Frame(plotgui.data_set_fitting_window)
    holder.pack(side=Tk.TOP)
    plotgui.fit_text = ScrolledText(holder, height=6, width=50, bd=1,
                                    relief=Tk.RIDGE, wrap=Tk.NONE)
    plotgui.fit_text.config(font=('courier', 12, 'bold'))
    plotgui.fit_text.pack()
    plotgui.set_fitting_list_area = tkinter.ttk.Combobox(holder, width=50,
                                                         height=10)
    plotgui.set_fitting_list_area.pack(side=Tk.TOP)
    setlist = []
    for loop in range(plotgui.nsets):
        label = (
            'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
            (loop+1, len(plotgui.xdata[loop]['values']),
             plotgui.xdata[loop]['minimum'], plotgui.xdata[loop]['maximum'])
            + 'y range %13.6g to %13.6g' %
            (plotgui.ydata[loop]['minimum'], plotgui.ydata[loop]['maximum']))
        setlist.append(label)
    plotgui.set_fitting_list_area['values'] = setlist
    plotgui.set_fitting_list_area.current(0)
    nsets_start = plotgui.nsets
    frame1 = Tk.Frame(holder)
    frame1.pack()
    labels = ['Function', 'Order/Smoothing']
    for loop in range(len(labels)):
        label = Tk.Label(frame1, text=labels[loop])
        label.grid(column=0, row=loop, sticky=Tk.W)
    plotgui.set_fitting_fields = []
    plotgui.set_fitting_fields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.set_fitting_fields[-1].grid(column=1, row=0, sticky=Tk.W)
    plotgui.set_fitting_fields[-1]['values'] = [
        'Polynomial', 'Legendre',
        'Chebyshev', 'Laguerre', 'Spline', 'Cubic Spline',
        'Least-Squares Spline', 'Internal Linear Fit']
    plotgui.set_fitting_fields[-1].current(1)
    plotgui.set_fitting_fields.append(Tk.Entry(frame1, width=15, text='4'))
    plotgui.set_fitting_fields[-1].grid(column=1, row=1, sticky=Tk.W)
    plotgui.set_fitting_fields[-1].delete(0, Tk.END)
    plotgui.set_fitting_fields[-1].insert(0, '4')
    frame2 = Tk.Frame(holder)
    frame2.pack()
    plotgui.fit_option = Tk.IntVar()
    label = Tk.Label(frame2, text='Fit Type: ')
    label.pack(side=Tk.LEFT)
    Tk.Radiobutton(
        frame2, text='y = f(x)',
        variable=plotgui.fit_option, value=0).pack(side=Tk.LEFT)
    Tk.Radiobutton(
        frame2, text='x = f(y)',
        variable=plotgui.fit_option, value=1).pack(side=Tk.LEFT)
    plotgui.fit_option.set(0)
    frame2 = Tk.Frame(holder)
    frame2.pack()
    cancel_button = Tk.Button(
        frame2, text="Cancel Fittings",
        command=lambda: cancel_fitting_fields(plotgui, nsets_start))
    cancel_button.pack(side=Tk.LEFT)
    set_button = Tk.Button(
        frame2, text="Fit Set",
        command=lambda: apply_fitting_fields(plotgui))
    set_button.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        frame2, text="Close",
        command=lambda: plotgui.close_data_window(
            plotgui.data_set_fitting_window))
    close_button.pack(side=Tk.LEFT)

def cancel_fitting_fields(plotgui, nsets):
    """
    Clear the set fitting results.

    This routine will "clear" the fitting data sets.  This is done by
    reverting the plotgui.nsets value to what it was before the fitting
    was started.  Hence all the fitting data sets are still present but
    are not plotted, and if more data sets are read in the fitting sets
    will be overwritten.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    plotgui.nsets = nsets
    make_plot.make_plot(plotgui)

def apply_fitting_fields(plotgui):
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

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    fit_type = plotgui.set_fitting_fields[0].get()
    if 'Cubic Spline' in fit_type:
        fit_order = float(plotgui.set_fitting_fields[1].get())
    else:
        try:
            fit_order = int(plotgui.set_fitting_fields[1].get())
        except ValueError:
            str1 = 'Error: bad fit order (%s).  Settng to 4.' % (
                plotgui.set_fitting_fields[1].get())
            plotgui.fit_text.insert(Tk.END, str1)
            plotgui.fit_text.see(Tk.END)
            fit_order = 4
    set_number = plotgui.set_fitting_list_area.current()
    fit_flag = plotgui.fit_option.get()
    if fit_flag == 0:
        xvalues = numpy.copy(plotgui.xdata[set_number]['values'])
        yvalues = numpy.copy(plotgui.ydata[set_number]['values'])
    else:
        xvalues = numpy.copy(plotgui.ydata[set_number]['values'])
        yvalues = numpy.copy(plotgui.xdata[set_number]['values'])
    inds = numpy.argsort(xvalues)
    xvalues = xvalues[inds]
    yvalues = yvalues[inds]
    npoints = len(xvalues)
    if ('Spline' not in fit_type) and ('Internal' not in fit_type):
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
            yerrors = (plotgui.ydata[set_number]['lowerror'] +
                       plotgui.ydata[set_number]['higherror'])/2.
        else:
            yerrors = (plotgui.xdata[set_number]['lowerror'] +
                       plotgui.xdata[set_number]['higherror'])/2.
        if (numpy.min(yerrors) == 0.) and (numpy.max(yerrors) == 0.):
            yerrors = yerrors + 1.
        slope, intercept, slope_error, intercept_error, covariance, \
            correlation = general_utilities.slope_calculation(
                xvalues, yvalues, yerrors)
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
        str1 = str1 + 'Slope: %g +/- %g\n' % (slope, slope_error)
        str1 = str1 + 'Intercept: %g +/- %g\n' % (
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
        general_utilities.list_polynomial_fitpars(fit_type, fit_order, fitpars)
    if fit_type == 'Legendre':
        fitpars = legendre.legfit(xvalues, yvalues, fit_order)
        yout = legendre.legval(xout, fitpars)
        yfit = legendre.legval(xvalues, fitpars)
        labelstring = 'Order %d Legendre polynomial fit' % (fit_order)
        general_utilities.list_polynomial_fitpars(fit_type, fit_order, fitpars)
    if fit_type == 'Laguerre':
        fitpars = laguerre.lagfit(xvalues, yvalues, fit_order)
        yout = laguerre.lagval(xout, fitpars)
        yfit = laguerre.lagval(xvalues, fitpars)
        labelstring = 'Order %d Laguerre polynomial fit' % (fit_order)
        general_utilities.list_polynomial_fitpars(fit_type, fit_order, fitpars)
    if fit_type == 'Chebyshev':
        fitpars = chebyshev.chebfit(xvalues, yvalues, fit_order)
        yout = chebyshev.chebval(xout, fitpars)
        yfit = chebyshev.chebval(xvalues, fitpars)
        labelstring = 'Order %d Chebyshev polynomial fit' % (fit_order)
        general_utilities.list_polynomial_fitpars(fit_type, fit_order, fitpars)
    if fit_type == 'Least-Squares Spline':
        if fit_flag == 0:
            yerrors = (self.ydata[set_number]['lowerror'] +
                       self.ydata[set_number]['higherror'])/2.
        else:
            yerrors = (self.xdata[set_number]['lowerror'] +
                       self.xdata[set_number]['higherror'])/2.
        if (numpy.min(yerrors) == 0.) and (numpy.max(yerrors) == 0.):
            yerrors = yerrors + 1.
        xmin1 = numpy.min(xvalues)
        xmax1 = numpy.max(xvalues)
        xrange = xmax1 - xmin1
        nknots = int(fit_order)
        if (nknots < 3) or (nknots > int(len(xvalues)/2)):
            nknots = 3
        xstep = xrange/(nknots-2)
        xknots = numpy.arange(numpy.min(xvalues)+xstep,
                              numpy.max(xvalues)*0.999999999, xstep)
        k = 3
        # Use cubic splines
        knotedges = numpy.r_[(xmin,)*(k+1), xknots, (xmax,)*(k+1)]
        weights = 1./yerrors
        weights[yerrors == 0.] = 0.
        fitobject = make_lsq_spline(xvalues, yvalues, knotedges, k,
                                    w=weights)
        yout = fitobject(xout)
        yfit = fitobject(xvalues)
        labelstring = 'Least squares spline fit, sections = %d' % (nknots)
    if fit_type == 'Spline':
        fitpars = UnivariateSpline(
            xvalues, yvalues, k=1, s=None, bbox=[xmin - delx, xmax + delx])
        labelstring = 'Default spline fit'
        yout = fitpars(xout)
        yfit = fitpars(xvalues)
    if fit_type == 'Cubic Spline':
        if fit_order < 0.:
            str1 = 'Error: smoothing value %f (< 0) is not allowed.'\
                   + '  Settng to 0.0' % (fit_order)
            plotgui.fit_text.insert(Tk.END, str1)
            plotgui.fit_text.see(Tk.END)
            fit_order = 0.0
        fitpars = UnivariateSpline(
            xvalues, yvalues, k=3, bbox=[xmin - delx,
                                         xmax + delx], s=fit_order)
        yout = fitpars(xout)
        yfit = fitpars(xvalues)
        labelstring = 'Cubic spline fit, smoothing = %f' % (fit_order)
    rms = fit_statistics(yvalues, yfit)
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
        plotgui.fit_text.insert(Tk.END, str1)
        plotgui.fit_text.see(Tk.END)
    if fit_flag == 0:
        plotgui.xdata[plotgui.nsets] = {'values': xout, 'lowerror': xlowerror,
                                        'higherror': xhigherror,
                                        'minimum': xmin, 'maximum': xmax,
                                        'errors': False, 'legend': True}
        plotgui.ydata[plotgui.nsets] = {'values': yout, 'lowerror': ylowerror,
                                        'higherror': yhigherror,
                                        'minimum': ymin, 'maximum': ymax,
                                        'errors': False, 'legend': True}
    else:
        plotgui.xdata[plotgui.nsets] = {'values': yout, 'lowerror': ylowerror,
                                        'higherror': yhigherror,
                                        'minimum': ymin, 'maximum': ymax,
                                        'errors': False, 'legend': True}
        plotgui.ydata[plotgui.nsets] = {'values': xout, 'lowerror': xlowerror,
                                        'higherror': xhigherror,
                                        'minimum': xmin, 'maximum': xmax,
                                        'errors': False, 'legend': True}
    m = plotgui.nsets % 10
    n = int(math.floor(plotgui.nsets / 10))
    plotgui.set_properties[plotgui.nsets]['symbol'] = None
    plotgui.set_properties[plotgui.nsets]['linestyle'] = '-'
    plotgui.set_properties[plotgui.nsets]['colour'] = plotgui.colourset[m]
    plotgui.set_properties[plotgui.nsets]['symbolsize'] = 4.0 + 0.3*n
    plotgui.set_properties[plotgui.nsets]['label'] = labelstring
    plotgui.nsets = plotgui.nsets + 1
    make_plot.make_plot(plotgui)

def fit_statistics(yvalues, yfit):
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
        rms :       the root mean square deviation between the two arrays,
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

def make_data_set_transformation_window(plotgui):
    """
    Create the data set transformation window.

    This routine makes a window within which a data set transformation
    can be defined.  The user can either transform an existing set or
    use the transformation to make a new set.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if plotgui.data_set_transformation_window is not None:
        return
    plotgui.data_set_transformation_window = Tk.Toplevel()
    plotgui.data_set_transformation_window.title(
        'Data Set Transformation Window')
    holder = Tk.Frame(plotgui.data_set_transformation_window)
    holder.pack(side=Tk.TOP)
    plotgui.set_list_area1 = tkinter.ttk.Combobox(holder, width=50, height=10)
    plotgui.set_list_area1.pack(side=Tk.TOP)
    if plotgui.nsets == 0:
        setlist = [' ']
    else:
        setlist = []
        for loop in range(plotgui.nsets):
            label = (
                'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                (loop+1, len(plotgui.xdata[loop]['values']),
                 plotgui.xdata[loop]['minimum'],
                 plotgui.xdata[loop]['maximum'])
                 + ' y range %13.6g to %13.6g' %
                 (plotgui.ydata[loop]['minimum'],
                  plotgui.ydata[loop]['maximum']))
            setlist.append(label)
    plotgui.set_list_area1['values'] = setlist
    plotgui.set_list_area1.current(0)
#    plotgui.set_list_area1.bind('<<ComboboxSelected>>', plotgui.set_fields1)
    frame1 = Tk.Frame(holder)
    frame1.pack()
    label1 = Tk.Label(frame1, text='Set Transformation:')
    label1.grid(row=0, column=0, columnspan=2)
    label2 = Tk.Label(frame1, text='x transformation')
    label2.grid(row=1, column=0)
    plotgui.set_x_transformation_entry_field = Tk.Entry(frame1, width=30)
    plotgui.set_x_transformation_entry_field.grid(row=1, column=1)
    label2 = Tk.Label(frame1, text='y transformation')
    label2.grid(row=2, column=0)
    plotgui.set_y_transformation_entry_field = Tk.Entry(frame1, width=30)
    plotgui.set_y_transformation_entry_field.grid(row=2, column=1)
    frame1 = Tk.Frame(holder)
    frame1.pack()
    label1 = Tk.Label(frame1, text='Create new set?')
    label1.grid(row=0, column=0)
    plotgui.new_set = Tk.IntVar()
    b1 = Tk.Frame(frame1)
    general_utilities.put_yes_no(b1, plotgui.new_set, ['Yes', 'No'], True)
    b1.grid(column=1, row=0, sticky=Tk.W)
    label_str = 'Enter the function you wish to apply.  Use $x for the x values of the \ncurrent set and $y for the y values of the current set selected above.  \nOne can also use $i for an index starting at 1.'
    label1 = Tk.Label(holder, text=label_str)
    label1.pack(side=Tk.TOP)
    frame2 = Tk.Frame(holder)
    frame2.pack()
    apply_button = Tk.Button(frame2, text="Apply",
                             command=lambda: apply_transformation(plotgui))
    apply_button.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        frame2, text="Close",
        command=lambda: plotgui.close_data_window(
            plotgui.data_set_transformation_window))
    close_button.pack(side=Tk.LEFT)

def apply_transformation(plotgui):
    """
    Apply the defined set transformation.

    This routine reads the strings defining a data transformation
    and attempts to apply them.  It uses the x and y values of the
    current set as inputs.  One can use either math or numpy functions
    in the expression that is defined, but nothing else aside from the
    builtin functions.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    try:
        ind1 = plotgui.set_list_area1.current()
        xvalues = numpy.copy(plotgui.xdata[ind1]['values'])
        yvalues = numpy.copy(plotgui.ydata[ind1]['values'])
        seq = numpy.arange(len(xvalues)) + 1.
        xstr = plotgui.set_x_transformation_entry_field.get()
        xstr = xstr.replace('$x', 'x')
        xstr = xstr.replace('$y', 'y')
        xstr = xstr.replace('$i', 'seq')
        ystr = plotgui.set_y_transformation_entry_field.get()
        ystr = ystr.replace('$x', 'x')
        ystr = ystr.replace('$y', 'y')
        ystr = ystr.replace('$i', 'seq')
        x1 = data_set_utilities.my_eval(xstr, seq, xvalues, yvalues)
        y1 = data_set_utilities.my_eval(ystr, seq, xvalues, yvalues)
        if (x1 is None) or (y1 is None):
            return
        newopt = plotgui.new_set.get()
        if newopt == 0:
            plotgui.xdata[ind1]['values'] = x1
            plotgui.ydata[ind1]['values'] = y1
            plotgui.xdata[ind1]['lowerror'] = x1 * 0.
            plotgui.xdata[ind1]['higherror'] = x1 * 0.
            plotgui.ydata[ind1]['lowerror'] = y1 * 0.
            plotgui.ydata[ind1]['higherror'] = y1 * 0.
        else:
            plotgui.add_set(x1, y1, current_plot=plotgui.current_plot)
        make_plot.make_plot(plotgui)
    except Exception:
        return

def make_data_set_delete_window(plotgui):
    """
    Create the window for deleting data sets.

    This routine makes a window within which a data set can be
    removed from the data list.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if plotgui.data_set_delete_window is not None:
        return
    plotgui.data_set_delete_window = Tk.Toplevel()
    plotgui.data_set_delete_window.title('Data Set Delete Window')
    holder = Tk.Frame(plotgui.data_set_delete_window)
    holder.pack(side=Tk.TOP)
    plotgui.set_list_area2 = tkinter.ttk.Combobox(holder, width=50, height=10)
    plotgui.set_list_area2.pack(side=Tk.TOP)
    if plotgui.nsets == 0:
        setlist = [' ']
    else:
        setlist = []
        for loop in range(plotgui.nsets):
            label = (
                'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                (loop+1, len(plotgui.xdata[loop]['values']),
                 plotgui.xdata[loop]['minimum'],
                 plotgui.xdata[loop]['maximum'])
                 + ' y range %13.6g to %13.6g' %
                 (plotgui.ydata[loop]['minimum'],
                  plotgui.ydata[loop]['maximum']))
            setlist.append(label)
    plotgui.set_list_area2['values'] = setlist
    plotgui.set_list_area2.current(0)
#    plotgui.set_list_area2.bind('<<ComboboxSelected>>', plotgui.set_fields1)
    label1 = Tk.Label(
        holder,
        text="Select a set above and remove with the 'delete' button.\n"
            + " Note that this action cannot be undone so be careful.")
    label1.pack(side=Tk.TOP)
    frame2 = Tk.Frame(holder)
    frame2.pack()
    apply_button = Tk.Button(frame2, text="Delete Set",
                             command=lambda: delete_set(plotgui))
    apply_button.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        frame2, text="Close",
        command=lambda: plotgui.close_data_window(
            plotgui.data_set_delete_window))
    close_button.pack(side=Tk.LEFT)

def delete_set(plotgui):
    """
    Carry out clearing a data set from the variables.

    This is a work routine that takes the input for deleting a set, checks
    with the user, and if requested removes a set from the list.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    ind1 = plotgui.set_list_area2.current()
    if (ind1 >= 0) and (ind1 < plotgui.nsets):
        response = tkinter.messagebox.askyesno(
            'Verify',
            'Do you really want to remove set %d?' % (ind1+1))
        if response:
            if plotgui.nsets > 0:
                plotgui.nsets = plotgui.nsets - 1
                for loop in range(ind1, plotgui.nsets):
                    plotgui.set_properties[loop] = deepcopy(
                        plotgui.set_properties[loop+1])
                    plotgui.xdata[loop] = deepcopy(plotgui.xdata[loop+1])
                    plotgui.ydata[loop] = deepcopy(plotgui.ydata[loop+1])
                plotgui.xdata[plotgui.nsets+1] = None
                plotgui.ydata[plotgui.nsets+1] = None
                plotgui.set_properties[plotgui.nsets+1]['symbol'] = None
                plotgui.set_properties[plotgui.nsets+1]['symbolsize'] = 4.0
                plotgui.set_properties[plotgui.nsets+1]['linestyle'] = 'None'
                plotgui.set_properties[plotgui.nsets+1]['linewidth'] = 1.0
                plotgui.set_properties[plotgui.nsets+1]['colour'] = 'black'
                plotgui.set_properties[plotgui.nsets+1]['label'] = ''
                plotgui.set_properties[plotgui.nsets+1]['xmin'] = 0.0
                plotgui.set_properties[plotgui.nsets+1]['xmax'] = 0.0
                plotgui.set_properties[plotgui.nsets+1]['ymin'] = 0.0
                plotgui.set_properties[plotgui.nsets+1]['ymax'] = 0.0
                plotgui.set_properties[plotgui.nsets+1]['display'] = True
                plotgui.set_properties[plotgui.nsets+1]['errors'] = False
                plotgui.set_properties[plotgui.nsets+1]['legend'] = True
                plotgui.set_properties[plotgui.nsets+1]['plot'] = 1
                make_plot.make_plot(plotgui)
                # revise the set list area
                setlist = []
                for loop in range(plotgui.nsets):
                    label = (
                        'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                        (loop+1, len(plotgui.xdata[loop]['values']),
                         plotgui.xdata[loop]['minimum'],
                         plotgui.xdata[loop]['maximum'])
                        + ' y range %13.6g to %13.6g' %
                        (plotgui.ydata[loop]['minimum'],
                         plotgui.ydata[loop]['maximum']))
                    setlist.append(label)
                plotgui.set_list_area2['values'] = setlist

def make_data_set_window(plotgui):
    """
    Create the data set properties window (symbol, colour, etc).

    The code makes a window within which one can set the display parameters
    for the data sets (symbol, colour, line properties, legend labels).  When
    one exits the window it is destroyed, so it gets regenerated each time
    the parent button is clicked.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    matplotlib_line_list = ['-', '--', '-.', ':', None]
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot', 'dotted',
                                 'None']
    matplotlib_symbol_list = [None, ".", ", ", "o", "v", "^", "<", ">", "1",
        "2", "3", "4", "8", "s", "p", "P", "*", "h", "H",
        "+", "x", "X", "D", "d", "|", "_", "histogram"]
    matplotlib_symbol_name_list = ['None', 'point', 'pixel', 'circle',
        'triangle down', 'triangle up', 'triangle left', 'triangle right',
        'tri_down', 'tri_up', 'tri_left', 'tri_right', 'octagon',
        'square', 'pentagon', 'plus (filled)', 'star', 'hexagon1',
        'hexagon2', 'plus', 'x', 'x (filled)', 'diamond', 'thin_diamond',
        'vline', 'hline', 'histogram']
    if plotgui.data_set_window is not None:
        return
    plotgui.data_set_window = Tk.Toplevel()
    plotgui.data_set_window.title('Data Sets Window')
    holder = Tk.Frame(plotgui.data_set_window)
    holder.pack(side=Tk.TOP)
    plotgui.set_list_area = tkinter.ttk.Combobox(holder, width=50, height=10)
    plotgui.set_list_area.pack(side=Tk.TOP)
    if plotgui.nsets == 0:
        setlist = [' ']
    else:
        setlist = []
        for loop in range(plotgui.nsets):
            label = (
                'Set %3d: %5d points, x range %13.6g to %13.6g, ' %
                (loop+1, len(plotgui.xdata[loop]['values']),
                 plotgui.xdata[loop]['minimum'], plotgui.xdata[loop]['maximum'])
                + ' y range %13.6g to %13.6g' %
                (plotgui.ydata[loop]['minimum'], plotgui.ydata[loop]['maximum'])
            )
            setlist.append(label)
    plotgui.set_list_area['values'] = setlist
    plotgui.set_list_area.current(0)
    plotgui.set_list_area.bind('<<ComboboxSelected>>', plotgui.set_fields)
    frame1 = Tk.Frame(holder)
    frame1.pack()
    plotgui.set_window_label = Tk.Label(frame1, text='Set Properties:')
    plotgui.set_window_label.grid(row=0, column=0, columnspan=2)
    labels = ['Symbol', 'Symbol Size', 'Line', 'Line Width', 'Colour',
              'Show', 'Error Bars', 'Label', 'Legend', 'Plot']
    for loop in range(len(labels)):
        label = Tk.Label(frame1, text=labels[loop])
        label.grid(column=0, row=loop+1, sticky=Tk.W)
    plotgui.set_entry_fields = []
    # item 0: menu for the symbol
    plotgui.set_entry_fields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.set_entry_fields[-1].grid(column=1, row=1, sticky=Tk.W)
    plotgui.set_entry_fields[-1]['values'] = matplotlib_symbol_name_list
    plotgui.set_entry_fields[-1].current(0)
    # item 1: entry field for the symbol size value
    plotgui.set_entry_fields.append(Tk.Entry(frame1, width=15))
    plotgui.set_entry_fields[-1].grid(column=1, row=2, sticky=Tk.W)
    plotgui.set_entry_fields[-1].insert(0, '1.0')
    # item 2: menu for the line style
    plotgui.set_entry_fields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.set_entry_fields[-1].grid(column=1, row=3, sticky=Tk.W)
    plotgui.set_entry_fields[-1]['values'] = matplotlib_line_name_list
    plotgui.set_entry_fields[-1].current(1)
    # item 3: entry field for the line width value
    plotgui.set_entry_fields.append(Tk.Entry(frame1, width=15))
    plotgui.set_entry_fields[-1].grid(column=1, row=4, sticky=Tk.W)
    # item 4: menu for the colour value
    plotgui.set_entry_fields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.set_entry_fields[-1].grid(column=1, row=5, sticky=Tk.W)
    plotgui.set_entry_fields[-1]['values'] = plotgui.colourset
    plotgui.set_entry_fields[-1].current(0)
    # item 5: variable for the show/hide set option
    plotgui.set_entry_fields.append(Tk.IntVar())
    b1 = Tk.Frame(frame1)
    general_utilities.put_yes_no(b1, plotgui.set_entry_fields[-1], ['Yes', 'No'], True)
    b1.grid(column=1, row=6, sticky=Tk.W)
    # item 6: variable for the error bar display option
    plotgui.set_entry_fields.append(Tk.IntVar())
    b1 = Tk.Frame(frame1)
    general_utilities.put_yes_no(b1, plotgui.set_entry_fields[-1], ['Yes', 'No'], False)
    b1.grid(column=1, row=7, sticky=Tk.W)
    # item 7: entry field for the set label (for the legend, if any)
    plotgui.set_entry_fields.append(Tk.Entry(frame1, width=35))
    plotgui.set_entry_fields[-1].grid(column=1, row=8, sticky=Tk.W)
    # item 8: flag for whether to include the set in the legend
    plotgui.set_entry_fields.append(Tk.IntVar())
    b1 = Tk.Frame(frame1)
    general_utilities.put_yes_no(b1, plotgui.set_entry_fields[-1], ['Yes', 'No'], True)
    b1.grid(column=1, row=9, sticky=Tk.W)
    # item 9: Entry field for which plot the set belongs to
    plotgui.set_entry_fields.append(Tk.Entry(frame1, width=5))
    plotgui.set_entry_fields[-1].grid(column=1, row=10, sticky=Tk.W)
    plotgui.set_entry_fields[-1].insert(0, str(plotgui.current_plot))
    frame2 = Tk.Frame(holder)
    frame2.pack()
    apply_button = Tk.Button(
        frame2, text="Apply",
        command=lambda: apply_data_set_fields(plotgui))
    apply_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(frame2, text="    ")
    label1.pack(side=Tk.LEFT)
    apply_all_button = Tk.Button(
        frame2, text="Apply to All",
        command=lambda: apply_data_set_fields_all(plotgui))
    apply_all_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(frame2, text="    ")
    label1.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        frame2, text="Close",
        command=lambda: plotgui.close_data_window(plotgui.data_set_window))
    close_button.pack(side=Tk.LEFT)
    if plotgui.nsets > 0:
        plotgui.set_property_fields(0)

def apply_data_set_fields_all(plotgui):
    """
    Apply the values in the data set window to all sets.

    This routine reads the values in the plotgui.data_set_window and applies
    these to all data sets on the list.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

    """
    index = plotgui.set_list_area.current()
    for loop in range(plotgui.nsets):
        plotgui.set_list_area.current(newindex=loop)
        if loop == index:
            apply_data_set_fields(plotgui, use_label=True)
        else:
            apply_data_set_fields(plotgui, use_label=False)
    plotgui.set_list_area.current(newindex=index)
    make_plot.make_plot(plotgui)

def apply_data_set_fields(plotgui, use_label=True):
    """
    Apply the values in the data set window.

    This routine reads the values in the plotgui.data_set_window and applies
    these to the selected data set on the list.  The only possible 
    exception is whether the label is changed.  

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

        use_label:  An optional boolean value; if True the label is 
                    applied, if False it is not; used when all sets 
                    are assigned the same parameters, as one would 
                    normally expect the labels to be distinct

    Returns
    -------

        None

    """
    matplotlib_symbol_list = [None, ".", ", ", "o", "v", "^", "<", ">", "1",
        "2", "3", "4", "8", "s", "p", "P", "*", "h", "H",
        "+", "x", "X", "D", "d", "|", "_", "histogram"]
    matplotlib_symbol_name_list = ['None', 'point', 'pixel', 'circle',
        'triangle down', 'triangle up', 'triangle left', 'triangle right',
        'tri_down', 'tri_up', 'tri_left', 'tri_right', 'octagon',
        'square', 'pentagon', 'plus (filled)', 'star', 'hexagon1',
        'hexagon2', 'plus', 'x', 'x (filled)', 'diamond', 'thin_diamond',
        'vline', 'hline', 'histogram']
    matplotlib_line_list = ['-', '--', '-.', ':', None]
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot', 'dotted',
                                 'None']
    try:
        index = plotgui.set_list_area.current()
        symbol_index = plotgui.set_entry_fields[0].current()
        line_index = plotgui.set_entry_fields[2].current()
        colour_index = plotgui.set_entry_fields[4].current()
        label = plotgui.set_entry_fields[7].get()
        symbol_size = float(plotgui.set_entry_fields[1].get())
        line_width = float(plotgui.set_entry_fields[3].get())
        show_flag = plotgui.set_entry_fields[5].get()
        error_bar_flag = plotgui.set_entry_fields[6].get()
        legend_flag = plotgui.set_entry_fields[8].get()
        target_plot = int(plotgui.set_entry_fields[9].get())
        if (target_plot < 1) or (target_plot > plotgui.number_of_plots):
            target_plot = plotgui.set_properties[index]['plot']
        plotgui.set_properties[index]['symbol'] = \
            matplotlib_symbol_list[symbol_index]
        plotgui.set_properties[index]['linestyle'] = \
            matplotlib_line_list[line_index]
        if plotgui.colourset[colour_index] == 'select':
            values = askcolor()
            plotgui.set_properties[index]['colour'] = values[1]
        else:
            plotgui.set_properties[index]['colour'] = \
                plotgui.colourset[colour_index]
        if use_label:
            plotgui.set_properties[index]['label'] = label
        plotgui.set_properties[index]['linewidth'] = line_width
        plotgui.set_properties[index]['symbolsize'] = symbol_size
        if show_flag == 1:
            plotgui.set_properties[index]['display'] = True
        else:
            plotgui.set_properties[index]['display'] = False
        if error_bar_flag == 1:
            plotgui.set_properties[index]['errors'] = True
        else:
            plotgui.set_properties[index]['errors'] = False
        if legend_flag == 1:
            plotgui.set_properties[index]['legend'] = True
        else:
            plotgui.set_properties[index]['legend'] = False
        plotgui.set_properties[index]['plot'] = target_plot
        make_plot.make_plot(plotgui)
    except Exception:
        tkinter.messagebox.showinfo(
            'Error',
            'There was some error applying the values.'
            + '  You may need to retry the settings.')

