import math
import numpy
import tkinter as Tk
from tkinter.scrolledtext import ScrolledText
import make_plot
import mpfitexpr

def activate_parameter(tkvars):
    """
    Process changes of state in the parameter selection

    Parameters
    ----------

    tkvars:  a list of tkinter variables for Entry and control of parameters
             in the non-linear fitting

    Returns
    -------

    None
    """
    flag = tkvars[0].get()
    if flag:
        tkvars[2].config(state='normal')
    else:
        tkvars[2].config(state='disabled')

def process_state(tkvars):
    """
    Process changes of state in the selection Checkbuttons

    Parameters
    ----------

    tkvars:  a list of tkinter variables for Entry and control of parameters
             in the non-linear fitting

    Returns
    -------

    None
    """
    flag = tkvars[3].get()
    if flag:
        tkvars[5].config(state='normal')
    else:
        tkvars[5].config(state='disabled')
    flag = tkvars[6].get()
    if flag:
        tkvars[8].config(state='normal')
    else:
        tkvars[8].config(state='disabled')

def fitting_parameter_group(outerframe, n, active):
    """
    Set a group of tkinter items for one non-linear fitting parameter.

    Parameters
    ----------

    outframe:  A tkinter frame variable that holds the set of widget items

    n:    An integer value for the parameter number (used in the labels)

    active:  A boolean value for whether the parameter is active or not

    Returns
    -------

    tkvars:  A list of tkinter variables

    Values in tkvars --

        0    parameter active IntVar
        1    parameter active Checkbutton
        2    parameter Entry field
        3    minimum range active IntVar
        4    minimum range active Checkbutton
        5    range minimum Entry field
        6    maximum range active IntVar
        7    maximum range active Checkbutton
        8    range maximum Entry field
        9    fixed IntVar
       10    fixed Checkbutton
    """
    tkvars=[]
    holder = Tk.Frame(outerframe)
    holder.pack(side=Tk.TOP)
    var1 = Tk.IntVar()
    tkvars.append(var1)
    control1 = Tk.Checkbutton(holder, variable=var1,
                              command=lambda: activate_parameter(tkvars))
    control1.pack(side=Tk.LEFT)
    tkvars.append(control1)
    str1 = 'Parameter %2d:' % n
    label = Tk.Label(holder, text=str1)
    label.pack(side=Tk.LEFT)
    entry1 = Tk.Entry(holder, width=10)
    entry1.pack(side=Tk.LEFT)
    entry1.insert(0, '1.0')
    tkvars.append(entry1)
    if active:
        var1.set(1)
    else:
        var1.set(0)
        entry1.config(state='disabled')
    label = Tk.Label(holder, text='Bounds:')
    label.pack(side=Tk.LEFT)
    var2 = Tk.IntVar()
    tkvars.append(var2)
    control2 = Tk.Checkbutton(holder, variable=var2,
                              command=lambda: process_state(tkvars))
    control2.pack(side=Tk.LEFT)
    tkvars.append(control2)
    var2.set(0)
    label = Tk.Label(holder, text='Minimum')
    label.pack(side=Tk.LEFT)
    entry2 = Tk.Entry(holder, width=10)
    entry2.pack(side=Tk.LEFT)
    entry2.insert(0, '1.0')
    entry2.config(state='disabled')
    tkvars.append(entry2)
    var3 = Tk.IntVar()
    tkvars.append(var3)
    control3 = Tk.Checkbutton(holder, variable=var3,
                              command=lambda: process_state(tkvars))
    control3.pack(side=Tk.LEFT)
    tkvars.append(control3)
    var3.set(0)
    label = Tk.Label(holder, text='Maximum')
    label.pack(side=Tk.LEFT)
    entry3 = Tk.Entry(holder, width=10)
    entry3.pack(side=Tk.LEFT)
    entry3.insert(0, '1.0')
    entry3.config(state='disabled')
    tkvars.append(entry3)
    label = Tk.Label(holder, text='Fixed')
    label.pack(side=Tk.LEFT)
    var4 = Tk.IntVar()
    tkvars.append(var4)
    control4 = Tk.Checkbutton(holder, variable=var4)
    control4.pack(side=Tk.LEFT)
    tkvars.append(control4)
    return tkvars

def make_fitting_window(plotgui):
    """
    Create the window for non-linear data set fitting.

    This routine creates a window for the data set fitting functions.
    One enters parameters and a function, and the code parses this into
    something that mpfit.py can work with for the fitting, similar to
    what is in xmgrace.

    If there are no sets to fit, the routine returns without doing
    anything.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if plotgui.non_linear_set_fitting_window is not None:
        plotgui.non_linear_set_fitting_window.deiconify()
        return
    if plotgui.nsets == 0:
        return
    plotgui.non_linear_set_fitting_window = Tk.Toplevel()
    plotgui.non_linear_set_fitting_window.title(
        'Non-linear Least Squares Fitting Window')
    holder = Tk.Frame(plotgui.non_linear_set_fitting_window)
    holder.pack(side=Tk.TOP)
    label1 = Tk.Label(holder, text='Levenberg-Marquardt Least Squares Fitting')
    label1.pack(side=Tk.TOP)
    plotgui.non_linear_set_fitting_text_field = ScrolledText(
        holder,height=20, width=80, wrap=Tk.NONE, relief='solid')
    plotgui.non_linear_set_fitting_text_field.config(font=('courier', 16))
    plotgui.non_linear_set_fitting_text_field.pack(side=Tk.TOP, padx=10,
                                                   pady=10)
    plotgui.tkcontrol = []
    for loop in range(10):
        active = bool(loop < 2)
        tkvals = fitting_parameter_group(holder, loop+1, active)
        plotgui.tkcontrol.append(tkvals)
    h1 = Tk.Frame(holder)
    h1.pack(side=Tk.TOP)
    label2 = Tk.Label(h1, text='Tolerance: ')
    label2.pack(side=Tk.LEFT)
    tolerance_entry = Tk.Entry(h1, width=15)
    tolerance_entry.pack(side=Tk.LEFT)
    tolerance_entry.insert(0, '1.0e-10')
    plotgui.tkcontrol.append(tolerance_entry)
    label3 = Tk.Label(h1, text='Set to Fit:')
    label3.pack(side=Tk.LEFT)
    set_entry = Tk.Entry(h1, width=5)
    set_entry.pack(side=Tk.LEFT)
    set_entry.insert(0, '1')
    plotgui.tkcontrol.append(set_entry)
    str1 = 'Function (use p[0]...p[9] for the parameters, and x for the input values)\n' + \
           'Use full names, e.g. numpy.sin(p[0]*x), for numpy functions.'
    label4 = Tk.Label(holder, text=str1, anchor='e')
    label4.pack(side=Tk.TOP)
    function_entry = Tk.Entry(holder, width=80)
    function_entry.pack(side=Tk.TOP)
    plotgui.tkcontrol.append(function_entry)
    bframe = Tk.Frame(holder)
    bframe.pack(side=Tk.TOP)
    select_button = Tk.Button(bframe, text="Do Fitting",
                              command=lambda: run_fitting(plotgui))
    select_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(bframe, text="    ")
    label1.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        bframe, text="Cancel/Close",
        command=lambda: plotgui.close_window(
            plotgui.non_linear_set_fitting_window,
            'non_linear_set_fitting_window'))
    close_button.pack(side=Tk.LEFT)

def run_fitting(plotgui):
    function_string = plotgui.tkcontrol[-1].get()
    try:
        tolerance = float(plotgui.tkcontrol[-3].get())
        if (tolerance <= 0.) or (tolerance > 0.0001):
            tolerance = 1.e-10
    except:
        tolerance = 1.e-10
    nset = int(plotgui.tkcontrol[-2].get())
    if (nset < 0) or (nset > plotgui.nsets):
        nset = 1
    params = []
    start = []
    lowbound = []
    highbound = []
    fixed = []
    lowflag = []
    highflag = []
    for loop in range(10):
        if plotgui.tkcontrol[loop][0].get():
            params.append(True)
            start.append(float(plotgui.tkcontrol[loop][2].get()))
            if plotgui.tkcontrol[loop][3].get():
                lowbound.append(float(plotgui.tkcontrol[loop][5].get()))
                lowflag.append(True)
            else:
                lowbound.append(0.)
                lowflag.append(False)
            if plotgui.tkcontrol[loop][6].get():
                highbound.append(float(plotgui.tkcontrol[loop][8].get()))
                highflag.append(True)
            else:
                highbound.append(0.)
                highflag.append(False)
            if plotgui.tkcontrol[loop][9].get():
                fixed.append(True)
            else:
                fixed.append(False)
        else:
            params.append(False)
            start.append(0.)
            lowflag.append(False)
            lowbound.append(0.)
            highflag.append(False)
            highbound.append(0.)
            fixed.append(True)
    params = numpy.asarray(params, dtype=bool)
    inds = numpy.where(params)
    start = numpy.asarray(start)
    lowbound = numpy.asarray(lowbound)
    highbound = numpy.asarray(highbound)
    fixed = numpy.asarray(fixed, dtype=bool)
    lowflag = numpy.asarray(lowflag, dtype=bool)
    highflag = numpy.asarray(highflag, dtype=bool)
    parinfo = []
    for loop in range(10):
        if lowflag[loop] and highflag[loop]:
            if (lowbound[loop] == highbound[loop]):
                str1 = '\nBad limits on parameter %d; check your inputs\n' % (
                    loop)
                plotgui.non_linear_set_fitting_text_field.insert(Tk.END, str1)
                return
        if params[loop]:
            parinfo1 = {'value': start[loop], 'fixed': fixed[loop],
                        'limited': [lowflag[loop], highflag[loop]],
                        'limits': [lowbound[loop], highbound[loop]]}
            parinfo.append(parinfo1)
    plotgui.non_linear_set_fitting_text_field.insert(
        Tk.END,'Function to fit:\n')
    plotgui.non_linear_set_fitting_text_field.insert(
        Tk.END, function_string+'\n')
    try:
        yerrors = (numpy.copy(plotgui.ydata[nset-1]['lowerror']) +
                   numpy.copy(plotgui.ydata[nset-1]['higherror']))/2.
    except:
        plotgui.non_linear_set_fitting_text_field.insert(
            'Constant uncertainties used.\n')
        yerrors = plotgui.xdata[nset-1]['values']*0.+0.01
    emin = numpy.min(yerrors)
    emax = numpy.max(yerrors)
    if (emin == 0.) and (emax == 0.):
        yerrors = yerrors+0.01
    else:
        inds = numpy.where(yerrors <= 0.)
        altinds = numpy.where(yerrors > 0.)
        if len(altinds[0]) == 0:
            yerrors = yerrors*0.+0.01
        else:
            meanerror = numpy.mean(yerrors[altinds])
            yerrors[inds] = meanerror
    start_values = start[params]
    fitparams, yfit = mpfitexpr.mpfitexpr(
        function_string, plotgui.xdata[nset-1]['values'],
        plotgui.ydata[nset-1]['values'], yerrors, start_values, check=True,
        full_output=True, parinfo=parinfo, ftol=tolerance)
    str1 = '\nMPFIT status: %d\n' % (fitparams.status)
    if fitparams.status < 0:
        str1 = 'Error in the fitting: '+str1
    else:
        cov = fitparams.covar
        pcov = cov*0.
        sh1 = cov.shape
        for i in range(sh1[0]):
            for j in range(sh1[1]):
                pcov[i, j] = cov[i,j]/math.sqrt(cov[i, i]*cov[j, j])
        for loop in range(len(fitparams.params)):
            str1 = str1 + 'parameter %2d : value %15.8g +/- %15.8g\n' % (
                loop+1, fitparams.params[loop], fitparams.perror[loop])
    plotgui.non_linear_set_fitting_text_field.insert(Tk.END, str1)
    xvalues = numpy.copy(plotgui.xdata[nset-1]['values'])
    yvalues = numpy.copy(plotgui.ydata[nset-1]['values'])
    plotgui.add_set(xvalues, yfit)
    oset = 1*plotgui.nsets
    mpfit_values = {
        'fit_parameters': fitparams, 'set_fit': nset, 'set_out': oset,
        'xvalues': xvalues, 'yvalues': yvalues, 'yerrors': yerrors
    }
    plotgui.mpfit_values.append(mpfit_values)
    make_plot.make_plot(plotgui)
    for loop in range(10):
        if params[loop]:
            plotgui.tkcontrol[loop][2].delete(0, Tk.END)
            plotgui.tkcontrol[loop][2].insert(0, str(fitparams.params[loop]))
