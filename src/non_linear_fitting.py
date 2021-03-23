import numpy
import tkinter as Tk
import make_plot
#import mpfit
#import mpfitexpr.py

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
        tkvars[6].config(state='normal')
    else:
        tkvars[5].config(state='disabled')
        tkvars[6].config(state='disabled')

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
        3    range active IntVar
        4    range active Checkbutton
        5    range minimum Entry field
        6    range maximum Entry field
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
    entry1 = Tk.Entry(holder, width=15)
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
    entry2 = Tk.Entry(holder, width=15, state='disabled')
    entry2.pack(side=Tk.LEFT)
    entry2.insert(0, '1.0')
    tkvars.append(entry2)
    label = Tk.Label(holder, text='Maximum')
    label.pack(side=Tk.LEFT)
    entry3 = Tk.Entry(holder, width=15, state='disabled')
    entry3.pack(side=Tk.LEFT)
    entry3.insert(0, '1.0')
    tkvars.append(entry3)
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
    plotgui.tkcontrol = []
    for loop in range(10):
        if loop < 2:
            active = True
        else:
            active = False
        tkvals = fitting_parameter_group(holder, loop+1, active)
        plotgui.tkcontrol.append(tkvals)
    h1 = Tk.Frame(holder)
    h1.pack(side=Tk.TOP)
    label2 = Tk.Label(h1, text='Tolerance: ')
    label2.pack(side=Tk.LEFT)
    tolerance_entry = Tk.Entry(h1, width=15)
    tolerance_entry.pack(side=Tk.LEFT)
    tolerance_entry.insert(0, '0.01')
    plotgui.tkcontrol.append(tolerance_entry)
    str1 = 'Function (use p[0]...p[9] for the parameters, and x)\n' + \
           'Use e.g. numpy.sin(p[0]*x) for numpy functions.'
    label3 = Tk.Label(holder, text=str1)
    label3.pack(side=Tk.TOP)
    function_entry = Tk.Entry(holder, width=60)
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
    tolerance = float(plotgui.tkcontrol[-2].get())
    if (tolerance > 0.1) or (tolerance <= 0.):
        tolerance = 0.01
    # It is not clear that the tolerance will be used....if not will take it
    # out later.
    params = []
    start = []
    lowbound = []
    highbound = []
    for loop in range(10):
        if plotgui.tkcontrol[loop][0]:
            params.append(True)
            start.append(float(plotgui.tkcontrol[loop][2].get()))
            if plotgui.tkcontrol[loop][3]:
                lowbound.append(float(plotgui.rkcontrol[loop][5].get()))
                highbound.append(float(plotgui.rkcontrol[loop][6].get()))
            else:
                lowbound.append(0.)
                highbound.append(0.)
        else:
            params.append(False)
            start.append(0.)
            lowbound.append(0.)
            highbound.append(0.)
    
