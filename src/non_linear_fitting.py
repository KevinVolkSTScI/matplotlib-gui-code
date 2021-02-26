import numpy
import tkinter as Tk
import make_plot
#import mpfit
#import mpfitexpr.py

def process_state(tkvars):
    flag = tkvars[1].get()
    if flag == 1:
        tkvars[3].config(state='normal')
        tkvars[4].config(state='normal')
    else:
        tkvars[3].config(state='disabled')
        tkvars[4].config(state='disabled')

def fitting_parameter_group(outerframe, n):
    tkvars=[]
    holder = Tk.Frame(outerframe)
    holder.pack(side=Tk.TOP)
    str1 = 'Parameter %2d:' % n
    label = Tk.Label(holder, text=str1)
    label.pack(side=Tk.LEFT)
    entry1 = Tk.Entry(holder, width=15)
    entry1.pack(side=Tk.LEFT)
    entry1.insert(0, '1.0')
    tkvars.append(entry1)
    label = Tk.Label(holder, text='Bounds:')
    label.pack(side=Tk.LEFT)
    var1 = Tk.IntVar()
    tkvars.append(var1)
    control1 = Tk.Checkbutton(holder, variable=var1,
                           command=lambda: process_state(tkvars))
    control1.pack(side=Tk.LEFT)
    tkvars.append(control1)
    var1.set(0)
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
        tkvals = fitting_parameter_group(holder, loop+1)
        plotgui.tkcontrol.append(tkvals)
    h1 = Tk.Frame(holder)
    h1.pack(side=Tk.TOP)
    label2 = Tk.Label(h1, text='Tolerance: ')
    label2.pack(side=Tk.LEFT)
    tolerance_entry = Tk.Entry(h1, width=15)
    tolerance_entry.pack(side=Tk.LEFT)
    tolerance_entry.insert(0, '0.01')
    plotgui.tkcontrol.append(tolerance_entry)
    label3 = Tk.Label(holder, text='Function (use A1...A10 for the parameters)')
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
    print(len(plotgui.tkcontrol))
