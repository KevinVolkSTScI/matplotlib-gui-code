import math
import numpy
import tkinter as Tk
import tkinter.ttk
import tkinter.filedialog
import tkinter.simpledialog
import tkinter.messagebox
from tkinter.colorchooser import askcolor
from tkinter.scrolledtext import ScrolledText
import general_utilities
import make_plot
import interpolation_utilities

def create_data_set(plotgui):
    """
    Open a window to define a function for making a data set.

    This routine makes the window that allows a limited capability to
    create data sets via a defined function.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

       None

    """
    function_window = Tk.Toplevel()
    function_window.title('Create Data Set')
    holder = Tk.Frame(function_window)
    holder.pack(side=Tk.TOP)
    h1 = Tk.Frame(holder)
    h1.pack(side=Tk.TOP)
    label1 = Tk.Label(h1, text='    Start values at ')
    label1.pack(side=Tk.LEFT)
    plotgui.start_value_field = Tk.Entry(h1, width=12)
    plotgui.start_value_field.pack(side=Tk.LEFT)
    label2 = Tk.Label(h1, text='    Stop values at ')
    label2.pack(side=Tk.LEFT)
    plotgui.stop_value_field = Tk.Entry(h1, width=12)
    plotgui.stop_value_field.pack(side=Tk.LEFT)
    label3 = Tk.Label(h1, text=' Number of values/step ')
    label3.pack(side=Tk.LEFT)
    plotgui.number_of_values_field = Tk.Entry(h1, width=6)
    plotgui.number_of_values_field.pack(side=Tk.LEFT)
    h2 = Tk.Frame(holder)
    h2.pack()
    plotgui.sequence_option = Tk.IntVar()
    lab1 = Tk.Label(
        h2,
        text='Spacing (logarithmic only works for postive range values): ')
    lab1.pack(side=Tk.LEFT)
    b1 = Tk.Radiobutton(
        h2, text='linear', variable=plotgui.sequence_option,
        value=0)
    b1.pack(side=Tk.LEFT)
    b2 = Tk.Radiobutton(
        h2, text='logarithmic', variable=plotgui.sequence_option,
        value=1)
    plotgui.sequence_option.set(0)
    b2.pack(side=Tk.LEFT)
    h3 = Tk.Frame(holder)
    h3.pack(side=Tk.TOP)
    label1 = Tk.Label(h3, text=' x function: ')
    label1.grid(column=0, row=0)
    plotgui.xfunction = Tk.Entry(h3, width=50)
    plotgui.xfunction.grid(column=1, row=0)
    label2 = Tk.Label(h3, text=' y function: ')
    label2.grid(column=0, row=1)
    plotgui.yfunction = Tk.Entry(h3, width=50)
    plotgui.yfunction.grid(column=1, row=1)
    plotgui.xfunction.insert(0, '$t')
    plotgui.yfunction.insert(0, '$x*$x')
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
    button1 = Tk.Button(h4, text='Apply',
                        command=lambda: parse_function(plotgui))
    button1.pack(side=Tk.LEFT)
    label1 = Tk.Label(h4, text='   ')
    label1.pack(side=Tk.LEFT)
    button2 = Tk.Button(h4, text='Close', command=function_window.destroy)
    button2.pack(side=Tk.LEFT)

def parse_function(plotgui):
    """
    Read parameters and try to parse to make a data set.

    This routine reads the parameters to make a new set and attempts to
    evaluate them using a parser function.  Whilst one could use the
    eval() function this is considered as too dangerous.

    """
    try:
        v1 = float(plotgui.start_value_field.get())
        v2 = float(plotgui.stop_value_field.get())
        if v1 > v2:
            temp = v1
            v1 = v2
            v2 = temp
        # check for a number of values or a step value
        try:
            n1 = int(plotgui.number_of_values_field.get())
        except ValueError:
            n1 = 0
            step = float(plotgui.number_of_values_field.get())
        # if the value n1 is not an integer, assume it is a step size
        option = plotgui.sequence_option.get()
        if (option == 0) or ((v1 <= 0.) or (v2 <= 0)):
            if n1 > 1:
                step = (v2-v1)/(n1-1)
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
        xstring = plotgui.xfunction.get()
        ystring = plotgui.yfunction.get()
        xstring = xstring.replace('$t', 'seq')
        ystring = ystring.replace('$t', 'seq')
        xstring = xstring.replace('$x', 'x')
        ystring = ystring.replace('$x', 'x')
        xstring = xstring.replace('$y', 'y')
        ystring = ystring.replace('$y', 'y')
        try:
            xvalues = my_eval(xstring, seq=seq, xvalues=None, yvalues=None)
            yvalues = my_eval(ystring, seq=seq, xvalues=xvalues, yvalues=None)
        except ValueError:
            yvalues = my_eval(ystring, seq=seq, xvalues=None, yvalues=None)
            xvalues = my_eval(xstring, seq=seq, xvalues=None, yvalues=yvalues)
        # deal with the case where one of x or y is entered as a constant
        try:
            n = len(xvalues)
        except ValueError:
            xvalues = numpy.asarray([xvalues, ])
        try:
            n = len(yvalues)
        except ValueError:
            yvalues = numpy.asarray([yvalues, ])
        if (len(xvalues) == 1) and (len(yvalues) > 1):
            xvalues = yvalues*0.+xvalues
        if (len(yvalues) == 1) and (len(xvalues) > 1):
            yvalues = xvalues*0.+yvalues
        if len(xvalues) != len(yvalues):
            tkinter.messagebox.showinfo(
                'Error',
                'There was some error trying to generate the sets. (1)')
            return
        if (xvalues is not None) and (yvalues is not None):
            plotgui.add_set(xvalues, yvalues, current_plot=plotgui.current_plot)
            make_plot.make_plot(plotgui)
        else:
            tkinter.messagebox.showinfo(
                'Error',
                'There was some error trying to generate the sets. (2)')
    except Exception:
        tkinter.messagebox.showinfo(
            'Error',
            'There was some error trying to generate the sets. (3)')

def my_eval(inputstring, seq, xvalues=None, yvalues=None):
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
    except AttributeError:
        sh2 = seq.shape
    try:
        sh3 = yvalues.shape
    except AttributeError:
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
    except Exception:
           return None

def create_data_set_by_editor(plotgui):
    """
    Allow the user to make a set by entering numbers in a window.


    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.data_entry_window is not None:
        return
    plotgui.data_entry_window = Tk.Toplevel()
    plotgui.data_entry_window.title('Enter Data Set Values')
    holder = Tk.Frame(plotgui.data_entry_window)
    holder.pack(side=Tk.TOP)
    holder.config(bg='black')
    plotgui.data_text = ScrolledText(holder, height=40, width=80,
                                     wrap=Tk.NONE, relief="solid")
    plotgui.data_text.config(font=('courier', 16))
    plotgui.data_text.pack(side=Tk.TOP, padx=10, pady=10)
    bframe = Tk.Frame(plotgui.data_entry_window)
    bframe.pack()
    set_button = Tk.Button(
        bframe, text="Apply",
        command=lambda: apply_data_input(plotgui))
    set_button.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        bframe, text="Close",
        command=lambda: plotgui.close_data_window(
            plotgui.data_entry_window))
    close_button.pack(side=Tk.LEFT)

def apply_data_input(plotgui):
    """
    Read the values from a data input text window and parse these to
    a data set for the plots.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    text = plotgui.data_text.get("1.0", Tk.END)
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
                'Unable to parse text from entry widget (2)')
            return
        plotgui.add_set(xvalues, yvalues, xlowerror=dxvalues1,
                        xhigherror=dxvalues2, ylowerror=dyvalues1,
                        yhigherror=dyvalues2, xerrorflag=errorflag,
                        yerrorflag=errorflag,
                        labelstring='Set from editor',
                        current_plot=plotgui.current_plot)
        make_plot.make_plot(plotgui)
    except Exception:
        tkinter.messagebox.showinfo(
            'error',
            'Unable to parse text from entry widget')
        return

def write_data_sets(plotgui):
    """
    Write the data values to an ascii output file.

    This routine writes the current set values (x, y) out to an ascii
    output file.  If no sets are defined, the routine simply returns.

    Parameters
    ----------
            
        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    The output is a simple (x, y) ascii dump with a blank line between
    data sets.
    """
    if plotgui.nsets == 0:
        return
    outfilename = tkinter.filedialog.asksaveasfilename()
    outfile = open(outfilename, "w")
    for loop in range(plotgui.nsets):
        xboth = False
        yboth = False
        for n in range(len(plotgui.xdata[loop]['values'])):
            if n == 0:
                if plotgui.xdata[loop]['errors']:
                    for n in range(len(plotgui.xdata[loop]['lowerror'])):
                        if plotgui.xdata[loop]['lowerror'][n] != \
                           plotgui.xdata[loop]['higherror'][n]:
                            xboth = True
                if plotgui.ydata[loop]['errors']:
                    for n in range(len(plotgui.ydata[loop]['lowerror'])):
                        if plotgui.ydata[loop]['lowerror'][n] != \
                           plotgui.ydata[loop]['higherror'][n]:
                            yboth = True
                headerstr = '# Set %d: %s' % (
                    loop+1, plotgui.set_properties[loop]['label'])
                print(headerstr, file=outfile)
                headerstr = '# X Value  |'
                if plotgui.xdata[loop]['errors']:
                    if xboth:
                        headerstr = headerstr \
                            + 'X Error minus   | X Error plus   |'
                    else:
                        headerstr = headerstr + 'X Error   |}'
                headerstr = headerstr + ' Y Value  |'
                if plotgui.ydata[loop]['errors']:
                    if yboth:
                        headerstr = headerstr \
                            + 'y Error minus   | X Error plus   |'
                    else:
                        headerstr = headerstr + 'Y Error   |'
                print(headerstr, file=outfile)
            str1 = ''
            str1 = str1 + my_format(plotgui.xdata[loop]['values'][n])
            if plotgui.xdata[loop]['errors']:
                str1 = str1 + my_format(plotgui.xdata[loop]['lowerror'][n])
                if xboth:
                    str1 = str1 + my_format(plotgui.xdata[loop]['higherror'][n])
            str1 = str1 + my_format(plotgui.ydata[loop]['values'][n])
            if plotgui.ydata[loop]['errors']:
                str1 = str1 + my_format(plotgui.ydata[loop]['lowerror'][n])
                if yboth:
                    str1 = str1 + my_format(plotgui.ydata[loop]['higherror'][n])
            print(str1, file=outfile)
        print(' ', file=outfile)
    outfile.close()

def my_format(value):
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

def set_statistics(plotgui):
    """
    Calculate set statistics.

    This routine prints some statistics about the different sets to
    a pop-up window.  If no sets are defined, then the routine just
    returns.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

    """
    if plotgui.nsets == 0:
        return
    outstr = 'Set  Number of points    range (minimum, maximum)' \
             + '       mean    standard deviation\n'
    for loop in range(plotgui.nsets):
        npoints = len(plotgui.xdata[loop]['values'])
        xdmin = numpy.min(plotgui.xdata[loop]['values'])
        xdmax = numpy.max(plotgui.xdata[loop]['values'])
        xdmean = numpy.mean(plotgui.xdata[loop]['values'])
        xdsigma = numpy.std(plotgui.xdata[loop]['values'])
        ydmin = numpy.min(plotgui.ydata[loop]['values'])
        ydmax = numpy.max(plotgui.ydata[loop]['values'])
        ydmean = numpy.mean(plotgui.ydata[loop]['values'])
        ydsigma = numpy.std(plotgui.ydata[loop]['values'])
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

def block_average(plotgui):
    """
    Make a new set by block averaging.

    A window is made in which to determine the parameters of the averaging, 
    whereupon a new set is made with the block averaged values.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

    """
    if plotgui.nsets == 0:
        return
    block_window = Tk.Toplevel()
    block_window.title('Block Average Data Set')
    holder = Tk.Frame(block_window)
    holder.pack(side=Tk.TOP)
    h1 = Tk.Frame(holder)
    h1.pack(side=Tk.TOP)
    label1 = Tk.Label(h1, text='Number of samples to average:')
    label1.pack(side=Tk.LEFT)
    plotgui.block_average_field = Tk.Entry(h1, width=12)
    plotgui.block_average_field.pack(side=Tk.LEFT)
    plotgui.block_average_field.insert(0, '10')
    h1 = Tk.Frame(holder)
    h1.pack(side=Tk.TOP)
    label2 = Tk.Label(h1, text='                   Clip set:')
    label2.pack(side=Tk.LEFT)
    plotgui.block_average_clip_variable = Tk.IntVar()
    bframe = Tk.Frame(h1)
    bframe.pack(side=Tk.LEFT)
    general_utilities.put_yes_no(bframe, plotgui.block_average_clip_variable,
                                 ['yes', 'no'], False)
    bholder = Tk.Frame(holder)
    bholder.pack(side=Tk.TOP)
    button1 = Tk.Button(bholder, text='Apply',
                        command=lambda: do_block_average(plotgui))
    button1.pack(side=Tk.LEFT)
    label1 = Tk.Label(bholder, text='   ')
    label1.pack(side=Tk.LEFT)
    button2 = Tk.Button(bholder, text='Close', command=block_window.destroy)
    button2.pack(side=Tk.LEFT)

def do_block_average(plotgui):
    """
    Take the block averaging parameters and make a new data set.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

    """
    if True:
#    try:
        nsample = int(plotgui.block_average_field.get())
        npoints = len(
            plotgui.xdata[plotgui.current_plot-1]['values'])
        xdata = numpy.copy(
            plotgui.xdata[plotgui.current_plot-1]['values'])
        ydata = numpy.copy(
            plotgui.ydata[plotgui.current_plot-1]['values'])
        if (nsample < 2) or (nsample > npoints//2):
            tkinter.messagebox.showinfo(
                'Error',
                'Bad Nsamples value (%d)' % (nsample))
            return
        flag = (plotgui.block_average_clip_variable.get() == 0)
        xnew, ynew = interpolation_utilities._smoother(
            xdata, ydata, nsample, clip=flag)
        str1 = '%d element block average of set %d' % (nsample,
                                                       plotgui.current_plot)
        plotgui.add_set(xnew, ynew, labelstring=str1,
                        current_plot=plotgui.current_plot)
        make_plot.make_plot(plotgui)
#    except:
#        tkinter.messagebox.showinfo(
#            'Error',
#            'Unable to do the block averaging.  Please check your inputs.')
