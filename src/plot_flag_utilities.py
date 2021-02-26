import numpy
import tkinter as Tk
import tkinter.messagebox
import make_plot

def toggle_equal_aspect(plotgui):
    """
    Toggle the equal aspect plot display flag.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    plotgui.equal_aspect[plotgui.current_plot-1] = not \
        plotgui.equal_aspect[plotgui.current_plot-1]
    make_plot.make_plot(plotgui)

def set_opposite_y(plotgui):
    """
    Arrange right side independent y axis on one of the plots.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    plot_number = plotgui.current_plot
    plotgui.subplot.append(plotgui.figure.add_subplot(
        plotgui.nxplots, plotgui.nyplots, plotgui.current_plot,
        sharex=plotgui.subplot[plotgui.current_plot-1], frameon=False))
#    plotgui.subplot.append(
#        plotgui.subplot[plotgui.current_plot-1].secondary_yaxis(
#        "right"))
    plotgui.share_axis.append(plot_number)
    plotgui.xparameters.append({
        'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
        'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
        'invert': 0, 'hide': 0, 'hideticks': 0,
        'hidelabels': 0, 'hybridlog': 0,
        'inverseticks': 0, 'ticklength': 6,
        'bothticks': 0, 'minorticks': 0, 'oppositeaxis': 0,
        'majorgridlines': 0, 'minorgridlines': 0})
    for key in plotgui.xparameters[-1].keys():
        plotgui.xparameters[-1][key] = plotgui.xparameters[
            plotgui.current_plot-1][key]
    plotgui.yparameters.append({
        'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
        'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
        'invert': 0, 'hide': 0, 'hideticks': 0,
        'hidelabels': 0, 'hybridlog': 0,
        'inverseticks': 0, 'ticklength': 6,
        'bothticks': 0, 'minorticks': 0,
        'oppositeaxis': 1-plotgui.yparameters[plotgui.current_plot-1][
            'oppositeaxis'],
        'majorgridlines': 0, 'minorgridlines': 0})
#    plotgui.yparameters[-1]['oppositeaxis'] = \
#        1-plotgui.yparameters[plotgui.current_plot-1]['oppositeaxis']
    plotgui.fontname.append(deepcopy(plotgui.fontname[plotgui.current_plot-1]))
    plotgui.fontsize.append(deepcopy(plotgui.fontsize[plotgui.current_plot-1]))
    plotgui.fontweight.append(deepcopy(plotgui.fontweight[plotgui.current_plot-1]))
    plotgui.legend_variable.append(None)
    plotgui.legend_frame.append(None)
    plotgui.legend_options.append(None)
    plotgui.legend_position.append(None)
    plotgui.legend_user_position.append(None)
    plotgui.plot_frame.append(deepcopy(plotgui.plot_frame[plotgui.current_plot-1]))
    plotgui.plot_range.append(deepcopy(plotgui.plot_range[plotgui.current_plot-1]))
    plotgui.bounding_box.append(deepcopy(plotgui.bounding_box[
        plotgui.current_plot-1]))
    plotgui.title.append('')
    plotgui.current_plot = len(plotgui.subplot)
    plotgui.number_of_plots = plotgui.number_of_plots+1
    make_plot.make_plot(plotgui)

def set_opposite_x(plotgui):
    """
    Arrange top independent x axis on one of the plots.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    plot_number = plotgui.current_plot
    plotgui.subplot.append(plotgui.figure.add_subplot(
        plotgui.nxplots, plotgui.nyplots, plotgui.current_plot,
        sharey=plotgui.subplot[plotgui.current_plot-1], frameon=False))
#    plotgui.subplot.append(
#        plotgui.subplot[plotgui.current_plot-1].secondary_xaxis(
#        "top"))
    plotgui.share_axis.append(-1*plot_number)
    plotgui.xparameters.append({
        'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
        'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
        'invert': 0, 'hide': 0, 'hideticks': 0,
        'hidelabels': 0, 'hybridlog': 0,
        'inverseticks': 0, 'ticklength': 6,
        'bothticks': 0, 'minorticks': 0,
        'oppositeaxis': 1-plotgui.xparameters[plotgui.current_plot-1][
            'oppositeaxis'],
        'majorgridlines': 0, 'minorgridlines': 0})
    plotgui.yparameters.append({
        'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
        'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
        'invert': 0, 'hide': 0, 'hideticks': 0,
        'hidelabels': 0, 'hybridlog': 0,
        'inverseticks': 0, 'ticklength': 6,
        'bothticks': 0, 'minorticks': 0,
        'oppositeaxis': 0,
        'majorgridlines': 0, 'minorgridlines': 0})
    for key in plotgui.yparameters[-1].keys():
        plotgui.yparameters[-1][key] = plotgui.yparameters[
            plotgui.current_plot-1][key]
    plotgui.xparameters[-1]['oppositeaxis'] = \
        1-plotgui.xparameters[plotgui.current_plot-1]['oppositeaxis']
    plotgui.fontname.append(deepcopy(plotgui.fontname[plotgui.current_plot-1]))
    plotgui.fontsize.append(deepcopy(plotgui.fontsize[plotgui.current_plot-1]))
    plotgui.fontweight.append(deepcopy(plotgui.fontweight[plotgui.current_plot-1]))
    plotgui.legend_variable.append(None)
    plotgui.legend_frame.append(None)
    plotgui.legend_options.append(None)
    plotgui.legend_position.append(None)
    plotgui.legend_user_position.append(None)
    plotgui.plot_frame.append(deepcopy(plotgui.plot_frame[plotgui.current_plot-1]))
    plotgui.plot_range.append(deepcopy(plotgui.plot_range[plotgui.current_plot-1]))
    plotgui.bounding_box.append(deepcopy(plotgui.bounding_box[
        plotgui.current_plot-1]))
    plotgui.title.append('')
    plotgui.grid_colour.append(None)
    plotgui.grid_linetype.append(None)
    plotgui.current_plot = len(plotgui.subplot)
    plotgui.number_of_plots = plotgui.number_of_plots+1
    make_plot.make_plot(plotgui)

def set_plot_number(plotgui):
    """
    Select the active plot, if there is more than one.

    This routine calls up a window so that the user can select the
    current active plot.  The routine is not useful unless more than
    one plot is defined.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if plotgui.setplot_window is not None:
        return
    plotgui.setplot_window = Tk.Toplevel()
    plotgui.setplot_window.title('Set Plot Number')
    holder = Tk.Frame(plotgui.setplot_window)
    holder.pack(side=Tk.TOP)
    label = Tk.Label(holder, text="Active plot:")
    label.grid(column=0, row=0)
    plotgui.currentplot_select = Tk.Entry(holder, width=10)
    plotgui.currentplot_select.grid(column=1, row=1)
    plotgui.currentplot_select.insert(0, str(plotgui.current_plot))
    buttonframe = Tk.Frame(plotgui.setplot_window)
    buttonframe.pack(side=Tk.TOP)
    apply_button = Tk.Button(buttonframe, text="Select",
                             command=lambda: set_current_plot(plotgui))
    apply_button.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        buttonframe, text="Close",
        command=lambda: plotgui.close_window(plotgui.setplot_window,
                                             'setplot_window'))
    close_button.pack(side=Tk.LEFT)

def set_plot_hide(plotgui):
    """
    Hide one of the plots, if there is more than one.

    This routine calls up a window so that the user can select the
    plot to hide.  The routine is not useful unless more than
    one plot is defined.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if plotgui.hideplot_window is not None:
        return
    plotgui.hideplot_window = Tk.Toplevel()
    plotgui.hideplot_window.title('Hide Plot Number')
    holder = Tk.Frame(plotgui.hideplot_window)
    holder.pack(side=Tk.TOP)
    label = Tk.Label(holder, text="Plot to Toggle Hide/Show:")
    label.grid(column=0, row=0)
    plotgui.hideplot_select = Tk.Entry(holder, width=10)
    plotgui.hideplot_select.grid(column=1, row=1)
    plotgui.hideplot_select.insert(0, str(plotgui.current_plot))
    buttonframe = Tk.Frame(plotgui.hideplot_window)
    buttonframe.pack(side=Tk.TOP)
    apply_button = Tk.Button(buttonframe, text="Select",
                             command=lambda: hide_selected_plot(plotgui))
    apply_button.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        buttonframe, text="Close",
        command=lambda: plotgui.close_window(plotgui.hideplot_window, 
                                             'hideplot_window'))
    close_button.pack(side=Tk.LEFT)

def tile_plots(plotgui):
    """
    Bring up a window to allow multiple plots in the main display.

    This subroutine creates a window for use in arranging multiple plots
    in a single display (i.e. canvas).

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if plotgui.tile_window is not None:
        return
    plotgui.tile_window = Tk.Toplevel()
    plotgui.tile_window.title('Set Number of Plots')
    holder = Tk.Frame(plotgui.tile_window)
    holder.pack(side=Tk.TOP)
    label = Tk.Label(holder, text="Number of plots in x direction:")
    label.grid(column=0, row=0)
    label = Tk.Label(holder, text="Number of plots in y direction:")
    label.grid(column=0, row=1)
    label = Tk.Label(holder, text="Active plot:")
    label.grid(column=0, row=2)
    plotgui.xplots_set = Tk.Entry(holder, width=10)
    plotgui.xplots_set.grid(column=1, row=0)
    plotgui.yplots_set = Tk.Entry(holder, width=10)
    plotgui.yplots_set.grid(column=1, row=1)
    plotgui.xplots_set.insert(0, '1')
    plotgui.yplots_set.insert(0, '1')
    plotgui.currentplot_set = Tk.Entry(holder, width=10)
    plotgui.currentplot_set.grid(column=1, row=2)
    plotgui.currentplot_set.insert(0, '1')
    buttonframe = Tk.Frame(plotgui.tile_window)
    buttonframe.pack(side=Tk.TOP)
    apply_button = Tk.Button(buttonframe, text="Apply Values",
                             command=plotgui.set_plot_layout)
    apply_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(buttonframe, text="     ")
    label1.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        buttonframe, text="Close",
        command=lambda: plotgui.close_window(plotgui.tile_window, 
                                             'tile_window'))
    close_button.pack(side=Tk.LEFT)

def set_current_plot(plotgui):
    """
    Set which plot is the current one, if there are multiple plots.

    The routine here brings up a window wherein the user can select
    the current plot number in the case where several plots are
    displayed.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    try:
        n1 = int(plotgui.currentplot_select.get())
        if (n1 < 1) | (n1 > plotgui.number_of_plots):
            raise ValueError
        if n1 != plotgui.current_plot:
            plotgui.current_plot = n1
    except ValueError:
        tkinter.messagebox.showinfo(
            "Error", "There was an error in the requested plot number.")

def hide_selected_plot(plotgui):
    """
    Set a plot to be hidden, if there are multiple plots.

    The routine here brings up a window wherein the user can select
    the current plot number in the case where several plots are
    displayed.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    try:
        n1 = int(plotgui.hideplot_select.get())
        if (n1 < 1) | (n1 > plotgui.number_of_plots):
            raise ValueError
        plotgui.hide_subplot[n1-1] = not plotgui.hide_subplot[n1-1]
        plotnumber = plotgui.current_plot
        plotgui.current_plot = n1
        make_plot.make_plot(plotgui)
        plotgui.current_plot = plotnumber
    except ValueError:
        tkinter.messagebox.showinfo(
            "Error", "There was an error in the requested plot number.")

def clear_plot(plotgui, query=True):
    """
    Clear the plot area.

    This routine clears the data sets and resets the various parameters to
    the initial values.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

        query :  An optional Boolean value for whether the user is queried
                 before the plot is cleared, which defaults to True

    Returns
    -------

        None

    """
    if query:
        response = tkinter.messagebox.askyesno(
            "Verify", "Do you want to abandon the plot?")
    else:
        response = True
    if not response:
        return
    plotgui.xdata = []
    plotgui.ydata = []
    plotgui.set_properties = []
    for loop in range(plotgui.max_sets):
        plotgui.set_properties.append({
            'symbol': None,
            'symbolsize': 4.0, 'linestyle': 'None',
            'linewidth': 1.0, 'colour': 'black',
            'label': '', 'xmin': 0.0, 'xmax': 1.0, 'ymin': 0.0,
            'ymax': 1.0, 'display': True, 'errors': False,
            'legend': True, 'plot': 1})
        plotgui.xdata.append(None)
        plotgui.ydata.append(None)
    plotgui.plot_range = [[0., 1., 0., 1.], ]
    plotgui.original_range = []
    for loop in range(plotgui.max_sets):
        plotgui.original_range.append(True)
    plotgui.nsets = 0
    plotgui.title = [' ', ]
    plotgui.xparameters = [{'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
                            'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
                            'invert': 0, 'hide': 0, 'hideticks': 0,
                            'hidelabels': 0, 'hybridlog': 0,
                            'inverseticks': 0, 'ticklength': 6,
                            'bothticks': 0, 'minorticks': 0,
                            'oppositeaxis': 0, 'majorgridlines': 0,
                            'minorgridlines': 0}, ]
    plotgui.yparameters = [{'label': ' ', 'minimum': 0.0, 'maximum': 1.0,
                            'major': 0.1, 'minor': 0.1, 'logarithmic': 0,
                            'invert': 0, 'hide': 0, 'hideticks': 0,
                            'hidelabels': 0, 'hybridlog': 0,
                            'inverseticks': 0, 'ticklength': 6,
                            'bothticks': 0, 'minorticks': 0,
                            'oppositeaxis': 0, 'majorgridlines': 0,
                            'minorgridlines': 0}, ]
    plotgui.position_stack = []
    plotgui.label_flag = False
    plotgui.line_flag = False
    plotgui.ellipse_flag = False
    plotgui.box_flag = False
    plotgui.vector_flag = False
    plotgui.number_of_lines = 0
    plotgui.number_of_labels = 0
    plotgui.number_of_ellipses = 0
    plotgui.number_of_boxes = 0
    plotgui.number_of_vectors = 0
    plotgui.positions = []
    for loop in range(len(plotgui.legend_variable)):
        try:
            plotgui.legend_variable[loop].set(0)
        except Exception:
            plotgui.legend_variable[loop] = None
    plotgui.legend_handles = [None, ]
    plotgui.legend_labels = [None, ]
    plotgui.legend_frame = [None, ]
    plotgui.legend_position = [None, ]
    plotgui.legend_user_position = [None, ]
    plotgui.plot_margin = 0.0
    plotgui.plot_frame = [0.0, ]
    plotgui.nxplots = 1
    plotgui.nyplots = 1
    plotgui.number_of_plots = plotgui.nxplots * plotgui.nyplots
    plotgui.current_plot = 1
    if len(plotgui.subplot) > 1:
        for loop in range(1, len(plotgui.subplot)):
            plotgui.subplot[loop].clear()
        plotgui.subplot = [plotgui.subplot[0], ]
        plotgui.hide_subplot = [plotgui.hide_subplot[0], ]
    make_plot.make_plot(plotgui)

def clear_current_plot(plotgui):
    """
    Clear the current plot if the user OKs this action.

    This routine clears the sets, parameters and plot items for the
    currently selected plot.  All values return to their initial values.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    if (plotgui.nxplots == 1) and (plotgui.nyplots == 1):
        clear_plot(plotgui)
    else:
        response = tkinter.messagebox.askyesno(
            "Verify",
            "Do you want to abandon the current plot?")
        if not response:
            return
        if plotgui.nsets > 0:
            for loop in range(plotgui.max_sets):
                if plotgui.set_properties[loop]['plot'] == plotgui.current_plot:
                    plotgui.set_properties[loop] = {
                        'symbol': None, 'symbolsize': 4.0,
                        'linestyle': 'None', 'linewidth': 1.0,
                        'colour': 'black', 'label': '', 'xmin': 0.0,
                        'xmax': 1.0, 'ymin': 0.0, 'ymax': 1.0,
                        'display': True, 'errors': False, 'legend': True,
                        'plot': 1}
                    plotgui.xdata[loop] = None
                    plotgui.ydata[loop] = None
                    plotgui.original_range[loop] = True
        plotgui.plot_range[plotgui.current_plot-1] = [0., 1., 0., 1.]
        plotgui.original_range[plotgui.current_plot-1] = True
        plotgui.title[plotgui.current_plot-1] = ' '
        plotgui.xparameters[plotgui.current_plot-1] = {
            'label': ' ', 'minimum': 0.0, 'maximum': 1.0, 'major': 0.1,
            'minor': 0.1, 'logarithmic': 0, 'invert': 0, 'hide': 0,
            'hideticks': 0, 'hidelabels': 0, 'hybridlog': 0,
            'inverseticks': 0, 'ticklength': 6, 'bothticks': 0,
            'minorticks': 0, 'oppositeaxis': 0,
            'majorgridlines': 0, 'minorgridlines': 0}
        plotgui.yparameters[plotgui.current_plot-1] = {
            'label': ' ', 'minimum': 0.0, 'maximum': 1.0, 'major': 0.1,
            'minor': 0.1, 'logarithmic': 0, 'invert': 0, 'hide': 0,
            'hideticks': 0, 'hidelabels': 0, 'hybridlog': 0,
            'inverseticks': 0, 'ticklength': 6, 'bothticks': 0,
            'minorticks': 0, 'oppositeaxis': 0,
            'majorgridlines': 0, 'minorgridlines': 0}
        try:
            plotgui.legend_variable[plotgui.current_plot-1].set(0)
        except Exception:
            plotgui.legend_variable[plotgui.current_plot-1] = None
    plotgui.legend_handles[plotgui.current_plot-1] = None
    plotgui.legend_labels[plotgui.current_plot-1] = None
    plotgui.legend_frame[plotgui.current_plot-1] = None
    plotgui.legend_position[plotgui.current_plot-1] = None
    plotgui.legend_user_position[plotgui.current_plot-1] = None
    plotgui.plot_frame[plotgui.current_plot-1] = 0.0
    object_utilities.clear_plot_objects(plotgui)
    clear_sets(plotgui)
    plotgui.current_plot = 1
    make_plot.make_plot(plotgui)

def clear_sets(plotgui):
    """
    Clear inactive sets from the set variable.

    Reset the set values to remove ones that have been removed with a plot.
    These are marked by plotgui.xdata[plotgui.current_set-1]['values'] set to
    None.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    nsets = 0
    properties = []
    xdata = []
    ydata = []
    original_range = []
    for loop in range(plotgui.nsets):
        if plotgui.xdata[loop]['values'] is not None:
            xdata.append(plotgui.xdata[loop])
            ydata.append(plotgui.ydata[loop])
            properties.append(plotgui.set_properties[loop])
            original_range.append(plotgui.original_range[loop])
            nsets = nsets + 1
    for loop in range(nsets, plotgui.max_sets):
        xdata.append(None)
        ydata.append(None)
        original_range.append(True)
        properties.append({
            'symbol': None, 'symbolsize': 4.0,
            'linestyle': 'None', 'linewidth': 1.0, 'colour': 'black',
            'label': '', 'xmin': 0.0, 'xmax': 1.0, 'ymin': 0.0,
            'ymax': 1.0, 'display': True, 'errors': False,
            'legend': True, 'plot': 1})
    plotgui.nsets = nsets
    plotgui.xdata = xdata
    plotgui.ydata = ydata
    plotgui.set_properties = properties
    plotgui.original_range = original_range

