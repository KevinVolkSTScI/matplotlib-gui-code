import tkinter as Tk
import tkinter.messagebox
import general_utilities
import make_plot

def make_plot_control_window(plotgui):
    """
    Create the plot control window (tick marks, axis labels, etc).

    This routine produces the plot control window wherein one sets the plot
    properties such as the axis labels and the title.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

    """
    matplotlib_line_list = ['-', '--', '-.', ':', None]
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    if plotgui.plot_control_window is not None:
        return
    plotgui.plot_control_window = Tk.Toplevel()
    plotgui.plot_control_window.title('Plot Parameters')
    outframe = Tk.Frame(plotgui.plot_control_window)
    outframe.pack(side=Tk.TOP)
    holder = Tk.Frame(outframe)
    holder.pack(side=Tk.LEFT)
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Plot Title:')
    label.pack(side=Tk.LEFT)
    plotgui.title_field = Tk.Entry(field1, width=20)
    plotgui.title_field.pack(side=Tk.LEFT)
    plotgui.title_field.insert(0, plotgui.title[plotgui.current_plot-1])
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='X Axis Label:')
    label.pack(side=Tk.LEFT)
    plotgui.xlabel_field = Tk.Entry(field1, width=20)
    plotgui.xlabel_field.pack(side=Tk.LEFT)
    plotgui.xlabel_field.insert(
        0, plotgui.xparameters[plotgui.current_plot-1]['label'])
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Y Axis Label:')
    label.pack(side=Tk.LEFT)
    plotgui.ylabel_field = Tk.Entry(field1, width=20)
    plotgui.ylabel_field.pack(side=Tk.LEFT)
    plotgui.ylabel_field.insert(
        0, plotgui.yparameters[plotgui.current_plot-1]['label'])
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='X Axis Minimum:')
    label.pack(side=Tk.LEFT)
    plotgui.xmin_field = Tk.Entry(field1, width=20)
    plotgui.xmin_field.pack(side=Tk.LEFT)
    if plotgui.xparameters[plotgui.current_plot-1]['hybridlog'] == 1:
        xlim = general_utilities.inverse_hybrid_transform(
            plotgui.xparameters[plotgui.current_plot-1]['minimum'])
        plotgui.xmin_field.insert(0, xlim)
    else:
        plotgui.xmin_field.insert(
            0, plotgui.xparameters[plotgui.current_plot-1]['minimum'])
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='X Axis Maximum:')
    label.pack(side=Tk.LEFT)
    plotgui.xmax_field = Tk.Entry(field1, width=20)
    plotgui.xmax_field.pack(side=Tk.LEFT)
    if plotgui.xparameters[plotgui.current_plot-1]['hybridlog'] == 1:
        xlim = general_utilities.inverse_hybrid_transform(
            plotgui.xparameters[plotgui.current_plot-1]['maximum'])
        plotgui.xmax_field.insert(0, xlim)
    else:
        plotgui.xmax_field.insert(
            0, plotgui.xparameters[plotgui.current_plot-1]['maximum'])
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Y Axis Minimum:')
    label.pack(side=Tk.LEFT)
    plotgui.ymin_field = Tk.Entry(field1, width=20)
    plotgui.ymin_field.pack(side=Tk.LEFT)
    if plotgui.yparameters[plotgui.current_plot-1]['hybridlog'] == 1:
        ylim = general_utilities.inverse_hybrid_transform(
            plotgui.yparameters[plotgui.current_plot-1]['minimum'])
        plotgui.ymin_field.insert(0, ylim)
    else:
        plotgui.ymin_field.insert(
            0, plotgui.yparameters[plotgui.current_plot-1]['minimum'])
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Y Axis Maximum:')
    label.pack(side=Tk.LEFT)
    plotgui.ymax_field = Tk.Entry(field1, width=20)
    plotgui.ymax_field.pack(side=Tk.LEFT)
    if plotgui.yparameters[plotgui.current_plot-1]['hybridlog'] == 1:
        ylim = general_utilities.inverse_hybrid_transform(
            plotgui.yparameters[plotgui.current_plot-1]['maximum'])
        plotgui.ymax_field.insert(0, ylim)
    else:
        plotgui.ymax_field.insert(
            0, plotgui.yparameters[plotgui.current_plot-1]['maximum'])
    holder1 = Tk.Frame(outframe)
    holder1.pack(side=Tk.LEFT)
    # Note it is not possible to put in the line height
    # automatically as the frames it is intended to match has
    # not been packed yet.
    sl = general_utilities.separator_line(holder1, 5, 550, 5, False)
    holder2 = Tk.Frame(outframe)
    holder2.pack(side=Tk.LEFT)
    flag = plotgui.xparameters[plotgui.current_plot-1]['majorgridlines'] == 1
    plotgui.majorxgrid_variable = general_utilities.add_yes_no_field(
        holder2, "major x grid", flag) 
    flag = plotgui.xparameters[plotgui.current_plot-1]['minorgridlines'] == 1
    plotgui.minorxgrid_variable = general_utilities.add_yes_no_field(
        holder2, "minor x grid", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['majorgridlines'] == 1
    plotgui.majorygrid_variable = general_utilities.add_yes_no_field(
        holder2, "major y grid", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['minorgridlines'] == 1
    plotgui.minorygrid_variable = general_utilities.add_yes_no_field(
        holder2, "minor y grid", flag) 
    sl = general_utilities.separator_line(holder2, 200, 5, 5, True)
    flag = plotgui.xparameters[plotgui.current_plot-1]['logarithmic'] == 1
    plotgui.logx_variable = general_utilities.add_yes_no_field(
        holder2, "logarithmic x axis", flag) 
    flag = plotgui.xparameters[plotgui.current_plot-1]['hybridlog'] == 1
    plotgui.hlogx_variable = general_utilities.add_yes_no_field(
        holder2, "hybrid log x axis", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['logarithmic'] == 1
    plotgui.logy_variable = general_utilities.add_yes_no_field(
        holder2, "logarithmic y axis", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['hybridlog'] == 1
    plotgui.hlogy_variable = general_utilities.add_yes_no_field(
        holder2, "hybrid log y axis", flag) 
    sl = general_utilities.separator_line(holder2, 200, 5, 5, True)
    flag = plotgui.xparameters[plotgui.current_plot-1]['invert'] == 1
    plotgui.invertx_variable = general_utilities.add_yes_no_field(
        holder2, "invert x axis", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['invert'] == 1
    plotgui.inverty_variable = general_utilities.add_yes_no_field(
        holder2, "invert y axis", flag) 
    flag = plotgui.xparameters[plotgui.current_plot-1]['hide'] == 1
    plotgui.hidex_variable = general_utilities.add_yes_no_field(
        holder2, "hide x axis", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['hide'] == 1
    plotgui.hidey_variable = general_utilities.add_yes_no_field(
        holder2, "hide y axis", flag) 
    flag = plotgui.xparameters[plotgui.current_plot-1]['oppositeaxis'] == 1
    plotgui.oppositex_variable = general_utilities.add_yes_no_field(
        holder2, "opposite x axis", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['oppositeaxis'] == 1
    plotgui.oppositey_variable = general_utilities.add_yes_no_field(
        holder2, "opposite y axis", flag) 
    sl = general_utilities.separator_line(holder2, 200, 5, 5, True)
    flag = plotgui.xparameters[plotgui.current_plot-1]['hidelabels'] == 1
    plotgui.hidexlabels_variable = general_utilities.add_yes_no_field(
        holder2, "hide x labels", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['hidelabels'] == 1
    plotgui.hideylabels_variable = general_utilities.add_yes_no_field(
        holder2, "hide y labels", flag) 
    flag = plotgui.xparameters[plotgui.current_plot-1]['hideticks'] == 1
    plotgui.hidexticks_variable = general_utilities.add_yes_no_field(
        holder2, "hide x ticks", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['hideticks'] == 1
    plotgui.hideyticks_variable = general_utilities.add_yes_no_field(
        holder2, "hide y ticks", flag) 
    flag = plotgui.xparameters[plotgui.current_plot-1]['bothticks'] == 1
    plotgui.bothxticks_variable = general_utilities.add_yes_no_field(
        holder2, "x ticks both sides", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['bothticks'] == 1
    plotgui.bothyticks_variable = general_utilities.add_yes_no_field(
        holder2, "y ticks both sides", flag) 
    flag = plotgui.xparameters[plotgui.current_plot-1]['inverseticks'] == 1
    plotgui.inversexticks_variable = general_utilities.add_yes_no_field(
        holder2, "invert x ticks", flag) 
    flag = plotgui.yparameters[plotgui.current_plot-1]['inverseticks'] == 1
    plotgui.inverseyticks_variable = general_utilities.add_yes_no_field(
        holder2, "invert y ticks", flag) 
    #
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Tick length:')
    label.pack(side=Tk.LEFT)
    plotgui.ticklengthfield = Tk.Entry(field1, width=3)
    plotgui.ticklengthfield.pack()
    try:
        plotgui.ticklengthfield.insert(
            0, str(plotgui.xparameters[plotgui.current_plot-1]['ticklength']))
    except Exception:
        plotgui.ticklengthfield.insert(0, '6')
        plotgui.xparameters[plotgui.current_plot-1]['ticklength'] = 6
        plotgui.yparameters[plotgui.current_plot-1]['ticklength'] = 6
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='X Minor tick interval:')
    label.pack(side=Tk.LEFT)
    plotgui.xminortickfield = Tk.Entry(field1, width=10)
    plotgui.xminortickfield.pack()
    plotgui.xminortickfield.insert(
        0, str(plotgui.xparameters[plotgui.current_plot-1]['minorticks']))
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Y Minor tick interval:')
    label.pack(side=Tk.LEFT)
    plotgui.yminortickfield = Tk.Entry(field1, width=10)
    plotgui.yminortickfield.pack()
    plotgui.yminortickfield.insert(
        0, str(plotgui.yparameters[plotgui.current_plot-1]['minorticks']))
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Plot Legend:')
    label.pack(side=Tk.LEFT)
    if plotgui.legend_variable[plotgui.current_plot-1] is None:
        plotgui.legend_variable[plotgui.current_plot-1] = Tk.IntVar()
        flag = False
    else:
        flag = plotgui.legend_variable[plotgui.current_plot-1].get()
    b1 = Tk.Frame(field1)
    general_utilities.put_yes_no(
        b1, plotgui.legend_variable[plotgui.current_plot-1],
        ['Yes', 'No'], flag)
    b1.pack()
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Legend type:')
    label.pack(side=Tk.LEFT)
    if plotgui.legend_options[plotgui.current_plot-1] is None:
        plotgui.legend_options[plotgui.current_plot-1] = Tk.StringVar()
        plotgui.legend_options[plotgui.current_plot-1].set('best')
        plotgui.legend_position[plotgui.current_plot-1] = 'best'
    else:
        plotgui.legend_options[plotgui.current_plot-1].set(
            plotgui.legend_position[plotgui.current_plot-1])
    plotgui.legend_option_list = ['user', 'best', 'upper right',
                                  'upper left', 'lower left',
                                  'lower right', 'right',
                                  'center left', 'center right',
                                  'lower center', 'upper center',
                                  'center']
    menu1 = Tk.OptionMenu(field1,
                          plotgui.legend_options[plotgui.current_plot-1],
                          *plotgui.legend_option_list,
                          command=plotgui.generate_legend)
    menu1.config(width=10)
    menu1.pack(side=Tk.LEFT)
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='User legend position:')
    label.pack(side=Tk.LEFT)
    plotgui.legend_position_field = Tk.Entry(field1, width=20)
    plotgui.legend_position_field.pack(side=Tk.LEFT)
    if plotgui.legend_user_position[plotgui.current_plot-1] is None:
        plotgui.legend_position_field.insert(0, '0.1 0.1')
    else:
        str1 = '%.3f %.3f' % (
            plotgui.legend_user_position[plotgui.current_plot-1][0],
            plotgui.legend_user_position[plotgui.current_plot-1][1])
        plotgui.legend_position_field.insert(0, str1)
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Legend frame:')
    label.pack(side=Tk.LEFT)
    if plotgui.legend_frame[plotgui.current_plot-1] is None:
        plotgui.legend_frame[plotgui.current_plot-1] = Tk.IntVar()
        plotgui.legend_frame[plotgui.current_plot-1].set(0)
        flag = False
    else:
        flag = plotgui.legend_frame[plotgui.current_plot-1].get()
    b1 = Tk.Frame(field1)
    general_utilities.put_yes_no(
        b1, plotgui.legend_frame[plotgui.current_plot-1],
        ['Yes', 'No'], flag)
    b1.pack()
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Plot margin:')
    label.pack(side=Tk.LEFT)
    plotgui.plot_margin_field = Tk.Entry(field1, width=10)
    plotgui.plot_margin_field.pack(side=Tk.LEFT)
    if plotgui.plot_margin is None:
        plotgui.plot_margin_field.insert(0, '0.0')
        plotgui.plot_margin = 0.0
    else:
        plotgui.plot_margin_field.insert(0, '%f' % (plotgui.plot_margin))
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Frame Width:')
    label.pack(side=Tk.LEFT)
    plotgui.plot_frame_field = Tk.Entry(field1, width=10)
    plotgui.plot_frame_field.pack(side=Tk.LEFT)
    if (plotgui.plot_frame[plotgui.current_plot-1] is None) or \
       (plotgui.plot_frame[plotgui.current_plot-1] <= 0.):
        plotgui.plot_frame_field.insert(0, '0.0')
        plotgui.plot_frame[plotgui.current_plot-1] = 0.0
    else:
        plotgui.plot_frame_field.insert(0, '%f' % (
            plotgui.plot_frame[plotgui.current_plot-1]))
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Grid Colour:')
    label.pack(side=Tk.LEFT)
    if plotgui.grid_colour_variable is None:
        plotgui.grid_colour_variable = Tk.StringVar()
        plotgui.grid_colour_variable.set('black')
    if plotgui.grid_colour[plotgui.current_plot-1] is None:
        plotgui.grid_colour[plotgui.current_plot-1] = 'black'
    try:
        plotgui.grid_colour_variable.set(
            plotgui.grid_colour[plotgui.current_plot-1])
    except:
        pass
    menu1 = Tk.OptionMenu(
        field1, plotgui.grid_colour_variable,
        *plotgui.colourset, command=plotgui.set_grid_colour)
    menu1.pack()
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Grid Line Type:')
    label.pack(side=Tk.LEFT)
    if plotgui.grid_linetype_variable is None:
        plotgui.grid_linetype_variable = Tk.StringVar()
        plotgui.grid_linetype_variable.set('solid')
    if plotgui.grid_linetype[plotgui.current_plot-1] is None:
        plotgui.grid_linetype[plotgui.current_plot-1] = '--'
    for loop in range(len(matplotlib_line_list)):
        if matplotlib_line_list[loop] == \
           plotgui.grid_linetype[plotgui.current_plot-1]:
            plotgui.grid_linetype_variable.set(matplotlib_line_name_list[loop])
    menu1 = Tk.OptionMenu(
        field1, plotgui.grid_linetype_variable,
        *matplotlib_line_name_list, command=plotgui.set_grid_linetype)
    menu1.pack()
    holder = Tk.Frame(plotgui.plot_control_window)
    holder.pack(side=Tk.TOP)
    # find the width for the separator line and apply it...
    outframe.update()
    sl = general_utilities.separator_line(
        holder, outframe.winfo_width(), 5, 5, True)
    field1 = Tk.Frame(plotgui.plot_control_window)
    field1.pack(side=Tk.TOP)
    apply_button = Tk.Button(field1, text="Apply",
                             command=lambda: apply_plot_parameters(plotgui))
    apply_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(field1, text="    ")
    label1.pack(side=Tk.LEFT)
    apply_all_button = Tk.Button(
        field1, text="Apply to All",
        command=lambda: apply_plot_parameters_all(plotgui))
    apply_all_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(field1, text="    ")
    label1.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        field1, text="Close",
        command=lambda: plotgui.close_window(plotgui.plot_control_window, 
                                             'plot_control_window'))
    close_button.pack(side=Tk.LEFT)

def apply_plot_parameters_all(plotgui):
    """
    Apply the plot parameters to all active plots.

    Parameters
    ----------



    Returns
    -------

        None

    """
    matplotlib_line_list = ['-', '--', '-.', ':', None]
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    plotnumber = plotgui.current_plot
    for loop in range(plotgui.number_of_plots):
        plotgui.current_plot = loop+1
        if plotgui.legend_options[plotgui.current_plot-1] is None:
            plotgui.legend_options[plotgui.current_plot-1] = Tk.StringVar()
            plotgui.legend_options[plotgui.current_plot-1].set('best')
            plotgui.legend_position[plotgui.current_plot-1] = 'best'
        plotgui.apply_plot_parameters()
    plotgui.current_plot = plotnumber

def apply_plot_parameters(plotgui):
    """
    Apply the values from the plot parameters window.

    This routine is called when the values in the plot control window
    are to be applied to the plot.  It reads the values and places them
    in the internal variables, and then replots the data.

    Parameters
    ----------


    Returns
    -------

        None

    """
    matplotlib_line_list = ['-', '--', '-.', ':', None]
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    try:
        plotgui.title[plotgui.current_plot-1] = plotgui.title_field.get()
        instring = plotgui.xlabel_field.get()
        # To do: test the string is OK as a label here....
        plotgui.xparameters[plotgui.current_plot-1]['label'] = instring
        instring = plotgui.ylabel_field.get()
        plotgui.yparameters[plotgui.current_plot-1]['label'] = instring
        plotgui.xparameters[plotgui.current_plot-1]['logarithmic'] = \
            plotgui.logx_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['logarithmic'] = \
            plotgui.logy_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['hybridlog'] = \
            plotgui.hlogx_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['hybridlog'] = \
            plotgui.hlogy_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['invert'] = \
            plotgui.invertx_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['invert'] = \
            plotgui.inverty_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['hide'] = \
            plotgui.hidex_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['hide'] = \
            plotgui.hidey_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['hideticks'] = \
            plotgui.hidexticks_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['hideticks'] = \
            plotgui.hideyticks_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['hidelabels'] = \
            plotgui.hidexlabels_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['hidelabels'] = \
            plotgui.hideylabels_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['inverseticks'] = \
            plotgui.inversexticks_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['inverseticks'] = \
            plotgui.inverseyticks_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['bothticks'] = \
            plotgui.bothxticks_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['bothticks'] = \
            plotgui.bothyticks_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['oppositeaxis'] = \
            plotgui.oppositex_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['oppositeaxis'] = \
            plotgui.oppositey_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['majorgridlines'] = \
            plotgui.majorxgrid_variable.get()
        plotgui.xparameters[plotgui.current_plot-1]['minorgridlines'] = \
            plotgui.minorxgrid_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['majorgridlines'] = \
            plotgui.majorygrid_variable.get()
        plotgui.yparameters[plotgui.current_plot-1]['minorgridlines'] = \
            plotgui.minorygrid_variable.get()
        try:
            plotgui.xparameters[plotgui.current_plot-1]['minorticks'] = \
                    float(plotgui.xminortickfield.get())
        except ValueError:
            plotgui.xparameters[plotgui.current_plot-1]['minorticks'] = 0.0
        if plotgui.xparameters[plotgui.current_plot-1]['minorticks'] < 0.0:
            plotgui.xparameters[plotgui.current_plot-1]['minorticks'] = 0.0
        try:
            plotgui.yparameters[plotgui.current_plot-1]['minorticks'] = \
                float(plotgui.yminortickfield.get())
        except ValueError:
            plotgui.yparameters[plotgui.current_plot-1]['minorticks'] = 0.0
        if plotgui.yparameters[plotgui.current_plot-1]['minorticks'] < 0.0:
            plotgui.yparameters[plotgui.current_plot-1]['minorticks'] = 0.0
        try:
            ticklength = int(plotgui.ticklengthfield.get())
            if ticklength < 1:
                ticklength = 1
        except ValueError:
            ticklength = 6
        plotgui.xparameters[plotgui.current_plot-1]['ticklength'] = ticklength
        plotgui.yparameters[plotgui.current_plot-1]['ticklength'] = ticklength
        try:
            frame = float(plotgui.plot_frame_field.get())
            if (frame <= 0.) or (frame > 5.):
                plotgui.plot_frame[plotgui.current_plot-1] = 0.
            else:
                plotgui.plot_frame[plotgui.current_plot-1] = frame
        except ValueError:
            plotgui.plot_frame[plotgui.current_plot-1] = 0.
        try:
            margin = float(plotgui.plot_margin_field.get())
            plotgui.plot_margin = margin
            if margin > 0.3:
                tkinter.messagebox.showinfo(
                    'Warning',
                    'Margins are limited to the range -0.1 to +0.3.')
                margin = 0.3
                plotgui.plot_margin_field.delete(0, Tk.END)
                plotgui.plot_margin_field.insert(0, '0.3')
            if margin < -0.1:
                tkinter.messagebox.showinfo(
                    'Warning',
                    'Margins are limited to the range -0.1 to +0.3.')
                margin = -0.1
                plotgui.plot_margin_field.delete(0, Tk.END)
                plotgui.plot_margin_field.insert, (0, '-0.1')
            plotgui.plot_margin = margin
        except Exception:
            tkinter.messagebox.showinfo(
                'Warning',
                'Margin value could not be read, setting to zero.')
            plotgui.plot_margin_field.delete(0, Tk.END)
            plotgui.plot_margin_field.insert, (0, '0.0')
        try:
            xmin = float(plotgui.xmin_field.get())
            xmax = float(plotgui.xmax_field.get())
            ymin = float(plotgui.ymin_field.get())
            ymax = float(plotgui.ymax_field.get())
            if (xmin >= xmax) | (ymin >= ymax):
                raise ValueError
            plotgui.plot_range[plotgui.current_plot-1][0] = xmin
            plotgui.plot_range[plotgui.current_plot-1][1] = xmax
            plotgui.plot_range[plotgui.current_plot-1][2] = ymin
            plotgui.plot_range[plotgui.current_plot-1][3] = ymax
        except Exception:
            tkinter.messagebox.showinfo(
                'Error',
                'There was some error in the axis range values.'
                + '  Please check your inputs.')
        plotgui.legend_position[plotgui.current_plot-1] = \
            plotgui.legend_options[plotgui.current_plot-1].get()
        str1 = plotgui.grid_linetype_variable.get()
        for loop in range(len(matplotlib_line_list)):
            if matplotlib_line_name_list[loop] == str1:
                plotgui.grid_linetype[plotgui.current_plot-1] = \
                    matplotlib_line_list[loop]
        for loop in range(len(plotgui.share_axis)):
            if plotgui.current_plot == abs(plotgui.share_axis[loop]):
                if plotgui.share_axis[loop] < 0:
                    for key in plotgui.yparameters[-1].keys():
                        plotgui.yparameters[loop][key] = plotgui.yparameters[
                            plotgui.current_plot-1][key]
                else:
                    for key in plotgui.xparameters[-1].keys():
                        plotgui.xparameters[loop][key] = plotgui.xparameters[
                            plotgui.current_plot-1][key]
        make_plot.make_plot(plotgui)
    except Exception:
        tkinter.messagebox.showinfo(
            'Error',
            'There was some error in the input values.  '
            + 'Please check your inputs.')

