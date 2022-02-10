"""
Routines to make the main window.  One need to pass in the
matplotlib_user_interface object to the routines.

Routines:

make_widget      Create the main GUI window

make_menus       Create the pull-down menus for the main window

set_font         Bring up a window to set the base font

set_font_values_all    Read the font window parameters and apply to all plots

set_fonr_values        Read the font window parameters and apply to the
                       current plot

make_controls     Create the plot controls area of the main GUI window

set_plot_layout   Read the plot layout parameters and apply them

make_plot_area    Make the main Figure area within the plot GUI

read_data_set     Bring up a window for reading in a data set from a file

select_file_window   Select a file and read data columns from it

check_columns     Check the columns of data in an ascii file for numerical
                  values

get_set           Read in data values and assign them to a data set


"""
import tkinter as Tk
import tkinter.filedialog
import tkinter.messagebox
import numpy
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import data_set_utilities
import data_set_operations
import general_utilities
import edit_objects
import make_plot
import label_utilities
import object_utilities
import plot_flag_utilities
import plot_controls
import histogram_utilities
import save_and_restore_plot
import non_linear_fitting

def make_widget(plotgui):
    """
    Create the main GUI window.

    This routine makes the main GUI window.  It uses the root variable
    passed to the class when it is initialized (although I do not know
    how this works as the value is not passed to this routine...but it
    does work).

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object


    Returns
    -------

        None

    """
    menuframe = Tk.Frame(plotgui.root)
    menuframe.pack(side=Tk.TOP, anchor=Tk.W)
    make_menus(plotgui, menuframe)
    controlframe = Tk.Frame(plotgui.root)
    controlframe.pack(side=Tk.LEFT, fill=Tk.Y, expand=1)
    make_controls(plotgui, controlframe)
    sl = general_utilities.separator_line(
        plotgui.root, 5, 750, 5, False, Tk.LEFT)
    plotgui.plotframe = Tk.Frame(plotgui.root)
    plotgui.plotframe.pack(side=Tk.LEFT, fill=Tk.Y, expand=1)
    make_plot_area(plotgui, plotgui.plotframe)

def make_menus(plotgui, parent):
    """
    Create pull-down menus for plot functionality.

    Given a Tk Frame variable "parent" this routine makes a pull-down
    menu area within this frame.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

        parent     A Tk.Frame variable, that holds the menus

    Returns
    -------

        None

    """
    menubutton1 = Tk.Menubutton(parent, text="Data")
    menubutton1.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
    menu1 = Tk.Menu(menubutton1)
    menubutton1['menu'] = menu1
    menu1.add_command(label='Read Data',
                      command=lambda: read_data_set(plotgui))
    menu1.add_command(
        label='Create Values by Formula',
        command=lambda: data_set_utilities.create_data_set(plotgui))
    menu1.add_command(
        label='Create Values in Widget',
        command=lambda: data_set_utilities.create_data_set_by_editor(plotgui))
    menu1.add_command(
        label='Write Data',
        command=lambda: data_set_utilities.write_data_sets(plotgui))
    #menu1.add_separator()
    #menu1.add_command(label='Read FITS Image', command=plotgui.read_image)
    label1 = Tk.Label(parent, text="    ")
    label1.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
    menubutton2 = Tk.Menubutton(parent, text="Sets")
    menubutton2.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
    menu2 = Tk.Menu(menubutton2)
    menubutton2['menu'] = menu2
    menu2.add_command(
        label='Set Properties',
        command=lambda: data_set_operations.make_data_set_window(plotgui))
    menu2.add_command(
        label='Fit Sets',
        command=\
        lambda: data_set_operations.make_data_set_fitting_window(plotgui))
    menu2.add_command(
        label='Non-linear Fitting',
        command=\
        lambda: non_linear_fitting.make_fitting_window(plotgui))
    menu2.add_command(
        label='Set Statistics',
        command=lambda: data_set_utilities.set_statistics(plotgui))
    menu2.add_command(
        label='Transform Set',
        command=lambda:\
        data_set_operations.make_data_set_transformation_window(plotgui))
    menu2.add_command(
        label='Edit Set in Widget',
        command=lambda: data_set_operations.make_data_set_edit_window(plotgui))
    menu2.add_command(
        label='Sort Set',
        command=lambda: \
        data_set_operations.make_data_set_sort_window(plotgui))
    menu2.add_command(
        label='Delete Set',
        command=lambda:\
        data_set_operations.make_data_set_delete_window(plotgui))
    label2 = Tk.Label(parent, text="    ")
    label2.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
    menubutton3 = Tk.Menubutton(parent, text="Plot")
    menubutton3.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
    menu3 = Tk.Menu(menubutton3)
    menubutton3['menu'] = menu3
    menu3.add_command(label='Plot Parameters',
        command=lambda: plot_controls.make_plot_control_window(plotgui))
    menu3.add_command(
        label='Clear Plot',
        command=lambda: plot_flag_utilities.clear_plot(plotgui))
    menu3.add_command(
        label='Clear Current Plot',
        command=lambda: plot_flag_utilities.clear_current_plot(plotgui))
    menu3.add_command(
        label='Opposite Y Axis',
        command=lambda: plot_flag_utilities.set_opposite_y(plotgui))
    menu3.add_command(
        label='Opposite X Axis',
        command=lambda: plot_flag_utilities.set_opposite_x(plotgui))
    menu3.add_command(
        label='Tile Plots',
        command=lambda: plot_flag_utilities.tile_plots(plotgui))
    menu3.add_command(
        label='Set Plot',
        command=lambda: plot_flag_utilities.set_plot_number(plotgui))
    menu3.add_command(
        label='Hide/Show Plot',
        command=lambda: plot_flag_utilities.set_plot_hide(plotgui))
    menu3.add_command(
        label='Toggle Equal Aspect',
        command=lambda: plot_flag_utilities.toggle_equal_aspect(plotgui))
    label3 = Tk.Label(parent, text="    ")
    label3.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
    menubutton4 = Tk.Menubutton(parent, text="Plot Items")
    menubutton4.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
    menu4 = Tk.Menu(menubutton4)
    menubutton4['menu'] = menu4
    menu4.add_command(label='Set Font', command=lambda: set_font(plotgui))
    menu4.add_separator()
    menu4.add_command(label='Add a label',
                      command=lambda: object_utilities.set_label(plotgui))
    menu4.add_command(label='Add a line',
                      command=lambda: object_utilities.add_line(plotgui))
    menu4.add_command(label='Add an ellipse',
                      command=lambda: object_utilities.add_ellipse(plotgui))
    menu4.add_command(label='Add a box',
                      command=lambda: object_utilities.add_box(plotgui))
    menu4.add_command(label='Add a vector',
                      command=lambda: object_utilities.add_vector(plotgui))
    menu4.add_separator()
    menu4.add_command(label='Edit labels',
                      command=lambda: label_utilities.edit_labels(plotgui))
    menu4.add_command(label='Edit lines',
                      command=lambda: edit_objects.edit_lines(plotgui))
    menu4.add_command(label='Edit ellipses',
                      command=lambda: edit_objects.edit_ellipses(plotgui))
    menu4.add_command(label='Edit boxes',
                      command=lambda: edit_objects.edit_boxes(plotgui))
    menu4.add_command(label='Edit vectors',
                      command=lambda: edit_objects.edit_vectors(plotgui))
    menu4.add_separator()
    menu4.add_command(label='Remove last line',
                      command=lambda: object_utilities.remove_line(plotgui))
    menu4.add_command(
        label='Remove last ellipse',
        command=lambda: object_utilities.remove_ellipse(plotgui))
    menu4.add_command(label='Remove last box',
                      command=lambda: object_utilities.remove_box(plotgui))
    menu4.add_command(label='Remove last vector',
                      command=lambda: object_utilities.remove_vector(plotgui))
    menu4.add_separator()
    menu4.add_command(label='Clear all labels',
                      command=lambda: label_utilities.clear_labels(plotgui))
    menu4.add_command(
        label='Remove all lines',
        command=lambda: object_utilities.remove_all_lines(plotgui))
    menu4.add_command(
        label='Remove all ellipses',
        command=lambda: object_utilities.remove_all_ellipses(plotgui))
    menu4.add_command(
        label='Remove all boxes',
        command=lambda: object_utilities.remove_all_boxes(plotgui))
    menu4.add_command(
        label='Remove all vectors',
        command=lambda: object_utilities.remove_all_vectors(plotgui))
    label4 = Tk.Label(parent, text="    ")
    label4.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
    menubutton5 = Tk.Menubutton(parent, text="Save/Restore Plot")
    menubutton5.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
    menu5 = Tk.Menu(menubutton5)
    menubutton5['menu'] = menu5
    menu5.add_command(
        label='Save configuration',
        command=lambda: save_and_restore_plot.save_plot(plotgui))
    menu5.add_command(
        label='Read configuration',
        command=lambda: save_and_restore_plot.load_plot(plotgui))
    menu5.add_command(
        label='Save as PNG',
        command=lambda: general_utilities.save_png_figure(plotgui.figure))
    menu5.add_command(
        label='Save as postscript',
        command=lambda: general_utilities.save_ps_figure(plotgui.figure))

def set_font(plotgui):
    """
    Bring up a window to set the plot font.

    This window allows one to set the font properties for the axis
    labels and titles.  The values are applied to all axis labels
    and the title in a given plot.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

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
    plotgui.font_name_list = Tk.ttk.Combobox(holder, width=20)
    plotgui.font_name_list.grid(column=1, row=0)
    plotgui.font_name_list['values'] = fontnames
    plotgui.font_name_list.set(plotgui.fontname[plotgui.current_plot-1])
    plotgui.font_size_list = Tk.ttk.Combobox(holder, width=20)
    plotgui.font_size_list.grid(column=1, row=1)
    plotgui.font_size_list['values'] = fontsizes
    plotgui.font_size_list.set(plotgui.fontsize[plotgui.current_plot-1])
    plotgui.font_weight_list = Tk.ttk.Combobox(holder, width=20)
    plotgui.font_weight_list.grid(column=1, row=2)
    plotgui.font_weight_list['values'] = fontweights
    plotgui.font_weight_list.set(plotgui.fontweight[plotgui.current_plot-1])
    bholder = Tk.Frame(font_window)
    bholder.pack(side=Tk.TOP)
    set_button = Tk.Button(bholder, text='Set',
                           command=lambda: set_font_values(plotgui))
    set_button.pack()
    set_all_button = Tk.Button(bholder, text='Set for All',
                               command=lambda: set_font_values_all(plotgui))
    set_all_button.pack()
    close_button = Tk.Button(bholder, text='Close Window',
                             command=font_window.destroy)
    close_button.pack()

def set_font_values_all(plotgui):
    """
    Read the font field and apply them to all plots.

    This routine reads the font fields and saves them to the internal
    variables.  It then calls the routine to re-plot.  This is done for
    all current plots (including hidden ones).

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

    """
    plotnumber = plotgui.current_plot
    for loop in range(plotgui.number_of_plots):
        plotgui.current_plot = loop + 1
        set_font_values(plotgui)
    plotgui.current_plot = plotnumber

def set_font_values(plotgui):
    """
    Read the font fields and apply them to the plot.

    This routine reads the font fields and saves them to the internal
    variables.  It then calls the routine to re-plot.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

        """
    plotgui.fontname[plotgui.current_plot-1] = plotgui.font_name_list.get()
    plotgui.fontsize[plotgui.current_plot-1] = plotgui.font_size_list.get()
    plotgui.fontweight[plotgui.current_plot-1] = plotgui.font_weight_list.get()
    make_plot.make_plot(plotgui)

def make_controls(plotgui, parent):
    """
    Make the control area within the main window.

    This routine makes a control area within the main window, under
    frame "parent".  The overall root value is also passed here, but
    not currently used.  It may be neeed for orderly closing of the
    window depending on what is done in the main window, hence it is
    included here.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

        parent :   A Tk.Frame variable for the holder of the controla

    Returns
    -------

        None

    """
    holder = Tk.Frame(parent)
    holder.pack(side=Tk.TOP)
    label1 = Tk.Label(holder, text=' ')
    label1.pack(side=Tk.TOP, fill=Tk.X)
    button1 = Tk.Button(holder, text='Plot',
                        command=lambda: make_plot.make_plot(plotgui))
    button1.pack(side=Tk.TOP, fill=Tk.X)
    button2 = Tk.Button(holder, text='Auto scale',
                        command=plotgui.autoscale_plot)
    button2.pack(side=Tk.TOP, fill=Tk.X)
    sl = general_utilities.separator_line(holder, 300, 25, 5, True)
    button3 = Tk.Button(
        holder, text='2-D Histogram',
        command=lambda: histogram_utilities.make_hess_plot(plotgui))
    button3.pack(side=Tk.TOP, fill=Tk.X)
    field = Tk.Frame(holder)
    field.pack(side=Tk.TOP)
    label1 = Tk.Label(field, text='Number of pixels: ')
    label1.pack(side=Tk.LEFT)
    plotgui.npixelfield = Tk.Entry(field, width=5)
    plotgui.npixelfield.pack(side=Tk.TOP)
    plotgui.npixelfield.insert(0, '500')
    plotgui.overplotflag = Tk.IntVar()
    h1 = Tk.Frame(holder)
    h1.pack(side=Tk.TOP)
    plotgui.hessindividualhistogramflag = Tk.IntVar()
    b1 = Tk.Frame(h1)
    general_utilities.put_yes_no(b1, plotgui.hessindividualhistogramflag,
                       ['all sets', 'one set'], True)
    b1.pack(side=Tk.LEFT)
    b1 = Tk.Frame(h1)
    plotgui.hess_set_field = Tk.Entry(b1, width=5)
    plotgui.hess_set_field.insert(0, '1')
    plotgui.hess_set_field.pack(side=Tk.LEFT)
    b1.pack(side=Tk.LEFT)
    h1 = Tk.Frame(holder)
    h1.pack(side=Tk.TOP)
    label1 = Tk.Label(h1, text='Overplot Sets: ')
    label1.pack(side=Tk.LEFT)
    b1 = Tk.Frame(h1)
    general_utilities.put_yes_no(b1, plotgui.overplotflag,
                                 ['Yes', 'No'], True, )
    b1.pack(side=Tk.LEFT)
    sl = general_utilities.separator_line(holder, 300, 25, 5, True)
    button4 = Tk.Button(
        holder, text='1-D Histogram',
        command=lambda: histogram_utilities.make_histogram(plotgui))
    button4.pack(side=Tk.TOP, fill=Tk.X)
    field = Tk.Frame(holder)
    field.pack(side=Tk.TOP)
    label1 = Tk.Label(field, text='Number of bins or bin size: ')
    label1.pack(side=Tk.LEFT)
    plotgui.nbinfield = Tk.Entry(field, width=10)
    plotgui.nbinfield.pack(side=Tk.TOP)
    plotgui.nbinfield.insert(0, '500')
    plotgui.histogramflag = Tk.IntVar()
    b1 = Tk.Frame(holder)
    general_utilities.put_yes_no(b1, plotgui.histogramflag,
                       ['x values', 'y values'], True)
    b1.pack(side=Tk.TOP)
    plotgui.histogramflag = Tk.IntVar()
    b1 = Tk.Frame(holder)
    plotgui.individualhistogramflag = Tk.IntVar()
    general_utilities.put_yes_no(b1, plotgui.individualhistogramflag,
                       ['all sets', 'individual sets'], True)
    b1.pack(side=Tk.TOP)
    sl = general_utilities.separator_line(holder, 300, 25, 5, True)
    plotgui.matplotlib_rounding = Tk.IntVar()
    b1 = Tk.Frame(holder)
    label1 = Tk.Label(b1, text='Axis limits rounding algorithm: ')
    label1.pack(side=Tk.TOP)
    b1.pack(side=Tk.TOP)
    b1 = Tk.Frame(holder)
    general_utilities.put_yes_no(b1, plotgui.matplotlib_rounding,
                       ['Matplotlib', 'Alternate'], True)
    b1.pack(side=Tk.TOP)
    sl = general_utilities.separator_line(holder, 300, 25, 5, True)
    button5 = Tk.Button(
        holder, text='Save as PNG',
        command=lambda: general_utilities.save_png_figure(plotgui.figure))
    button5.pack(side=Tk.TOP, fill=Tk.X)
    button6 = Tk.Button(
        holder, text='Save as PS',
        command=lambda: general_utilities.save_ps_figure(plotgui.figure))
    button6.pack(side=Tk.TOP, fill=Tk.X)
    button7 = Tk.Button(
        holder, text='Clear Current Plot',
        command=lambda: plot_flag_utilities.clear_current_plot(plotgui))
    button7.pack(side=Tk.TOP, fill=Tk.X)
    button8 = Tk.Button(
        holder, text='Tile Plots',
        command=lambda: plot_flag_utilities.tile_plots(plotgui))
    button8.pack(side=Tk.TOP, fill=Tk.X)
    button9 = Tk.Button(
        holder, text='Set Plot',
        command=lambda: plot_flag_utilities.set_plot_number(plotgui))
    button9.pack(side=Tk.TOP, fill=Tk.X)
    button10 = Tk.Button(holder, text='Close Window',
                         command=plotgui.root.destroy)
    button10.pack(side=Tk.TOP, fill=Tk.X)

def set_plot_layout(plotgui):
    """
    Read the plot layout parmeters and apply them.

    This sub-routine reads the plot layout parameters and attempts
    to apply them.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

    """
    try:
        nx1 = int(plotgui.xplots_set.get())
        ny1 = int(plotgui.yplots_set.get())
        n1 = int(plotgui.currentplot_set.get())
        if (nx1 < 1) | (ny1 < 1) | (nx1 > 5) | (ny1 > 5) | (n1 < 1) \
           | (n1 > nx1*ny1):
            raise ValueError
        plotgui.make_plot_layout(nx1, ny1, n1)
    except ValueError:
        tkinter.messagebox.showinfo(
            "Error",
            "There was some issue with the plot layout parameters.")

def make_plot_area(plotgui, parent):
    """
    Set up the main figure area for the plot.

    This routine makes the figure area and the sub-plot, and then
    sets up some event call-backs.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

        parent :   A Tk.Frame variable that holds the plot

    Returns
    -------

        None

    """
    # The size of the plot area is determined here (values in inches).
    # The x/y ratio is 10 to 7 as in xmgrace.  One could also use the
    # golden ratio 1.618 to 1.  Take whatever values seem to be the
    # most aesthetically pleasing.  The DPI value should probably not
    # be set below 100.
    if plotgui.plot_area_flag:
        plotgui.figure = Figure(figsize=(7.142857, 5), dpi=150)
    else:
        plotgui.figure.clf()
    # The code allows a grid of plotgui.nxplots by plotgui.nyplots panels
    plotgui.bounding_box = []
    plotgui.subplot = []
    plotgui.hide_subplot = []
    plotgui.share_axis = []
    plot_number = 1
    for xloop in range(plotgui.nxplots):
        for yloop in range(plotgui.nyplots):
            plotgui.subplot.append(
                plotgui.figure.add_subplot(plotgui.nyplots,
                                           plotgui.nxplots, plot_number))
            plotgui.share_axis.append(0)
            plotgui.hide_subplot.append(False)
            bbox = plotgui.subplot[-1].get_position()
            bound_values = bbox.bounds
            plotgui.bounding_box.append(bound_values)
            plot_number = plot_number + 1
    if plotgui.plot_area_flag:
        # The following sets up a label above the plot, for the plot
        # position.
        plotgui.position_label_text = Tk.StringVar()
        plotgui.position_label = Tk.Label(
            parent, textvariable=plotgui.position_label_text)
        plotgui.position_label.pack(side=Tk.TOP)
        plotgui.position_label_text.set("Position:\n")
        plotgui.canvas = FigureCanvasTkAgg(plotgui.figure, master=parent)
        # Here are defined the events that the program responds to for
        # the plot area.
        plotgui.canvas.mpl_connect("motion_notify_event", plotgui.plot_position)
        plotgui.canvas.mpl_connect("button_press_event", plotgui.plot_marker_set)
        plotgui.canvas.mpl_connect("button_release_event",
                                   plotgui.plot_marker_release)
        plotgui.canvas.mpl_connect("key_press_event", plotgui.key_press_event)
        plotgui.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH,
                                            expand=Tk.YES)
        plotgui.plot_area_flag = False
    plotgui.canvas.draw()

def read_data_set(plotgui):
    """
    Make a window for reading in data sets.

    This routine produces a window to read in ascii data sets from a file.
    The window stays until one clicks on the "Close" button.  The window
    is regenerated each time the parent button is clicked.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------
        None

    """
    if plotgui.read_window is not None:
        return
    plotgui.labelstring = ' '
    plotgui.read_window = Tk.Toplevel()
    plotgui.read_window.title('Read in Data')
    holder = Tk.Frame(plotgui.read_window)
    holder.pack(side=Tk.TOP)
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Data File Filter')
    label.pack(side=Tk.LEFT)
    plotgui.file_template_field = Tk.Entry(field1, width=20)
    plotgui.file_template_field.pack(side=Tk.LEFT)
    plotgui.file_template_field.insert(0, '*')
    field2 = Tk.Frame(holder)
    field2.pack(side=Tk.TOP)
    label = Tk.Label(field2, text='Set Type')
    label.pack(side=Tk.LEFT)
    plotgui.set_options = Tk.StringVar()
    plotgui.set_options.set('XY')
    plotgui.set_option_list = ['XY', 'XYdY', 'XdXY', 'XdXYdY',
                               'XYdYdY', 'XdXdXY', 'XdXdXYdYdY']
    menu1 = Tk.OptionMenu(field2, plotgui.set_options,
                          *plotgui.set_option_list,
                          command=plotgui.error_fields)
    menu1.config(width=10)
    menu1.pack(side=Tk.LEFT)
    field3 = Tk.Frame(holder)
    field3.pack(side=Tk.TOP)
    label = Tk.Label(field3, text='Autoscale')
    label.pack(side=Tk.LEFT)
    plotgui.autoscale_options = Tk.StringVar()
    plotgui.autoscale_options.set('XY')
    plotgui.autoscale_option_list = ['XY', 'X', 'Y', 'None']
    menu1 = Tk.OptionMenu(field3, plotgui.autoscale_options,
                          *plotgui.autoscale_option_list)
    menu1.config(width=10)
    menu1.pack(side=Tk.LEFT)
    holder = Tk.Frame(plotgui.read_window)
    holder.pack(side=Tk.TOP)
    plotgui.read_window_label = Tk.Label(holder, text=plotgui.labelstring)
    plotgui.read_window_label.pack(side=Tk.TOP)
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='X data from column (0 for index): ')
    label.pack(side=Tk.LEFT)
    plotgui.xdata_field = Tk.Entry(field1, width=10)
    plotgui.xdata_field.pack(side=Tk.LEFT)
    plotgui.xdata_field.insert(0, '1')
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Y data from column (0 for index): ')
    label.pack(side=Tk.LEFT)
    plotgui.ydata_field = Tk.Entry(field1, width=10)
    plotgui.ydata_field.pack(side=Tk.LEFT)
    plotgui.ydata_field.insert(0, '2')
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='X uncertainties from column(s): ')
    label.pack(side=Tk.LEFT)
    plotgui.dxdata_field = Tk.Entry(field1, width=10)
    plotgui.dxdata_field.pack(side=Tk.LEFT)
    plotgui.dxdata_field.insert(0, '3')
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Y uncertainties from column(s): ')
    label.pack(side=Tk.LEFT)
    plotgui.dydata_field = Tk.Entry(field1, width=10)
    plotgui.dydata_field.pack(side=Tk.LEFT)
    plotgui.dydata_field.insert(0, '4')
    plotgui.error_fields(None)
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='All other columns as Y: ')
    label.pack(side=Tk.LEFT)
    plotgui.ally = Tk.IntVar()
    Tk.Radiobutton(field1, text='Yes', variable=plotgui.ally, value=1).pack()
    Tk.Radiobutton(field1, text='No', variable=plotgui.ally, value=0).pack()
    plotgui.ally.set(0)
    field1 = Tk.Frame(holder)
    label = Tk.Label(field1, text=' ')
    label.pack(side=Tk.TOP)
    field1.pack(side=Tk.TOP)
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text='Sort data points: ')
    label.pack(side=Tk.LEFT)
    plotgui.sortdata = Tk.IntVar()
    Tk.Radiobutton(field1, text='Sort on x',
                   variable=plotgui.sortdata, value=2).pack()
    Tk.Radiobutton(field1, text='Sort on y',
                   variable=plotgui.sortdata, value=1).pack()
    Tk.Radiobutton(field1, text='Do not sort',
                   variable=plotgui.sortdata, value=0).pack()
    plotgui.sortdata.set(0)
    field1 = Tk.Frame(holder)
    field1.pack(side=Tk.TOP)
    plotgui.labelsstring = ' '
    field4 = Tk.Frame(holder)
    field4.pack(side=Tk.TOP)
    select_button = Tk.Button(field4, text="Select File",
                              command=lambda: select_file_window(plotgui))
    select_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(field4, text="    ")
    label1.pack(side=Tk.LEFT)
    select_button = Tk.Button(
        field4, text="Get Values",
        command=lambda: \
        get_set(plotgui, plotgui.datavalues, plotgui.labelstring))
    select_button.pack(side=Tk.LEFT)
    label2 = Tk.Label(field4, text="    ")
    label2.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        field4, text="Cancel/Close",
        command=lambda: plotgui.close_window(plotgui.read_window,
                                             'read_window', True))
    close_button.pack(side=Tk.LEFT)

def select_file_window(plotgui):
    """
    Select a file and read data columns from it for input.

    This window queries the user for an input file name and tries to
    read the data therein.  If this is successful, some of the fields
    in the "Read Data" window are updated.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

        None

    """
    try:
        pattern = plotgui.file_template_field.get()
        values = pattern.split(', ')
        filetypes = []
        for loop in range(len(values)):
            label1 = values[loop].replace('*.', '')
            label2 = values[loop].replace('*', '')
            label1 = label1.upper()
            template1 = (label1, label2)
            filetypes.append(template1)
    except ValueError:
        tkinter.messagebox.showinfo(
            'Error', 'There was '
            + 'some error reading the file pattern.  Using the default.')
        pattern = '*'
    if (pattern == '*') | (pattern == ''):
        plotgui.filename = tkinter.filedialog.askopenfilename()
    else:
        try:
            plotgui.filename = tkinter.filedialog.askopenfilename(
                filetypes=filetypes)
        except Exception:
            plotgui.filename = tkinter.filedialog.askopenfilename()
    try:
        plotgui.shortname = None
        inds, ncols = check_columns(plotgui.filename)
        if inds is None:
            tkinter.messagebox.showinfo(
                'Error', 'There was '
                + 'some error reading file %s' % (plotgui.filename))
        else:
            plotgui.datavalues = numpy.loadtxt(
                plotgui.filename,
                comments=['#', '\\', '|', 'x_or_RA'], usecols=inds)
            values = plotgui.filename.split('/')
            plotgui.shortname = values[-1]
            plotgui.datashape = plotgui.datavalues.shape
            if len(plotgui.datashape) > 1:
                newlabel = '%d values ' % (plotgui.datashape[0])
                newlabel = newlabel + '(%s) ' % (plotgui.shortname)
                plotgui.labelstring = newlabel
                plotgui.ncolumns = plotgui.datashape[1]
            else:
                newlabel = '1 value '
                newlabel = newlabel + '(%s) ' % (plotgui.shortname)
                plotgui.labelstring = '(%s) ' % (plotgui.shortname)
                plotgui.ncolumns = plotgui.datashape[0]
            labeltext = '%d data columns, %d values ' % (
                plotgui.datashape[1], plotgui.datashape[0])
            plotgui.read_window_label['text'] = labeltext
            plotgui.data_indexes = inds
            plotgui.ndatacolumns = ncols
    except Exception:
        try:
            if plotgui.shortname is None:
                tkinter.messagebox.showinfo(
                    'Error',
                    'There was some error reading file %s' % (
                        plotgui.filename))
            else:
                tkinter.messagebox.showinfo(
                    'Error',
                    'There was some error reading file %s' % (
                        plotgui.shortname))
        except Exception:
            tkinter.messagebox.showinfo(
                'Error',
                'There was some error trying to read a file')

def check_columns(filename):
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
    except Exception:
        return None, None

def get_set(plotgui, datavalues, labelstring):
    """
    Extract data values to a new data set.

    This routine tries to make a data set (x, y) pair according to the
    specified values from the read sets window.

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

        datavalues :  This is a numpy array from the "loadtxt" function.

        labelstring : This is a string variable, used to set the labels
                      for the set(s)

    Returns
    -------

       None

    """
    if datavalues is None:
        tkinter.messagebox.showinfo(
            'Error',
            'You must specify a file before reading data.')
        return
    if len(datavalues.shape) == 1:
        datavalues = numpy.expand_dims(datavalues, 0)
    if True:
#    try:
        set_option_string = plotgui.set_options.get()
        if set_option_string is None:
            set_option_string = 'None'
        option1 = plotgui.set_option_list.index(set_option_string)
        autoscale_option_string = plotgui.autoscale_options.get()
        option2 = plotgui.autoscale_option_list.index(
            autoscale_option_string)
        option3 = plotgui.ally.get()
        option4 = plotgui.sortdata.get()
        shape1 = datavalues.shape
        xerrorinds = [2, 3, 5, 6]
        yerrorinds = [1, 3, 4, 6]
        if option3 == 0:
            try:
                inds = numpy.zeros((6), dtype=numpy.int8) - 1
                str1 = plotgui.xdata_field.get()
                n1 = int(str1) - 1
                inds[0] = n1
                str1 = plotgui.ydata_field.get()
                n1 = int(str1) - 1
                inds[1] = n1
                if option1 in xerrorinds:
                    str1 = plotgui.dxdata_field.get()
                    values = str1.split(', ')
                    n1 = int(values[0]) - 1
                    inds[2] = n1
                    if len(values) > 1:
                        n1 = int(values[1])-1
                        inds[3] = n1
                    else:
                        inds[3] = -1
                if option1 in yerrorinds:
                    str1 = plotgui.dydata_field.get()
                    values = str1.split(', ')
                    n1 = int(values[0]) - 1
                    inds[4] = n1
                    if len(values) > 1:
                        n1 = int(values[1])-1
                        inds[5] = n1
                    else:
                        inds[5] = -1
            # exit here if there is an exception setting values
            except Exception:
                raise ValueError
            newinds = inds*0 - 1
            for loop in range(len(inds)):
                for n in range(len(plotgui.data_indexes)):
                    if inds[loop] == plotgui.data_indexes[n]:
                        newinds[loop] = n
            if max(newinds) > plotgui.ndatacolumns:
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
            plotgui.add_set(xvalues, yvalues, xlowerror, xhigherror,
                            ylowerror, yhigherror, xerrorflag, yerrorflag,
                            option2, newlabelstring, plotgui.current_plot)
        else:
            str1 = plotgui.xdata_field.get()
            n1 = int(str1) - 1
            for loop in range(len(plotgui.data_indexes)):
                if n1 == plotgui.data_indexes[loop]:
                    xind = loop
            xvalues = numpy.squeeze(datavalues[:, xind])
            xlowerror = xvalues * 0.
            xhigherror = xvalues * 0.
            xerrorflag = False
            yerrorflag = False
            for loop in range(len(plotgui.data_indexes)):
                if loop == n1:
                    pass
                else:
                    yvalues = numpy.squeeze(datavalues[:, loop])
                    ylowerror = yvalues * 0.
                    yhigherror = yvalues * 0.
                    newlabelstring = 'columns %d and %d from file %s' % (
                        n1+1, loop+1, labelstring)
                    plotgui.add_set(xvalues, yvalues, xlowerror, xhigherror,
                                    ylowerror, yhigherror, xerrorflag,
                                    yerrorflag, option2, labelstring,
                                    plotgui.current_plot)
        plotgui.filename = None
        make_plot.make_plot(plotgui)
        return
#    except Exception:
#        plotgui.datavalues = None
#        plotgui.filename = None
#        tkinter.messagebox.showinfo(
#            'Error',
#            'There was some error in finding the data values.')
#        return
