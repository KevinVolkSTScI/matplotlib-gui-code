import math
import tkinter as Tk
import tkinter.messagebox
from tkinter.colorchooser import askcolor
from tkinter.scrolledtext import ScrolledText
import matplotlib
import matplotlib.lines as mlines
from matplotlib.patches import Rectangle, Ellipse, FancyArrow
import general_utilities
import make_plot

def initialize_plot_objects(plotgui):
    """
    Initialize the plot objects (lines, labels, and so on).

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    # set the number of labels, lines, and so on
    plotgui.max_labels = 100
    plotgui.max_ellipses = 100
    plotgui.max_boxes = 100
    plotgui.max_lines = 100
    plotgui.max_vectors = 100
    # For labels, define a variable to hold the position and the labels.
    plotgui.plot_labels = []
    for loop in range(plotgui.max_labels):
        plotgui.plot_labels.append({'xposition': None, 'yposition': None,
                                 'labelstring': '', 'plot': 1,
                                 'colour': 'black', 'size': 12,
                                 'font': 'times new roman',
                                 # 'font': 'sans-serif',
                                 'fontweight': 'normal'})
    # Similar thing for other drawing objects
    plotgui.plot_lines = []
    for loop in range(plotgui.max_lines):
        plotgui.plot_lines.append({'xstart': None, 'ystart': None,
                                'xend': None, 'yend': None, 'plot': 1,
                                'line_type': 'solid', 'line_colour':
                                'black', 'line_thickness': 1.0})
    plotgui.plot_vectors = []
    for loop in range(plotgui.max_vectors):
        plotgui.plot_vectors.append({'xstart': None, 'ystart': None,
                                  'xend': None, 'yend': None,
                                  'delx': None, 'dely': None,
                                  'plot': 1, 'line_type': 'solid',
                                  'line_colour': 'black',
                                  'line_thickness': 1.0, 'fill': True,
                                  'fill_colour': 'black'})
    plotgui.plot_ellipses = []
    for loop in range(plotgui.max_ellipses):
        plotgui.plot_ellipses.append({'xposition': None, 'yposition': None,
                                   'major': None, 'minor': None,
                                   'rotation': 0.0, 'plot': 1,
                                   'line_type': 'solid', 'line_colour':
                                   'black', 'line_thickness': 1.0,
                                   'fill_colour': 'none'})
    plotgui.plot_boxes = []
    for loop in range(plotgui.max_boxes):
        plotgui.plot_boxes.append({'xstart': None, 'ystart': None,
                                'xend': None, 'yend': None,
                                'rotation': 0.0, 'plot': 1,
                                'line_type': 'solid', 'line_colour':
                                'black', 'line_thickness': 1.0,
                                'fill_colour': 'none'})
    # flags for setting the different objects, variables for the number
    # of objects of each type
    plotgui.label_flag = False
    plotgui.box_flag = False
    plotgui.ellipse_flag = False
    plotgui.line_flag = False
    plotgui.vector_flag = False
    plotgui.number_of_labels = 0
    plotgui.number_of_lines = 0
    plotgui.number_of_boxes = 0
    plotgui.number_of_ellipses = 0
    plotgui.number_of_vectors = 0

def clear_plot_objects(plotgui):
    """
    Clear any labels, lines, and so on for the current plot

    Parameters
    ----------

        plotgui:   by assumption a matplotlib_user_interface object

    Returns
    -------

       None

    """
    for loop in range(plotgui.number_of_lines):
        if plotgui.plot_lines[loop]['plot'] == plotgui.current_plot:
            plotgui.plot_lines[loop] = {
                'xstart': None, 'ystart': None, 'xend': None,
                'yend': None, 'plot': 1, 'line_type': 'solid',
                'line_colour': 'black', 'line_thickness': 1.0}
    for loop in range(plotgui.number_of_labels):
        if plotgui.plot_labels[loop]['plot'] == plotgui.current_plot:
            plotgui.plot_labels = {
                'xposition': None, 'yposition': None, 'labelstring': '',
                'plot': 1, 'colour': 'black', 'size': 12,
                'font': 'times new roman', 'fontweight': 'normal'}
                # 'font': 'sans-serif', 'fontweight': 'normal'}
    for loop in range(plotgui.number_of_vectors):
        if plotgui.plot_vectors[loop]['plot'] == plotgui.current_plot:
            plotgui.plot_vectors[loop] = {
                'xstart': None, 'ystart': None, 'xend': None,
                'yend': None, 'delx': None, 'dely': None, 'plot': 1,
                'line_type': 'solid', 'line_colour': 'black',
                'line_thickness': 1.0, 'fill': True,
                'fill_colour': 'black'}
    for loop in range(plotgui.number_of_boxes):
        if plotgui.plot_boxes[loop]['plot'] == plotgui.current_plot:
            plotgui.plot_boxes[loop] = {
                'xstart': None, 'ystart': None, 'xend': None,
                'yend': None, 'rotation': 0.0, 'plot': 1,
                'line_type': 'solid', 'line_colour': 'black',
                'line_thickness': 1.0, 'fill_colour': 'none'}
    for loop in range(plotgui.number_of_ellipses):
        if plotgui.plot_ellipses[loop]['plot'] == plotgui.current_plot:
            plotgui.plot_ellipses[loop] = {
                'xposition': None, 'yposition': None, 'major': None,
                'minor': None, 'rotation': 0.0, 'plot': 1,
                'line_type': 'solid', 'line_colour': 'black',
                'line_thickness': 1.0, 'fill_colour': 'none'}

def remove_line(plotgui):
    """
    Remove the last line from the list of plot elements.

    This routine sets the line counter down by 1 so the last defined
    line is no longer plotted by the code.  The next new line
    subsequently defined will then overwrite the last current one.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_lines > 0:
        plotgui.number_of_lines = plotgui.number_of_lines - 1
        make_plot.make_plot(plotgui)

def remove_vector(plotgui):
    """
    Remove the last vector from the list of plot elements.

    This routine sets the vector counter down by 1 so the last defined
    vector is no longer plotted by the code.  The next new vector
    subsequently defined will then overwrite the last current one.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_vectors > 0:
        plotgui.number_of_vectors = plotgui.number_of_vectors - 1
        make_plot.make_plot(plotgui)

def remove_ellipse(plotgui):
    """
    Remove the last ellipse from the list of plot elements.

    This routine sets the ellipse counter down by 1 so the last defined
    ellipse is no longer plotted by the code.  Any new ellipses
    subsquently defined will then overwrite the last current one.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None


    """
    if plotgui.number_of_ellipses > 0:
        plotgui.number_of_ellipses = plotgui.number_of_ellipses - 1
        make_plot.make_plot(plotgui)

def remove_box(plotgui):
    """
    Remove the last box from the list of plot elements.

    This routine sets the box counter down by 1 so the last box is no
    longer plotted by the code.  Any new box subsequently defined will
    then overwrite the last existing one.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_boxes > 0:
        plotgui.number_of_boxes = plotgui.number_of_boxes - 1
        make_plot.make_plot(plotgui)

def remove_all_lines(plotgui):
    """
    Remove all lines from the plot elements.

    This routine sets the line counter variable to zero so that any
    defined lines are not plotted, and subsequently defined lines will
    overwrite any that are currently in the plotgui.plot_lines variable.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_lines > 0:
        plotgui.number_of_lines = 0
        plotgui.positions = []
        make_plot.make_plot(plotgui)

def remove_all_vectors(plotgui):
    """
    Remove all vectors from the plot elements.

    This routine sets the vector counter variable to zero so that any
    defined vectors are not plotted, and subsequently defined vectors
    will overwrite any that are currently in the plotgui.plot_vectors
    variable.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_vectors > 0:
        plotgui.number_of_vectors = 0
        plotgui.positions = []
        make_plot.make_plot(plotgui)

def remove_all_ellipses(plotgui):
    """
    Remove all ellipses from the plot elements.

    This routine sets the ellipse counter variable to zero so that any
    defined ellipses are not plotted, and subsequently defined ellipses
    will overwrite any that are currently in the plotgui.plot_ellipses
    variable.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_ellipses > 0:
        plotgui.number_of_ellipses = 0
        plotgui.positions = []
        make_plot.make_plot(plotgui)

def remove_all_boxes(plotgui):
    """
    Remove all boxes from the plot elements.

    This routine sets the box counter variable to zero so that any
    defined boxes are not plotted, and subsequently defined boxes
    will overwrite any that are currently in the plotgui.plot_boxes variable.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_boxes > 0:
        plotgui.number_of_boxes = 0
        plotgui.positions = []
        make_plot.make_plot(plotgui)

def add_line(plotgui):
    """
    Set a flag for marking line points on the plot.

    This routine sets a flag so that mouse clicks are used for defining
    a line.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    plotgui.line_flag = True

def add_vector(plotgui):
    """
    Set a flag for marking vector points on the plot.

     This routine sets a flag so that mouse clicks are used for defining
     a vector.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    plotgui.vector_flag = True

def add_ellipse(plotgui):
    """
    Set a flag for marking ellipse points on the plot.

    This routine sets a flag so that mouse clicks are used for defining
    an ellipse.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    plotgui.ellipse_flag = True

def add_box(plotgui):
    """
    Set a flag for marking box points on the plot.

    This routine sets a flag so that mouse clicks are used for defining
    a box.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    plotgui.box_flag = True

def set_label(plotgui):
    """
    Set the flag so that a mouse click is used to position a label.

    This routine sets the label flag in response to the "Put Label"
    button.  When the flag is set, any key pressed leads to putting
    a label on the plot.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    plotgui.label_flag = True


def add_box_values(plotgui):
    """
    Create a box on the plot.

    This code is activated when the box definition option is selected.
    When a button press event and then a button release event are
    received then the positions are recorded in plotgui.positions.  This
    routine reads these positions and presents a window with the box
    parameters for the user to change as they wish.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    BGCOL = '#F8F8FF'
    try:
        plotgui.positions[-2][0]
        plotgui.positions[-2][1]
        plotgui.positions[-1][0]
        plotgui.positions[-1][1]
    except TypeError:
        tkinter.messagebox.showinfo(
            "Error",
            "The required start and stop positions are not available"
             + " to make a box.")
        return
    except ValueError:
        tkinter.messagebox.showinfo(
            "Error",
            "The required start and stop positions are not available"
             + " to make a box.")
        return
    plotgui.box_flag = False
    plotgui.box_window = Tk.Toplevel()
    plotgui.box_window.title('Set Box Properties')
    plotgui.box_window.config(bg=BGCOL)
    frame1 = Tk.Frame(plotgui.box_window)
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
    plotgui.boxfields = []
    plotgui.boxfields.append(Tk.Entry(frame1, width=20))
    plotgui.boxfields[-1].grid(column=1, row=0, sticky=Tk.W)
    plotgui.boxfields.append(Tk.Entry(frame1, width=20))
    plotgui.boxfields[-1].grid(column=1, row=1, sticky=Tk.W)
    plotgui.boxfields.append(Tk.Entry(frame1, width=20))
    plotgui.boxfields[-1].grid(column=1, row=2, sticky=Tk.W)
    plotgui.boxfields.append(Tk.Entry(frame1, width=20))
    plotgui.boxfields[-1].grid(column=1, row=3, sticky=Tk.W)
    plotgui.boxfields.append(Tk.Entry(frame1, width=20))
    plotgui.boxfields[-1].grid(column=1, row=4, sticky=Tk.W)
    plotgui.boxfields[0].insert(0, str(plotgui.positions[-2][0]))
    plotgui.boxfields[1].insert(0, str(plotgui.positions[-2][1]))
    plotgui.boxfields[2].insert(0, str(plotgui.positions[-1][0]))
    plotgui.boxfields[3].insert(0, str(plotgui.positions[-1][1]))
    plotgui.boxfields[4].insert(0, '0.0')
    plotgui.boxfields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.boxfields[-1].grid(column=1, row=5, sticky=Tk.W)
    plotgui.boxfields[-1]['values'] = matplotlib_line_name_list
    plotgui.boxfields[-1].current(0)
    plotgui.boxfields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.boxfields[-1].grid(column=1, row=6, sticky=Tk.W)
    plotgui.boxfields[-1]['values'] = plotgui.colourset
    plotgui.boxfields[-1].current(0)
    plotgui.boxfields.append(Tk.Entry(frame1, width=15))
    plotgui.boxfields[-1].grid(column=1, row=7, sticky=Tk.W)
    plotgui.boxfields[-1].insert(0, '1.0')
    plotgui.boxfields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.boxfields[-1].grid(column=1, row=8, sticky=Tk.W)
    plotgui.boxfields[-1]['values'] = plotgui.altcolourset
    plotgui.boxfields[-1].current(0)
    frame2 = Tk.Frame(plotgui.box_window)
    frame2.pack(side=Tk.TOP)
    apply_button = Tk.Button(frame2, text="Apply",
                             command=lambda: apply_box_values(plotgui))
    apply_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(frame2, text="    ")
    label1.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        frame2, text="Close Window",
        command=lambda: plotgui.close_window(plotgui.box_window))
    close_button.pack(side=Tk.LEFT)

def apply_line_values(plotgui):
    """
    Create a line on the plot.

    This code reads the values in the line properties defintion window
    and applies them to the next available line.  The plot is then
    redone and the line properties window is removed.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    try:
        x1 = float(plotgui.linefields[0].get())
        y1 = float(plotgui.linefields[1].get())
        x2 = float(plotgui.linefields[2].get())
        y2 = float(plotgui.linefields[3].get())
        t1 = float(plotgui.linefields[6].get())
        plotgui.plot_lines[plotgui.number_of_lines]['xstart'] = x1
        plotgui.plot_lines[plotgui.number_of_lines]['ystart'] = y1
        plotgui.plot_lines[plotgui.number_of_lines]['xend'] = x2
        plotgui.plot_lines[plotgui.number_of_lines]['yend'] = y2
        line_index = plotgui.linefields[4].current()
        colour_index = plotgui.linefields[5].current()
        plotgui.plot_lines[plotgui.number_of_lines]['line_thickness'] = t1
        if plotgui.colourset[colour_index] == 'select':
            values = askcolor()
            plotgui.plot_lines[plotgui.number_of_lines]['line_colour'] = \
                values[1]
        else:
            plotgui.plot_lines[plotgui.number_of_lines]['line_colour'] = \
                plotgui.colourset[colour_index]
        plotgui.plot_lines[plotgui.number_of_lines]['line_type'] = \
            matplotlib_line_name_list[line_index]
        plotgui.plot_lines[plotgui.number_of_lines]['plot'] = plotgui.current_plot
        plotgui.number_of_lines = plotgui.number_of_lines + 1
        make_plot.make_plot(plotgui)
        plotgui.close_data_window(plotgui.line_window)
    except Exception:
        return

def apply_vector_values(plotgui):
    """
    Create a vector for the plot.

    This code reads the values in the vector properties defintion window
    and applies them to the next available vector.  The plot is then
    redone and the vector properties window is removed.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    try:
        x1 = float(plotgui.vectorfields[0].get())
        y1 = float(plotgui.vectorfields[1].get())
        x2 = float(plotgui.vectorfields[2].get())
        y2 = float(plotgui.vectorfields[3].get())
        delx = float(plotgui.vectorfields[4].get())
        dely = float(plotgui.vectorfields[5].get())
        t1 = float(plotgui.vectorfields[8].get())
        plotgui.plot_vectors[plotgui.number_of_vectors]['xstart'] = x1
        plotgui.plot_vectors[plotgui.number_of_vectors]['ystart'] = y1
        plotgui.plot_vectors[plotgui.number_of_vectors]['xend'] = x2
        plotgui.plot_vectors[plotgui.number_of_vectors]['yend'] = y2
        plotgui.plot_vectors[plotgui.number_of_vectors]['delx'] = delx
        plotgui.plot_vectors[plotgui.number_of_vectors]['dely'] = dely
        plotgui.plot_vectors[plotgui.number_of_vectors]['line_thickness'] = t1
        line_index = plotgui.vectorfields[6].current()
        colour_index = plotgui.vectorfields[7].current()
        if plotgui.colourset[colour_index] == 'select':
            values = askcolor()
            plotgui.plot_vectors[plotgui.number_of_vectors]['line_colour'] = \
                values[1]
        else:
            plotgui.plot_vectors[plotgui.number_of_vectors]['line_colour'] = \
                plotgui.colourset[colour_index]
        plotgui.plot_vectors[plotgui.number_of_vectors]['line_type'] = \
            matplotlib_line_name_list[line_index]
        plotgui.plot_vectors[plotgui.number_of_vectors]['plot'] = \
            plotgui.current_plot
        flag = plotgui.vectorfields[9].get()
        if flag == 1:
            plotgui.plot_vectors[plotgui.number_of_vectors]['fill'] = True
        else:
            plotgui.plot_vectors[plotgui.number_of_vectors]['fill'] = False
        colour_index = plotgui.vectorfields[10].current()
        if plotgui.colourset[colour_index] == 'select':
            values = askcolor()
            plotgui.plot_vectors[plotgui.number_of_vectors]['fill_colour'] = \
                values[1]
        else:
            plotgui.plot_vectors[plotgui.number_of_vectors]['fill_colour'] = \
                plotgui.colourset[colour_index]
        plotgui.number_of_vectors = plotgui.number_of_vectors + 1
        make_plot.make_plot(plotgui)
        plotgui.close_data_window(plotgui.vector_window)
    except Exception:
        return

def apply_box_values(plotgui):
    """
    Create a box for the plot.

    This code reads the values in the box properties defintion window and
    applies them to the next available box.  The plot is then redone and
    the box properties window is removed.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    try:
        x1 = float(plotgui.boxfields[0].get())
        y1 = float(plotgui.boxfields[1].get())
        x2 = float(plotgui.boxfields[2].get())
        y2 = float(plotgui.boxfields[3].get())
        angle = float(plotgui.boxfields[4].get())
        angle = angle % 360.
        t1 = float(plotgui.boxfields[7].get())
        plotgui.plot_boxes[plotgui.number_of_boxes]['xstart'] = x1
        plotgui.plot_boxes[plotgui.number_of_boxes]['ystart'] = y1
        plotgui.plot_boxes[plotgui.number_of_boxes]['xend'] = x2
        plotgui.plot_boxes[plotgui.number_of_boxes]['yend'] = y2
        plotgui.plot_boxes[plotgui.number_of_boxes]['rotation'] = angle
        line_index = plotgui.boxfields[5].current()
        colour_index1 = plotgui.boxfields[6].current()
        colour_index2 = plotgui.boxfields[8].current()
        plotgui.plot_boxes[plotgui.number_of_boxes]['line_thickness'] = t1
        if plotgui.colourset[colour_index1] == 'select':
            values = askcolor()
            plotgui.plot_boxes[plotgui.number_of_boxes]['line_colour'] = \
                values[1]
        else:
            plotgui.plot_boxes[plotgui.number_of_boxes]['line_colour'] = \
                plotgui.colourset[colour_index1]
        plotgui.plot_boxes[plotgui.number_of_boxes]['line_type'] = \
            matplotlib_line_name_list[line_index]
        plotgui.plot_boxes[plotgui.number_of_boxes]['plot'] = plotgui.current_plot
        plotgui.plot_boxes[plotgui.number_of_boxes]['fill_colour'] = \
            plotgui.altcolourset[colour_index2]
        plotgui.number_of_boxes = plotgui.number_of_boxes + 1
        make_plot.make_plot(plotgui)
        plotgui.close_window(plotgui.box_window)
    except Exception:
        return

def add_ellipse_values(plotgui):
    """
    Create an ellipse for the plot.

    This code is activated when the ellipse definition option is selected.
    When a button press event and then a button release event are received
    then the positions are recorded in plotgui.positions.  This routine reads
    these positions and presents a window with the ellipse parameters for
    the user to change as they wish.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    BGCOL = '#F8F8FF'
    try:
        plotgui.positions[-2][0]
        plotgui.positions[-2][1]
        plotgui.positions[-1][0]
        plotgui.positions[-1][1]
    except ValueError:
        tkinter.messagebox.showinfo(
            "Error", "The required start and"
             + " stop positions are not available to make an ellipse.")
        return
    plotgui.ellipse_flag = False
    plotgui.ellipse_window = Tk.Toplevel()
    plotgui.ellipse_window.title('Set Ellipse Properties')
    plotgui.ellipse_window.config(bg=BGCOL)
    frame1 = Tk.Frame(plotgui.ellipse_window)
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
    plotgui.ellipsefields = []
    plotgui.ellipsefields.append(Tk.Entry(frame1, width=20))
    plotgui.ellipsefields[-1].grid(column=1, row=0, sticky=Tk.W)
    plotgui.ellipsefields.append(Tk.Entry(frame1, width=20))
    plotgui.ellipsefields[-1].grid(column=1, row=1, sticky=Tk.W)
    plotgui.ellipsefields.append(Tk.Entry(frame1, width=20))
    plotgui.ellipsefields[-1].grid(column=1, row=2, sticky=Tk.W)
    plotgui.ellipsefields.append(Tk.Entry(frame1, width=20))
    plotgui.ellipsefields[-1].grid(column=1, row=3, sticky=Tk.W)
    plotgui.ellipsefields.append(Tk.Entry(frame1, width=20))
    plotgui.ellipsefields[-1].grid(column=1, row=4, sticky=Tk.W)
    xcenter = (plotgui.positions[-2][0]+plotgui.positions[-1][0])/2.
    ycenter = (plotgui.positions[-2][1]+plotgui.positions[-1][1])/2.
    xwidth = abs(plotgui.positions[-2][0]-plotgui.positions[-1][0])
    ywidth = abs(plotgui.positions[-2][1]-plotgui.positions[-1][1])
    plotgui.ellipsefields[0].insert(0, str(xcenter))
    plotgui.ellipsefields[1].insert(0, str(ycenter))
    plotgui.ellipsefields[2].insert(0, str(xwidth))
    plotgui.ellipsefields[3].insert(0, str(ywidth))
    plotgui.ellipsefields[4].insert(0, '0.0')
    plotgui.ellipsefields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.ellipsefields[-1].grid(column=1, row=5, sticky=Tk.W)
    plotgui.ellipsefields[-1]['values'] = matplotlib_line_name_list
    plotgui.ellipsefields[-1].current(0)
    plotgui.ellipsefields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.ellipsefields[-1].grid(column=1, row=6, sticky=Tk.W)
    plotgui.ellipsefields[-1]['values'] = plotgui.colourset
    plotgui.ellipsefields[-1].current(0)
    plotgui.ellipsefields.append(Tk.Entry(frame1, width=15))
    plotgui.ellipsefields[-1].grid(column=1, row=7, sticky=Tk.W)
    plotgui.ellipsefields[-1].insert(0, '1.0')
    plotgui.ellipsefields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.ellipsefields[-1].grid(column=1, row=8, sticky=Tk.W)
    plotgui.ellipsefields[-1]['values'] = plotgui.altcolourset
    plotgui.ellipsefields[-1].current(0)
    frame2 = Tk.Frame(plotgui.ellipse_window)
    frame2.pack(side=Tk.TOP)
    apply_button = Tk.Button(frame2, text="Apply",
                             command=lambda: apply_ellipse_values(plotgui))
    apply_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(frame2, text="    ")
    label1.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        frame2, text="Close Window",
        command=lambda: plotgui.close_window(plotgui.ellipse_window))
    close_button.pack(side=Tk.LEFT)

def apply_ellipse_values(plotgui):
    """
    Create an ellipse for the plot.

    This code reads the values in the ellipse properties defintion window
    and applies them to the next available ellipse.  The plot is then
    redone and the ellipse properties window is removed.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    try:
        x1 = float(plotgui.ellipsefields[0].get())
        y1 = float(plotgui.ellipsefields[1].get())
        x2 = float(plotgui.ellipsefields[2].get())
        y2 = float(plotgui.ellipsefields[3].get())
        angle = float(plotgui.ellipsefields[4].get())
        angle = angle % 360.
        t1 = float(plotgui.ellipsefields[7].get())
        plotgui.plot_ellipses[plotgui.number_of_ellipses]['xposition'] = x1
        plotgui.plot_ellipses[plotgui.number_of_ellipses]['yposition'] = y1
        plotgui.plot_ellipses[plotgui.number_of_ellipses]['major'] = x2
        plotgui.plot_ellipses[plotgui.number_of_ellipses]['minor'] = y2
        plotgui.plot_ellipses[plotgui.number_of_ellipses]['rotation'] = angle
        line_index = plotgui.ellipsefields[5].current()
        colour_index1 = plotgui.ellipsefields[6].current()
        colour_index2 = plotgui.ellipsefields[8].current()
        plotgui.plot_ellipses[plotgui.number_of_ellipses]['line_thickness'] = t1
        if plotgui.colourset[colour_index1] == 'select':
            values = askcolor()
            plotgui.plot_ellipses[plotgui.number_of_ellipses]['line_colour'] = \
                values[1]
        else:
            plotgui.plot_ellipses[plotgui.number_of_ellipses]['line_colour'] = \
                plotgui.colourset[colour_index1]
        plotgui.plot_ellipses[plotgui.number_of_ellipses]['line_type'] = \
            matplotlib_line_name_list[line_index]
        plotgui.plot_ellipses[plotgui.number_of_ellipses]['plot'] = \
            plotgui.current_plot
        plotgui.plot_ellipses[plotgui.number_of_ellipses]['fill_colour'] = \
            plotgui.altcolourset[colour_index2]
        plotgui.number_of_ellipses = plotgui.number_of_ellipses + 1
        make_plot.make_plot(plotgui)
        plotgui.close_window(plotgui.ellipse_window)
    except Exception:
        return

def add_line_values(plotgui):
    """
    Create a line for the plot.

    This code is activated when the line definition option is selected.
    When a button press event and then a button release event are
    received then the positions are recorded in plotgui.positions.  This
    routine reads these positions and presents a window with the line
    parameters for the user to change as they wish.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    BGCOL = '#F8F8FF'
    try:
        plotgui.positions[-2][0]
        plotgui.positions[-2][1]
        plotgui.positions[-1][0]
        plotgui.positions[-1][1]
    except ValueError:
        tkinter.messagebox.showinfo(
            "Error", "The required start and "
            + "stop positions are not available to make a line.")
        return
    plotgui.line_flag = False
    plotgui.line_window = Tk.Toplevel()
    plotgui.line_window.title('Set Line Properties')
    plotgui.line_window.config(bg=BGCOL)
    frame1 = Tk.Frame(plotgui.line_window)
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
    plotgui.linefields = []
    plotgui.linefields.append(Tk.Entry(frame1, width=20))
    plotgui.linefields[-1].grid(column=1, row=0, sticky=Tk.W)
    plotgui.linefields.append(Tk.Entry(frame1, width=20))
    plotgui.linefields[-1].grid(column=1, row=1, sticky=Tk.W)
    plotgui.linefields.append(Tk.Entry(frame1, width=20))
    plotgui.linefields[-1].grid(column=1, row=2, sticky=Tk.W)
    plotgui.linefields.append(Tk.Entry(frame1, width=20))
    plotgui.linefields[-1].grid(column=1, row=3, sticky=Tk.W)
    plotgui.linefields[0].insert(0, str(plotgui.positions[-2][0]))
    plotgui.linefields[1].insert(0, str(plotgui.positions[-2][1]))
    plotgui.linefields[2].insert(0, str(plotgui.positions[-1][0]))
    plotgui.linefields[3].insert(0, str(plotgui.positions[-1][1]))
    plotgui.linefields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.linefields[-1].grid(column=1, row=4, sticky=Tk.W)
    plotgui.linefields[-1]['values'] = matplotlib_line_name_list[0:-1]
    plotgui.linefields[-1].current(0)
    plotgui.linefields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.linefields[-1].grid(column=1, row=5, sticky=Tk.W)
    plotgui.linefields[-1]['values'] = plotgui.colourset
    plotgui.linefields[-1].current(0)
    plotgui.linefields.append(Tk.Entry(frame1, width=15))
    plotgui.linefields[-1].grid(column=1, row=6, sticky=Tk.W)
    plotgui.linefields[-1].insert(0, '1.0')
    frame2 = Tk.Frame(plotgui.line_window)
    frame2.pack(side=Tk.TOP)
    apply_button = Tk.Button(frame2, text="Apply",
                             command=lambda: apply_line_values(plotgui))
    apply_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(frame2, text="    ")
    label1.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        frame2, text="Close",
        command=lambda: plotgui.close_window(plotgui.line_window))
    close_button.pack(side=Tk.LEFT)

def add_vector_values(plotgui):
    """
    Create a vector for the plot.

    This code is activated when the vector definition option is selected.
    When a button press event and then a button release event are received
    then the positions are recorded in plotgui.positions.  This routine reads
    these positions and presents a window with the vector parameters for
    the user to change as they wish.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    matplotlib_line_name_list = ['solid', 'dashed', 'dashdot',
                                 'dotted', 'None']
    BGCOL = '#F8F8FF'
    try:
        plotgui.positions[-2][0]
        plotgui.positions[-2][1]
        plotgui.positions[-1][0]
        plotgui.positions[-1][1]
    except ValueError:
        tkinter.messagebox.showinfo(
            "Error", "The required start "
            + "and stop positions are not available to make a vector.")
        return
    plotgui.vector_flag = False
    plotgui.vector_window = Tk.Toplevel()
    plotgui.vector_window.title('Set Vector Properties')
    plotgui.vector_window.config(bg=BGCOL)
    frame1 = Tk.Frame(plotgui.vector_window)
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
    plotgui.vectorfields = []
    plotgui.vectorfields.append(Tk.Entry(frame1, width=20))
    plotgui.vectorfields[-1].grid(column=1, row=0, sticky=Tk.W)
    plotgui.vectorfields.append(Tk.Entry(frame1, width=20))
    plotgui.vectorfields[-1].grid(column=1, row=1, sticky=Tk.W)
    plotgui.vectorfields.append(Tk.Entry(frame1, width=20))
    plotgui.vectorfields[-1].grid(column=1, row=2, sticky=Tk.W)
    plotgui.vectorfields.append(Tk.Entry(frame1, width=20))
    plotgui.vectorfields[-1].grid(column=1, row=3, sticky=Tk.W)
    plotgui.vectorfields.append(Tk.Entry(frame1, width=20))
    plotgui.vectorfields[-1].grid(column=1, row=4, sticky=Tk.W)
    plotgui.vectorfields.append(Tk.Entry(frame1, width=20))
    plotgui.vectorfields[-1].grid(column=1, row=5, sticky=Tk.W)
    plotgui.vectorfields[0].insert(0, str(plotgui.positions[-2][0]))
    plotgui.vectorfields[1].insert(0, str(plotgui.positions[-2][1]))
    plotgui.vectorfields[2].insert(0, str(plotgui.positions[-1][0]))
    plotgui.vectorfields[3].insert(0, str(plotgui.positions[-1][1]))
    plotgui.vectorfields[4].insert(0, str(0.1))
    plotgui.vectorfields[5].insert(0, str(0.1))
    plotgui.vectorfields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.vectorfields[-1].grid(column=1, row=6, sticky=Tk.W)
    plotgui.vectorfields[-1]['values'] = matplotlib_line_name_list[0:-1]
    plotgui.vectorfields[-1].current(0)
    plotgui.vectorfields.append(tkinter.ttk.Combobox(frame1, width=15))
    plotgui.vectorfields[-1].grid(column=1, row=7, sticky=Tk.W)
    plotgui.vectorfields[-1]['values'] = plotgui.colourset
    plotgui.vectorfields[-1].current(0)
    plotgui.vectorfields.append(Tk.Entry(frame1, width=15))
    plotgui.vectorfields[-1].grid(column=1, row=8, sticky=Tk.W)
    plotgui.vectorfields[-1].insert(0, '1.0')
    plotgui.vector_head_fill_flag = Tk.IntVar()
    plotgui.vectorfields.append(plotgui.vector_head_fill_flag)
    plotgui.vectorfields.append(tkinter.ttk.Combobox(frame1, width=15))
    b1 = Tk.Frame(frame1, )
    b1.grid(column=1, row=9, sticky=Tk.W)
    general_utilities.put_yes_no(
        b1, plotgui.vector_head_fill_flag, ['yes', 'no'], True)
    plotgui.vectorfields[-1].grid(column=1, row=10, sticky=Tk.W)
    plotgui.vectorfields[-1]['values'] = plotgui.colourset
    plotgui.vectorfields[-1].current(0)
    frame2 = Tk.Frame(plotgui.vector_window)
    frame2.pack(side=Tk.TOP)
    apply_button = Tk.Button(frame2, text="Apply",
                             command=lambda: apply_vector_values(plotgui))
    apply_button.pack(side=Tk.LEFT)
    label1 = Tk.Label(frame2, text="    ")
    label1.pack(side=Tk.LEFT)
    close_button = Tk.Button(
        frame2, text="Close",
        command=lambda: plotgui.close_window(plotgui.vector_window))
    close_button.pack(side=Tk.LEFT)

