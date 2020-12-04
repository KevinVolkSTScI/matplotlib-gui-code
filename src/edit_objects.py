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

def edit_lines(plotgui):
    """
    Make window for editing the line values.

    This routine produces a text box in a window, within which one can
    edit the line values.  If no lines are defined the routine just
    exits with no action.

    Line values are presented one per text line, with the start x
    position, the start y position, the end x position, the end y
    position, the plot number, the line type, the line colour, and the
    line thickness separated by tab symbols.  One can edit the values
    within the text window and then these are applied when one clicks
    on the "Close Window" button.

    Parameters
    ----------

        plotgui : the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_lines == 0:
        return
    str1 = 'Edit values below: fields are separated by tab characters.\n'\
           + '     start x       start y       end x         end y   ' \
           + 'plot        colour   line type  thickness\n------------' \
           + '-------------------------------------------------------' \
           + '-----------------------------\n'
    for loop in range(plotgui.number_of_lines):
        str1 = str1 + '%12.6g\t%12.6g\t%12.6g\t%12.6g\t%6d\t' % (
            plotgui.plot_lines[loop]['xstart'],
            plotgui.plot_lines[loop]['ystart'],
            plotgui.plot_lines[loop]['xend'],
            plotgui.plot_lines[loop]['yend'],
            plotgui.plot_lines[loop]['plot'])
        str1 = str1 + '%15s\t%8s\t%7.3f\n' % (
            plotgui.plot_lines[loop]['line_colour'],
            plotgui.plot_lines[loop]['line_type'],
            plotgui.plot_lines[loop]['line_thickness'])
    line_window = Tk.Toplevel()
    line_window.title('Lines:')
    holder = Tk.Frame(line_window)
    holder.pack(side=Tk.TOP)
    line_message_text = ScrolledText(holder, height=40, width=100,
                                     wrap=Tk.NONE)
    line_message_text.config(font=('courier', 16, 'bold'))
    line_message_text.pack(side=Tk.TOP)
    line_message_text.insert(0.0, str1)
    bholder = Tk.Frame(line_window)
    bholder.pack(side=Tk.TOP)
    close_button = Tk.Button(
        bholder, text='Close Window',
        command=lambda: read_lines(plotgui, line_message_text, line_window))
    close_button.pack()

def read_lines(plotgui, line_message_text, line_window):
    """
    Read and apply the lines text field.

    This routine reads the line text field and makes the new set of lines
    and line positions.  It then applies these and closes the line window.

    Parameters
    ----------

        plotgui: the matplotlib_user_interface object holding the plot

        line_message_text:  a tkinter text field variable

        line_window:  a tkinter Toplevel or Tk variable that holds the 
                      text field

    Returns
    -------
        None

    The code does, however, change the plotgui.plot_lines values as needed
    to match what is in the line text field.

    """
    linetext = line_message_text.get(0.0, Tk.END)
    lines = linetext.split('\n')
    newlines = []
    nlines = 0
    for line in lines:
        values = line.split('\t')
        if len(values) == 8:
            try:
                x1 = float(values[0])
                y1 = float(values[1])
                x2 = float(values[2])
                y2 = float(values[3])
                nplot = int(values[4])
                colour = values[5].strip(' ')
                linetype = values[6].strip(' ')
                thickness = float(values[7])
                newlines.append({'xstart': x1, 'ystart': y1, 'xend': x2,
                                 'yend': y2, 'plot': nplot,
                                 'line_type': linetype,
                                 'line_colour': colour,
                                 'line_thickness': thickness})
                nlines = nlines + 1
            except ValueError:
                pass
    line_window.destroy()
    if nlines > plotgui.max_lines:
        plotgui.max_lines = nlines
    else:
        for loop in range(nlines, plotgui.max_lines):
            newlines.append({'xstart': None, 'ystart': None,
                             'xend': None, 'yend': None,
                             'plot': 1, 'line_type': 'solid',
                             'line_colour': 'black',
                             'line_thickness': 1.0})
    plotgui.plot_lines = newlines
    plotgui.number_of_lines = nlines
    make_plot.make_plot(plotgui)

def edit_boxes(plotgui):
    """
    Make a window for editing the box values.

    This routine produces a text box in a window, within which one can
    edit the box values.  If no boxes are defined the routine just exits
    with no action.

    Box values are presented one per line, with the start x position, the
    start y position, the end x position, the end y position, the
    orientation, the plot number, the line type, the line colour, the
    line thickness, and the fill colour separated by tab symbols.
    One can edit the values within the text window and then these are
    applied when one clicks on the "Close Window" button.

    Parameters
    ----------

        plotgui: the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_boxes == 0:
        return
    str1 = 'Edit values below: fields are separated by tab characters.\n' \
           + '     start x       start y       end x      end y      ' \
           + 'orient.   plot        colour   type  thickness  fill \n' \
           + '-------------------------------------------------------' \
           + '--------------------------------------------------------\n'
    for loop in range(plotgui.number_of_boxes):
        str1 = str1 + '%12.6g\t%12.6g\t%12.6g\t%12.6g\t%8.3f\t%6d' % (
            plotgui.plot_boxes[loop]['xstart'],
            plotgui.plot_boxes[loop]['ystart'],
            plotgui.plot_boxes[loop]['xend'],
            plotgui.plot_boxes[loop]['yend'],
            plotgui.plot_boxes[loop]['rotation'],
            plotgui.plot_boxes[loop]['plot'])
        str1 = str1 + '\t%15s\t%8s\t%7.3f\t%15s\n' % (
            plotgui.plot_boxes[loop]['line_colour'],
            plotgui.plot_boxes[loop]['line_type'],
            plotgui.plot_boxes[loop]['line_thickness'],
            plotgui.plot_boxes[loop]['fill_colour'])
    box_window = Tk.Toplevel()
    box_window.title('Boxes:')
    holder = Tk.Frame(box_window)
    holder.pack(side=Tk.TOP)
    box_message_text = ScrolledText(holder, height=40, width=115,
                                    wrap=Tk.NONE)
    box_message_text.config(font=('courier', 16, 'bold'))
    box_message_text.pack(side=Tk.TOP)
    box_message_text.insert(0.0, str1)
    bholder = Tk.Frame(box_window)
    bholder.pack(side=Tk.TOP)
    close_button = Tk.Button(
        bholder, text='Close Window',
        command=lambda: read_boxes(plotgui, box_message_text, box_window))
    close_button.pack()

def read_boxes(plotgui, box_message_text, box_window):
    """
    Read and apply the box text field.

    This routine reads the box text field and makes the new set of boxes
    and box positions.  It then applies these and closes the box window.

    Parameters
    ----------

        plotgui: the matplotlib_user_interface object holding the plot

        box_message_text:  a tkinter text field variable

        box_window:  a tkinter Toplevel or Tk variable that holds the 
                      text field

    Returns
    -------

        None

    The code does, however, change the plotgui.plot_boxes values as needed
    to match what is in the box text field.

    """
    boxtext = box_message_text.get(0.0, Tk.END)
    lines = boxtext.split('\n')
    newboxes = []
    nboxes = 0
    for line in lines:
        values = line.split('\t')
        if len(values) == 10:
            try:
                x1 = float(values[0])
                y1 = float(values[1])
                x2 = float(values[2])
                y2 = float(values[3])
                theta = float(values[4])
                nplot = int(values[5])
                colour = values[6].strip(' ')
                linetype = values[7].strip(' ')
                thickness = float(values[8])
                fillcolour = values[9].strip(' ')
                newboxes.append({'xstart': x1, 'ystart': y1,
                                 'xend': x2, 'yend': y2,
                                 'rotation': theta, 'plot': nplot,
                                 'line_type': linetype,
                                 'line_colour': colour,
                                 'line_thickness': thickness,
                                 'fill_colour': fillcolour})
                nboxes = nboxes + 1
            except ValueError:
                pass
    box_window.destroy()
    if nboxes > plotgui.max_boxes:
        plotgui.max_boxes = nboxes
    else:
        for loop in range(nboxes, plotgui.max_boxes):
            newboxes.append({'xstart': None, 'ystart': None,
                             'xend': None, 'yend': None,
                             'rotation': 0.0, 'plot': 1,
                             'line_type': 'solid',
                             'line_colour': 'black',
                             'line_thickness': 1.0,
                             'fill_colour': 'none'})
    plotgui.plot_boxes = newboxes
    plotgui.number_of_boxes = nboxes
    make_plot.make_plot(plotgui)

def edit_ellipses(plotgui):
    """
    Make a window for editing the ellipse values.

    This routine produces a text box in a window, within which one can
    edit the ellipse values.  If no ellipses are defined the routine
    just exits with no action.

    Ellipse values are presented one per line, with the center x position,
    the center y position, the major axis length (x), the minor axis
    length (y), the orientation, the plot number, the line type, the
    line colour, the line thickness, and the fill colour separated by
    tab symbols.  One can edit the values within the text window and
    then these are applied when one clicks on the "Close Window" button.

    Parameters
    ----------

        plotgui : the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_ellipses == 0:
        return
    str1 = 'Edit values below: fields are separated by tab characters.\n' \
           + '     center x       center y       major      minor      ' \
           + 'orient.   plot        colour   type  thickness  fill \n--' \
           + '---------------------------------------------------------' \
           + '----------------------------------------------------\n'
    for loop in range(plotgui.number_of_ellipses):
        str1 = str1 + '%12.6g\t%12.6g\t%12.6g\t%12.6g\t%8.3f\t%6d' % (
            plotgui.plot_ellipses[loop]['xposition'],
            plotgui.plot_ellipses[loop]['yposition'],
            plotgui.plot_ellipses[loop]['major'],
            plotgui.plot_ellipses[loop]['minor'],
            plotgui.plot_ellipses[loop]['rotation'],
            plotgui.plot_ellipses[loop]['plot'])
        str1 = str1 + '\t%15s\t%8s\t%7.3f\t%15s\n' % (
            plotgui.plot_ellipses[loop]['line_colour'],
            plotgui.plot_ellipses[loop]['line_type'],
            plotgui.plot_ellipses[loop]['line_thickness'],
            plotgui.plot_ellipses[loop]['fill_colour'])
    ellipse_window = Tk.Toplevel()
    ellipse_window.title('Ellipses:')
    holder = Tk.Frame(ellipse_window)
    holder.pack(side=Tk.TOP)
    ellipse_message_text = ScrolledText(holder, height=40, width=115,
                                        wrap=Tk.NONE)
    ellipse_message_text.config(font=('courier', 16, 'bold'))
    ellipse_message_text.pack(side=Tk.TOP)
    ellipse_message_text.insert(0.0, str1)
    bholder = Tk.Frame(ellipse_window)
    bholder.pack(side=Tk.TOP)
    close_button = Tk.Button(
        bholder, text='Close Window',
        command=lambda: read_ellipses(plotgui, ellipse_message_text,
                                      ellipse_window))
    close_button.pack()

def read_ellipses(plotgui, ellipse_message_text, ellipse_window):
    """
    Read and apply the ellipse text field.

    This routine reads the ellipse text field and makes the new set of
    ellipses and ellipse positions.  It then applies these and closes
    the ellipse window.

    Parameters
    ----------

        plotgui: the matplotlib_user_interface object holding the plot

        ellipse_message_text:  a tkinter text field variable

        ellipse_window:  a tkinter Toplevel or Tk variable that holds the 
                         text field

    Returns
    -------

        None

    The code does, however, change the plotgui.plot_ellipses values as
    needed to match what is in the ellipse text field.

    """
    ellipsetext = ellipse_message_text.get(0.0, Tk.END)
    lines = ellipsetext.split('\n')
    newellipses = []
    nellipses = 0
    for line in lines:
        values = line.split('\t')
        if len(values) == 10:
            try:
                x1 = float(values[0])
                y1 = float(values[1])
                x2 = float(values[2])
                y2 = float(values[3])
                theta = float(values[4])
                nplot = int(values[5])
                colour = values[6].strip(' ')
                linetype = values[7].strip(' ')
                thickness = float(values[8])
                fillcolour = values[9].strip(' ')
                newellipses.append({'xposition': x1, 'yposition': y1,
                                    'major': x2, 'minor': y2,
                                    'rotation': theta, 'plot': nplot,
                                    'line_type': linetype,
                                    'line_colour': colour,
                                    'line_thickness': thickness,
                                    'fill_colour': fillcolour})
                nellipses = nellipses + 1
            except ValueError:
                pass
    ellipse_window.destroy()
    if nellipses > plotgui.max_ellipses:
        plotgui.max_ellipses = nellipses
    else:
        for loop in range(nellipses, plotgui.max_ellipses):
            newellipses.append({'xposition': None, 'yposition': None,
                                'major': None, 'minor': None,
                                'rotation': 0.0, 'plot': 1,
                                'line_type': 'solid',
                                'line_colour': 'black',
                                'line_thickness': 1.0,
                                'fill_colour': 'none'})
    plotgui.plot_ellipses = newellipses
    plotgui.number_of_ellipses = nellipses
    make_plot.make_plot(plotgui)

def edit_vectors(plotgui):
    """
    Make a window for editing the vector values.

    This routine produces a text box in a window, within which one
    can edit the vector values.  If no vectors are defined the routine
    just exits with no action.

    Vectors are presented one per line, with the start x position, the
    start y position, the end x position, the end y position, the head
    width, the head length, the plot number, the line type, the line
    colour, and the line thickness separated by tab symbols.  One can
    edit the values within the text window and then these are applied
    when one clicks on the "Close Window" button.

    Parameters
    ----------

        plotgui : the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    if plotgui.number_of_vectors == 0:
        return
    str1 = 'Edit values below: fields are separated by tab characters.\n' \
           + '     start x       start y       end x         end y   ' \
           + 'head width head length   plot     line colour   line type' \
           + '  thickness head fill head colour\n----------------------' \
           + '---------------------------------------------------------' \
           + '---------------------------------------------------------' \
           + '------\n'
    for loop in range(plotgui.number_of_vectors):
        if plotgui.plot_vectors[loop]['fill']:
            flag1 = 'True'
        else:
            flag1 = 'False'
        str1 = str1 + '%12.6g\t%12.6g\t%12.6g\t%12.6g\t%12.6g' % (
            plotgui.plot_vectors[loop]['xstart'],
            plotgui.plot_vectors[loop]['ystart'],
            plotgui.plot_vectors[loop]['xend'],
            plotgui.plot_vectors[loop]['yend'],
            plotgui.plot_vectors[loop]['delx'])
        str1 = str1 + '\t%12.6g\t%6d\t%15s\t%8s\t%10.3f' % (
            plotgui.plot_vectors[loop]['dely'],
            plotgui.plot_vectors[loop]['plot'],
            plotgui.plot_vectors[loop]['line_colour'],
            plotgui.plot_vectors[loop]['line_type'],
            plotgui.plot_vectors[loop]['line_thickness'])
        str1 = str1 + '\t     %s\t    %s\n' % (
            flag1,
            plotgui.plot_vectors[loop]['fill_colour'])
    vector_window = Tk.Toplevel()
    vector_window.title('Vectors:')
    holder = Tk.Frame(vector_window)
    holder.pack(side=Tk.TOP)
    vector_message_text = ScrolledText(holder, height=40, width=145,
                                       wrap=Tk.NONE)
    vector_message_text.config(font=('courier', 16, 'bold'))
    vector_message_text.pack(side=Tk.TOP)
    vector_message_text.insert(0.0, str1)
    bholder = Tk.Frame(vector_window)
    bholder.pack(side=Tk.TOP)
    close_button = Tk.Button(
        bholder, text='Close Window',
        command=lambda: read_vectors(plotgui, vector_message_text,
                                     vector_window))
    close_button.pack()

def read_vectors(plotgui, vector_message_text, vector_window):
    """
    Read and apply the vectors text field.

    This routine reads the vector text field and makes the new set of
    vectors and vector positions.  It then applies these and closes
    the vector window.

    Parameters
    ----------

        plotgui: the matplotlib_user_interface object holding the plot

        vector_message_text:  a tkinter text field variable

        vector_window:  a tkinter Toplevel or Tk variable that holds the 
                        text field

    Returns
    -------

        None

    The code does, however, change the plotgui.plot_vectors values as
    needed to match what is in the vector text field.

    """
    vectortext = vector_message_text.get(0.0, Tk.END)
    vectors = vectortext.split('\n')
    newvectors = []
    nvectors = 0
    for vector in vectors:
        values = vector.split('\t')
        if len(values) == 12:
            try:
                x1 = float(values[0])
                y1 = float(values[1])
                x2 = float(values[2])
                y2 = float(values[3])
                delx = float(values[4])
                dely = float(values[5])
                nplot = int(values[6])
                colour = values[7].strip(' ')
                linetype = values[8].strip(' ')
                thickness = float(values[9])
                if 'true' in values[10].lower():
                    flag = True
                else:
                    flag = False
                hcolour = values[11].strip(' ')
                newvectors.append({'xstart': x1, 'ystart': y1,
                                   'xend': x2, 'yend': y2,
                                   'delx': delx, 'dely': dely,
                                   'plot': nplot, 'line_type': linetype,
                                   'line_colour': colour,
                                   'line_thickness': thickness,
                                   'fill': flag, 'fill_colour': hcolour})
                nvectors = nvectors + 1
            except ValueError:
                pass
    vector_window.destroy()
    if nvectors > plotgui.max_vectors:
        plotgui.max_vectors = nvectors
    else:
        for loop in range(nvectors, plotgui.max_vectors):
            newvectors.append({'xstart': None, 'ystart': None,
                               'xend': None, 'yend': None,
                               'plot': 1, 'line_type': 'solid',
                               'line_colour': 'black',
                               'line_thickness': 1.0, 'fill': True,
                               'fill_colour': 'black'})
    plotgui.plot_vectors = newvectors
    plotgui.number_of_vectors = nvectors
    make_plot.make_plot(plotgui)


