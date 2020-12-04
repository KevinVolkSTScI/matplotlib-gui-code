import numpy
import tkinter as Tk
import tkinter.messagebox
from tkinter.colorchooser import askcolor
from tkinter.scrolledtext import ScrolledText
import general_utilities
import make_plot

def edit_labels(plotgui):
    """
    Bring up a window to edit the labels in the plot.

    This routine produces a text box in a window, within which one can
    edit the label values.  If no labels are defined the routine just
    exits with no action.

    Labels are presented one per line, with the parameter values
    separated by tab symbols.  One can edit the values within the text
    window and then these are applied when one clicks on the "Close
    Window" button.

    Parameters
    ----------
        plotgui:  A matplotlib_user_interface object, for the plot holding 
                  the labels

    Returns
    -------
        None

    """
    if plotgui.number_of_labels == 0:
        return
    str1 = 'Edit values below: fields are separated by tab ' \
           + 'characters\n x position    y position  plot   ' \
           + 'label   colour size     font    fontweight' \
           + '\n----------------------------------------------' \
           + '-------------------------------------\n'
    for loop in range(plotgui.number_of_labels):
        if plotgui.xparameters[plotgui.plot_labels[loop]['plot']-1]['hybridlog']:
            xpos1 = general_utilities.inverse_hybrid_transform(
                plotgui.plot_labels[loop]['xposition'])
        else:
            xpos1 = plotgui.plot_labels[loop]['xposition']
        if plotgui.yparameters[plotgui.plot_labels[loop]['plot']-1]['hybridlog']:
            ypos1 = general_utilities.inverse_hybrid_transform(
                plotgui.plot_labels[loop]['yposition'])
        else:
            ypos1 = plotgui.plot_labels[loop]['yposition']
        str1 = str1 + '%12.6g\t%12.6g\t%3d\t%s\t%s\t%d\t%s\t%s\n' % (
            xpos1,
            ypos1,
            plotgui.plot_labels[loop]['plot'],
            plotgui.plot_labels[loop]['labelstring'],
            plotgui.plot_labels[loop]['colour'],
            plotgui.plot_labels[loop]['size'],
            plotgui.plot_labels[loop]['font'],
            plotgui.plot_labels[loop]['fontweight'])
    label_window = Tk.Toplevel()
    label_window.title('Labels:')
    holder = Tk.Frame(label_window)
    holder.pack(side=Tk.TOP)
    label_message_text = ScrolledText(holder, height=40, width=90,
                                      wrap=Tk.NONE)
    label_message_text.config(font=('courier', 16, 'bold'))
    label_message_text.pack(side=Tk.TOP)
    label_message_text.insert(0.0, str1)
    bholder = Tk.Frame(label_window)
    bholder.pack(side=Tk.TOP)
    close_button = Tk.Button(
        bholder, text='Close Window',
        command=lambda: read_labels(plotgui, label_message_text,
                                    label_window))
    close_button.pack()

def read_labels(plotgui, label_message_text, label_window):
    """
    Read and parse the label text field.

    This routine reads the label text field and makes the new set of
    labels and label positions.  It then applies these and closes the
    label window.

    Parameters
    ----------
        plotgui:   A matplotlib_user_interface object

        label_message_text:   A tkinter text field variable

        label_window:   A tkinter top level variable for the label entry

    Returns
    -------
        Nothing

    The code does, however, change the plotgui.plot_labels values as
    needed to match what is in the label text field.

    """
    labeltext = label_message_text.get(0.0, Tk.END)
    lines = labeltext.split('\n')
    newlabels = []
    nlabels = 0
    for line in lines:
        values = line.split('\t')
        if len(values) == 8:
            try:
                x1 = float(values[0])
                y1 = float(values[1])
                nplot = int(values[2])
                label = values[3]
                colour = values[4]
                size = int(values[5])
                font = values[6]
                fontweight = values[7]
                label = label.strip('\n')
                if plotgui.xparameters[nplot-1]['hybridlog']:
                    x1 = general_utilities.hybrid_transform(x1)
                if plotgui.yparameters[nplot-1]['hybridlog']:
                    y1 = general_utilities.hybrid_transform(y1)
                newlabels.append({'xposition': x1, 'yposition': y1,
                                  'labelstring': label, 'plot': nplot,
                                  'colour': colour, 'size': size,
                                  'font': font,
                                  'fontweight': fontweight})
                nlabels = nlabels + 1
            except ValueError:
                pass
    label_window.destroy()
    if nlabels > plotgui.max_labels:
        plotgui.max_labels = nlabels
    else:
        for loop in range(nlabels, plotgui.max_labels):
            newlabels.append({'xposition': None, 'yposition': None,
                              'labelstring': '', 'plot': 1,
                              'colour': 'black', 'size': 12,
                              # 'font': 'sans-serif',
                              'font': 'times new roman',
                              'fontweight': 'normal'})
        plotgui.plot_labels = newlabels
        plotgui.number_of_labels = nlabels
    make_plot.make_plot(plotgui)

def clear_labels(plotgui):
    """
    Clear all the labels on the plot.

    This routine is called when the "Clear Labels" button is pressed.
    It asks whether the labels should be cleared, and if so all labels
    are removed.

    Parameters
    ----------
        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------
        None

    """
    response = tkinter.messagebox.askyesno("Verify", "Delete all labels.")
    if response:
        plotgui.label_flag = False
        plotgui.plot_labels = []
        plotgui.max_labels = 100
        for loop in range(plotgui.max_labels):
            plotgui.plot_labels.append({'xposition': None,
                                     'yposition': None,
                                     'labelstring': '', 'plot': 1,
                                     'colour': 'black', 'size': 12,
                                     # 'font': 'sans-serif',
                                     'font': 'times new roman',
                                    'fontweight': 'normal'})
        plotgui.positions = []
        make_plot.make_plot(plotgui)

def set_label_properties(plotgui, ind1):
    """
    Make a window in which the label properties are assigned.

    Parameters
    ----------
        plotgui:   the matplotlib_user_interface object holding the plot

        ind1:  an integer value >= 0, the index of the label that is
                having the properties set

    Returns
    -------
        Nothing

    The routine gets input from the window each time the "set properties"
    button is activated, and applies these to label ind1.
    """
    label_property_window = Tk.Toplevel()
    label_property_window.title('Label properties')
    holder = Tk.Frame(label_property_window)
    holder.pack(side=Tk.TOP)
    fontnames = ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace',
                 'times new roman']
    fontsizes = ['8', '9', '10', '11', '12', '13', '14', '16', '18', '20',
                 '24', '30']
    fontweights = ['ultralight', 'light', 'normal', 'regular', 'book',
                   'medium', 'roman', 'semibold', 'demibold', 'demi',
                   'bold', 'heavy', 'extra bold', 'black']
    label = Tk.Label(holder, text='Font Name:')
    label.grid(column=0, row=0)
    label = Tk.Label(holder, text='Font Size:')
    label.grid(column=0, row=1)
    label = Tk.Label(holder, text='Font Weight:')
    label.grid(column=0, row=2)
    label = Tk.Label(holder, text='Font Colour:')
    label.grid(column=0, row=3)
    label = Tk.Label(holder, text='Label Text:')
    label.grid(column=0, row=4)
    label = Tk.Label(holder, text='x position:')
    label.grid(column=0, row=5)
    label = Tk.Label(holder, text='y position:')
    label.grid(column=0, row=6)
    plotgui.label_font_name_list = Tk.ttk.Combobox(holder, width=20)
    plotgui.label_font_name_list.grid(column=1, row=0)
    plotgui.label_font_name_list['values'] = fontnames
#    plotgui.label_font_name_list.set('sans-serif')
    plotgui.label_font_name_list.set('times new roman')
    for loop in range(len(fontnames)):
        if plotgui.plot_labels[ind1]['font'] == fontnames[loop]:
            plotgui.label_font_name_list.set(fontnames[loop])
    plotgui.label_font_size_list = Tk.ttk.Combobox(holder, width=20)
    plotgui.label_font_size_list.grid(column=1, row=1)
    plotgui.label_font_size_list['values'] = fontsizes
    plotgui.label_font_size_list.set('12')
    for loop in range(len(fontsizes)):
        if plotgui.plot_labels[ind1]['size'] == int(fontsizes[loop]):
            plotgui.label_font_size_list.set(fontsizes[loop])
    plotgui.label_font_weight_list = Tk.ttk.Combobox(holder, width=20)
    plotgui.label_font_weight_list.grid(column=1, row=2)
    plotgui.label_font_weight_list['values'] = fontweights
    for loop in range(len(fontweights)):
        if plotgui.plot_labels[ind1]['fontweight'] == fontweights[loop]:
            plotgui.label_font_weight_list.set(fontweights[loop])
    plotgui.label_font_colour_list = Tk.ttk.Combobox(holder, width=20)
    plotgui.label_font_colour_list.grid(column=1, row=3)
    plotgui.label_font_colour_list['values'] = plotgui.colourset
    plotgui.label_font_colour_list.set('black')
    for loop in range(len(plotgui.colourset)):
        if plotgui.plot_labels[ind1]['colour'] == plotgui.colourset[loop]:
            plotgui.label_font_colour_list.set(plotgui.colourset[loop])
    plotgui.label_text_entry = Tk.Entry(holder, width=20)
    plotgui.label_text_entry.grid(column=1, row=4)
    plotgui.label_text_entry.insert(0, plotgui.plot_labels[ind1]['labelstring'])
    plotgui.label_x_position = Tk.Entry(holder, width=20)
    plotgui.label_x_position.grid(column=1, row=5)
    if plotgui.xparameters[plotgui.plot_labels[ind1]['plot']-1]['hybridlog']:
        xpos1 = general_utilities.inverse_hybrid_transform(
            plotgui.plot_labels[ind1]['xposition'])
    else:
        xpos1 = plotgui.plot_labels[ind1]['xposition']
    plotgui.label_x_position.insert(0, str(xpos1))
    plotgui.label_y_position = Tk.Entry(holder, width=20)
    plotgui.label_y_position.grid(column=1, row=6)
    if plotgui.yparameters[plotgui.plot_labels[ind1]['plot']-1]['hybridlog']:
        xpos1 = general_utilities.inverse_hybrid_transform(
            plotgui.plot_labels[ind1]['yposition'])
    else:
        ypos1 = plotgui.plot_labels[ind1]['yposition']
    plotgui.label_y_position.insert(0, str(ypos1))
    bholder = Tk.Frame(label_property_window)
    bholder.pack(side=Tk.TOP)
    set_button = Tk.Button(bholder, text='Set Properties',
                           command=lambda: set_label_values(plotgui, ind1))
    set_button.pack()
    close_button = Tk.Button(bholder, text='Close Window',
                             command=label_property_window.destroy)
    close_button.pack()

def set_label_values(plotgui, ind1):
    """
    Read parameters from the label input fields and apply them.

    Parameters
    ----------
        plotgui :   the matplotlib_user_interface object holding the plot

        ind1 : an integer value >= 0, the index of the label that is
               having the properties set

    Returns
    -------
        Nothing

    """
    try:
        font = plotgui.label_font_name_list.get()
        fontweight = plotgui.label_font_weight_list.get()
        fontsize = int(plotgui.label_font_size_list.get())
        fontcolour = plotgui.label_font_colour_list.get()
        labelstring = plotgui.label_text_entry.get()
        xpos1 = float(plotgui.label_x_position.get())
        if plotgui.xparameters[plotgui.plot_labels[ind1]['plot']-1]['hybridlog']:
            xpos1 = general_utilities.hybrid_transform(xpos1)
        ypos1 = float(plotgui.label_y_position.get())
        if plotgui.yparameters[plotgui.plot_labels[ind1]['plot']-1]['hybridlog']:
            ypos1 = general_utilities.hybrid_transform(ypos1)
        plotgui.plot_labels[ind1]['xposition'] = xpos1
        plotgui.plot_labels[ind1]['yposition'] = ypos1
        plotgui.plot_labels[ind1]['labelstring'] = labelstring
        plotgui.plot_labels[ind1]['colour'] = fontcolour
        plotgui.plot_labels[ind1]['size'] = fontsize
        plotgui.plot_labels[ind1]['font'] = font
        plotgui.plot_labels[ind1]['fontweight'] = fontweight
        make_plot.make_plot(plotgui)
    except Exception:
        pass

