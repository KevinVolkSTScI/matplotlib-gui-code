"""
Routines to save the current plot state to an ascii file, and to read back
such a file to make a plot.  All the routines need to be passed a
matplotlib_user_interface variable ("plotgui").

Routines:

save_plot     save the current plot information to an ascii file

load_plot     read an ascii file to make a plot

parse_save_file    attempt to read an ascii file for plot values, either to
                   test the file or to actually load in all the values


"""
import os
from copy import deepcopy
import tkinter as Tk
import tkinter.filedialog
import numpy
import general_utilities
import make_plot
import plot_flag_utilities
import window_creation

def save_plot(plotgui):
    """
    Save the plot state to an ascii file that can be loaded later.

    This routine writes the current sets and parameters to an ascii output
    file in a set format that can be read back in later.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None
    """
    outfilename = tkinter.filedialog.asksaveasfilename()
    if isinstance(outfilename, type('string')):
        outfile = open(outfilename, 'w')
        print('# matplotlib_user_interface.py save file version 1.1', file=outfile)
        print('# number of sets: %d maximum: %d' %
              (plotgui.nsets, plotgui.max_sets), file=outfile)
        for loop in range(plotgui.nsets):
            print('# set number %d: %d points' %
                  (loop+1, len(plotgui.xdata[loop]['values'])), file=outfile)
            for n1 in range(len(plotgui.xdata[loop]['values'])):
                print(plotgui.xdata[loop]['values'][n1],
                      plotgui.ydata[loop]['values'][n1],
                      plotgui.xdata[loop]['lowerror'][n1],
                      plotgui.ydata[loop]['lowerror'][n1],
                      plotgui.xdata[loop]['higherror'][n1],
                      plotgui.ydata[loop]['higherror'][n1], file=outfile)
            print('# set parameters', file=outfile)
            print(plotgui.xdata[loop]['minimum'],
                  plotgui.ydata[loop]['minimum'],
                  plotgui.xdata[loop]['maximum'],
                  plotgui.ydata[loop]['maximum'], file=outfile)
            print(plotgui.xdata[loop]['errors'], plotgui.ydata[loop]['errors'],
                  plotgui.xdata[loop]['legend'], plotgui.ydata[loop]['legend'],
                  file=outfile)
        print('# set properties', file=outfile)
        for loop in range(plotgui.nsets):
            str1 = '%s\t%g\t%s\t%g\t%s\t%s\t%g\t%g\t%g\t%g' % (
                plotgui.set_properties[loop]['symbol'],
                plotgui.set_properties[loop]['symbolsize'],
                plotgui.set_properties[loop]['linestyle'],
                plotgui.set_properties[loop]['linewidth'],
                plotgui.set_properties[loop]['colour'],
                plotgui.set_properties[loop]['label'],
                plotgui.set_properties[loop]['xmin'],
                plotgui.set_properties[loop]['xmax'],
                plotgui.set_properties[loop]['ymin'],
                plotgui.set_properties[loop]['ymax'])
            str1 = str1 + '\t' \
                + str(plotgui.set_properties[loop]['display']) \
                + '\t' + str(plotgui.set_properties[loop]['errors']) \
                + '\t' + str(plotgui.set_properties[loop]['legend']) \
                + '\t%d' % (plotgui.set_properties[loop]['plot'])
            print(str1, file=outfile)
        print('# plot properties', file=outfile)
        print('# nxplots, nyplots, total: %d %d %d ' %
              (plotgui.nxplots, plotgui.nyplots, plotgui.number_of_plots),
              file=outfile)
        str1 = '# hide plots: '
        for loop in range(plotgui.number_of_plots):
            str1 = str1 + str(plotgui.hide_subplot[loop]) + '\t'
        str1 = str1.rstrip('\t')
        print(str1, file=outfile)
        print('# margin: %f ' % (plotgui.plot_margin), file=outfile)
        nplot = 0
        for n1 in range(plotgui.nxplots):
            for n2 in range(plotgui.nyplots):
                print('# plot %d index values: %d %d' %
                      (nplot+1, n1, n2), file=outfile)
                print('# title: %s ' % (plotgui.title[nplot]), file=outfile)
                print('# frame: %f ' %
                      (plotgui.plot_frame[nplot]), file=outfile)
                print('# font name: %s ' %
                      (plotgui.fontname[nplot]), file=outfile)
                print('# font size: %s ' %
                      (plotgui.fontsize[nplot]), file=outfile)
                print('# font weight: %s ' %
                      (plotgui.fontweight[nplot]), file=outfile)
                strformat = '# x parameters: \t%s\t%g\t%g\t%g\t%g\t%d' \
                            + '\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%d' \
                            + '\t%d\t%d'
                print(strformat % (
                    plotgui.xparameters[nplot]['label'],
                    plotgui.xparameters[nplot]['minimum'],
                    plotgui.xparameters[nplot]['maximum'],
                    plotgui.xparameters[nplot]['major'],
                    plotgui.xparameters[nplot]['minor'],
                    plotgui.xparameters[nplot]['logarithmic'],
                    plotgui.xparameters[nplot]['invert'],
                    plotgui.xparameters[nplot]['hide'],
                    plotgui.xparameters[nplot]['hideticks'],
                    plotgui.xparameters[nplot]['hidelabels'],
                    plotgui.xparameters[nplot]['hybridlog'],
                    plotgui.xparameters[nplot]['inverseticks'],
                    plotgui.xparameters[nplot]['ticklength'],
                    plotgui.xparameters[nplot]['bothticks'],
                    plotgui.xparameters[nplot]['minorticks'],
                    plotgui.xparameters[nplot]['oppositeaxis'],
                    plotgui.xparameters[nplot]['majorgridlines'],
                    plotgui.xparameters[nplot]['minorgridlines']),
                      file=outfile)
                strformat = '# y parameters: \t%s\t%g\t%g\t%g\t%g\t%d' \
                    + '\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%d\t%d\t%d'
                print(strformat % (
                    plotgui.yparameters[nplot]['label'],
                    plotgui.yparameters[nplot]['minimum'],
                    plotgui.yparameters[nplot]['maximum'],
                    plotgui.yparameters[nplot]['major'],
                    plotgui.yparameters[nplot]['minor'],
                    plotgui.yparameters[nplot]['logarithmic'],
                    plotgui.yparameters[nplot]['invert'],
                    plotgui.yparameters[nplot]['hide'],
                    plotgui.yparameters[nplot]['hideticks'],
                    plotgui.yparameters[nplot]['hidelabels'],
                    plotgui.yparameters[nplot]['hybridlog'],
                    plotgui.yparameters[nplot]['inverseticks'],
                    plotgui.yparameters[nplot]['ticklength'],
                    plotgui.yparameters[nplot]['bothticks'],
                    plotgui.yparameters[nplot]['minorticks'],
                    plotgui.yparameters[nplot]['oppositeaxis'],
                    plotgui.yparameters[nplot]['majorgridlines'],
                    plotgui.yparameters[nplot]['minorgridlines']),
                      file=outfile)
                print('# plot range: %g %g %g %g %s' % (
                    plotgui.plot_range[nplot][0],
                    plotgui.plot_range[nplot][1],
                    plotgui.plot_range[nplot][2],
                    plotgui.plot_range[nplot][3],
                    plotgui.original_range[nplot]), file=outfile)
                nplot = nplot + 1
        if len(plotgui.xparameters) > plotgui.nxplots*plotgui.nyplots:
            for n1 in range(nplot, len(plotgui.xparameters)):
                print('# plot %d index values: %d %d' %
                      (nplot, -1, -1), file=outfile)
                print('# title: %s ' % (plotgui.title[nplot]), file=outfile)
                print('# frame: %f ' %
                      (plotgui.plot_frame[nplot]), file=outfile)
                print('# font name: %s ' %
                      (plotgui.fontname[nplot]), file=outfile)
                print('# font size: %s ' %
                      (plotgui.fontsize[nplot]), file=outfile)
                print('# font weight: %s ' %
                      (plotgui.fontweight[nplot]), file=outfile)
                strformat = '# x parameters: \t%s\t%g\t%g\t%g\t%g\t%d' \
                    + '\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%d\t%d\t%d'
                print(strformat % (
                    plotgui.xparameters[nplot]['label'],
                    plotgui.xparameters[nplot]['minimum'],
                    plotgui.xparameters[nplot]['maximum'],
                    plotgui.xparameters[nplot]['major'],
                    plotgui.xparameters[nplot]['minor'],
                    plotgui.xparameters[nplot]['logarithmic'],
                    plotgui.xparameters[nplot]['invert'],
                    plotgui.xparameters[nplot]['hide'],
                    plotgui.xparameters[nplot]['hideticks'],
                    plotgui.xparameters[nplot]['hidelabels'],
                    plotgui.xparameters[nplot]['hybridlog'],
                    plotgui.xparameters[nplot]['inverseticks'],
                    plotgui.xparameters[nplot]['ticklength'],
                    plotgui.xparameters[nplot]['bothticks'],
                    plotgui.xparameters[nplot]['minorticks'],
                    plotgui.xparameters[nplot]['oppositeaxis'],
                    plotgui.xparameters[nplot]['majorgridlines'],
                    plotgui.xparameters[nplot]['minorgridlines']),
                      file=outfile)
                strformat = '# y parameters: \t%s\t%g\t%g\t%g\t%g\t%d' \
                    + '\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%d\t%d\t%d'
                print(strformat % (
                    plotgui.yparameters[nplot]['label'],
                    plotgui.yparameters[nplot]['minimum'],
                    plotgui.yparameters[nplot]['maximum'],
                    plotgui.yparameters[nplot]['major'],
                    plotgui.yparameters[nplot]['minor'],
                    plotgui.yparameters[nplot]['logarithmic'],
                    plotgui.yparameters[nplot]['invert'],
                    plotgui.yparameters[nplot]['hide'],
                    plotgui.yparameters[nplot]['hideticks'],
                    plotgui.yparameters[nplot]['hidelabels'],
                    plotgui.yparameters[nplot]['hybridlog'],
                    plotgui.yparameters[nplot]['inverseticks'],
                    plotgui.yparameters[nplot]['ticklength'],
                    plotgui.yparameters[nplot]['bothticks'],
                    plotgui.yparameters[nplot]['minorticks'],
                    plotgui.yparameters[nplot]['oppositeaxis'],
                    plotgui.yparameters[nplot]['majorgridlines'],
                    plotgui.yparameters[nplot]['minorgridlines']),
                      file=outfile)
                print('# plot range: %g %g %g %g %s' % (
                    plotgui.plot_range[nplot][0],
                    plotgui.plot_range[nplot][1],
                    plotgui.plot_range[nplot][2],
                    plotgui.plot_range[nplot][3],
                    plotgui.original_range[nplot]), file=outfile)
                nplot = nplot + 1
        print('# number of labels: %d maximum: %d' %
              (plotgui.number_of_labels, plotgui.max_labels), file=outfile)
        for loop in range(plotgui.number_of_labels):
            print('%g\t%g\t%d\t%s\t%s\t%d\t%s\t%s' % (
                plotgui.plot_labels[loop]['xposition'],
                plotgui.plot_labels[loop]['yposition'],
                plotgui.plot_labels[loop]['plot'],
                plotgui.plot_labels[loop]['labelstring'],
                plotgui.plot_labels[loop]['colour'],
                plotgui.plot_labels[loop]['size'],
                plotgui.plot_labels[loop]['font'],
                plotgui.plot_labels[loop]['fontweight']), file=outfile)
        print('# number of lines: %d maximum: %d ' %
              (plotgui.number_of_lines, plotgui.max_lines), file=outfile)
        for loop in range(plotgui.number_of_lines):
            print('%g\t%g\t%g\t%g\t%d\t%s\t%s\t%g' % (
                plotgui.plot_lines[loop]['xstart'],
                plotgui.plot_lines[loop]['ystart'],
                plotgui.plot_lines[loop]['xend'],
                plotgui.plot_lines[loop]['yend'],
                plotgui.plot_lines[loop]['plot'],
                plotgui.plot_lines[loop]['line_type'],
                plotgui.plot_lines[loop]['line_colour'],
                plotgui.plot_lines[loop]['line_thickness']), file=outfile)
        print('# number of vectors: %d maximum: %d ' %
              (plotgui.number_of_vectors, plotgui.max_vectors), file=outfile)
        for loop in range(plotgui.number_of_vectors):
            print('%g\t%g\t%g\t%g\t%g\t%g\t%d\t%s\t%s\t%g\t%5.5s\t%s' % (
                plotgui.plot_vectors[loop]['xstart'],
                plotgui.plot_vectors[loop]['ystart'],
                plotgui.plot_vectors[loop]['xend'],
                plotgui.plot_vectors[loop]['yend'],
                plotgui.plot_vectors[loop]['delx'],
                plotgui.plot_vectors[loop]['dely'],
                plotgui.plot_vectors[loop]['plot'],
                plotgui.plot_vectors[loop]['line_type'],
                plotgui.plot_vectors[loop]['line_colour'],
                plotgui.plot_vectors[loop]['line_thickness'],
                str(plotgui.plot_vectors[loop]['fill']),
                plotgui.plot_vectors[loop]['fill_colour']), file=outfile)
        print('# number of ellipses: %d maximum: %d ' %
              (plotgui.number_of_ellipses, plotgui.max_ellipses), file=outfile)
        for loop in range(plotgui.number_of_ellipses):
            print('%g\t%g\t%g\t%g\t%g\t%d\t%s\t%s\t%g\t%s' % (
                plotgui.plot_ellipses[loop]['xposition'],
                plotgui.plot_ellipses[loop]['yposition'],
                plotgui.plot_ellipses[loop]['major'],
                plotgui.plot_ellipses[loop]['minor'],
                plotgui.plot_ellipses[loop]['rotation'],
                plotgui.plot_ellipses[loop]['plot'],
                plotgui.plot_ellipses[loop]['line_colour'],
                plotgui.plot_ellipses[loop]['line_type'],
                plotgui.plot_ellipses[loop]['line_thickness'],
                plotgui.plot_ellipses[loop]['fill_colour']), file=outfile)
        print('# number of boxes: %d maximum: %d ' %
              (plotgui.number_of_boxes, plotgui.max_boxes), file=outfile)
        for loop in range(plotgui.number_of_boxes):
            print('%g\t%g\t%g\t%g\t%g\t%d\t%s\t%s\t%g\t%s' % (
                plotgui.plot_boxes[loop]['xstart'],
                plotgui.plot_boxes[loop]['ystart'],
                plotgui.plot_boxes[loop]['xend'],
                plotgui.plot_boxes[loop]['yend'],
                plotgui.plot_boxes[loop]['rotation'],
                plotgui.plot_boxes[loop]['plot'],
                plotgui.plot_boxes[loop]['line_colour'],
                plotgui.plot_boxes[loop]['line_type'],
                plotgui.plot_boxes[loop]['line_thickness'],
                plotgui.plot_boxes[loop]['fill_colour']), file=outfile)
        print('# legend values:', file=outfile)
        nplot = 0
        for n1 in range(plotgui.nxplots):
            for n2 in range(plotgui.nyplots):
                str1 = '# plot %d index values %d %d: \t' % (nplot, n1, n2)
                if plotgui.legend_variable[nplot] is None:
                    str1 = str1 + 'None\tNone\tNone\tNone\tNone\t'
                else:
                    str1 = str1 + '%d\t%d\t%s\t%s\t' % (
                        plotgui.legend_variable[nplot].get(),
                        plotgui.legend_frame[nplot].get(),
                        plotgui.legend_position[nplot],
                        plotgui.legend_options[nplot].get())
                    try:
                        str1 = str1 + '%g\t%g' % (
                            plotgui.legend_user_position[nplot][0],
                            plotgui.legend_user_position[nplot][1])
                    except Exception:
                        str1 = str1 + 'None\tNone'
                print(str1, file=outfile)
                nplot = nplot + 1
        if len(plotgui.legend_variable) > plotgui.nxplots*plotgui.nyplots:
            for n1 in range(nplot, len(plotgui.xparameters)):
                str1 = '# plot %d index values %d %d: \t' % (nplot, -1, -1)
                if plotgui.legend_variable[nplot] is None:
                    str1 = str1 + 'None\tNone\tNone\tNone\tNone\t'
                else:
                    str1 = str1 + '%d\t%d\t%s\t%s\t' % (
                        plotgui.legend_variable[nplot].get(),
                        plotgui.legend_frame[nplot].get(),
                        plotgui.legend_position[nplot],
                        plotgui.legend_options[nplot].get())
                    try:
                        str1 = str1 + '%g\t%g' % (
                            plotgui.legend_user_position[nplot][0],
                            plotgui.legend_user_position[nplot][1])
                    except Exception:
                        str1 = str1 + 'None\tNone'
                print(str1, file=outfile)
                nplot = nplot + 1
        for n1 in range(len(plotgui.share_axis)):
            if n1 == 0:
                str1 = ' %d ' % (plotgui.share_axis[n1])
            else:
                str1 = str1 + '\t %d ' % (plotgui.share_axis[n1])
        print('# share_axis', file=outfile)
        print(str1, file=outfile)
        str1 = ''
        for n1 in range(len(plotgui.equal_aspect)):
            if plotgui.equal_aspect[n1]:
                str1 = str1 + ' True' + '\t'
            else:
                str1 = str1 + 'False' + '\t'
        str1 = str1.rstrip('\t')
        print('# equal_aspect', file=outfile)
        print(str1, file=outfile)
        print('# end', file=outfile)
        outfile.close()
    else:
        tkinter.messagebox.showinfo(
            "Error",
            "File "+outfilename+" was not written properly.")

def load_plot(plotgui):
    """
    Read an ascii file of plot parameters to make a plot.

    This routine asks for the ascii file to load and calls the routine
    to parse the file and load in the data and parameters.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    savefilename = tkinter.filedialog.askopenfilename()
    if savefilename is None:
        return
    if not os.path.isfile(savefilename):
        tkinter.messagebox.showinfo(
            "Error",
            "The file %s was not found." % (savefilename))
        return
    savefile = open(savefilename, 'r')
    lines = savefile.readlines()
    savefile.close()
    flag = parse_save_file(plotgui, lines, False)
    if not flag:
        tkinter.messagebox.showinfo(
            "Error",
            "The file %s is not formatted properly for a save file." %
            (savefilename))
        return
    response = tkinter.messagebox.askyesno(
        "Verify",
        "Do you want to abandon the current plot for the saved one?")
    if response:
        plot_flag_utilities.clear_plot(plotgui, False)
        flag = parse_save_file(plotgui, lines, True)
        save_plot_range = deepcopy(plotgui.plot_range)
        save_xpars = deepcopy(plotgui.xparameters)
        share_axis = deepcopy(plotgui.share_axis)
        save_ypars = deepcopy(plotgui.yparameters)
        nx1 = 1*plotgui.nxplots
        ny1 = 1*plotgui.nyplots
        nplots = 1*plotgui.number_of_plots
        plotgui.nxplots = 0
        plotgui.nyplots = 0
        n1 = plotgui.current_plot
        plotgui.subplot = []
        plotgui.hide_subplot = [False, ]
        plotgui.make_plot_layout(nx1, ny1, 1)
        plotgui.share_axis = deepcopy(share_axis)
        if plotgui.number_of_plots < nplots:
            for loop in range(len(plotgui.share_axis)):
                if plotgui.share_axis[loop] != 0:
                    n2 = abs(plotgui.share_axis[loop])
                    plotgui.current_plot = n2
                    if plotgui.share_axis[loop] < 0:
                        plotgui.subplot.append(plotgui.figure.add_subplot(
                            plotgui.nxplots, plotgui.nyplots, n2,
                            sharey=plotgui.subplot[n2-1],
                            frameon=False))
                        plotgui.bounding_box.append(plotgui.bounding_box[n2-1])
                        plotgui.current_axis = len(plotgui.subplot)
                    else:
                        plotgui.subplot.append(plotgui.figure.add_subplot(
                            plotgui.nxplots, plotgui.nyplots, n2,
                            sharex=plotgui.subplot[n2-1],
                            frameon=False))
                        plotgui.bounding_box.append(plotgui.bounding_box[n2-1])
                        plotgui.current_axis = len(plotgui.subplot)
            plotgui.number_of_plots = nplots
        plotgui.xparameters = deepcopy(save_xpars)
        plotgui.yparameters = deepcopy(save_ypars)
        for loop in range(len(plotgui.subplot)):
            plotgui.current_plot = loop+1
            make_plot.make_plot(plotgui)
        plotgui.current_plot = n1
        plotgui.plot_range = deepcopy(save_plot_range)
        make_plot.make_plot(plotgui)

def parse_save_file(plotgui, lines, flag):
    """
    Parse an ascii save file and optionally load the values.

    This program reads the lines from a matplotlib_user_interface.py save
    file.  It determines whether the file is structured properly.  If the
    flag value is True it also tries to set the parameters for the plot.

    Parameters
    ----------

        plotgui : the matplotlib_user_interface object holding the plot

        lines : A set of limes (assumed to be from the .readlines()
                function) from a matplotlib_user_interface.py save file

        flag :  A boolean value, if True the code tries to assign the
                values for the plot; if False it only tries to parse
                the file

    Returns
    -------

        goodfile :  A boolean value for whether the lines are of the
                    expected structure for a matplotlib_user_interface.py
                    save file

    """
    goodfile = False
    equal_aspect = None
    ndatasets = 0
    xdata = []
    ydata = []
    set_properties = []
    set_parameters = []
    plot_parameters = []
    labels = []
    plines = []
    vectors = []
    ellipses = []
    boxes = []
    legend = []
    for ind1 in range(len(lines)):
        line = lines[ind1]
        line = line.strip('\n')
        if '# matplotlib_user_interface.py save file' in line:
            values = line.split()
            try:
                version = float(values[-1])
            except ValueError:
                return goodfile
            if version != plotgui.save_version:
                return goodfile
        if '# number of sets:' in line:
            values = line.split()
            if len(values) != 7:
                return goodfile
            try:
                nsets = int(values[4])
                maxsets = int(values[6])
            except ValueError:
                return goodfile
            if (nsets < 1) or (nsets > maxsets):
                return goodfile
        if '# set number' in line:
            ind2 = general_utilities.line_range(lines, ind1)
            npoints = ind2 - ind1 - 1
            xset = numpy.zeros((npoints), dtype=numpy.float32)
            yset = numpy.zeros((npoints), dtype=numpy.float32)
            xseterr1 = numpy.zeros((npoints), dtype=numpy.float32)
            yseterr1 = numpy.zeros((npoints), dtype=numpy.float32)
            xseterr2 = numpy.zeros((npoints), dtype=numpy.float32)
            yseterr2 = numpy.zeros((npoints), dtype=numpy.float32)
            for loop in range(ind1+1, ind2):
                n1 = loop - ind1 - 1
                line1 = lines[loop].strip('\n')
                values = line1.split()
                if len(values) == 6:
                    xset[n1] = values[0]
                    yset[n1] = values[1]
                    xseterr1[n1] = values[2]
                    yseterr1[n1] = values[3]
                    xseterr2[n1] = values[4]
                    yseterr2[n1] = values[5]
                else:
                    return goodfile
            xdata.append(xset)
            xdata.append(xseterr1)
            xdata.append(xseterr2)
            ydata.append(yset)
            ydata.append(yseterr1)
            ydata.append(yseterr2)
            ndatasets = ndatasets + 1
        if '# set parameters' in line:
            p1 = {}
            line1 = lines[ind1+1].strip('\n')
            values = line1.split()
            if len(values) == 4:
                p1['xminimum'] = float(values[0])
                p1['yminimum'] = float(values[1])
                p1['xmaximum'] = float(values[2])
                p1['ymaximum'] = float(values[3])
            line2 = lines[ind1+2].strip('\n')
            values = line2.split()
            if len(values) == 4:
                if values[0] == 'True':
                    p1['xerrors'] = True
                elif values[0] == 'False':
                    p1['xerrors'] = False
                else:
                    return goodfile
                if values[1] == 'True':
                    p1['yerrors'] = True
                elif values[1] == 'False':
                    p1['yerrors'] = False
                else:
                    return goodfile
                if values[2] == 'True':
                    p1['xlegend'] = True
                elif values[2] == 'False':
                    p1['xlegend'] = False
                else:
                    return goodfile
                if values[3] == 'True':
                    p1['ylegend'] = True
                elif values[3] == 'False':
                    p1['ylegend'] = False
                else:
                    return goodfile
            set_parameters.append(p1)
        if '# set properties' in line:
            if ndatasets != nsets:
                return goodfile
            ind2 = general_utilities.line_range(lines, ind1)
            npoints = ind2 - ind1 - 1
            if npoints == nsets:
                for loop in range(ind1+1, ind2):
                    p2 = {}
                    line = lines[loop].strip('\n')
                    values = line.split('\t')
                    if len(values) != 14:
                        return goodfile
                    p2['symbol'] = values[0]
                    if p2['symbol'] == 'None':
                        p2['symbol'] = None
                    p2['symbolsize'] = float(values[1])
                    p2['linestyle'] = values[2]
                    if p2['linestyle'] == 'None':
                        p2['linestyle'] = None
                    p2['linewidth'] = float(values[3])
                    p2['colour'] = values[4]
                    p2['label'] = values[5]
                    p2['xmin'] = float(values[6])
                    p2['xmax'] = float(values[7])
                    p2['ymin'] = float(values[8])
                    p2['ymax'] = float(values[9])
                    if values[10] == 'True':
                        p2['display'] = True
                    elif values[10] == 'False':
                        p2['display'] = False
                    else:
                        return goodfile
                    if values[11] == 'True':
                        p2['errors'] = True
                    elif values[11] == 'False':
                        p2['errors'] = False
                    else:
                        return goodfile
                    if values[12] == 'True':
                        p2['legend'] = True
                    elif values[12] == 'False':
                        p2['legend'] = False
                    else:
                        return goodfile
                    p2['plot'] = int(values[13])
                    set_properties.append(p2)
            else:
                return goodfile
        if ('# plot ' in line) and ('index values:' in line):
            pp = {}
        if '# nxplots, nyplots, total:' in line:
            values = line.split()
            nxplots = int(values[-3])
            nyplots = int(values[-2])
            nplots = int(values[-1])
            if (nxplots < 1) or (nyplots < 1):
                return goodfile
            hide_subplot = []
            for loop in range(nplots):
                hide_subplot.append(False)
        if '# hide plots' in line:
            line1 = line.replace('# hide plots: ', '')
            values = line.split('\t')
            if len(values) == nplots:
                for loop in range(nplots):
                    if 'True' in values[loop]:
                        hide_subplot[loop] = True
                    else:
                        hide_subplot[loop] = False
        if '# margin' in line:
            values = line.split()
            if len(values) != 3:
                return goodfile
            plot_margin = float(values[2])
            if plot_margin < 0.:
                return goodfile
        if '# title:' in line:
            frag = line.replace('# title:', '')
            frag = frag.lstrip()
            frag = frag.rstrip()
            pp['title'] = frag
        if '# frame' in line:
            values = line.split()
            if len(values) != 3:
                return goodfile
            pp['frame'] = float(values[2])
        if '# font name:' in line:
            frag = line.replace('# font name:', '')
            frag = frag.lstrip()
            frag = frag.rstrip()
            pp['fontname'] = frag
        if '# font size' in line:
            values = line.split()
            if len(values) != 4:
                return goodfile
            pp['fontsize'] = int(values[3])
        if '# font weight:' in line:
            frag = line.replace('# font weight:', '')
            frag = frag.lstrip()
            frag = frag.rstrip()
            pp['fontweight'] = frag
        if '# x parameters:' in line:
            line = line.strip('\n')
            values = line.split('\t')
            if len(values) != 19:
                return goodfile
            xp = {}
            xp['label'] = values[1]
            xp['minimum'] = values[2]
            xp['maximum'] = values[3]
            xp['major'] = float(values[4])
            xp['minor'] = float(values[5])
            xp['logarithmic'] = int(values[6])
            xp['invert'] = int(values[7])
            xp['hide'] = int(values[8])
            xp['hideticks'] = int(values[9])
            xp['hidelabels'] = int(values[10])
            xp['hybridlog'] = int(values[11])
            xp['inverseticks'] = int(values[12])
            xp['ticklength'] = int(values[13])
            xp['bothticks'] = int(values[14])
            xp['minorticks'] = float(values[15])
            xp['oppositeaxis'] = int(values[16])
            xp['majorgridlines'] = int(values[17])
            xp['minorgridlines'] = int(values[18])
            pp['xparameters'] = xp
        if '# y parameters:' in line:
            values = line.split('\t')
            if len(values) != 19:
                return goodfile
            yp = {}
            yp['label'] = values[1]
            yp['minimum'] = values[2]
            yp['maximum'] = values[3]
            yp['major'] = float(values[4])
            yp['minor'] = float(values[5])
            yp['logarithmic'] = int(values[6])
            yp['invert'] = int(values[7])
            yp['hide'] = int(values[8])
            yp['hideticks'] = int(values[9])
            yp['hidelabels'] = int(values[10])
            yp['hybridlog'] = int(values[11])
            yp['inverseticks'] = int(values[12])
            yp['ticklength'] = int(values[13])
            yp['bothticks'] = int(values[14])
            yp['minorticks'] = float(values[15])
            yp['oppositeaxis'] = int(values[16])
            yp['majorgridlines'] = int(values[17])
            yp['minorgridlines'] = int(values[18])
            pp['yparameters'] = yp
        if '# plot range:' in line:
            values = line.split()
            if len(values) != 8:
                return goodfile
            try:
                xmin = float(values[3])
                xmax = float(values[4])
                ymin = float(values[5])
                ymax = float(values[6])
            except ValueError:
                return goodfile
            pp['plot_range'] = [0., 1., 0., 1.]
            pp['plot_range'][0] = 1.*xmin
            pp['plot_range'][1] = 1.*xmax
            pp['plot_range'][2] = 1.*ymin
            pp['plot_range'][3] = 1.*ymax
            if 'True' in values[7]:
                pp['original_range'] = True
            elif 'False' in values[7]:
                pp['original_range'] = False
            else:
                return goodfile
            plot_parameters.append(pp)
        if '# number of labels:' in line:
            values = line.split()
            if len(values) != 7:
                return goodfile
            if int(values[4]) > 0:
                ind2 = general_utilities.line_range(lines, ind1)
                npoints = ind2 - ind1 - 1
                for loop in range(ind1+1, ind2):
                    lv = {}
                    line = lines[loop].strip('\n')
                    values = line.split('\t')
                    if len(values) != 8:
                        return goodfile
                    lv['xposition'] = float(values[0])
                    lv['yposition'] = float(values[1])
                    lv['plot'] = int(values[2])
                    lv['labelstring'] = values[3]
                    lv['colour'] = values[4]
                    lv['size'] = int(values[5])
                    lv['font'] = values[6]
                    lv['fontweight'] = values[7]
                    labels.append(lv)
        if '# number of lines:' in line:
            values = line.split()
            if len(values) != 7:
                return goodfile
            if int(values[4]) > 0:
                ind2 = general_utilities.line_range(lines, ind1)
                npoints = ind2 - ind1 - 1
                for loop in range(ind1+1, ind2):
                    lv = {}
                    line = lines[loop].strip('\n')
                    values = line.split('\t')
                    if len(values) != 8:
                        return goodfile
                    lv['xstart'] = float(values[0])
                    lv['ystart'] = float(values[1])
                    lv['xend'] = float(values[2])
                    lv['yend'] = float(values[3])
                    lv['plot'] = int(values[4])
                    lv['line_type'] = values[5]
                    lv['line_colour'] = values[6]
                    lv['line_width'] = float(values[7])
                    plines.append(lv)
        if '# number of vectors:' in line:
            values = line.split()
            if len(values) != 7:
                return goodfile
            if int(values[4]) > 0:
                ind2 = general_utilities.line_range(lines, ind1)
                npoints = ind2 - ind1 - 1
                for loop in range(ind1+1, ind2):
                    v1 = {}
                    line = lines[loop].strip('\n')
                    values = line.split('\t')
                    if len(values) != 12:
                        return goodfile
                    v1['xstart'] = float(values[0])
                    v1['ystart'] = float(values[1])
                    v1['xend'] = float(values[2])
                    v1['yend'] = float(values[3])
                    v1['delx'] = float(values[4])
                    v1['dely'] = float(values[5])
                    v1['plot'] = int(values[6])
                    v1['line_type'] = values[7]
                    v1['line_colour'] = values[8]
                    v1['line_thickness'] = float(values[9])
                    if values[10] == ' True':
                        v1['fill'] = True
                    elif values[10] == 'False':
                        v1['fill'] = False
                    else:
                        return goodfile
                    v1['fill_colour'] = values[11]
                    vectors.append(v1)
        if '# number of ellipses:' in line:
            values = line.split()
            if len(values) != 7:
                return goodfile
            if int(values[4]) > 0:
                ind2 = general_utilities.line_range(lines, ind1)
                npoints = ind2 - ind1 - 1
                for loop in range(ind1+1, ind2):
                    e1 = {}
                    line = lines[loop].strip('\n')
                    values = line.split('\t')
                    if len(values) != 10:
                        return goodfile
                    e1['xstart'] = float(values[0])
                    e1['ystart'] = float(values[1])
                    e1['major'] = float(values[2])
                    e1['minor'] = float(values[3])
                    e1['rotation'] = float(values[4])
                    e1['plot'] = int(values[5])
                    e1['line_colour'] = values[6]
                    e1['line_type'] = values[7]
                    e1['line_thickness'] = float(values[8])
                    e1['fill_colour'] = values[9]
                    ellipses.append(e1)
        if '# number of boxes:' in line:
            values = line.split()
            if len(values) != 7:
                return goodfile
            if int(values[4]) > 0:
                ind2 = general_utilities.line_range(lines, ind1)
                npoints = ind2 - ind1 - 1
                for loop in range(ind1+1, ind2):
                    b1 = {}
                    line = lines[loop].strip('\n')
                    values = line.split('\t')
                    if len(values) != 10:
                        return goodfile
                    b1['xstart'] = float(values[0])
                    b1['ystart'] = float(values[1])
                    b1['xend'] = float(values[2])
                    b1['yend'] = float(values[3])
                    b1['rotation'] = float(values[4])
                    b1['plot'] = int(values[5])
                    b1['line_colour'] = values[6]
                    b1['line_type'] = values[7]
                    b1['line_thickness'] = float(values[8])
                    b1['fill_colour'] = values[9]
                    boxes.append(b1)
        if '# legend values:' in line:
            legend = []
            for loop in range(nplots):
                line1 = lines[ind1+1+loop].strip('\n')
                values = line1.split('\t')
                if len(values) != 7:
                    return goodfile
                lg = {}
                if values[1] == 'None':
                    lg['legend_variable_value'] = None
                    lg['legend_frame_value'] = None
                    lg['legend_position'] = None
                    lg['legend_option_value'] = None
                    lg['legend_user_position'] = None
                else:
                    lg['legend_variable_value'] = int(values[1])
                    lg['legend_frame_value'] = int(values[2])
                    lg['legend_option_value'] = values[3]
                    lg['legend_position'] = values[4]
                    if values[5] == 'None':
                        lg['legend_user_position'] = None
                    else:
                        lg['legend_user_xposition'] = [float(values[5]),
                                                       float(values[6])]
                legend.append(lg)
        if '# share_axis' in line:
            share_axis = []
            line1 = lines[ind1+1].strip('\n')
            values = line1.split('\t')
            for loop in range(len(values)):
                share_axis.append(int(values[loop]))
        if '# equal_aspect' in line:
            equal_aspect = []
            line1 = lines[ind1+1].strip('\n')
            values = line1.split('\t')
            for loop in range(len(values)):
                if 'True' in values[loop]:
                    equal_aspect.append(True)
                else:
                    equal_aspect.append(False)
    goodfile = True
    if flag:
        plotgui.maxsets = maxsets
        for loop in range(nsets):
            ind1 = 3 * loop
            plotgui.add_set(xdata[ind1], ydata[ind1],
                xlowerror=xdata[ind1+1],
                xhigherror=xdata[ind1+2],
                ylowerror=ydata[ind1+1],
                yhigherror=ydata[ind1+2])
            plotgui.xdata[loop]['minimum'] = set_parameters[loop]['xminimum']
            plotgui.ydata[loop]['minimum'] = set_parameters[loop]['yminimum']
            plotgui.xdata[loop]['maximum'] = set_parameters[loop]['xmaximum']
            plotgui.ydata[loop]['maximum'] = set_parameters[loop]['ymaximum']
            plotgui.xdata[loop]['errors'] = set_parameters[loop]['xerrors']
            plotgui.ydata[loop]['errors'] = set_parameters[loop]['yerrors']
            plotgui.xdata[loop]['legend'] = set_parameters[loop]['xlegend']
            plotgui.ydata[loop]['legend'] = set_parameters[loop]['ylegend']
            plotgui.set_properties[loop]['symbol'] = \
                set_properties[loop]['symbol']
            plotgui.set_properties[loop]['symbolsize'] = \
                set_properties[loop]['symbolsize']
            plotgui.set_properties[loop]['linestyle'] = \
                set_properties[loop]['linestyle']
            plotgui.set_properties[loop]['linewidth'] = \
                set_properties[loop]['linewidth']
            plotgui.set_properties[loop]['colour'] = \
                set_properties[loop]['colour']
            plotgui.set_properties[loop]['label'] = \
                set_properties[loop]['label']
            plotgui.set_properties[loop]['xmin'] = \
                set_properties[loop]['xmin']
            plotgui.set_properties[loop]['xmax'] = \
                set_properties[loop]['xmax']
            plotgui.set_properties[loop]['ymin'] = \
                set_properties[loop]['ymin']
            plotgui.set_properties[loop]['ymax'] = \
                set_properties[loop]['ymax']
            plotgui.set_properties[loop]['display'] = \
                set_properties[loop]['display']
            plotgui.set_properties[loop]['errors'] = \
                set_properties[loop]['errors']
            plotgui.set_properties[loop]['legend'] = \
                set_properties[loop]['legend']
            plotgui.set_properties[loop]['plot'] = \
                set_properties[loop]['plot']
            plotgui.plot_margin = plot_margin
        plotgui.nsets = nsets
        plotgui.title = []
        plotgui.plot_frame = []
        plotgui.plot_range = []
        plotgui.fontname = []
        plotgui.fontsize = []
        plotgui.fontweight = []
        plotgui.xparameters = []
        plotgui.yparameters = []
        for loop in range(nplots):
            plotgui.title.append(plot_parameters[loop]['title'])
            plotgui.plot_frame.append(plot_parameters[loop]['frame'])
            plotgui.fontname.append(plot_parameters[loop]['fontname'])
            plotgui.fontsize.append(plot_parameters[loop]['fontsize'])
            plotgui.fontweight.append(plot_parameters[loop]['fontweight'])
            plotgui.xparameters.append({'label': ' ', 'minimum': 0.0,
                                        'maximum': 1.0, 'major': 0.1,
                                        'minor': 0.1, 'logarithmic': 0,
                                        'invert': 0, 'hide': 0,
                                        'hideticks': 0, 'hidelabels': 0,
                                        'hybridlog': 0, 'inverseticks': 0,
                                        'ticklength': 6, 'bothticks': 0,
                                        'minorticks': 0,
                                        'oppositeaxis': 0,
                                        'majorgridlines': 0,
                                        'minorgridlines': 0})
            plotgui.yparameters.append({'label': ' ', 'minimum': 0.0,
                                        'maximum': 1.0, 'major': 0.1,
                                        'minor': 0.1, 'logarithmic': 0,
                                        'invert': 0, 'hide': 0,
                                        'hideticks': 0, 'hidelabels': 0,
                                        'hybridlog': 0, 'inverseticks': 0,
                                        'ticklength': 6, 'bothticks': 0,
                                        'minorticks': 0,
                                        'oppositeaxis': 0,
                                        'majorgridlines': 0,
                                        'minorgridlines': 0})
            plotgui.plot_range.append(plot_parameters[loop]['plot_range'])
            plotgui.original_range.append(True)
            plotgui.xparameters[loop]['label'] = \
                plot_parameters[loop]['xparameters']['label']
            plotgui.xparameters[loop]['minimum'] = \
                plot_parameters[loop]['xparameters']['minimum']
            plotgui.xparameters[loop]['maximum'] = \
                plot_parameters[loop]['xparameters']['maximum']
            plotgui.xparameters[loop]['major'] = \
                plot_parameters[loop]['xparameters']['major']
            plotgui.xparameters[loop]['minor'] = \
                plot_parameters[loop]['xparameters']['minor']
            plotgui.xparameters[loop]['logarithmic'] = \
                plot_parameters[loop]['xparameters']['logarithmic']
            plotgui.xparameters[loop]['invert'] = \
                plot_parameters[loop]['xparameters']['invert']
            plotgui.xparameters[loop]['hide'] = \
                plot_parameters[loop]['xparameters']['hide']
            plotgui.xparameters[loop]['hideticks'] = \
                plot_parameters[loop]['xparameters']['hideticks']
            plotgui.xparameters[loop]['hidelabels'] = \
                plot_parameters[loop]['xparameters']['hidelabels']
            plotgui.xparameters[loop]['hybridlog'] = \
                plot_parameters[loop]['xparameters']['hybridlog']
            plotgui.xparameters[loop]['inverseticks'] = \
                plot_parameters[loop]['xparameters']['inverseticks']
            plotgui.xparameters[loop]['ticklength'] = \
                plot_parameters[loop]['xparameters']['ticklength']
            plotgui.xparameters[loop]['bothticks'] = \
                plot_parameters[loop]['xparameters']['bothticks']
            plotgui.xparameters[loop]['minorticks'] = \
                plot_parameters[loop]['xparameters']['minorticks']
            plotgui.xparameters[loop]['oppositeaxis'] = \
                plot_parameters[loop]['xparameters']['oppositeaxis']
            plotgui.xparameters[loop]['majorgridlines'] = \
                plot_parameters[loop]['xparameters']['majorgridlines']
            plotgui.xparameters[loop]['minorgridlines'] = \
                plot_parameters[loop]['xparameters']['minorgridlines']
            plotgui.yparameters[loop]['label'] = \
                plot_parameters[loop]['yparameters']['label']
            plotgui.yparameters[loop]['minimum'] = \
                plot_parameters[loop]['yparameters']['minimum']
            plotgui.yparameters[loop]['maximum'] = \
                plot_parameters[loop]['yparameters']['maximum']
            plotgui.yparameters[loop]['major'] = \
                plot_parameters[loop]['yparameters']['major']
            plotgui.yparameters[loop]['minor'] = \
                plot_parameters[loop]['yparameters']['minor']
            plotgui.yparameters[loop]['logarithmic'] = \
                plot_parameters[loop]['yparameters']['logarithmic']
            plotgui.yparameters[loop]['invert'] = \
                plot_parameters[loop]['yparameters']['invert']
            plotgui.yparameters[loop]['hide'] = \
                plot_parameters[loop]['yparameters']['hide']
            plotgui.yparameters[loop]['hideticks'] = \
                plot_parameters[loop]['yparameters']['hideticks']
            plotgui.yparameters[loop]['hidelabels'] = \
                plot_parameters[loop]['yparameters']['hidelabels']
            plotgui.yparameters[loop]['hybridlog'] = \
                plot_parameters[loop]['yparameters']['hybridlog']
            plotgui.yparameters[loop]['inverseticks'] = \
                plot_parameters[loop]['yparameters']['inverseticks']
            plotgui.yparameters[loop]['ticklength'] = \
                plot_parameters[loop]['yparameters']['ticklength']
            plotgui.yparameters[loop]['bothticks'] = \
                plot_parameters[loop]['yparameters']['bothticks']
            plotgui.yparameters[loop]['minorticks'] = \
                plot_parameters[loop]['yparameters']['minorticks']
            plotgui.yparameters[loop]['oppositeaxis'] = \
                plot_parameters[loop]['yparameters']['oppositeaxis']
            plotgui.original_range[loop] = \
                plot_parameters[loop]['original_range']
            plotgui.yparameters[loop]['majorgridlines'] = \
                plot_parameters[loop]['yparameters']['majorgridlines']
            plotgui.yparameters[loop]['minorgridlines'] = \
                plot_parameters[loop]['yparameters']['minorgridlines']
        if len(labels) > 0:
            plotgui.number_of_labels = len(labels)
            for loop in range(len(labels)):
                plotgui.plot_labels[loop]['xposition'] = \
                    labels[loop]['xposition']
                plotgui.plot_labels[loop]['yposition'] = \
                    labels[loop]['yposition']
                plotgui.plot_labels[loop]['plot'] = labels[loop]['plot']
                plotgui.plot_labels[loop]['labelstring'] = \
                    labels[loop]['labelstring']
                plotgui.plot_labels[loop]['colour'] = labels[loop]['colour']
                plotgui.plot_labels[loop]['size'] = labels[loop]['size']
                plotgui.plot_labels[loop]['font'] = labels[loop]['font']
                plotgui.plot_labels[loop]['fontweight'] = \
                    labels[loop]['fontweight']
        if len(plines) > 0:
            plotgui.number_of_lines = len(plines)
            for loop in range(len(plines)):
                plotgui.plot_lines[loop]['xstart'] = plines[loop]['xstart']
                plotgui.plot_lines[loop]['ystart'] = plines[loop]['ystart']
                plotgui.plot_lines[loop]['xend'] = plines[loop]['xend']
                plotgui.plot_lines[loop]['yend'] = plines[loop]['yend']
                plotgui.plot_lines[loop]['plot'] = plines[loop]['plot']
                plotgui.plot_lines[loop]['line_type'] = \
                    plines[loop]['line_type']
                plotgui.plot_lines[loop]['line_colour'] = \
                    plines[loop]['line_colour']
                plotgui.plot_lines[loop]['line_width'] = \
                    plines[loop]['line_width']
        if len(vectors) > 0:
            plotgui.number_of_vectors = len(vectors)
            for loop in range(len(vectors)):
                plotgui.plot_vectors[loop]['xstart'] = vectors[loop]['xstart']
                plotgui.plot_vectors[loop]['ystart'] = vectors[loop]['ystart']
                plotgui.plot_vectors[loop]['xend'] = vectors[loop]['xend']
                plotgui.plot_vectors[loop]['yend'] = vectors[loop]['yend']
                plotgui.plot_vectors[loop]['delx'] = vectors[loop]['delx']
                plotgui.plot_vectors[loop]['dely'] = vectors[loop]['dely']
                plotgui.plot_vectors[loop]['plot'] = vectors[loop]['plot']
                plotgui.plot_vectors[loop]['line_type'] = \
                    vectors[loop]['line_type']
                plotgui.plot_vectors[loop]['line_colour'] = \
                    vectors[loop]['line_colour']
                plotgui.plot_vectors[loop]['line_thickness'] = \
                    vectors[loop]['line_thickness']
                plotgui.plot_vectors[loop]['fill'] = vectors[loop]['fill']
                plotgui.plot_vectors[loop]['fill_colour'] = vectors[loop]['fill_colour']
        if len(ellipses) > 0:
            plotgui.number_of_ellipses = len(ellipses)
            for loop in range(len(ellipses)):
                plotgui.plot_ellipses[loop]['xstart'] = \
                    ellipses[loop]['xstart']
                plotgui.plot_ellipses[loop]['ystart'] = \
                    ellipses[loop]['ystart']
                plotgui.plot_ellipses[loop]['major'] = ellipses[loop]['major']
                plotgui.plot_ellipses[loop]['minor'] = ellipses[loop]['minor']
                plotgui.plot_ellipses[loop]['rotation'] = \
                    ellipses[loop]['rotation']
                plotgui.plot_ellipses[loop]['plot'] = ellipses[loop]['plot']
                plotgui.plot_ellipses[loop]['line_type'] = \
                    ellipses[loop]['line_type']
                plotgui.plot_ellipses[loop]['line_colour'] = \
                    ellipses[loop]['line_colour']
                plotgui.plot_ellipses[loop]['line_thickness'] = \
                    ellipses[loop]['line_thickness']
                plotgui.plot_ellipses[loop]['fill_colour'] = \
                    ellipses[loop]['fill_colour']
        if len(boxes) > 0:
            plotgui.number_of_boxes = len(boxes)
            for loop in range(len(boxes)):
                plotgui.plot_boxes[loop]['xstart'] = boxes[loop]['xstart']
                plotgui.plot_boxes[loop]['ystart'] = boxes[loop]['ystart']
                plotgui.plot_boxes[loop]['xend'] = boxes[loop]['xend']
                plotgui.plot_boxes[loop]['yend'] = boxes[loop]['yend']
                plotgui.plot_boxes[loop]['rotation'] = boxes[loop]['rotation']
                plotgui.plot_boxes[loop]['plot'] = boxes[loop]['plot']
                plotgui.plot_boxes[loop]['line_type'] = \
                    boxes[loop]['line_type']
                plotgui.plot_boxes[loop]['line_colour'] = \
                    boxes[loop]['line_colour']
                plotgui.plot_boxes[loop]['line_thickness'] = \
                    boxes[loop]['line_thickness']
                plotgui.plot_boxes[loop]['fill_colour'] = \
                    boxes[loop]['fill_colour']
#        nplots = plotgui.nxplots * plotgui.nyplots
        plotgui.legend_position = []
        plotgui.legend_user_position = []
        for loop in range(nplots):
            plotgui.legend_position.append(legend[loop]['legend_position'])
            plotgui.legend_user_position.append(legend[loop][
                'legend_user_position'])
            if legend[loop]['legend_variable_value'] is None:
                if loop < len(plotgui.legend_variable):
                    try:
                        plotgui.legend_variable[loop].set(0)
                    except Exception:
                        plotgui.legend_variable.append(None)
                else:
                    plotgui.legend_variable.append(None)
            else:
                if loop < len(plotgui.legend_variable):
                    try:
                        # check if the legend_variable exists, if not
                        # create it
                        plotgui.legend_variable[loop].set(
                            legend[loop]['legend_variable_value'])
                    except Exception:
                        plotgui.legend_variable[loop] = Tk.IntVar()
                        plotgui.legend_variable[loop].set(
                            legend[loop]['legend_variable_value'])
                else:
                    plotgui.legend_variable.append(Tk.IntVar())
                    plotgui.legend_variable[loop].set(
                        legend[loop]['legend_variable_value'])
            if legend[loop]['legend_frame_value'] is None:
                if loop < len(plotgui.legend_frame):
                    # Set the legend_frame value, if possible
                    try:
                        plotgui.legend_frame[loop].set(0)
                    except Exception:
                        plotgui.legend_frame.append(None)
                else:
                    plotgui.legend_frame.append(None)
            else:
                if loop < len(plotgui.legend_frame):
                    # set the legend_frame value, if possible
                    try:
                        plotgui.legend_frame[loop].set(
                            legend[loop]['legend_frame_value'])
                    except Exception:
                        plotgui.legend_frame[loop] = Tk.IntVar()
                        plotgui.legend_frame[loop].set(
                            legend[loop]['legend_frame_value'])
                else:
                    plotgui.legend_frame.append(Tk.IntVar())
                    plotgui.legend_frame[loop].set(
                        legend[loop]['legend_frame_value'])
            if legend[loop]['legend_option_value'] is None:
                if loop < len(plotgui.legend_options):
                    # get the legend_option, if possible
                    try:
                        plotgui.legend_options[loop].get()
                    except Exception:
                        plotgui.legend_options.append(None)
                else:
                    plotgui.legend_options.append(None)
            else:
                if loop < len(plotgui.legend_options):
                    try:
                        # set the legend_option, if possible, or create
                        # the variable
                        plotgui.legend_options[loop].set(
                            legend[loop]['legend_option_value'])
                    except Exception:
                        plotgui.legend_options[loop] = Tk.StringVar()
                        plotgui.legend_options[loop].set(
                            legend[loop]['legend_option_value'])
                else:
                    plotgui.legend_options.append(Tk.StringVar())
                    plotgui.legend_options[loop].set(
                        legend[loop]['legend_option_value'])
        plotgui.share_axis = deepcopy(share_axis)
        if equal_aspect is None:
            plotgui.equal_aspect = []
            for loop in range(nplots):
                plotgui.equal_aspect.append(False)
        else:
            plotgui.equal_aspect = deepcopy(equal_aspect)
        plotgui.number_of_plots = 1*nplots
        plotgui.hide_subplot = deepcopy(hide_subplot)
    return goodfile
