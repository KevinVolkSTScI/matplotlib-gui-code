"""
This file contains code for making histogram plots.

Routines:

make_histogram     Code to make a histogram plot in a new window

write_histogram    Write the histogram values out to an ascii file

make_hess_plot     Code to make a two-dimensional "histogram" or density of
                   points plot for a plot (called a Hess plot when used
                   for a colour-magnitude diagram)

makeFits           Write out the two-dimensional histogram values to a
                   FITS image

All these routines need to get data from the matplotlib_user_interface
object ("plotgui").

"""
import tkinter as Tk
import tkinter.ttk
import tkinter.filedialog
import tkinter.simpledialog
import tkinter.messagebox
import numpy
from astropy.io import fits
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
import general_utilities


def make_histogram(plotgui):
    """
    Create a histogram plot in a new window.

    This routine creates a new plot window within which a histogram
    plot is made for the current active main window plot.

    The new window has options for output of the histogram plot to a
    file.  The colours of the bars are the same as the colours of the
    points in the main plot.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    BGCOL = '#F8F8FF'
    try:
        histogramwindow = Tk.Toplevel()
        histogramwindow.config(bg=BGCOL)
        optionflag = plotgui.histogramflag.get()
        setoptionflag = plotgui.individualhistogramflag.get()
        try:
            value = float(plotgui.nbinfield.get())
            if value == 0:
                tkinter.messagebox.showinfo(
                    "Error",
                    "Zero value for the number of bins/bin size."
                    + "  Check your inputs.")
                return
            if value < 1.:
                delx = value
                nbins = 0
            else:
                nbins = int(value+0.001)
                delx = 0.
        except ValueError:
            tkinter.messagebox.showinfo(
                "Error",
                "Unable to read the number of bins/bin size.  "
                + "Check your inputs.")
            return
        if optionflag == 1:
            xmin = plotgui.plot_range[plotgui.current_plot-1][0]
            xmax = plotgui.plot_range[plotgui.current_plot-1][1]
        else:
            xmin = plotgui.plot_range[plotgui.current_plot-1][2]
            xmax = plotgui.plot_range[plotgui.current_plot-1][3]
        xp = None
        for loop in range(plotgui.nsets):
            if (plotgui.set_properties[loop]['display']) and \
               (plotgui.set_properties[loop]['plot'] == plotgui.current_plot):
                mycolour = plotgui.set_properties[loop]['colour']
                if optionflag == 1:
                    values = numpy.copy(plotgui.xdata[loop]['values'])
                else:
                    values = numpy.copy(plotgui.ydata[loop]['values'])
                if xp is None:
                    xp = [values, ]
                    histcolours = [mycolour, ]
                else:
                    if setoptionflag == 1:
                        oldvalues = numpy.copy(xp[0])
                        newvalues = numpy.append(oldvalues, values)
                        xp[0] = newvalues
                    else:
                        xp.append(values)
                        histcolours.append(mycolour)
        plotgui.histogramLabelText = Tk.StringVar()
        plotgui.histogramLabel = Tk.Label(
            histogramwindow,
            textvariable=plotgui.histogramLabelText, anchor=Tk.N, width=70)
        plotgui.histogramLabel.pack()
        plotgui.histogramLabelText.set("Value:")
        plotgui.p2 = Figure(figsize=(6, 6), dpi=100)
        sp1 = plotgui.p2.add_subplot(1, 1, 1)
        c1 = FigureCanvasTkAgg(plotgui.p2, master=histogramwindow)
        c1.mpl_connect("motion_notify_event", plotgui.histogram_position)
        try:
            if delx > 0.:
                npixels = int(abs((xmax-xmin)/delx))
                delx = abs(delx)
            else:
                npixels = nbins
                delx = (xmax-xmin)/npixels
        except ValueError:
            nbins = 500
        histx = []
        histy = []
        for loop in range(len(xp)):
            histogramy, hxedges = numpy.histogram(
                xp[loop], npixels, range=[xmin, xmax])
            histogramx = (hxedges[1:]+hxedges[0:-1])/2.
            sp1.bar(histogramx, histogramy,
                    color=histcolours[loop],
                    width=delx*0.9)
            histx.append(histogramx)
            histy.append(histogramy)
        if optionflag == 1:
            sp1.set_xlabel(plotgui.xparameters[plotgui.current_plot-1]['label'],
                           family=plotgui.fontname[plotgui.current_plot-1],
                           size=plotgui.fontsize[plotgui.current_plot-1],
                           weight=plotgui.fontweight[plotgui.current_plot-1])
        else:
            sp1.set_xlabel(plotgui.yparameters[plotgui.current_plot-1]['label'],
                           family=plotgui.fontname[plotgui.current_plot-1],
                           size=plotgui.fontsize[plotgui.current_plot-1],
                           weight=plotgui.fontweight[plotgui.current_plot-1])
        sp1.set_ylabel('Number of points per bin')
        invertxflag = plotgui.xparameters[plotgui.current_plot-1]['invert']
        invertyflag = plotgui.yparameters[plotgui.current_plot-1]['invert']
        if (invertxflag == 1) and (optionflag == 1):
            sp1.invert_xaxis()
        if (invertyflag == 1) and (optionflag == 0):
            sp1.invert_xaxis()
        c1.draw()
        c1.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES)
        h1 = Tk.Frame(histogramwindow)
        h1.pack(side=Tk.TOP)
        h1.config(bg=BGCOL)
        button = Tk.Button(
            h1, text="Save as PS",
            command=lambda: general_utilities.save_ps_figure(plotgui.p2))
        button.pack(side=Tk.LEFT)
        button.config(bg=BGCOL)
        button = Tk.Button(
            h1, text="Save as PNG",
            command=lambda: general_utilities.save_png_figure(plotgui.p2))
        button.pack(side=Tk.LEFT)
        button.config(bg=BGCOL)
        button = Tk.Button(
            h1, text="Write out values",
            command=lambda: writeHistogram(histx, histy))
        button.pack(side=Tk.LEFT)
        button.config(bg=BGCOL)
        button = Tk.Button(
            h1, text="Close", command=histogramwindow.destroy)
        button.pack()
        button.config(bg=BGCOL)
    except Exception:
        pass

def writeHistogram(xvalues, yvalues):
    """
    Write out the histgram values to a selected output file.

    Parameters
    ----------
        xvalues :  a list of vectors of x values for the histogram, each a
                   numpy float array

        yvalues :  a list of vector of y values for the histogram, each a
                       numpy int array

    Returns
    -------
        None

    """
    outfilename = tkinter.filedialog.asksaveasfilename()
    if isinstance(outfilename, type('string')):
        outfile = open(outfilename, 'w')
        for n1 in range(len(xvalues)):
            for loop in range(len(xvalues[n1])):
                print('%13.6g %10d' % (xvalues[n1][loop],
                                       yvalues[n1][loop]),
                      file=outfile)
            print(' ', file=outfile)
        outfile.close()

def make_hess_plot(plotgui):
    """
    Make a two-dimensional histogram plot.

    This routine creates a new plot window within which a Hess plot (i.e.
    a two-dimensional histogram) is made for the current active plot.

    The new window has options for control of the two-dimensional
    histogram.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    BGCOL = '#F8F8FF'
    try:
        hesswindow = Tk.Toplevel()
        hesswindow.config(bg=BGCOL)
        xmin = plotgui.plot_range[plotgui.current_plot-1][0]
        xmax = plotgui.plot_range[plotgui.current_plot-1][1]
        ymin = plotgui.plot_range[plotgui.current_plot-1][2]
        ymax = plotgui.plot_range[plotgui.current_plot-1][3]
        allflag = plotgui.hessindividualhistogramflag.get()
        overplot = plotgui.overplotflag.get()
        try:
            set_number = int(plotgui.hess_set_field.get())
            if (set_number < 1) or (set_number < plotgui.nsets):
                set_number = 1
        except:
            set_number = 1
        xp = None
        yp = None
        if allflag == 1:
            for loop in range(plotgui.nsets):
                if (plotgui.set_properties[loop]['display']) and \
                   (plotgui.set_properties[loop]['plot'] == plotgui.current_plot):
                    if xp is None:
                        xp = numpy.copy(plotgui.xdata[loop]['values'])
                        yp = numpy.copy(plotgui.ydata[loop]['values'])
                    else:
                        xp = numpy.append(xp, plotgui.xdata[loop]['values'])
                        yp = numpy.append(yp, plotgui.ydata[loop]['values'])
        else:
            xp = numpy.copy(plotgui.xdata[set_number-1]['values'])
            yp = numpy.copy(plotgui.ydata[set_number-1]['values'])
        plotgui.hessLabelText = Tk.StringVar()
        plotgui.hessLabel = Tk.Label(
            hesswindow,
            textvariable=plotgui.hessLabelText, anchor=Tk.N, width=70)
        plotgui.hessLabel.pack()
        plotgui.hessLabelText.set("Value:")
        plotgui.p1 = Figure(figsize=(6, 6), dpi=100)
        sp1 = plotgui.p1.add_subplot(1, 1, 1)
        c1 = FigureCanvasTkAgg(plotgui.p1, master=hesswindow)
        c1.mpl_connect("motion_notify_event", plotgui.hess_position)
        try:
            npixels = int(plotgui.npixelfield.get())
            if (npixels < 50) or (npixels > 5000):
                print('Error: too many or too few pixels requested '
                      + '(%d).  Using the default value of 500.' %
                      (npixels))
                npixels = 500
                plotgui.npixelfield.delete(0, Tk.END)
                plotgui.npixelfield.insert(0, '500')
        except ValueError:
            npixels = 500
        plotgui.hist2d, plotgui.xedges, plotgui.yedges, image = sp1.hist2d(
            xp, yp, npixels, range=[[xmin, xmax], [ymin, ymax]],
            norm=LogNorm())
        if (overplot == 1) and (allflag == 0):
            for loop in range(plotgui.nsets):
                if loop != set_number-1:
                    if (plotgui.set_properties[loop]['display']) and \
                       (plotgui.set_properties[loop]['plot'] == \
                        plotgui.current_plot):
                        if plotgui.set_properties[loop]['symbol'] is None:
                            sp1.plot(
                                plotgui.xdata[loop]['values'],
                                plotgui.ydata[loop]['values'],
                                color=plotgui.set_properties[loop]['colour'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
                        elif plotgui.set_properties[loop]['linestyle'] is None:
                            sp1.plot(
                                plotgui.xdata[loop]['values'],
                                plotgui.ydata[loop]['values'],
                                color=plotgui.set_properties[loop]['colour'],
                                marker=plotgui.set_properties[loop]['symbol'],
                                linestyle='none', markersize=
                                plotgui.set_properties[loop]['symbolsize'])
                        else:
                            sp1.plot(
                                plotgui.xdata[loop]['values'],
                                plotgui.ydata[loop]['values'],
                                color=plotgui.set_properties[loop]['colour'],
                                marker=plotgui.set_properties[loop]['symbol'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                markersize=
                                plotgui.set_properties[loop]['symbolsize'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
        sp1.set_xlabel(plotgui.xparameters[plotgui.current_plot-1]['label'],
                       family=plotgui.fontname[plotgui.current_plot-1],
                       size=plotgui.fontsize[plotgui.current_plot-1],
                       weight=plotgui.fontweight[plotgui.current_plot-1])
        sp1.set_ylabel(plotgui.yparameters[plotgui.current_plot-1]['label'],
                       family=plotgui.fontname[plotgui.current_plot-1],
                       size=plotgui.fontsize[plotgui.current_plot-1],
                       weight=plotgui.fontweight[plotgui.current_plot-1])
        invertxflag = plotgui.xparameters[plotgui.current_plot-1]['invert']
        invertyflag = plotgui.yparameters[plotgui.current_plot-1]['invert']
        hidexflag = plotgui.xparameters[plotgui.current_plot-1]['hide']
        hideyflag = plotgui.yparameters[plotgui.current_plot-1]['hide']
        hidexticksflag = plotgui.xparameters[plotgui.current_plot-1]['hideticks']
        hideyticksflag = plotgui.yparameters[plotgui.current_plot-1]['hideticks']
        hidexlabelsflag = plotgui.xparameters[
            plotgui.current_plot-1]['hidelabels']
        hideylabelsflag = plotgui.yparameters[
            plotgui.current_plot-1]['hidelabels']
        inversexticksflag = plotgui.xparameters[
            plotgui.current_plot-1]['inverseticks']
        inverseyticksflag = plotgui.yparameters[
            plotgui.current_plot-1]['inverseticks']
        bothxticksflag = plotgui.xparameters[plotgui.current_plot-1]['bothticks']
        bothyticksflag = plotgui.yparameters[plotgui.current_plot-1]['bothticks']
        try:
            xminorticks = float(plotgui.xparameters[
                plotgui.current_plot-1]['minorticks'])
        except ValueError:
            xminorticks = 0.0
        try:
            yminorticks = float(plotgui.yparameters[
                plotgui.current_plot-1]['minorticks'])
        except ValueError:
            yminorticks = 0.0
        if invertxflag == 1:
            sp1.invert_xaxis()
        if invertyflag == 1:
            sp1.invert_yaxis()
        if bothyticksflag == 1:
            sp1.tick_params(left=True, right=True, which='both')
        else:
            sp1.tick_params(left=True, right=False, which='both')
        if bothxticksflag == 1:
            sp1.tick_params(bottom=True, top=True, which='both')
        else:
            sp1.tick_params(bottom=True, top=False, which='both')
        if xminorticks > 0:
            sp1.xaxis.set_minor_locator(MultipleLocator(xminorticks))
            xtl = int(plotgui.xparameters[plotgui.current_plot-1]['ticklength']
                      // 2)
            sp1.tick_params(axis='x', which='minor', length=xtl)
        if yminorticks > 0:
            sp1.yaxis.set_minor_locator(MultipleLocator(yminorticks))
            ytl = int(plotgui.yparameters[plotgui.current_plot-1]['ticklength']
                      // 2)
            sp1.tick_params(axis='y', which='minor', length=ytl)
        if plotgui.plot_frame[plotgui.current_plot-1] > 0.:
            sp1.spines['bottom'].set_linewidth(
                plotgui.plot_frame[plotgui.current_plot-1])
            sp1.spines['top'].set_linewidth(
                plotgui.plot_frame[plotgui.current_plot-1])
            sp1.spines['left'].set_linewidth(
                plotgui.plot_frame[plotgui.current_plot-1])
            sp1.spines['right'].set_linewidth(
                plotgui.plot_frame[plotgui.current_plot-1])
        else:
            sp1.spines['bottom'].set_linewidth(0.5)
            sp1.spines['top'].set_linewidth(0.5)
            sp1.spines['left'].set_linewidth(0.5)
            sp1.spines['right'].set_linewidth(0.5)
        if inversexticksflag == 1:
            sp1.tick_params(
                axis='x', direction='in',
                length=
                plotgui.xparameters[plotgui.current_plot-1]['ticklength'])
            sp1.tick_params(axis='x', direction='in', which='minor')
        else:
            sp1.tick_params(
                axis='x', direction='out',
                length=plotgui.xparameters[plotgui.current_plot-1]['ticklength'])
            sp1.tick_params(axis='x', direction='out', which='minor')
        if inverseyticksflag == 1:
            sp1.tick_params(
                axis='y', direction='in',
                length=plotgui.yparameters[plotgui.current_plot-1]['ticklength'])
            sp1.tick_params(axis='y', direction='in', which='minor')
        else:
            sp1.tick_params(
                axis='y', direction='out',
                length=plotgui.yparameters[plotgui.current_plot-1]['ticklength'])
            sp1.tick_params(axis='y', direction='out', which='minor')
        if hidexticksflag == 1:
            sp1.get_xaxis().set_ticks([])
        if hideyticksflag == 1:
            sp1.get_yaxis().set_ticks([])
        if hidexlabelsflag == 1:
            sp1.get_xaxis().set_ticklabels([])
        if hideylabelsflag == 1:
            sp1.get_yaxis().set_ticklabels([])
        if hidexflag == 0:
            sp1.get_xaxis().set_visible(True)
        else:
            sp1.get_xaxis().set_visible(False)
        if hideyflag == 0:
            sp1.get_yaxis().set_visible(True)
        else:
            sp1.get_yaxis().set_visible(False)
        c1.draw()
        c1.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES)
        h1 = Tk.Frame(hesswindow)
        h1.pack(side=Tk.TOP)
        h1.config(bg=BGCOL)
        button = Tk.Button(
            h1, text="Save as PS",
            command=lambda: general_utilities.save_ps_figure(plotgui.p1))
        button.pack(side=Tk.LEFT)
        button.config(bg=BGCOL)
        button = Tk.Button(
            h1, text="Save as PNG",
            command=lambda: general_utilities.save_png_figure(plotgui.p1))
        button.pack(side=Tk.LEFT)
        button.config(bg=BGCOL)
        button = Tk.Button(
            h1, text="Save as FITS", command=lambda: makeFits(plotgui))
        button.pack(side=Tk.LEFT)
        button.config(bg=BGCOL)
        button = Tk.Button(h1, text="Close", command=hesswindow.destroy)
        button.pack()
        button.config(bg=BGCOL)
    except Exception:
        pass

def makeFits(plotgui):
    """
    Write out the two-dimensional histogram image as a FITS file.

    This routine makes a FITS output file of the two-dimensional histogram
    values in the plot, for use by other codes.

    Parameters
    ----------

        plotgui:  the matplotlib_user_interface object holding the plot

    Returns
    -------

        None

    """
    filename = tkinter.filedialog.asksaveasfilename(
        filetypes=[('FITS', '*.fits')])
    st1 = filename.split('.')
    if 'fits' not in st1[-1]:
        filename = filename + '.fits'
    newimage = numpy.transpose(plotgui.hist2d)
    hdu = fits.PrimaryHDU(newimage)
    hdulist = fits.HDUList(hdu)
    primary = hdulist[0].header
    primary['CRPIX1'] = (1.0, 'Axis 1 reference pixel')
    primary['CRPIX2'] = (1.0, 'Axis 2 reference pixel')
    x1 = (plotgui.xedges[1]+plotgui.xedges[0])/2.
    delx = (plotgui.xedges[1]-plotgui.xedges[0])
    y1 = (plotgui.yedges[1]+plotgui.yedges[0])/2.
    dely = (plotgui.yedges[1]-plotgui.yedges[0])
    primary['CRVAL1'] = (x1, 'mean value at reference pixel')
    primary['CRVAL2'] = (y1, 'sigma value at reference pixel')
    primary['CDELT1'] = (delx, 'change in mean value per pixel')
    primary['CDELT2'] = (dely, 'change in sigma value per pixel')
    primary['CTYPE1'] = (' ', 'axis 1 type')
    primary['CTYPE2'] = (' ', 'axis 2 type')
    primary['CUNIT1'] = (' ', 'axis 1 unit')
    primary['CUNIT2'] = (' ', 'axis 2 unit')
    hdulist.writeto(filename, overwrite=True)

