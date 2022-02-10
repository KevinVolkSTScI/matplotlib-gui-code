#! /usr/bin/env python
#
"""
This code uses matplotlib and numpy to produce a window within which a FITS 
image can be displayed.  The reason for having this and not using the usual 
packages already in existence is that I will want specific functions on the 
image for data reduction.

Usage:

fits_image_display.py imagename.fits

or just

fits_image_display.py

In the first case the image name given is loaded (if possible) and displayed.

In the second case the widget comes up and one can read in an image.

Note that if the image is of dimension larger than 2 then the first "plane" 
is used.  There is no mechanism here for using other planes.

"""
import math
import sys
import numpy
import astropy.io.fits as fits
from astropy.modeling import models, fitting
import tkinter as Tk
import tkinter.ttk
import tkinter.filedialog
import tkinter.simpledialog
import tkinter.messagebox
from tkinter.scrolledtext import ScrolledText
import matplotlib
import matplotlib.lines as mlines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from photutils import aperture
import general_utilities

class ImageGUI(Tk.Frame):
    """
    This class brings up a separate image display window.

    Parameters
    ----------

    Tk.Frame:   The base class of the object, matching a Tkinter root or
                Toplevel variable

    Returns
    -------

    The class variable is returned, effectively.
    """
    # The following section of code concerns the image display functionality.
    #

    def __init__(self, parent=None, **args):
        self.image = None
        self.imagefilename = None
        self.namestring = None
        self.zscale_flag = False
        self.root = None
        self.indpi = 100
        self.zoom = [1, 0, 0]
        self.segment = 10
        self.cross = 2
        self.xposition = None
        self.yposition = None
        self.catalogue = None
        self.fit_image = None
        if parent is not None:
            # initialize the window and make the plot area.
            Tk.Frame.__init__(self, parent, args)
            self.root = parent

    def show_plot(self):
        """
        This is a wrapper routine to show the plot.

        Parameters
        ----------

           None

        Returns
        -------

           None
        """
        make_plot.make_plot(self)
            
    def make_image_window(self):
        """
        Make the main image display window.

        Returns
        -------
        None.

        """
        # make the window
        BGCOL = '#F8F8FF'
        if self.root is not None:
            imagewindow = self.root
        else:
            imagewindow = Tk.Toplevel()
        imagewindow.config(bg=BGCOL)
        self.showImageAxes = True
        imageLabelFrame = Tk.Frame(imagewindow)
        imageLabelFrame.pack(side=Tk.TOP)
        self.imagePosLabelText = Tk.StringVar()
        self.imagePosLabel = Tk.Label(imageLabelFrame,
                                      textvariable=self.imagePosLabelText,
                                      anchor=Tk.N, width=70)
        self.imagePosLabel.pack(side=Tk.LEFT)
        self.imagePosLabelText.set("Position:  Value:")
        controlFrame = Tk.Frame(imagewindow)
        controlFrame.pack(side=Tk.LEFT, fill=Tk.Y, expand=1)
        self.plotFrame = Tk.Frame(imagewindow)
        self.plotFrame.pack()
        self.mplfig1 = Figure(figsize=(6, 6), dpi=self.indpi)
        self.mplsubplot1 = self.mplfig1.add_subplot(1, 1, 1)
        self.canvas1 = FigureCanvasTkAgg(self.mplfig1, master=self.plotFrame)
        self.canvas1.draw()
        self.canvas1.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH,
                                          expand=Tk.YES)
        self.canvas1.mpl_connect("motion_notify_event", self.setPlotPosition)
        self.canvas1.mpl_connect("button_press_event", self.buttonPress)
        self.canvas1.mpl_connect("button_release_event", self.buttonRelease)
        self.canvas1.mpl_connect("key_press_event", self.keyPress)
        newframe = Tk.Frame(controlFrame)
        newframe.pack(side=Tk.TOP)
        lb = Tk.Label(newframe, text='Colour Scheme')
        lb.pack(side=Tk.TOP)
        self.colourScheme = tkinter.ttk.Combobox(newframe, width=15)
        self.colourLabels = ['jet', 'rainbow', 'gist_ncar', 'viridis',
                             'gnuplot', 'gist_gray', 'nipy_spectral']
        self.colourScheme['values'] = self.colourLabels
        self.colourScheme.pack()
        self.colourScheme.current(0)
        #
        lb = Tk.Label(newframe, text='Show Colour Bar')
        lb.pack()
        selectFrame = Tk.Frame(newframe)
        selectFrame.pack()
        self.colourBar = Tk.IntVar()
        t1 = Tk.Radiobutton(selectFrame, text='vertical',
                            variable=self.colourBar, value=0,
                            command=self.displayImage)
        t1.pack(side=Tk.LEFT)
        t2 = Tk.Radiobutton(selectFrame, text='horizontal',
                            variable=self.colourBar, value=1,
                            command=self.displayImage)
        t2.pack(side=Tk.LEFT)
        t3 = Tk.Radiobutton(selectFrame, text='none', variable=self.colourBar,
                            value=2, command=self.displayImage)
        t3.pack(side=Tk.LEFT)
        self.colourBar.set(2)
        lb = Tk.Label(newframe, text='Colour Bar Label')
        lb.pack()
        self.barLabel = Tk.Entry(newframe, width=30)
        self.barLabel.pack()
        rangeframe = Tk.Frame(newframe)
        rangeframe.pack()
        fr1 = Tk.Frame(rangeframe)
        fr1.pack(side=Tk.LEFT)
        lb = Tk.Label(fr1, text='Display Minimum')
        lb.pack(side=Tk.TOP)
        self.minField = Tk.Entry(fr1, width=10)
        self.minField.pack()
        fr1 = Tk.Frame(rangeframe)
        fr1.pack(side=Tk.LEFT)
        Tk.Label(fr1, text=' ').pack()
        fr1 = Tk.Frame(rangeframe)
        fr1.pack(side=Tk.LEFT)
        lb = Tk.Label(fr1, text='Display Maximum')
        lb.pack(side=Tk.TOP)
        self.maxField = Tk.Entry(fr1, width=10)
        self.maxField.pack()
        zmin = numpy.min(self.image)
        zmax = numpy.max(self.image)
        general_utilities.put_value(zmin, self.minField)
        general_utilities.put_value(zmax, self.maxField)
        rangeframe = Tk.Frame(newframe)
        rangeframe.pack()
        fr1 = Tk.Frame(rangeframe)
        fr1.pack(side=Tk.LEFT)
        lb = Tk.Label(fr1, text='Zscale Minimum')
        lb.pack(side=Tk.TOP)
        self.zsminField = Tk.Entry(fr1, width=10)
        self.zsminField.pack()
        fr1 = Tk.Frame(rangeframe)
        fr1.pack(side=Tk.LEFT)
        Tk.Label(fr1, text=' ').pack()
        fr1 = Tk.Frame(rangeframe)
        fr1.pack(side=Tk.LEFT)
        lb = Tk.Label(fr1, text='Zscale Maximum')
        lb.pack(side=Tk.TOP)
        self.zsmaxField = Tk.Entry(fr1, width=10)
        self.zsmaxField.pack()
        try:
            zmin1, zmax1 = self.get_limits(self.image)
        except:
            zmin1 = 0.
            zmax1 = 1.
        general_utilities.put_value(zmin1, self.zsminField)
        general_utilities.put_value(zmax1, self.zsmaxField)
        lb = Tk.Label(newframe, text='Image Scaling')
        lb.pack()
        selectFrame = Tk.Frame(newframe)
        selectFrame.pack()
        self.scaleType = Tk.IntVar()
        t1 = Tk.Radiobutton(selectFrame, text='linear',
                            variable=self.scaleType, value=0,
                            command=self.displayImage)
        t1.pack(side=Tk.LEFT)
        t2 = Tk.Radiobutton(selectFrame, text='log', variable=self.scaleType,
                            value=1, command=self.displayImage)
        t2.pack(side=Tk.LEFT)
        t3 = Tk.Radiobutton(selectFrame, text='sqrt',
                            variable=self.scaleType, value=2,
                            command=self.displayImage)
        t3.pack(side=Tk.LEFT)
        self.scaleType.set(0)
        lb = Tk.Label(newframe, text='Image Range')
        lb.pack()
        selectFrame = Tk.Frame(newframe)
        selectFrame.pack()
        self.rangeType = Tk.IntVar()
        t1 = Tk.Radiobutton(
            selectFrame, text='full', variable=self.rangeType,
            value=0, command=self.toggle_zscale)
        t1.pack(side=Tk.LEFT)
        t2 = Tk.Radiobutton(
            selectFrame, text='zscale', variable=self.rangeType,
            value=1, command=self.toggle_zscale)
        t2.pack(side=Tk.LEFT)
        self.rangeType.set(0)
        buttonFrame = Tk.Frame(controlFrame)
        buttonFrame.pack(side=Tk.TOP)
        subFrame = Tk.Frame(buttonFrame)
        subFrame.pack(side=Tk.TOP)
        side1 = Tk.Frame(subFrame)
        side1.pack(side=Tk.LEFT)
        b1 = Tk.Button(side1, text='Toggle Axes',
                       command=self.toggleAxes)
        b1.pack(side=Tk.TOP)
        b1 = Tk.Button(side1, text='Auto Scale',
                       command=self.imageAutoscale)
        b1.pack(side=Tk.TOP)
        side2 = Tk.Frame(subFrame)
        side2.pack(side=Tk.LEFT)
        b1 = Tk.Button(side2, text='Image Histogram',
                       command=self.imageHistogram)
        b1.pack(side=Tk.TOP)
        b1 = Tk.Button(side2, text='Set Zoom',
                       command=self.set_zoom)
        b1.pack(side=Tk.TOP)
        bin_frame = Tk.Frame(buttonFrame)
        bin_frame.pack(side=Tk.TOP)
        label = Tk.Label(bin_frame, text='bin size/number')
        label.grid(row=0, column=0)
        self.bin_field = Tk.Entry(bin_frame, width=10)
        self.bin_field.grid(row=0, column=1)
        self.bin_field.insert(0, '100')
        label = Tk.Label(
            bin_frame, text='Positive for bin number, negative for \nbin size')
        label.grid(row=1, column=0, columnspan=2)
        label = Tk.Label(buttonFrame, text='Histogram y scaling:')
        label.pack()
        yscaleFrame = Tk.Frame(buttonFrame)
        yscaleFrame.pack(side=Tk.TOP)
        self.yscaleType = Tk.IntVar()
        t1 = Tk.Radiobutton(
            yscaleFrame, text='linear', variable=self.yscaleType,
            value=0)
        t1.pack(side=Tk.LEFT)
        t2 = Tk.Radiobutton(
            yscaleFrame, text='hybrid log', variable=self.yscaleType,
            value=1)
        t2.pack(side=Tk.LEFT)
        self.rangeType.set(0)
        b1 = Tk.Button(buttonFrame, text='Read New Image',
                       command=self.readNewImage)
        b1.pack(side=Tk.TOP)
        b1 = Tk.Button(buttonFrame, text='Save as PNG',
                       command=lambda: general_utilities.save_png_figure(
                           self.mplfig1))
        b1.pack(side=Tk.TOP)
        b1 = Tk.Button(buttonFrame, text='Save as PS',
                       command=lambda: general_utilities.save_ps_figure(
                           self.mplfig1))
        b1.pack(side=Tk.TOP)
        b1 = Tk.Button(buttonFrame, text='Redisplay',
                       command=self.displayImage)
        b1.pack(side=Tk.TOP)
        b1 = Tk.Button(buttonFrame, text='Close',
                       command=lambda: self.imageExit(imagewindow))
        b1.pack(side=Tk.TOP)
        self.displayImage()

    def zoom_corner(self, sh1, zoom, x1, y1):
        """
        Given the zoom parameters find the array lower left corner.

        Parameters
        ----------

        sh1:  A two-element list of the shape of the input image, values being
              integers

        zoom:  A positive integer zoom function to be applied to the image

        x1:    The x pixel value for the centre of the field to display
               (float or integer)

        y1:    The y pixel value for the centre of the field to display
               (float or integer)

        Returns
        -------

        xmin:  An integer value for the lower left corner x pixel index

        ymin:  An integer value for the lower left corner y pixel index


        """
        nxpixel = sh1[1] // zoom
        nypixel = sh1[0] // zoom
        xmin = x1 - nxpixel/2.
        ymin = y1 - nypixel/2.
        xmin = int(xmin)
        ymin = int(ymin)
        if xmin < 0:
            xmin = 0
        if ymin < 0:
            ymin = 0
        xmax = xmin + nxpixel
        ymax = ymin + nypixel
        if ymax > sh1[0]:
            ymax = sh1[0]
            ymin = ymax - nypixel
        if xmax > sh1[1]:
            xmax = sh1[1]
            xmin = xmax - nxpixel
        return xmin, ymin

    def set_zoom(self):
        """
        Bring up a window to set the zoom parameter.

        No values are passed to this routine or returned from it.  The
        self.zoom variable is changed by the routine.
        """
        sh1 = self.image.shape
        npixel = min(sh1[0], sh1[1])
        zoommax = int(npixel/64.)
        if zoommax <= 1:
            tkinter.messagebox.showinfo(
                "Error",
                "Zoom is disabled for minimum image size < 128 pixels.")
            return
        if self.xposition is None:
            x1 = sh1[1]/2.
            y1 = sh1[0]/2.
        else:
            x1 = self.xposition
            y1 = self.yposition
        zoom = tkinter.simpledialog.askinteger(
            'Input',
            'Set the integer zoom value (1 to %d)' % (zoommax))
        if zoom is None:
            return
        else:
            xmin, ymin = self.zoom_corner(sh1, zoom, x1, y1)
            self.zoom[0] = zoom
            self.zoom[1] = int(xmin)
            self.zoom[2] = int(ymin)
            self.displayImage()

    def toggle_zscale(self):
        """
        Toggle the zscale option in the image display

        This routine is called in response to the "Image Range" radio button.
        It turns the zscale display option on or off via the self.zscale_flag
        boolean variable.

        No values are passed to this routine or returned form the routine.
        """
        ind = self.rangeType.get()
        if ind == 1:
            self.zscale_flag = True
        else:
            self.zscale_flag = False
        self.displayImage()

    def readNewImage(self):
        """
        Read a FITS image from a file and display it.

        Routine to read a FITS files and extract a two-dimensional image if
        possible.  The image is then displayed.  This routine will only work
        if the image display window exists.

        No parameters are passed to this routine or returned from this routine.
        """
        try:
            filename = tkinter.filedialog.askopenfilename(
                filetypes=[('FITS', '*.fits')])
            if filename is not None:
                self.imagefilename = filename
                values = filename.split('/')
                self.namestring = values[-1]
                values = filename.split('/')
                self.namestring = values[-1]
                self.image = self.get_image()
                if self.image is None:
                    self.imagefilename = None
                    return
                sh1 = self.image.shape
                self.xposition = sh1[1] // 2
                self.yposition = sh1[0] // 2
                print('centre position: ', self.xposition, self.yposition)
                self.displayImage()
                self.canvas1.draw()
                values = filename.split('/')
                self.namestring = values[-1]
        except Exception:
            self.namestring = ' '

    def get_limits(self, values, nsamples=1000, contrast=0.25, max_reject=0.5,
                   min_npixels=5, krej=2.5, max_iterations=5):
        """
        Find the IRAF-like "zscale" signal limits for an image.

        This routine is copied from astropy.visualization.

        Aside from a change to the passing of the arguments the code has
        not been changed.  The original code is part of ZScaleInterval.
        It is a recoding of the IRAF zscale algorithm in python.

        All parameters except the input image array are optional.

        Parameters
        ----------
        values :   a two-dimensional numpy array for which the zscale limit
                   values are to be calculated.  Can be float or integer values.

        nsamples : the number of pixels to use to estimate the median and the
                   range (integer).

        contrast : The constrast parameter from IRAF imexam which controls the
                   range of values considered to estimate the minimum and
                   maximum values to use in the display, a real value between
                   0.0 and 1.0.

        max_reject : Parameter for the maximum fraction of rejected pixels,
                     a real values between 0.0 and 1.0; if more than this
                     fraction of pixels are rejected then the full range
                     of the data values is returned.

        min_npixels : An integer value for the minimum number of pixels that
                      are rejected by the iterative algorithm; if less than
                      this number of pixels is rejected the full data range is
                      returned.

        krej :  A float value, The number of standard deviations used for
                rejection.  It must be positive.

        max_iterations : An integer value giving the maximum number of
                         rejection iterations to use.

        Returns
        -------
        vmin :  the minimum value for the zscale range, a real number

        vmax :  the maximum value for the zscale range, a real number

        """
        # Sample the image
        values = numpy.asarray(values)
        values = values[numpy.isfinite(values)]
        stride = int(max(1.0, values.size / nsamples))
        samples = values[::stride][:nsamples]
        samples.sort()

        npix = len(samples)
        vmin = samples[0]
        vmax = samples[-1]

        # Fit a line to the sorted array of samples
        minpix = max(min_npixels, int(npix * max_reject))
        xvalues = numpy.arange(npix)
        ngoodpix = npix
        last_ngoodpix = npix + 1

        # Bad pixels mask used in k-sigma clipping
        badpix = numpy.zeros(npix, dtype=bool)

        # Kernel used to dilate the bad pixels mask
        ngrow = max(1, int(npix * 0.01))
        kernel = numpy.ones(ngrow, dtype=bool)

        for niter in range(max_iterations):
            if ngoodpix >= last_ngoodpix or ngoodpix < minpix:
                break

            fit = numpy.polyfit(xvalues, samples, deg=1,
                                w=(~badpix).astype(int))
            fitted = numpy.poly1d(fit)(xvalues)

            # Subtract fitted line from the data array
            flat = samples - fitted

            # Compute the k-sigma rejection threshold
            threshold = krej * flat[~badpix].std()

            # Detect and reject pixels further than k*sigma from the
            # fitted line
            badpix[(flat < - threshold) | (flat > threshold)] = True

            # Convolve with a kernel of length ngrow
            badpix = numpy.convolve(badpix, kernel, mode='same')

            last_ngoodpix = ngoodpix
            ngoodpix = numpy.sum(~badpix)

        slope, intercept = fit

        if ngoodpix >= minpix:
            if contrast > 0:
                slope = slope / contrast
            center_pixel = (npix - 1) // 2
            median = numpy.median(samples)
            vmin = max(vmin, median - (center_pixel - 1) * slope)
            vmax = min(vmax, median + (npix - center_pixel) * slope)

        return vmin, vmax

    def get_image(self):
        """
        Read a FITS image from the 0th or 1st extension.

        This routine tries to read a FITS file and returns the image, or None
        if there is an issue:

        Parameters
        ----------
            None

        Returns
        -------
            image :    a numpy two-dimensional array of image values, or None
                       if there is an issue.

        """
        try:
            image = fits.getdata(self.imagefilename)
        except (ValueError, IndexError):
            image = fits.getdata(self.imagefilename, ext=1)
        sh1 = image.shape
        if len(sh1) < 2:
            print('Bad image dimensions in file %s.' %
                  (self.imagefilename))
            return None
        if len(sh1) == 3:
            image = numpy.squeeze(image[0, :, :])
        if len(sh1) == 4:
            image = numpy.squeeze(image[0, 0, :, :])
        if len(sh1) == 5:
            image = numpy.squeeze(image[0, 0, 0, :, :])
        if len(sh1) == 6:
            image = numpy.squeeze(image[0, 0, 0, 0, :, :])
        zmin = numpy.min(image)
        zmax = numpy.max(image)
        general_utilities.put_value(zmin, self.minField)
        general_utilities.put_value(zmax, self.maxField)
        return image

    def imageHistogram(self):
        """
        Plot an IRAF-like image histogram for the current image.

        This routine plots a histogram of the image pixel values in
        a new window.  No values are passed to this routine or returned from
        this routine.
        """
        if self.image is None:
            return
        BGCOL = '#F8F8FF'
        try:
            histogramwindow = Tk.Toplevel()
            histogramwindow.config(bg=BGCOL)
            if self.zscale_flag:
                xmin = float(self.zsminField.get())
                xmax = float(self.zsmaxField.get())
            else:
                xmin = float(self.minField.get())
                xmax = float(self.maxField.get())
            yscale_option = self.yscaleType.get()
            try:
                value = float(self.bin_field.get())
                if value == 0:
                    nbins = 100
                if value < 0.:
                    xstep = abs(value)
                    xmin = xmin - xstep
                    xmax = xmax + 2.0*xstep
                    nbins = int((xmax - xmin)/xstep)
                    xmax = xmin + nbins*xstep
                else:
                    nbins = int(value)
                    nbins = max(nbins, 10)
            except ValueError:
                nbins = 100
            xstep = (xmax - xmin)/nbins
            xmin = xmin - xstep
            xmax = xmax + 2.0*xstep
            nbins = int((xmax - xmin)/xstep)
            xmax = xmin + nbins*xstep
            self.imageHistogramLabelText = Tk.StringVar()
            self.imageHistogramLabel = Tk.Label(
                histogramwindow, textvariable=self.imageHistogramLabelText,
                anchor=Tk.N, width=70)
            self.imageHistogramLabel.pack()
            self.imageHistogramLabelText.set("Value:")
            self.p3 = Figure(figsize=(6, 6), dpi=100)
            sp1 = self.p3.add_subplot(1, 1, 1)
            c1 = FigureCanvasTkAgg(self.p3, master=histogramwindow)
            c1.mpl_connect("motion_notify_event", self.imageHistogramPosition)
            histogramy, hxedges = numpy.histogram(
                self.image.flatten(), nbins, range=[xmin, xmax])
            histogramx = (hxedges[1:]+hxedges[0:-1])/2.
            if yscale_option == 1:
                newyvalues = general_utilities.hybrid_transform(histogramy)
                sp1.plot(histogramx, newyvalues, color='blue')
            else:
                sp1.plot(histogramx, histogramy, color='blue')
            sp1.set_xlabel('Signal')
            sp1.set_ylabel('Number of points per bin')
            if yscale_option == 1:
                tickmarks, ticklabels = general_utilities.hybrid_labels(
                    newyvalues)
                sp1.set_yticks(tickmarks)
                sp1.set_yticklabels(ticklabels)
            label = 'Bin size: %.5g\nNumber of Bins: %d' % (xstep, nbins)
            xpos = xmin + 0.01*(xmax - xmin)
            ymin, ymax = sp1.get_ybound()
            ypos = ymax + (ymax - ymin)*0.02
            if self.imagefilename is None:
                outstring = None
            else:
                outstring = '# Histogram from file ' + self.imagefilename
            sp1.text(xpos, ypos, label)
            c1.draw()
            c1.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES)
            h1 = Tk.Frame(histogramwindow)
            h1.pack(side=Tk.TOP)
            h1.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save values",
                command=lambda: general_utilities.save_data_set_values(
                    histogramx, histogramy, outstring))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save as PS",
                command=lambda: general_utilities.save_ps_figure(self.p3))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save as PNG",
                command=lambda: general_utilities.save_png_figure(self.p3))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(h1, text="Close",
                               command=histogramwindow.destroy)
            button.pack()
            button.config(bg=BGCOL)
        except Exception:
            pass

    def imageHistogramPosition(self, event):
        """
        Post mouse position on histogram plot to the status line.

        When a normal histogram plot exists, this routine takes the mouse
        position events and updates the position values at the top of the
        window.

        Parameters
        ----------
            event     a standard Tkinter event variable.

        Returns
        -------
            No values are returned by this routine.

        """
        try:
            xpos = float(event.xdata)
            ypos = float(event.ydata)
            if self.yscaleType.get() == 1:
                ypos = general_utilities.inverse_hybrid_transform(ypos)
            s1 = 'Value: [%g, %g]' % (xpos, ypos)
            self.imageHistogramLabelText.set(s1)
        except Exception:
            pass

    def plotPosition(self, event):
        """
        Post mouse position on (x, y) plot to the status line.

        When a normal plot exists, this routine takes the mouse
        position events and updates the position values at the top of the
        window.

        Parameters
        ----------
            event     a standard Tkinter event variable.

        Returns
        -------
            No values are returned by this routine.

        """
        try:
            xpos = float(event.xdata)
            ypos = float(event.ydata)
            s1 = 'Value: [%g, %g]' % (xpos, ypos)
            self.plotLabelText.set(s1)
        except Exception:
            pass

    def profilePosition(self, event):
        """
        Post mouse position on the profile plot to the status line.

        When a peofile plot exists, this routine takes the mouse
        position events and updates the position values at the top of the
        window.

        Parameters
        ----------
            event     a standard Tkinter event variable.

        Returns
        -------
            No values are returned by this routine.

        """
        try:
            xpos = float(event.xdata)
            ypos = float(event.ydata)
            s1 = 'Value: [%g, %g]' % (xpos, ypos)
            self.profileLabelText.set(s1)
        except Exception:
            pass

    def put_value(self, value, field):
        """
        Place a value in a widgit text field.

        Any current contents of the field are deleted.

        Parameters
        ----------
            value :  the string value to be placed in the text field

            field :  the tkinter text field variable where the string is to
                     be put

        No values are returned from this routine.

        """
        try:
            s1 = field.get()
            field.delete(0, last=len(s1))
            field.insert(0, str(value))
        except Exception:
            pass

    def toggleAxes(self):
        """
        Toggle the axis display variable.

        Each call to this routine toggles the logical variable determining
        whether the axes are plotted with the image.  No values are passed
        to this routine or returned from it.
        """
        self.showImageAxes = not self.showImageAxes
        self.displayImage()

    def imageAutoscale(self):
        """
        Autoscale the image display.

        This routine resets the minimum and maximum image display values to
        the full range of the current image.

        No values are passed to this routine or returned from this routine.
        """
        zmin = numpy.min(self.image)
        zmax = numpy.max(self.image)
        general_utilities.put_value(zmin, self.minField)
        general_utilities.put_value(zmax, self.maxField)
        zmin1, zmax1 = self.get_limits(self.image)
        general_utilities.put_value(zmin1, self.zsminField)
        general_utilities.put_value(zmax1, self.zsmaxField)
        self.displayImage()

    def imageExit(self, window):
        """
        Close a Tkinter window.

        This routine closes the window for the image display (or
        whichever top level window variable is passed into the routine).

        Parameters
        ----------
            window :  A tkinter Toplevel variable (or equivalent), the window
                      to be closed.

        No values are returned by this routine.

        """
        window.destroy()

    def keyPress(self, event):
        """
        Routine for applying imaging key press events.

        Holder routine for key press events in the image window.  Sets the
        image position.
        """
        if (event.xdata is None) or (event.ydata is None):
            return
        xpixel = int(self.zoom[1]+event.xdata+0.5)
        ypixel = int(self.zoom[2]+event.ydata+0.5)
        if (xpixel is None) or (ypixel is None):
            return
        self.xposition = self.zoom[1]+event.xdata
        self.yposition = self.zoom[2]+event.ydata
        sh1 = self.image.shape
        if event.key == 'l':
            yvalues = numpy.squeeze(self.image[ypixel, :])
            xvalues = numpy.arange(sh1[0])+1
            self.plotxy(xvalues, yvalues, colour='blue', symb=None,
                        xlabel='Column (Pixels)', ylabel='Pixel Value',
                        title='Line %d' % (ypixel))
            return
        if event.key == 'c':
            yvalues = numpy.squeeze(self.image[:, xpixel])
            xvalues = numpy.arange(sh1[1])+1
            self.plotxy(xvalues, yvalues, colour='blue', symb=None,
                        xlabel='Line (Pixels)', ylabel='Pixel Value',
                        title='Column %d' % (xpixel))
            return
        if event.key == 'j':
            xstart = xpixel - self.segment
            if xstart < 0:
                xstart = 0
            xend = xstart + self.segment+self.segment + 2
            if xend > sh1[1]:
                xend = sh1[1]
                xstart = xend - self.segment-self.segment - 2
            ystart = ypixel - self.cross
            if ystart < 0:
                ystart = 0
            yend = ystart + self.cross + self.cross + 2
            if yend > sh1[0]:
                yend = sh1[0]
                ystart = yend - self.cross - self.cross - 2
            subim = numpy.copy(self.image[ystart:yend, xstart:xend])
            yvalues = numpy.mean(subim, axis=0)
            xvalues = numpy.arange(len(yvalues))+xstart
            tstring = 'Mean of columns (y): %d:%d' % (ystart, yend)
            self.plotxy(xvalues, yvalues, symbol=None, colour='blue',
                        xlabel='x pixel position', ylabel='Mean Signal',
                        title=tstring)
            return
        if event.key == 'k':
            ystart = ypixel-self.segment
            if ystart < 0:
                ystart = 0
            yend = ystart + self.segment + self.segment + 2
            if yend > sh1[0]:
                yend = sh1[0]
                ystart = yend - self.segment - self.segment - 2
            xstart = xpixel-2
            if xstart < 0:
                xstart = 0
            xend = xstart + self.cross + self.cross + 2
            if xend >= sh1[1]:
                xend = sh1[1]
                xstart = xend - self.cross - self.cross - 2
            subim = numpy.copy(self.image[ystart:yend, xstart:xend])
            yvalues = numpy.mean(subim, axis=1)
            xvalues = numpy.arange(len(yvalues))+ystart
            tstring = 'Mean of rows (x) %d:%d' % (ystart, yend)
            self.plotxy(xvalues, yvalues, symbol='-', colour='blue',
                        xlabel='y pixel position', ylabel='Signal (ADU/s)',
                        title=tstring)
            return
        if event.key == 'r':
            tstring = 'Radial profile at (%.3f %.3f)' % (
                event.xdata+self.zoom[1], event.ydata+self.zoom[2])
            self.plot_radial_profile(event.xdata+self.zoom[1], 
                                     event.ydata+self.zoom[2],
                                     xlabel='Radius (pixels)',
                                     ylabel='Signal', title=tstring)
            return
        if event.key == 's':
            self.plot_surface_fit()
            return
        if event.key == 'h':
            self.show_key_help()
            return
        # all other keys move the zoom window to be centred on the position
        xmin, ymin = self.zoom_corner(sh1, self.zoom[0], self.xposition,
                                      self.yposition)
        self.zoom[1] = xmin
        self.zoom[2] = ymin
        self.displayImage()
        return

    def buttonPress(self, event):
        """
        Routine for applying imaging button press events.

        Holder routine for button press events in the image window.
        Not currently active.
        """
        return

    def buttonRelease(self, event):
        """
        Routine for applying imaging button release events.

        Holder routine for button release events in the image window.

        """
        if (event.xdata is None) or (event.ydata is None):
            return
        sh1 = self.image.shape
        xpixel = int(self.zoom[1]+event.xdata+0.5)
        ypixel = int(self.zoom[2]+event.ydata+0.5)
        if (xpixel is None) or (ypixel is None):
            return
        self.xposition = self.zoom[1]+event.xdata
        self.yposition = self.zoom[2]+event.ydata
        xmin, ymin = self.zoom_corner(sh1, self.zoom[0], self.xposition,
                                      self.yposition)
        self.zoom[1] = xmin
        self.zoom[2] = ymin
        self.displayImage()
        return

    def setPlotPosition(self, event):
        """
        Post the image position to the information line on the image display.

        Routine to post the image position and the image value (if possible)
        to the text area above the image display.

        Parameters
        ----------
            event :   a motion-notify event from the image display window

        Returns
        -------
            No values are returned by this routine.

        """
        try:
            event.canvas.get_tk_widget().focus_set()
            x1 = int(self.zoom[1]+event.xdata+0.5)
            y1 = int(self.zoom[2]+event.ydata+0.5)
            try:
                value = '%.6g' % (self.image[y1, x1])
            except ValueError:
                value = ' '
            try:
                s1 = self.namestring+"\n"
            except:
                s1 = ''
            s1 = s1+"Position: x = %.2f y = %.2f Value: %s" % (x1, y1, value)
            self.imagePosLabelText.set(s1)
            self.imagexpos = event.xdata
            self.imageypos = event.ydata
        except Exception:
            pass

    def plotxy(self, xvalues, yvalues, **parameters):
        BGCOL = '#F8F8FF'
        if self.image is None:
            return
        try:
            plotwindow = Tk.Toplevel()
            plotwindow.config(bg=BGCOL)
            self.plotLabelText = Tk.StringVar()
            self.plotLabel = Tk.Label(
                plotwindow, textvariable=self.plotLabelText,
                anchor=Tk.N, width=70)
            self.plotLabel.pack()
            self.plotLabelText.set("Value:")
            self.p4 = Figure(figsize=(6, 6), dpi=100)
            sp1 = self.p4.add_subplot(1, 1, 1)
            c1 = FigureCanvasTkAgg(self.p4, master=plotwindow)
            c1.mpl_connect("motion_notify_event", self.plotPosition)
            symbol = parameters.get('symb')
            if symbol is None:
                symbol='-'
            colour = parameters.get('colour')
            if colour is None:
                colour = 'blue'
            sp1.plot(xvalues, yvalues, symbol, color=colour)
            sp1.set_xlabel(parameters.get('xlabel'))
            sp1.set_ylabel(parameters.get('ylabel'))
            sp1.set_title(parameters.get('title'))
            c1.draw()
            c1.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES)
            h1 = Tk.Frame(plotwindow)
            h1.pack(side=Tk.TOP)
            h1.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save values",
                command=lambda: general_utilities.save_data_set_values(
                    xvalues, yvalues, outstring))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save as PS",
                command=lambda: general_utilities.save_ps_figure(self.p4))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save as PNG",
                command=lambda: general_utilities.save_png_figure(self.p4))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(h1, text="Close",
                               command=plotwindow.destroy)
            button.pack()
            button.config(bg=BGCOL)
        except Exception:
            pass

    def plot_radial_profile(self, xposition, yposition, **parameters):
        BGCOL = '#F8F8FF'
        if self.image is None:
            return
        try:
            xvalues, yvalues, yerror = self.radial_profile(
                self.image, 1.0, 10., centre=[xposition, yposition])
            profilewindow = Tk.Toplevel()
            profilewindow.config(bg=BGCOL)
            self.profileLabelText = Tk.StringVar()
            self.profileLabel = Tk.Label(
                profilewindow, textvariable=self.profileLabelText,
                anchor=Tk.N, width=70)
            self.profileLabel.pack()
            self.profileLabelText.set("Value:")
            self.p5 = Figure(figsize=(6, 6), dpi=100)
            sp1 = self.p5.add_subplot(1, 1, 1)
            c1 = FigureCanvasTkAgg(self.p5, master=profilewindow)
            c1.mpl_connect("motion_notify_event", self.profilePosition)
            symbol = parameters.get('symb')
            if symbol is None:
                symbol='-'
            colour = parameters.get('colour')
            if colour is None:
                colour = 'blue'
            sp1.plot(xvalues, yvalues, symbol, color=colour)
            sp1.set_xlabel(parameters.get('xlabel'))
            sp1.set_ylabel(parameters.get('ylabel'))
            sp1.set_title(parameters.get('title'))
            c1.draw()
            c1.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES)
            h1 = Tk.Frame(profilewindow)
            h1.pack(side=Tk.TOP)
            h1.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save values",
                command=lambda: general_utilities.save_data_set_values(
                    xvalues, yvalues, outstring))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save as PS",
                command=lambda: general_utilities.save_ps_figure(self.p5))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save as PNG",
                command=lambda: general_utilities.save_png_figure(self.p5))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(h1, text="Close",
                               command=profilewindow.destroy)
            button.pack()
            button.config(bg=BGCOL)
        except Exception:
            pass

    def radial_profile(self, array, rstep=0.5, rmax=0.0, centre=None, error=None):
        """
        This routine calculates the radial profile of values for an image.  The
        image is passed as the first argument.  The centre pixel is used for the
        centre point of the calculations unless a value is passed to the routine.

        Parameters
        ----------
        array :   A two-dimensional numpy image (assumed to be float values).

        rstep :   Optional floating point step value in pixels for the encircled
                  energy function.

        rmax :   Optional floating point value, the maximum radius in pixels for
                 the encircled energy function. If no value is given, the
                 distance from the centre position to the nearest image edge
                 is used.
 
        centre : Optional two-element list of the [x, y] values of the centre from
                 which the radius is calculated.  It is assumed to be two float
                 values, x and then y.  If not provided, the image centre is used.

        error :  An optional two-dimensional numpy image (assumed to be float 
                 values) of the uncertainties per pixel; if nothing is passed, 
                 the uncertainties are all set to zero.  The errors must all be 
                 positive.

        Returns
        -------

        radius :  A one-dimensional numpy float array of radius values in pixels.

        signal :  A one-dimensional numpy float array of the signal values

        signal_error : A one-dimensional numpy float array of the uncertainties 
                       in the signal values

        """
        shape = array.shape
        if len(shape) != 2:
            print('Error: input image needs to be two dimensional.')
            return None, None, None
        if error is None:
            uncertainties = array*0.
        else:
            uncertainties = numpy.copy(error)
            if uncertainites.shape != shape:
                print('Error: the uncertainty array is not the same shape as the ')
                print('  signal array, setting to zero.')
                uncertainties = array*0.
        # The following assumes that an integer value corresponds to the pixel
        # centre, as opposed to the IRAF convention that a value of 0.5 denotes
        # the pixel centre
        if centre is None:
            centre = ((shape[1]-1)/2., (shape[0]-1)/2.)
        if rmax <= 0.:
            rmax = max(centre[0], abs(centre[0] - shape[1]),
                       centre[1], abs(centre[1] - shape[0]))
        if rmax <= rstep:
            print('Error: maxumum radius value less than the step value.')
            return None, None, None
        nrad = int(rmax/rstep) + 2
        rout = numpy.zeros((nrad), dtype=numpy.float32)
        signal = numpy.zeros((nrad), dtype=numpy.float32)
        signal_error = numpy.zeros((nrad), dtype=numpy.float32)
        nout = nrad
        for loop in range(nrad):
            rinner = loop * rstep
            router = (loop+1) * rstep
            rout[loop] = (router+rinner)/2.
            if router <= rmax:
                if rinner > 0.:
                    aper = aperture.CircularAnnulus(centre, rinner, router)
                else:
                    aper = aperture.CircularAperture(centre, router)
                values = aperture.aperture_photometry(array, aper,
                                                      method='exact')
                signal[loop] = values['aperture_sum']
                evalues = aperture.aperture_photometry(uncertainties, aper,
                                                       method='exact')
                signal_error[loop] = evalues['aperture_sum']
            else:
                nout = loop
        return rout[0:nout], signal[0:nout], signal_error[0:nout]
        
    def displayImage(self):
        """
        Display the current image in the display area.

        This routine has no parameters and no return value.  It is a
        service routine to display the current image stored in self.image
        into the image display window.
        """
        if self.image is not None:
            self.mplsubplot1.clear()
            s1 = self.minField.get()
            zmin = float(s1)
            s1 = self.maxField.get()
            zmax = float(s1)
            cind = self.colourScheme.current()
            scaleOption = self.scaleType.get()
            try:
                # if the colourBarVariable exists, remove it
                self.colourBarVariable.remove()
            except Exception:
                pass
            startimage = general_utilities.get_subimage(self.image, self.zoom)
            if self.zscale_flag:
                zmin, zmax = self.get_limits(startimage)
                scaleOption = 0
                self.scaleType.set(0)
                general_utilities.put_value(zmin, self.zsminField)
                general_utilities.put_value(zmax, self.zsmaxField)
            else:
                s1 = self.minField.get()
                zmin = float(s1)
                s1 = self.maxField.get()
                zmax = float(s1)
            if (scaleOption == 0) or self.zscale_flag:
                newimage = numpy.copy(startimage)
                im1 = self.mplsubplot1.imshow(
                    newimage, cmap=self.colourLabels[cind],
                    origin='lower', vmin=zmin, vmax=zmax)
            elif scaleOption == 1:
                newimage = self.logTransform(startimage, zmin, zmax)
                zmin1 = numpy.min(newimage)
                zmax1 = numpy.max(newimage)
                im1 = self.mplsubplot1.imshow(
                    newimage, cmap=self.colourLabels[cind], origin='lower',
                    vmin=zmin1, vmax=zmax1)
            else:
                newimage = self.sqrtTransform(startimage, zmin, zmax)
                zmin1 = numpy.min(newimage)
                zmax1 = numpy.max(newimage)
                im1 = self.mplsubplot1.imshow(
                    newimage, cmap=self.colourLabels[cind], origin='lower',
                    vmin=zmin1, vmax=zmax1)
            self.mplsubplot1.get_xaxis().set_visible(self.showImageAxes)
            self.mplsubplot1.get_yaxis().set_visible(self.showImageAxes)
            if self.showImageAxes:
                self.mplsubplot1.set_xlabel('x Pixel Position')
                self.mplsubplot1.set_ylabel('y Pixel Position')
            cbflag = self.colourBar.get()
            cblabel = self.barLabel.get()
            if scaleOption == 1:
                ticklist = numpy.zeros((11), dtype=numpy.float32)
                label1 = numpy.zeros((11), dtype=numpy.float32)
                for lv in range(11):
                    ticklist[lv] = 0.3*lv
                    vout = self.invLogTransform(ticklist[lv], zmin, zmax)
                    label1[lv] = '%.3g' % (vout)
            if scaleOption == 2:
                ticklist = numpy.zeros((11), dtype=numpy.float32)
                label1 = numpy.zeros((11), dtype=numpy.float32)
                for lv in range(11):
                    zrange = numpy.max(newimage) - numpy.min(newimage)
                    ticklist[lv] = numpy.min(newimage) + (zrange * lv / 10.)
                    vout = self.invSqrtTransform(ticklist[lv], zmin, zmax)
                    label1[lv] = '%.3g' % (vout)
            if cbflag == 0:
                if scaleOption > 0:
                    self.colourBarVariable = self.mplfig1.colorbar(
                        im1, cmap=self.colourLabels[cind],
                        orientation='vertical', ticks=ticklist)
                    self.colourBarVariable.ax.set_yticklabels(label1)
                else:
                    self.colourBarVariable = self.mplfig1.colorbar(
                        im1, cmap=self.colourLabels[cind],
                        orientation='vertical')
                self.colourBarVariable.ax.get_yaxis().labelpad = 15
                self.colourBarVariable.ax.set_ylabel(cblabel, rotation=90)
            if cbflag == 1:
                if scaleOption > 0:
                    self.colourBarVariable = self.mplfig1.colorbar(
                        im1, cmap=self.colourLabels[cind],
                        orientation='horizontal', ticks=ticklist)
                    self.colourBarVariable.ax.set_xticklabels(label1)
                else:
                    self.colourBarVariable = self.mplfig1.colorbar(
                        im1, cmap=self.colourLabels[cind],
                        orientation='horizontal')
                self.colourBarVariable.ax.set_xlabel(cblabel, rotation=0)
            if self.catalogue is not None:
                try:
                    starpos = numpy.transpose((self.catalogue['xcentroid'],
                        self.catalogue['ycentroid']))
                except:
                    try:
                        starpos = numpy.transpose((
                            self.catalogue['x_centroid'],
                            self.catalogue['y_centroid']))
                    except:
                        starpos = None
                if not starpos is None:
                    if self.zoom[0] > 1:
                        starpos[:,0] = starpos[:,0]-self.zoom[1]
                        starpos[:,1] = starpos[:,1]-self.zoom[2]
                    rad1 = newimage.shape[0]/50.
                    ap1 = aperture.CircularAperture(starpos, r=rad1)
                    ap1.plot(axes=self.mplsubplot1, color='black', lw=1.0)
                    for loop in range(len(self.catalogue['label'])):
                        radx = rad1
                        rady = rad1
                        if starpos[loop, 0]+radx >= newimage.shape[1]:
                            radx = -radx
                        if starpos[loop, 1]+rady >= newimage.shape[0]:
                            rady = -rady
                        self.mplsubplot1.annotate(
                            self.catalogue['label'][loop],
                            [starpos[loop, 0]+radx, starpos[loop, 1]+rady],
                            color='black', size=10)
            self.canvas1.draw()
#            self.canvas1.get_tk_widget().focus_set()

    def invLogTransform(self, value, zmin, zmax):
        """
        Transform a log value back to the original value in an image.

        This routine returns the original value corresponding to a given
        logarithmic display value in the range from 1 to 1000.

        Parameters
        ----------
           value : a real value, by assumption, in the range from 1 to 1000

            zmin :  a real value, the signal minimum for the logarithmic
                    mapping

            zmax :  a real value, the signal maximum for the logarithmic
                    mapping

        Returns
        -------
            vout   The original image value corresponding to the new image
                   value, a real number

        """
        newvalue = value*1
        if value < 0.:
            newvalue = 0.
        if value > 3.:
            newvalue = 3.
        v1 = math.pow(10., newvalue)
        if v1 < 1.:
            v1 = 1.
        if v1 > 1000.:
            v1 = 1000.
        v2 = (v1 - 1.)/999.9
        vout = zmin + (zmax - zmin) * v2
        return vout

    def logTransform(self, image, zmin, zmax):
        """
        Apply an IRAF-style logarithmic transformation to an image.

        This routine applies an iraf-style logaritmic transform to an
        image; the requested range is mapped to logarithmic values from
        0.0 to 3.0.  The range can be negative since the mapping is for the
        signal values with respect to the defined range.

        Parameters
        ----------
            image :  a numpy array of (assumed) floating point or integer
                     values

            zmin :   a real value, the minimum of the range for the
                     transformation

            zmax :   a real value, the maximum of the range for the
                     transformation

        Returns
        -------
            newimage : a numpy image of the same dimensions as the input
                       image, with the transformation applied; floating
                       point values in the range between 0.0 and 3.0 are
                       contained in the new image
        """
        newimage = numpy.copy(image)
        newimage[newimage < zmin] = zmin
        newimage[newimage > zmax] = zmax
        zrange = zmax - zmin
        newimage = 1. + 999.*(newimage - zmin)/zrange
        newimage = numpy.log10(newimage)
        self.transvalues = [zmin, zmax]
        return newimage

    def invSqrtTransform(self, value, zmin, zmax):
        """
        Transform a sqaure-root scaled image value back to the original value.

        Routine to map the data values from the sqrt transform back to the
        original range of values.

        Parameters
        ----------
            value : a real value, by assumption

            zmin :  a real value, the signal minimum for the sqrt mapping
                    (not used, present for uniformity with the log call)

            zmax :  a real value, the signal maximum for the sqrt mapping
                   (not used, present for uniformity with the log call)

        Returns
        -------
            v1 :  The original image value corresponding to the new image
                  value, a real number

        """
        newvalue = abs(value)
        v1 = newvalue * newvalue
        if value < 0.:
            v1 = -v1
        return v1

    def sqrtTransform(self, image, zmin, zmax):
        """
        Apply a square-root scaling to an image.

        Given an image, this routine applies a square-root scaling of the
        absolute value, preserving the original sign, and returns the
        transformed image.

        Parameters
        ----------
            image :  a numpy array of (assumed) floating point or integer
                     values

            zmin :   a real value, the minimum of the range for the
                     transformation (not currently used, present so the
                     form of the call matches the other transformations)

            zmax  :   a real value, the maximum of the range for the
                      transformation (not currently used)

        Returns
        -------
            newimage  a numpy image of the same dimensions as the input
                       image, with the transformation applied; all values
                       in the image are replaced by the square-root of the
                       absolute value times the original sign
        """
        newimage = numpy.sqrt(numpy.abs(image))
        newimage[image < 0.] = -1. * newimage[image < 0.]
        self.transvalues = [1.]
        return newimage

    def plot_surface_fit(self):
        if self.image is None:
            return
        ypix, xpix = numpy.mgrid[:self.image.shape[0], :self.image.shape[1]]
        nfit = tkinter.simpledialog.askinteger(
            'Input',
            'Set the fitting order in x and y')
        if nfit is None:
            nfit = 4
        if nfit < 1:
            nfit = 4
        polynomial_init = models.Polynomial2D(degree=nfit)
        fit_polynomial = fitting.LevMarLSQFitter()
        fit_image = fit_polynomial(polynomial_init, xpix, ypix, self.image)
        self.fit_image = fit_image(xpix, ypix)
        BGCOL = '#F8F8FF'
        if True:
#        try:
            surfacewindow = Tk.Toplevel()
            surfacewindow.config(bg=BGCOL)
            self.surfaceLabelText = Tk.StringVar()
            self.surfaceLabel = Tk.Label(
                surfacewindow, textvariable=self.surfaceLabelText,
                anchor=Tk.N, width=70)
            self.surfaceLabel.pack()
            self.surfaceLabelText.set("Value:")
            self.p6 = Figure(figsize=(6, 6), dpi=100)
            sp1 = self.p6.add_subplot(1, 1, 1)
            cind = self.colourScheme.current()
            im1 = sp1.imshow(
                self.fit_image, cmap=self.colourLabels[cind],
                origin='lower')
            bar = self.p6.colorbar(
                im1, cmap=self.colourLabels[cind],
                orientation='vertical')
            c1 = FigureCanvasTkAgg(self.p6, master=surfacewindow)
            c1.mpl_connect("motion_notify_event", self.surfacePosition)
            xmin, xmax = sp1.get_xbound()
            xpos = xmin + (xmax - xmin)*0.02
            ymin, ymax = sp1.get_ybound()
            ypos = ymax + (ymax - ymin)*0.02
            outstring = '# Order %d surface fit ' % (nfit)
            if self.imagefilename is not None:
                outstring = outstring+'from file\n'+self.imagefilename
            sp1.text(xpos, ypos, outstring)
            c1.draw()
            c1.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=Tk.YES)
            h1 = Tk.Frame(surfacewindow)
            h1.pack(side=Tk.TOP)
            h1.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save image",
                command=lambda: general_utilities.save_fits_image(self.image))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save as PS",
                command=lambda: general_utilities.save_ps_figure(self.p6))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(
                h1, text="Save as PNG",
                command=lambda: general_utilities.save_png_figure(self.p6))
            button.pack(side=Tk.LEFT)
            button.config(bg=BGCOL)
            button = Tk.Button(h1, text="Close",
                               command=surfacewindow.destroy)
            button.pack()
            button.config(bg=BGCOL)
            residuals = self.image - self.fit_image
            stats1 = general_utilities.image_stats(residuals)
            stats2 = general_utilities.image_stats(self.image)
            str1 = 'Image statistics:\n'
            str1 = str1+'   mean:   '+str(stats2[0])+'\n'
            str1 = str1+'   sigma:  '+str(stats2[1])+'\n'
            str1 = str1+'   median: '+str(stats2[2])+'\n'
            str1 = str1+'   min:    '+str(stats2[3])+'\n'
            str1 = str1+'   max:    '+str(stats2[4])+'\n'
            str1 = str1+'   clipped mean:   '+str(stats2[5])+'\n'
            str1 = str1+'   clipped sigma:  '+str(stats2[6])+'\n'
            str1 = str1+'   clipped median: '+str(stats2[7])+'\n'
            str1 = str1+'   clipped min:    '+str(stats2[8])+'\n'
            str1 = str1+'   clipped max:    '+str(stats2[9])+'\n'
            str1 = str1+'\nResidual image statistics:\n'
            str1 = str1+'   mean:   '+str(stats1[0])+'\n'
            str1 = str1+'   sigma:  '+str(stats1[1])+'\n'
            str1 = str1+'   median: '+str(stats1[2])+'\n'
            str1 = str1+'   min:    '+str(stats1[3])+'\n'
            str1 = str1+'   max:    '+str(stats1[4])+'\n'
            str1 = str1+'   clipped mean:   '+str(stats1[5])+'\n'
            str1 = str1+'   clipped sigma:  '+str(stats1[6])+'\n'
            str1 = str1+'   clipped median: '+str(stats1[7])+'\n'
            str1 = str1+'   clipped min:    '+str(stats1[8])+'\n'
            str1 = str1+'   clipped max:    '+str(stats1[9])+'\n'
            general_utilities.show_text(str1)
#        except Exception:
#            pass
        return

    def show_key_help(self):
        """
        Show the key commands in a scrolled text window.

        Parameters
        ----------

        None

        Returns
        -------

        None

        """
        str1 = """
The following key commands are available when the cursor is in the plot area:

(c)  plot the column at the current position (that is, all y at the current x 
     pixel position) in a new window

(l)  plot the line at the current position (that is, all x at the current y 
     pixel position) in a new window

(j)  take a slice in x of +/- 10 pixels from the current position over a 
     y range of +/-2 pixels, plot the mean profile and attempt to fit this 
     with a Gaussian function plus a baseline; if the fitting works, also 
     plot the fitting function

(k)  take a slice in y of +/- 10 pixels from the current position over a 
     x range of +/-2 pixels, plot the mean profile and attempt to fit this 
     with a Gaussian function plus a baseline; if the fitting works, also 
     plot the fitting function

(r)  calculate the radial profile from the current cursor position and plot 
     the resulting function in a new window

(s)  attempt to make a surface fit of the currently displayed (sub-)image, 
     and if successful show the fit image in a new window; some statistics of 
     the images are also calculated and shown in a scrolled text window

(h)  show the key command help text

Any other key command is used to set the centre position for a zoomed image, 
if any zoom factor is defined.

            """
        general_utilities.show_text(str1, "Key Commands")

    def surfacePosition(self, event):
        """
        Post the image position to the information line on the image display.

        Routine to post the image position and the image value (if possible)
        to the text area above the image display in the surface image window.

        Parameters
        ----------
            event :   a motion-notify event from the image display window

        Returns
        -------
            No values are returned by this routine.

        """
        try:
            event.canvas.get_tk_widget().focus_set()
            x1 = int(self.zoom[1]+event.xdata+0.5)
            y1 = int(self.zoom[2]+event.ydata+0.5)
            try:
                value = '%.6g' % (self.fit_image[y1, x1])
            except ValueError:
                value = ' '
            s1 = "Position: x = %.2f y = %.2f Value: %s" % (x1, y1, value)
            self.surfaceLabelText.set(s1)
        except Exception:
            pass

if __name__ == "__main__":
    # create the window
    root = Tk.Tk()
    root.title('Image Display Widget')
    imdisp = ImageGUI(root)
    imdisp.make_image_window()
    if '.fits' in sys.argv[-1]:
        imdisp.imagefilename = sys.argv[-1]
        imdisp.image = imdisp.get_image()
        imdisp.displayImage()
    root.mainloop()

