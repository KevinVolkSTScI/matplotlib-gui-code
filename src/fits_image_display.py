import math
import sys
import numpy
import astropy.io.fits as fits
import tkinter as Tk
import tkinter.ttk
import tkinter.filedialog
import tkinter.simpledialog
import tkinter.messagebox
import matplotlib
import matplotlib.lines as mlines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
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
        self.zscale_flag = False
        self.root = None
        self.indpi = 100
        self.zoom = [1, 0, 0]
        self.xposition = None
        self.yposition = None
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
        zmin1, zmax1 = self.get_limits(self.image)
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
        except Exception:
            pass

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
        except IndexError:
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
        Post mouse position on image to the status line.

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
            s1 = "Position: x = %.2f y = %.2f Value: %s" % (x1, y1, value)
            self.imagePosLabelText.set(s1)
            self.imagexpos = event.xdata
            self.imageypos = event.ydata
        except Exception:
            pass

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
            self.canvas1.draw()

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

