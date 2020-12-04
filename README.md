# matplotlib_user_interface

The code here is a general driver code for the Matplotlib package, similar in
intent to the xmgrace code that the author has used for a long time (the
xmgrace web site is at https://plasma-gate.weizmann.ac.il/Grace/, but since it
depends on X11 and Motif it does not seem to be a viable option for the future).

The code is a package that can be installed and imported in the usual way.  It
is also possible to run the code from the command line if ones $PATH environment
variable is set properly (the author does not use Microsoft Windows and so does
not know what the equivalent is in such computers).  The code uses matplotlib
and tkinter as well as other standard Python packages.

## Note: MacOS 10.14 Mojave Issue With Python version 3.7

The code here uses the Tkinter package that is part of the standard Python 
distribution.  In the MacOS version 10.14 (Majave) there is an issue that
causes the use of the Tkinter package functions with Python 3.7 to crash
the user's session.  This issue does not exist in Python versions 3.6 or 
3.8 with this version of MacOS, and it did not exist in the previous MacOS 
version as far as the author is aware.    There is no solution for this issue
within Python.  Hence, do not use this code with that particular combination
of MacOS 10.14 and Python 3.7.

## Installation and Set-Up

To install the code one needs to copy the repository to some directory, which
will be called /your/path/to/the/files/ in the documentation.  Substitute the
actual path on your system for this dummy path.  It is best to put the package
in a new directory.

When the code is unpacked in the target directory you should see the following
files in the directory:

LICENSE
README.md
matplotlib_driver_code_documentation.docx
matplotlib_driver_code_documentation.pdf
setup.py
src
tests

The tests directory is a placeholder and is currently empty.  The code itself
is in the src directory.

To install the package run the normal command

python setup.py install

in the target directory /your /path/to/the/files.  This will add build and dist
directories.  Do not issue this command in the src directory, it must be done
in the upper level directory.

If one wishes to be able to start the interface from the command line one
needs to add the target directory to the $PATH environment variable.  This
is done on Unix-type systems via the command

setenv PATH /your/path/to/the/files:${PATH}

in tcsh/csh or

export PATH='/your/path/to/the/files:$PATH'

in bash.  This assumes that the $PATH variable is already defined.

## Starting the code

Assuming that the installation worked as expected, one can then import
the plot interface and start it via Python commands such as

\>\>\> import matplotlib_user_interface as mui
\>\>\> root, myplot = mui.startup()

This will cause a new window to open with the interface controls.  See the
main documentation in the matplotlib_driver_code_documentation.pdf or
matplotlib_driver_code_documentation.docx files for more information about
using the interface and adding data sets from the Python interpreter.  The
"myplot" variable is an object that runs the interface.  The "root" variable
is a Tkinter variable for the plot window.

## Use of the code

The code assumes that one will want to read in data values, possibly with
associated uncertainties, from an ascii file.  One can have one or more "sets"
of data values in one or more plots in the display area.  The code also has
some capability to create a data set via a formula or in a text area window
created by the code.  One can add sets from the Python interpreter as well if
the code was started that way.  The latter process is to do the following,
once root and myplot are created as noted above.  First create the x and y
values in Python as one-dimensional numpy arrays, for example

\>\>\> import numpy
\>\>\> xvalues = numpy.arange(0, 100, 1)
\>\>\> yvalues = numpy.sqrt(xvalues)

Once the x and y values are defined in the interpreter, one can load them into
the plot object

>>> myplot.add_set(xvalues, yvalues)

and then tell the code to replot the plot(s) associated with the myplot object

>>> mui.make_plot.make_plot(myplot)

One can add multiple sets to the interface by making additional calls to
the add_set method of the object.

The code also has the machinery to add vectors/boxes/lines/ellipses/text to
a given plot.

Once a plot has been made, one can save the state of the code variables to an
ascii file for later re-loading.  This is analogous to the save/open plot
function in xmgrace.

One has control of the symbols/lines/legend entries/error bars for each set,
and control of the axis limits/labels/tick marks for each axis of each plot
in the display area.

The code is written in Python 3 and will not function in Python 2.  The code
has been tested on a MacOS system and on a Linux system, but not on a
Microsoft Windows system.

