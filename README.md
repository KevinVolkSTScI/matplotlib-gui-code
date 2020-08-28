# matplotlib-gui-code

The code here is a general driver code for the Matplotlib package, similar in
intent to the xmgrace code that the author has used for a long time (the
xmgrace web site is at https://plasma-gate.weizmann.ac.il/Grace/, but since it
depends on X11 and Motif it does not seem to be a viable option for the future).

This code is a single monolithic code intended to be run from the command line.
It uses matplotlib and tkinter as well as other standard Python packages.  One
can also use the interface from the Python interpreter.

** Note: MacOS 10.14 Mojave Issue With Python version 3.7

The code here uses the Tkinter package that is part of the standard Python 
distribution.  In the MacOS version 10.14 (Majave) there is an issue that
causes the use of the Tkinter package functions with Python 3.7 to crash
the user's session.  This issue does not exist in Python versions 3.6 or 
3.8 with this version of MacOS, and it did not exist in the previous MacOS 
version as far as the author is aware.    There is no solution for this issue
within Python.  Hence, do not use this code with that particular combination
of MacOS 10.14 and Python 3.7.

** Use of the code

The code assumes that one will want to read in data values, possibly with
associated uncertainties, from an ascii file.  One can have one or more "sets"
of data values in one or more plots in the display area.  The code also has
some capability to create a data set via a formula or in a text area window
created by the code.

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

