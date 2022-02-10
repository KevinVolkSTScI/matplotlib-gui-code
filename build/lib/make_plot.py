"""
This routine is the main plotting routine.  The routine needs to be passed
the matplotlib_user_interface variable ("plotgui") which which the data
values and options (plus the plot variables) are stored.

Routines:

make_plot    The main plottine routine, which reads many options and applies
             them as needed.

"""
import numpy
import matplotlib
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Ellipse, FancyArrow
from matplotlib.ticker import MultipleLocator
import general_utilities

def make_plot(plotgui, current=None):
    """
    Read the parameters and the data sets, and produce the plot.

    This is the main plotting routine where the graphs are produced.  It
    looks at the (many) parameters to determine how the plot should be
    made.

    Parameters
    ----------
        plotgui:  a matplotlib_user_interface object for the plot

        current:  an optional integer value to override plotgui.current_plot

    Returns
    -------
        None

    """
    if current is None:
        current = plotgui.current_plot
    myfont = {'family': plotgui.fontname[current-1],
              'weight': plotgui.fontweight[current-1],
              'size': plotgui.fontsize[current-1]}
    matplotlib.rc('font', **myfont)
    plotgui.subplot[current-1].clear()
    if plotgui.hide_subplot[current-1]:
        plotgui.subplot[current-1].axis('off')
    logxflag = plotgui.xparameters[current-1]['logarithmic']
    logyflag = plotgui.yparameters[current-1]['logarithmic']
    hlogxflag = plotgui.xparameters[current-1]['hybridlog']
    hlogyflag = plotgui.yparameters[current-1]['hybridlog']
    invertxflag = plotgui.xparameters[current-1]['invert']
    invertyflag = plotgui.yparameters[current-1]['invert']
    hidexflag = plotgui.xparameters[current-1]['hide']
    hideyflag = plotgui.yparameters[current-1]['hide']
    hidexticksflag = plotgui.xparameters[current-1]['hideticks']
    hideyticksflag = plotgui.yparameters[current-1]['hideticks']
    hidexlabelsflag = plotgui.xparameters[current-1]['hidelabels']
    hideylabelsflag = plotgui.yparameters[current-1]['hidelabels']
    inversexticksflag = \
        plotgui.xparameters[current-1]['inverseticks']
    inverseyticksflag = \
        plotgui.yparameters[current-1]['inverseticks']
    bothxticksflag = plotgui.xparameters[current-1]['bothticks']
    bothyticksflag = plotgui.yparameters[current-1]['bothticks']
    oppositexflag = plotgui.xparameters[current-1]['oppositeaxis']
    oppositeyflag = plotgui.yparameters[current-1]['oppositeaxis']
    majorxgridlinesflag = plotgui.xparameters[ current-1][
        'majorgridlines']
    minorxgridlinesflag = plotgui.xparameters[ current-1][
        'minorgridlines']
    majorygridlinesflag = plotgui.yparameters[ current-1][
        'majorgridlines']
    minorygridlinesflag = plotgui.yparameters[ current-1][
        'minorgridlines']
    gridoptionx = minorxgridlinesflag + 2*majorxgridlinesflag
    gridoptiony = minorygridlinesflag + 2*majorygridlinesflag
    gridflags1 = [[False, 'both', 'x'], [True, 'minor', 'x'],
                  [True, 'major', 'x'], [True, 'both', 'x']]
    gridflags2 = [[False, 'both', 'y'], [True, 'minor', 'y'],
                  [True, 'major', 'y'], [True, 'both', 'y']]
    try:
        xminorticks = \
            float(plotgui.xparameters[current-1]['minorticks'])
    except ValueError:
        xminorticks = 0.0
    try:
        yminorticks = \
            float(plotgui.yparameters[current-1]['minorticks'])
    except ValueError:
        yminorticks = 0.0
    for loop in range(plotgui.nsets):
        if (plotgui.set_properties[loop]['display']) \
           and (plotgui.set_properties[loop]['plot'] == current):
            if plotgui.set_properties[loop]['errors']:
                xerrors = numpy.zeros((2, len(plotgui.xdata[loop]['values'])),
                                      dtype=numpy.float32)
                xerrors[0, :] = numpy.copy(plotgui.xdata[loop]['lowerror'])
                xerrors[1, :] = numpy.copy(plotgui.xdata[loop]['higherror'])
                yerrors = numpy.zeros((2, len(plotgui.xdata[loop]['values'])),
                                      dtype=numpy.float32)
                yerrors[0, :] = numpy.copy(plotgui.ydata[loop]['lowerror'])
                yerrors[1, :] = numpy.copy(plotgui.ydata[loop]['higherror'])
            try:
                if plotgui.set_properties[loop]['symbol'] == 'histogram':
                    barwidth = plotgui.xdata[loop]['values'][1] - \
                        plotgui.xdata[loop]['values'][0]
                    xoff = plotgui.xdata[loop]['values'][1:] - \
                           plotgui.xdata[loop]['values'][0:-1]
                    xoff = numpy.append(xoff, xoff[-1])
                    plotgui.subplot[current-1].bar(
                        plotgui.xdata[loop]['values'] - xoff/2.,
                        plotgui.ydata[loop]['values'],
                        width=barwidth,
                        color='white',
                        edgecolor=plotgui.set_properties[loop]['colour'],
                        linewidth=1)
                elif logyflag == 0 and logxflag == 0 and hlogxflag == 0 \
                     and hlogyflag == 0:
                    if plotgui.set_properties[loop]['symbol'] is None:
                        plotgui.subplot[current-1].plot(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=plotgui.set_properties[loop]['colour'],
                            linestyle=
                            plotgui.set_properties[loop]['linestyle'],
                            linewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    elif plotgui.set_properties[loop]['linestyle'] is None:
                        plotgui.subplot[current-1].plot(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=plotgui.set_properties[loop]['colour'],
                            marker=plotgui.set_properties[loop]['symbol'],
                            linestyle='none', markersize=
                            plotgui.set_properties[loop]['symbolsize'])
                    else:
                        plotgui.subplot[current-1].plot(
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
                    if plotgui.set_properties[loop]['errors']:
                        plotgui.subplot[current-1].errorbar(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'], yerrors, xerrors,
                            fmt='None',
                            ecolor=plotgui.set_properties[loop]['colour'],
                            elinewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    if xminorticks > 0:
                        plotgui.subplot[current-1].xaxis.\
                            set_minor_locator(MultipleLocator(xminorticks))
                        xtl = int(
                            plotgui.xparameters[
                                current-1]['ticklength'] // 2)
                        plotgui.subplot[current-1].tick_params(
                            axis='x', which='minor', length=xtl)
                    if yminorticks > 0:
                        plotgui.subplot[current-1].yaxis.\
                            set_minor_locator(MultipleLocator(yminorticks))
                        ytl = int(
                            plotgui.yparameters[
                                current-1]['ticklength'] // 2)
                        plotgui.subplot[current-1].tick_params(
                            axis='x', which='minor', length=ytl)
                elif hlogyflag == 0 and hlogxflag == 1:
                    newxvalues = general_utilities.hybrid_transform(
                        plotgui.xdata[loop]['values'])
                    if logyflag == 0:
                        if plotgui.set_properties[loop]['symbol'] is None:
                            plotgui.subplot[current-1].plot(
                                newxvalues, plotgui.ydata[loop]['values'],
                                color=plotgui.set_properties[loop]['colour'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
                        elif plotgui.set_properties[loop]['linestyle'] is \
                             None:
                            plotgui.subplot[current-1].plot(
                                newxvalues, plotgui.ydata[loop]['values'],
                                color=plotgui.set_properties[loop]['colour'],
                                marker=
                                plotgui.set_properties[loop]['symbol'],
                                linestyle='none',
                                markersize=
                                plotgui.set_properties[loop]['symbolsize'])
                        else:
                            plotgui.subplot[current-1].plot(
                                newxvalues, plotgui.ydata[loop]['values'],
                                color=plotgui.set_properties[loop]['colour'],
                                marker=
                                plotgui.set_properties[loop]['symbol'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                markersize=
                                plotgui.set_properties[loop]['symbolsize'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
                        if yminorticks > 0:
                            plotgui.subplot[current-1].yaxis.\
                                set_minor_locator(
                                    MultipleLocator(yminorticks))
                            ytl = int(
                                plotgui.xparameters[current-1]
                                ['ticklength'] // 2)
                            plotgui.subplot[current-1].\
                                tick_params(axis='y', which='minor',
                                            length=ytl)
                    else:
                        if plotgui.set_properties[loop]['symbol'] is None:
                            plotgui.subplot[current-1].semilogy(
                                newxvalues, plotgui.ydata[loop]['values'],
                                color=plotgui.set_properties[loop]['colour'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
                        elif plotgui.set_properties[loop]['linestyle'] is \
                             None:
                            plotgui.subplot[
                                current-1].semilogy(
                                    newxvalues, plotgui.ydata[loop]['values'],
                                    color=
                                    plotgui.set_properties[loop]['colour'],
                                    marker=
                                    plotgui.set_properties[loop]['symbol'],
                                    linestyle='none', markersize=
                                    plotgui.set_properties[loop]['symbolsize']
                                )
                        else:
                            plotgui.subplot[current-1].semilogy(
                                newxvalues, plotgui.ydata[loop]['values'],
                                color=plotgui.set_properties[loop]['colour'],
                                marker=
                                plotgui.set_properties[loop]['symbol'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                markersize=
                                plotgui.set_properties[loop]['symbolsize'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
                        if xminorticks > 0:
                            plotgui.subplot[current-1].xaxis.\
                                set_minor_locator(MultipleLocator(xminorticks))
                            xtl = int(
                                plotgui.xparameters[
                                    current-1]['ticklength'] // 2)
                            plotgui.subplot[current-1].tick_params(
                                axis='x', which='minor', length=xtl)
                elif hlogyflag == 1 and hlogxflag == 0:
                    newyvalues = general_utilities.hybrid_transform(
                        plotgui.ydata[loop]['values'])
                    if logxflag == 0:
                        if plotgui.set_properties[loop]['symbol'] is None:
                            plotgui.subplot[current-1].plot(
                                plotgui.xdata[loop]['values'], newyvalues,
                                color=plotgui.set_properties[loop]['colour'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
                        elif plotgui.set_properties[loop]['linestyle'] is \
                             None:
                            plotgui.subplot[current-1].plot(
                                plotgui.xdata[loop]['values'], newyvalues,
                                color=plotgui.set_properties[loop]['colour'],
                                marker=
                                plotgui.set_properties[loop]['symbol'],
                                linestyle='none',
                                markersize=
                                plotgui.set_properties[loop]['symbolsize'])
                        else:
                            plotgui.subplot[current-1].plot(
                                plotgui.xdata[loop]['values'], newyvalues,
                                color=plotgui.set_properties[loop]['colour'],
                                marker=
                                plotgui.set_properties[loop]['symbol'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                markersize=
                                plotgui.set_properties[loop]['symbolsize'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
                        if xminorticks > 0:
                            plotgui.subplot[current-1].xaxis.\
                                set_minor_locator(MultipleLocator(
                                    xminorticks))
                            xtl = int(
                                plotgui.xparameters[
                                    current-1]['ticklength']
                                // 2)
                            plotgui.subplot[current-1].tick_params(
                                axis='x', which='minor', length=xtl)
                    else:
                        if plotgui.set_properties[loop]['symbol'] is None:
                            plotgui.subplot[current-1].semilogx(
                                plotgui.xdata[loop]['values'],
                                newyvalues,
                                color=plotgui.set_properties[loop]['colour'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
                        elif plotgui.set_properties[loop]['linestyle'] is \
                             None:
                            plotgui.subplot[current-1].semilogx(
                                plotgui.xdata[loop]['values'],
                                newyvalues,
                                color=plotgui.set_properties[loop]['colour'],
                                marker=plotgui.set_properties[loop]['symbol'],
                                linestyle='none',
                                markersize=
                                plotgui.set_properties[loop]['symbolsize'])
                        else:
                            plotgui.subplot[current-1].semilogx(
                                plotgui.xdata[loop]['values'],
                                newyvalues,
                                color=plotgui.set_properties[loop]['colour'],
                                marker=
                                plotgui.set_properties[loop]['symbol'],
                                linestyle=
                                plotgui.set_properties[loop]['linestyle'],
                                markersize=
                                plotgui.set_properties[loop]['symbolsize'],
                                linewidth=
                                plotgui.set_properties[loop]['linewidth'])
                        if xminorticks > 0:
                            plotgui.subplot[current-1].xaxis.\
                                set_minor_locator(MultipleLocator(
                                    xminorticks))
                            xtl = int(
                                plotgui.xparameters[
                                    current-1]['ticklength']
                                // 2)
                            plotgui.subplot[current-1].tick_params(
                                axis='x', which='minor', length=xtl)
                elif hlogyflag == 1 and hlogxflag == 1:
                    newxvalues = general_utilities.hybrid_transform(
                        plotgui.xdata[loop]['values'])
                    newyvalues = general_utilities.hybrid_transform(
                        plotgui.ydata[loop]['values'])
                    if plotgui.set_properties[loop]['symbol'] is None:
                        plotgui.subplot[current-1].plot(
                            newxvalues, newyvalues,
                            color=plotgui.set_properties[loop]['colour'],
                            linestyle=
                            plotgui.set_properties[loop]['linestyle'],
                            linewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    elif plotgui.set_properties[loop]['linestyle'] is None:
                        plotgui.subplot[current-1].plot(
                            newxvalues, newyvalues,
                            color=plotgui.set_properties[loop]['colour'],
                            marker=plotgui.set_properties[loop]['symbol'],
                            linestyle='none',
                            markersize=
                            plotgui.set_properties[loop]['symbolsize'])
                    else:
                        plotgui.subplot[current-1].plot(
                            newxvalues, newyvalues,
                            color=plotgui.set_properties[loop]['colour'],
                            marker=plotgui.set_properties[loop]['symbol'],
                            linestyle=
                            plotgui.set_properties[loop]['linestyle'],
                            markersize=
                            plotgui.set_properties[loop]['symbolsize'],
                            linewidth=
                            plotgui.set_properties[loop]['linewidth'])
                elif logyflag == 0 and logxflag == 1 and hlogxflag == 0:
                    if plotgui.set_properties[loop]['symbol'] is None:
                        plotgui.subplot[current-1].semilogx(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=plotgui.set_properties[loop]['colour'],
                            linestyle=
                            plotgui.set_properties[loop]['linestyle'],
                            linewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    elif plotgui.set_properties[loop]['linestyle'] is None:
                        plotgui.subplot[current-1].semilogx(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=plotgui.set_properties[loop]['colour'],
                            marker=plotgui.set_properties[loop]['symbol'],
                            linestyle='None',
                            markersize=
                            plotgui.set_properties[loop]['symbolsize'])
                    else:
                        plotgui.subplot[current-1].semilogx(
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
                    if plotgui.set_properties[loop]['errors']:
                        plotgui.subplot[current-1].errorbar(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'], xerrors, yerrors,
                            fmt=None, ecolor=
                            plotgui.set_properties[loop]['colour'],
                            elinewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    if yminorticks > 0:
                        plotgui.subplot[current-1].\
                            yaxis.set_minor_locator(
                                MultipleLocator(yminorticks))
                        ytl = int(
                            plotgui.xparameters[current-1][
                                'ticklength'] // 2)
                        plotgui.subplot[current-1].tick_params(
                            axis='y', which='minor', length=ytl)
                elif logxflag == 0 and logyflag == 1 and hlogyflag == 0:
                    if plotgui.set_properties[loop]['symbol'] is None:
                        plotgui.subplot[current-1].semilogy(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=plotgui.set_properties[loop]['colour'],
                            linestyle=
                            plotgui.set_properties[loop]['linestyle'],
                            linewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    elif plotgui.set_properties[loop]['linestyle'] is None:
                        plotgui.subplot[current-1].semilogy(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=plotgui.set_properties[loop]['colour'],
                            marker=plotgui.set_properties[loop]['symbol'],
                            linestyle='None',
                            markersize=
                            plotgui.set_properties[loop]['symbolsize'])
                    else:
                        plotgui.subplot[current-1].semilogy(
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
                    if plotgui.set_properties[loop]['errors']:
                        plotgui.subplot[current-1].errorbar(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'], xerrors, yerrors,
                            fmt=None,
                            ecolor=plotgui.set_properties[loop]['colour'],
                            elinewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    if xminorticks > 0:
                        plotgui.subplot[current-1].xaxis.\
                            set_minor_locator(MultipleLocator(xminorticks))
                        xtl = int(
                            plotgui.xparameters[
                                current-1]['ticklength'] // 2)
                        plotgui.subplot[current-1].tick_params(
                            axis='x', which='minor', length=xtl)
                elif logxflag == 1 and logyflag == 1:
                    if plotgui.set_properties[loop]['symbol'] is None:
                        plotgui.subplot[current-1].loglog(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=
                            plotgui.set_properties[loop]['colour'],
                            linestyle=
                            plotgui.set_properties[loop]['linestyle'],
                            linewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    elif plotgui.set_properties[loop]['linestyle'] is None:
                        plotgui.subplot[current-1].loglog(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=plotgui.set_properties[loop]['colour'],
                            marker=plotgui.set_properties[loop]['symbol'],
                            linestyle='None',
                            markersize=
                            plotgui.set_properties[loop]['symbolsize'])
                    else:
                        plotgui.subplot[current-1].loglog(
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
                    if plotgui.set_properties[loop]['errors']:
                        plotgui.subplot[current-1].errorbar(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            xerrors, yerrors, fmt=None,
                            ecolor=plotgui.set_properties[loop]['colour'],
                            elinewidth=
                            plotgui.set_properties[loop]['linewidth'])
                else:
                    # One should not get here, but if the flags do not
                    # correspond to the expected options use a linear
                    # plot as the default
                    if plotgui.set_properties[loop]['symbol'] is None:
                        plotgui.subplot[current-1].plot(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=plotgui.set_properties[loop]['colour'],
                            linestyle=
                            plotgui.set_properties[loop]['linestyle'],
                            linewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    elif plotgui.set_properties[loop]['linestyle'] is None:
                        plotgui.subplot[current-1].plot(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'],
                            color=plotgui.set_properties[loop]['colour'],
                            marker=plotgui.set_properties[loop]['symbol'],
                            linestyle='none',
                            markersize=
                            plotgui.set_properties[loop]['symbolsize'])
                    else:
                        plotgui.subplot[current-1].plot(
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
                    if plotgui.set_properties[loop]['errors']:
                        plotgui.subplot[current-1].errorbar(
                            plotgui.xdata[loop]['values'],
                            plotgui.ydata[loop]['values'], yerrors, xerrors,
                            fmt='None',
                            ecolor=plotgui.set_properties[loop]['colour'],
                            elinewidth=
                            plotgui.set_properties[loop]['linewidth'])
                    if xminorticks > 0:
                        plotgui.subplot[current-1].xaxis.\
                            set_minor_locator(MultipleLocator(xminorticks))
                        xtl = int(plotgui.xparameters[
                            current-1]['ticklength'] // 2)
                        plotgui.subplot[current-1].tick_params(
                            axis='x', which='minor', length=xtl)
                    if yminorticks > 0:
                        plotgui.subplot[
                            current-1].yaxis.set_minor_locator(
                                MultipleLocator(yminorticks))
                        ytl = int(plotgui.yparameters[
                            current-1]['ticklength'] // 2)
                        plotgui.subplot[current-1].tick_params(
                            axis='y', which='minor', length=ytl)
            except Exception:
                pass
    if bothyticksflag == 1:
        plotgui.subplot[current-1].tick_params(
            left=True, right=True, which='both')
    else:
        plotgui.subplot[current-1].tick_params(
            left=True, right=False, which='both')
    if bothxticksflag == 1:
        plotgui.subplot[current-1].tick_params(
            bottom=True, top=True, which='both')
    else:
        plotgui.subplot[current-1].tick_params(
            bottom=True, top=False, which='both')
    for n1 in range(plotgui.number_of_lines):
        if plotgui.plot_lines[n1]['plot'] == current:
            xlvalues = numpy.asarray([plotgui.plot_lines[n1]['xstart'],
                                      plotgui.plot_lines[n1]['xend']])
            if hlogxflag == 1:
                xlvalues[0] = general_utilities.hybrid_transform(
                    xlvalues[0])
                xlvalues[1] = general_utilities.hybrid_transform(
                    xlvalues[1])
            ylvalues = numpy.asarray([plotgui.plot_lines[n1]['ystart'],
                                      plotgui.plot_lines[n1]['yend']])
            if hlogyflag == 1:
                ylvalues[0] = general_utilities.hybrid_transform(
                    ylvalues[0])
                ylvalues[1] = general_utilities.hybrid_transform(
                    ylvalues[1])
            plotgui.subplot[current-1].plot(
                xlvalues, ylvalues,
                color=plotgui.plot_lines[n1]['line_colour'],
                linestyle=plotgui.plot_lines[n1]['line_type'],
                linewidth=plotgui.plot_lines[n1]['line_thickness'])
    patches = []
    for n1 in range(plotgui.number_of_boxes):
        if plotgui.plot_boxes[n1]['plot'] == current:
            xcorner = plotgui.plot_boxes[n1]['xstart']
            ycorner = plotgui.plot_boxes[n1]['ystart']
            delx = plotgui.plot_boxes[n1]['xend'] \
                - plotgui.plot_boxes[n1]['xstart']
            dely = plotgui.plot_boxes[n1]['yend'] \
                - plotgui.plot_boxes[n1]['ystart']
            if delx < 0.:
                xcorner = xcorner + delx
                delx = -delx
            if dely < 0.:
                ycorner = ycorner + dely
                dely = -dely
            if hlogxflag == 1:
                xc1 = xcorner
                xc2 = xcorner+delx
                xcorner = general_utilities.hybrid_transform(xc1)
                delx = general_utilities.hybrid_transform(xc2) - xcorner
            if hlogyflag == 1:
                yc1 = ycorner
                yc2 = ycorner+dely
                ycorner = general_utilities.hybrid_transform(yc1)
                dely = general_utilities.hybrid_transform(yc2) - ycorner
            rect = Rectangle(
                (xcorner, ycorner), delx, dely,
                angle=plotgui.plot_boxes[n1]['rotation'],
                edgecolor=plotgui.plot_boxes[n1]['line_colour'],
                linestyle=plotgui.plot_boxes[n1]['line_type'],
                linewidth=plotgui.plot_boxes[n1]['line_thickness'],
                facecolor=plotgui.plot_boxes[n1]['fill_colour'])
            plotgui.subplot[current-1].add_artist(rect)
            patches.append(rect)
    if len(patches) > 0:
        pc = PatchCollection(patches)
        plotgui.subplot[current-1].add_collection(pc)
    patches = []
    for n1 in range(plotgui.number_of_ellipses):
        if plotgui.plot_ellipses[n1]['plot'] == current:
            xcenter = plotgui.plot_ellipses[n1]['xposition']
            ycenter = plotgui.plot_ellipses[n1]['yposition']
            delx = plotgui.plot_ellipses[n1]['major']
            dely = plotgui.plot_ellipses[n1]['minor']
            if hlogxflag == 1:
                xc1 = xcorner
                xc2 = xcorner+delx
                xcorner = general_utilities.hybrid_transform(xc1)
                delx = general_utilities.hybrid_transform(xc2) - xcorner
            if hlogyflag == 1:
                yc1 = ycorner
                yc2 = ycorner+dely
                ycorner = general_utilities.hybrid_transform(yc1)
                dely = general_utilities.hybrid_transform(yc2) - ycorner
            ellip = Ellipse(
                (xcenter, ycenter), delx, dely,
                angle=plotgui.plot_ellipses[n1]['rotation'],
                edgecolor=plotgui.plot_ellipses[n1]['line_colour'],
                linestyle=plotgui.plot_ellipses[n1]['line_type'],
                linewidth=plotgui.plot_ellipses[n1]['line_thickness'],
                facecolor=plotgui.plot_ellipses[n1]['fill_colour'])
            plotgui.subplot[current-1].add_artist(ellip)
            patches.append(ellip)
    if len(patches) > 0:
        pc = PatchCollection(patches)
        plotgui.subplot[current-1].add_collection(pc)
    patches = []
    for n1 in range(plotgui.number_of_vectors):
        if plotgui.plot_vectors[n1]['plot'] == current:
            x1 = plotgui.plot_vectors[n1]['xstart']
            y1 = plotgui.plot_vectors[n1]['ystart']
            xlength = plotgui.plot_vectors[n1]['xend'] \
                - plotgui.plot_vectors[n1]['xstart']
            ylength = plotgui.plot_vectors[n1]['yend'] \
                - plotgui.plot_vectors[n1]['ystart']
            delx = plotgui.plot_vectors[n1]['delx']
            dely = plotgui.plot_vectors[n1]['dely']
            if hlogxflag == 1:
                xs1 = general_utilities.hybrid_transform(x1)
                xs2 = general_utilities.hybrid_transform(
                    plotgui.plot_vectors[n1]['xend'])
                xlength = xs2 - xs1
                delx = xs2 - general_utilities.hybrid_transform(
                    plotgui.plot_vectors[n1]['xend']-delx)
            if hlogyflag == 1:
                ys1 = general_utilities.hybrid_transform(y1)
                ys2 = general_utilities.hybrid_transform(
                    plotgui.plot_vectors[n1]['yend'])
                ylength = ys2 - ys1
                delx = ys2 - general_utilities.hybrid_transform(
                    plotgui.plot_vectors[n1]['yend']-dely)
            arrow = FancyArrow(
                x1, y1, xlength, ylength, head_width=delx,
                head_length=dely, fill=plotgui.plot_vectors[n1]['fill'],
                edgecolor=plotgui.plot_vectors[n1]['line_colour'],
                linestyle=plotgui.plot_vectors[n1]['line_type'],
                linewidth=plotgui.plot_vectors[n1]['line_thickness'],
                facecolor=plotgui.plot_vectors[n1]['fill_colour'])
            plotgui.subplot[current-1].add_artist(arrow)
            patches.append(arrow)
    if len(patches) > 0:
        pc = PatchCollection(patches)
        plotgui.subplot[current-1].add_collection(pc)
    if plotgui.plot_frame[current-1] > 0.:
        plotgui.subplot[current-1].spines['bottom'].set_linewidth(
            plotgui.plot_frame[current-1])
        plotgui.subplot[current-1].spines['top'].set_linewidth(
            plotgui.plot_frame[current-1])
        plotgui.subplot[current-1].spines['left'].set_linewidth(
            plotgui.plot_frame[current-1])
        plotgui.subplot[current-1].spines['right'].set_linewidth(
            plotgui.plot_frame[current-1])
    else:
        plotgui.subplot[current-1].spines['bottom'].set_linewidth(
            0.5)
        plotgui.subplot[current-1].spines['top'].set_linewidth(
            0.5)
        plotgui.subplot[current-1].spines['left'].set_linewidth(
            0.5)
        plotgui.subplot[current-1].spines['right'].set_linewidth(
            0.5)
    xisinverted = plotgui.subplot[current-1].xaxis_inverted()
    yisinverted = plotgui.subplot[current-1].yaxis_inverted()
    if (invertxflag == 1 and (not xisinverted)) or \
       (invertxflag == 0 and xisinverted):
        plotgui.subplot[current-1].invert_xaxis()
    if (invertyflag == 1 and (not yisinverted)) or \
       (invertyflag == 0 and yisinverted):
        plotgui.subplot[current-1].invert_yaxis()
    if plotgui.matplotlib_rounding.get():
        xmin1, xmax1 = plotgui.subplot[current-1].get_xbound()
        ymin1, ymax1 = plotgui.subplot[current-1].get_ybound()
        if plotgui.original_range[current-1]:
            plotgui.plot_range[current-1][0] = xmin1
            plotgui.plot_range[current-1][1] = xmax1
            plotgui.plot_range[current-1][2] = ymin1
            plotgui.plot_range[current-1][3] = ymax1
            plotgui.original_range[current-1] = False
    if hlogxflag == 0:
        plotgui.subplot[current-1].set_xbound(
            plotgui.plot_range[current-1][0],
            plotgui.plot_range[current-1][1])
    else:
        xrange = numpy.asarray([plotgui.plot_range[current-1][0],
                                plotgui.plot_range[current-1][1]],
                               dtype=numpy.float32)
        xrange1 = general_utilities.hybrid_transform(xrange)
        plotgui.subplot[current-1].set_xbound(xrange1[0],
                                                     xrange1[1])
        newxvalues = general_utilities.hybrid_transform(
            plotgui.xdata[loop]['values'])
        tickmarks, ticklabels = general_utilities.hybrid_labels(xrange1)
        plotgui.subplot[current-1].set_xticks(tickmarks)
        plotgui.subplot[current-1].set_xticklabels(ticklabels)
    if hlogyflag == 0:
        plotgui.subplot[current-1].set_ybound(
            plotgui.plot_range[current-1][2],
            plotgui.plot_range[current-1][3])
    else:
        yrange = numpy.asarray([plotgui.plot_range[current-1][2],
                                plotgui.plot_range[current-1][3]],
                               dtype=numpy.float32)
        yrange1 = general_utilities.hybrid_transform(yrange)
        plotgui.subplot[current-1].set_ybound(yrange1[0],
                                                     yrange1[1])
        newyvalues = general_utilities.hybrid_transform(
            plotgui.ydata[loop]['values'])
        tickmarks, ticklabels = general_utilities.hybrid_labels(yrange1)
        plotgui.subplot[current-1].set_yticks(tickmarks)
        plotgui.subplot[current-1].set_yticklabels(ticklabels)
    plotgui.subplot[current-1].set_xlabel(
        plotgui.xparameters[current-1]['label'],
        family=plotgui.fontname[current-1],
        size=plotgui.fontsize[current-1],
        weight=plotgui.fontweight[current-1])
    plotgui.subplot[current-1].set_ylabel(
        plotgui.yparameters[current-1]['label'],
        family=plotgui.fontname[current-1],
        size=plotgui.fontsize[current-1],
        weight=plotgui.fontweight[current-1])
    try:
        plotgui.subplot[current-1].set_title(
            plotgui.title[current-1],
            family=plotgui.fontname[current-1],
            size=plotgui.fontsize[current-1],
            weight=plotgui.fontweight[current-1])
    except Exception:
        pass
    # Adjust the margins if the plot_margin value is set.
    left = plotgui.bounding_box[current-1][0] + plotgui.plot_margin
    right = plotgui.bounding_box[current-1][0] \
        + plotgui.bounding_box[current-1][2] - plotgui.plot_margin
    bottom = plotgui.bounding_box[current-1][1] + plotgui.plot_margin
    top = plotgui.bounding_box[current-1][1] \
          + plotgui.bounding_box[current-1][3] - plotgui.plot_margin
    plotgui.subplot[current-1].set_position(
        [left, bottom, right-left, top-bottom], which='both')
    if plotgui.number_of_labels > 0:
        for loop in range(plotgui.number_of_labels):
            if current == plotgui.plot_labels[loop]['plot']:
                xpos1 = plotgui.plot_labels[loop]['xposition']
                if plotgui.xparameters[current-1]['hybridlog'] == 1:
                    xpos1 = general_utilities.hybrid_transform(xpos1)
                ypos1 = plotgui.plot_labels[loop]['yposition']
                if plotgui.xparameters[current-1]['hybridlog'] == 1:
                    ypos1 = general_utilities.hybrid_transform(ypos1)
                plotgui.subplot[current-1].text(
                    xpos1, ypos1,
                    plotgui.plot_labels[loop]['labelstring'],
                    {'color': plotgui.plot_labels[loop]['colour'],
                     'fontsize': plotgui.plot_labels[loop]['size'],
                     'fontfamily': plotgui.plot_labels[loop]['font'],
                     'fontweight': plotgui.plot_labels[loop]['fontweight']})
    try:
        legend_flag = plotgui.legend_variable[current-1].get()
    except Exception:
        legend_flag = 0
    if legend_flag:
        plotgui.generate_legend(None)
        legend_option = plotgui.legend_options[current-1].get()
        legend_position = None
        if legend_option == 'user':
            try:
                str1 = plotgui.legend_position_field.get()
                values = str1.split()
                if len(values) == 2:
                    xlpos = float(values[0])
                    ylpos = float(values[1])
                    legend_position = [xlpos, ylpos]
                    plotgui.legend_user_position[current-1] = \
                        legend_position
            except ValueError:
                legend_option = 'best'
        else:
            legend_option = plotgui.legend_position[current-1]
            if legend_option is None:
                legend_position = 'best'
        try:
            legend_frame = plotgui.legend_frame[current-1].get()
        except Exception:
            legend_frame = 0
        if legend_position is None:
            plotgui.subplot[current-1].legend(
                handles=plotgui.legend_handles[current-1],
                labels=plotgui.legend_labels[current-1],
                loc=legend_option, frameon=legend_frame)
        else:
            plotgui.subplot[current-1].legend(
                handles=plotgui.legend_handles[current-1],
                labels=plotgui.legend_labels[current-1],
                loc=legend_position, frameon=legend_frame)
    if oppositexflag == 1:
        plotgui.subplot[current-1].get_xaxis().set_ticks_position(
            "top")
        plotgui.subplot[current-1].get_xaxis().set_label_position(
            "top")
    if oppositeyflag == 1:
        plotgui.subplot[current-1].get_yaxis().set_ticks_position(
            "right")
        plotgui.subplot[current-1].get_yaxis().set_label_position(
            "right")
    if inversexticksflag == 1:
        plotgui.subplot[current-1].tick_params(
            axis='x', direction='in', length=
            plotgui.xparameters[current-1]['ticklength'])
        plotgui.subplot[current-1].tick_params(
            axis='x', direction='in', which='minor')
    else:
        plotgui.subplot[current-1].tick_params(
            axis='x', direction='out', length=
            plotgui.xparameters[current-1]['ticklength'])
        plotgui.subplot[current-1].tick_params(
            axis='x', direction='out', which='minor')
    if inverseyticksflag == 1:
        plotgui.subplot[current-1].tick_params(
            axis='y', direction='in', length=
            plotgui.yparameters[current-1]['ticklength'])
        plotgui.subplot[current-1].tick_params(axis='y',
                                                      direction='in',
                                                      which='minor')
    else:
        plotgui.subplot[current-1].tick_params(
            axis='y', direction='out', length=
            plotgui.yparameters[current-1]['ticklength'])
        plotgui.subplot[current-1].tick_params(axis='y',
                                                      direction='out',
                                                      which='minor')
    if plotgui.grid_linetype[current-1] is None:
        stylestring = 'None'
    else:
        stylestring = plotgui.grid_linetype[current-1]
    if plotgui.grid_colour[current-1] is None:
        colourvalue = 'black'
    else:
        colourvalue = plotgui.grid_colour[current-1]
    if (gridoptionx == 0) and (gridoptiony == 0):
        plotgui.subplot[current-1].grid(b=False, axis='both',
                                               which='both')
    else:
        if gridoptionx == 0:
            plotgui.subplot[current-1].grid(b=False, axis='x',
                                                   which='both')
        else:
            plotgui.subplot[current-1].grid(
                b=gridflags1[gridoptionx][0],
                which=gridflags1[gridoptionx][1],
                axis=gridflags1[gridoptionx][2], linestyle=stylestring,
                color=colourvalue)
        if gridoptiony == 0:
            plotgui.subplot[current-1].grid(b=False, axis='y',
                                                   which='both')
        else:
            plotgui.subplot[current-1].grid(
                b=gridflags2[gridoptiony][0],
                which=gridflags2[gridoptiony][1],
                axis=gridflags2[gridoptiony][2], linestyle=stylestring,
                color=colourvalue)
    if hidexticksflag == 1:
        plotgui.subplot[current-1].get_xaxis().set_ticks([])
    if hideyticksflag == 1:
        plotgui.subplot[current-1].get_yaxis().set_ticks([])
    if hidexlabelsflag == 1:
        plotgui.subplot[current-1].get_xaxis().set_ticklabels([])
    if hideylabelsflag == 1:
        plotgui.subplot[current-1].get_yaxis().set_ticklabels([])
    if hidexflag == 1:
        plotgui.subplot[current-1].get_xaxis().set_visible(False)
    else:
        plotgui.subplot[current-1].get_xaxis().set_visible(True)
    if hideyflag == 1:
        plotgui.subplot[current-1].get_yaxis().set_visible(False)
    else:
        plotgui.subplot[current-1].get_yaxis().set_visible(True)
    try:
        if plotgui.equal_aspect[current-1]:
            plotgui.figure.gca().set_aspect('equal', adjustable='box')
        else:
            plotgui.figure.gca().set_aspect('auto')
    except:
        print(plotgui.equal_aspect)
        plotgui.figure.gca().set_aspect('auto')
    plotgui.canvas.draw()
    # The follow is duplicate information, but may be used to set tick
    # intervals....
    if not plotgui.matplotlib_rounding.get():
        plotgui.xparameters[current-1]['minimum'] = \
            plotgui.plot_range[current-1][0]
        plotgui.xparameters[current-1]['maximum'] = \
            plotgui.plot_range[current-1][1]
        plotgui.yparameters[current-1]['minimum'] = \
            plotgui.plot_range[current-1][2]
        plotgui.yparameters[current-1]['maximum'] = \
            plotgui.plot_range[current-1][3]
    else:
        xr1, xr2 = plotgui.subplot[current-1].get_xbound()
        yr1, yr2 = plotgui.subplot[current-1].get_ybound()
        plotgui.xparameters[current-1]['minimum'] = xr1
        plotgui.xparameters[current-1]['maximum'] = xr2
        plotgui.yparameters[current-1]['minimum'] = yr1
        plotgui.yparameters[current-1]['maximum'] = yr2
