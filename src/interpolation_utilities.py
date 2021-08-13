#! /usr/bin/env python
#
"""
This file contains a couple of resampling routines for one-dimensional
functions.

The calls to the functions assume the inputs are first an array of x values
and then an array of corresponding y values.  The x values are assumed to be
strictly increasing.

Routines defined here:

_blocksum(xvalues, yvalues)

  This routine carries out trapezoidal averaging over the array of x values
and returns float values for the mean x and y values.  The mean x value is
taken to be the midpoint of the range defined by the xvalues array.  The
mean y values is the trapezoidal integration over the interval divided by the
the x range.

_smother(xvalues, yvalues, nsmooth, clip=True)

  This routine does a block averaging over some number of input points given
by the nsmooth value.  If nsmooth is 10 for example, each group of 10 points
in the input arrays are averaged over to produce new output x and y arrays.
It is assumed that the length of the input arrays xvalues and yvalues is
larger than nsmooth.  If clip is True any remaining points at the end of the
input arrays that form a group of length less than nsmooth points are ignored,
and if this value is False the remaining values are block averaged to make a
final point in each of the output arrays.

_resample(xnew, xvalues, yvalues)

  The input array yvalues at points xvalues is resampled to the points xnew.
It is assumed that the xnew array overlaps at least partially with the xvalues
array.  There are options for extrapolation and how the resampling is done.
Resampling is not simply interpolation, it involves integration over the input
function in the case where the new sampling is coarser than the original
sampling.  If the xnew array represents finer sampling than the original
xvalues array then the output here is very similar to what interpolation
using numpy.interp would give.  The difference here is that the resampling
checks the range between a given point in the xnew array and the previous and
subsequent points, and in effect the code integrates over this interval in
the case where the interval includes one of the input xvalues points.  It the
interval does not include such a point, one gets an interpolated value that
is the mean over the range.  If the range is symemtric about the new x value
then this is the same as using numpy.interp, but if the xnew points are
asymemtrically spaced the results will differ slightly from the simple
interpolation case.  This is a matter of concern when the xnew points have
varying resolution in different ranges when the transition from one
resolution to another takes place.  The resampling will also differ a bit
from the simple interpolation when the slopes between the xvalues points are
varying significantly from one point to the next.

"""
import numpy
from scipy import interpolate

def _blocksum(xvalues, yvalues):
    """
    Carry out trapezoidal averaging over a range of points.

    Parameters
    ----------

    xvalues:   a numpy 1-D float array of values, assumed to be sorted in
               ascending order

    yvalues:  a numpy 1-D float array of values, of the same length as xvalues,
              containing the corresponding values to be averaged over

    Returns
    -------

    xout:   the mean x value over the xvalues array, a float value

    yout:   the integrated mean y value over the yvalues array, a float value
    """
    xout = (xvalues[-1]+xvalues[0])/2.
    xrange = (xvalues[-1]-xvalues[0])
    if xrange <= 0.:
        raise ValueError
    delx = xvalues[1:] - xvalues[0:-1]
    if numpy.min(delx) <= 0.:
        raise ValueError
    ypair = (yvalues[0:-1] + yvalues[1:])/2.
    yout = numpy.sum(ypair*delx)/xrange
    return xout, yout

def _smoother(xvalues, yvalues, nsmooth, clip=True):
    """
    Code to block average the x/y values over nsmooth samples.

    Parameters
    ----------

    xvalues:  by assumption a numpy 1-D array of float values

    yvalues:  by assumption a numpy 1-D array of float values the same size as
              xvalues

    nsmooth:  an integer value > 1, the number of values to block average;
              should be less than len(xvalues) // 2 and greator than 1

    clip:     an optional boolean value, if True any "extra" points at the
              end of the input arrays are ignored

    Returns
    -------

    newxvalues:  a 1-D numpy float array of the block averaged x values;

    newyvalues:  a 1-D numpy float array of the block averaged y values;

    The code raises ValueError if the input values are not as expected.

    It is assumed that the input sampling may irregular, so the returned value
    is the trapezoidal integral over the range in x divided by this range.

    If there are N values in the xvalues array and the nsmooth value is M,
    the number of returned values is N // M if N % M is zero or clip is True,
    or (N // M)+1 if N % M is not zero and clip is False.  The last point
    carries out the block averaging over the N % M values when clip is False.
    The clip flag determines whether the entire input x range is used or not.

    """
    xwork = numpy.copy(xvalues)
    ywork = numpy.copy(yvalues)
    if (len(xwork.shape) != 1) or (xwork.shape != ywork.shape) or \
       (nsmooth > len(xvalues) // 2) or (nsmooth < 2):
        raise ValueError
    nin = len(xwork)
    nout = nin // nsmooth
    nmain = nout
    nclip = nin % nsmooth
    if (nclip > 0) and (not clip):
        nout = nout + 1
    newxvalues = numpy.zeros((nout), dtype=numpy.float32)
    newyvalues = numpy.zeros((nout), dtype=numpy.float32)
    for loop in range(nmain):
        n1 = loop*nsmooth
        n2 = (loop+1)*nsmooth
        subx = numpy.copy(xvalues[n1:n2])
        suby = numpy.copy(yvalues[n1:n2])
        newxvalues[loop], newyvalues[loop] = _blocksum(subx, suby)
    if nout > nmain:
        subx = numpy.copy(xvalues[n2:])
        suby = numpy.copy(yvalues[n2:])
        newxvalues[-1], newyvalues[-1] = _blocksum(subx, suby)
    return newxvalues, newyvalues

def _resample(xnew, xvalues, yvalues, extrapolate=False, ends=False,
             cubic=False, noversample=10):
    """
    Routine to do general resampling at any resolution.  Can use either a
    cubic spline interpolation or linear interpolation.

    Parameters
    ----------

    xnew:      a 1-D numpy float array of the new sample points (assumed to
               be sorted in ascending order)

    xvalues:   a 1-D numpy float array of the current sample points (assumed
               to be sorted in ascending order) which must have 5 or more
               values...normally many more than 5 values are expected

    yvalues:   a 1-D numpy float array of the current function values at the
               xvalues points

    extrapolate:   a boolean flag for whether to extrapolate the interpolation
                   function or return a set value (see "ends"); defaults to
                   False

    ends:     a boolean flag for returning the nearest end value if True or
              zero otherwise for points outside the original range, if the
              extrapolate parameter is False; defaults to False

    cubic:   a boolean flag for whether to use cubic spline interpolation
             between points, the alternative being linear interpolation; in
             the linear interpolation case trapezoidal integration is used
             with the input data values where needed, otherwise an oversampled
             version of the cubic spline fit is calculated and this is put to
             the trapezoidal integration step; defaults to False

    noversample:  an integer value larger than 5 for how much to oversample
                  between the input points in the case where the cubic spline
                  option is used; defaults to 10

    Returns
    -------

    ynew:  a numpy 1-D float array of the new y values at the xnew points

    Given a data point in xnew a range around the point is defined, with the
    minimum and maximum from the mid-points between the new values, with the
    intervals on the ends taken to be the same as for the next points in to
    the xvalues array.  The input yvalues are evaluated within this range in
    the following manner:

    (1) if the cubic spline option is used the ranges between the xvalues
        points are over-sampled by the number given in the noversample
        parameter and the values are fed into the numpy.trapz function, then
        the result is divided by the x range;

    (2) if the linear option is used and if the current data range includes
        1 or more of the input xvalues, the end points are determined by
        linear interpolation--subject to the end points option--and these
        along with the xvalues/yvalues points within the interval are given
        to numpy.trapz amd the result is divided by the x range;

    (3) if the linear option is used and if the x range is entirely within the
        interval between points in the original x values array, then the
        interpolated value at the mid-point of the range is used.

    This should work for arbitrary input x points.  Note that the cubic
    spline option will not work with too few input data points.

    The routine assumes some overlap between the xnew array and the xvalues
    array; if there is none, an ValueError exception is raised.
    """
    delx = xvalues[1:] - xvalues[0:-1]
    if numpy.min(delx) <= 0.:
        raise ValueError
    if (xvalues.shape != yvalues.shape) or (len(xvalues.shape) > 1) or \
       (len(xvalues) < 5) or ((noversample < 5) and cubic):
        raise ValueError
    if (numpy.max(xnew) < numpy.min(xvalues)) or (numpy.min(xnew)
                                                  > numpy.max(xvalues)):
        raise ValueError
    ynew = xnew*0.
    ncalc = len(xnew)
    xmid = numpy.zeros((ncalc+1), dtype=numpy.float32)
    xmid[1:-1] = (xnew[1:] + xnew[0:-1])/2.
    xmid[0] = xmid[1] - (xnew[1]-xnew[0])
    xmid[-1] = xmid[-2] + (xnew[-1] - xnew[-2])
    ymid = numpy.interp(xmid, xvalues, yvalues)
    if not extrapolate:
        inds = numpy.where(xmid < xvalues[0])
        try:
            if ends:
                ymid[inds] = yvalues[0]
            else:
                ymid[inds] = 0.
        except:
            pass
        inds = numpy.where(xmid > xvalues[-1])
        try:
            if ends:
                ymid[inds] = yvalues[-1]
            else:
                ymid[inds] = 0.
        except:
            pass
    if cubic:
        csfit = interpolate.CubicSpline(xvalues, yvalues)
    for loop in range(ncalc):
        inds1 = numpy.where(xvalues <= xmid[loop])
        nstart = inds1[0][-1]
        inds2 = numpy.where(xvalues >= xmid[loop+1])
        nend = inds2[0][0]
        if nstart == nend-1:
            xcalc = (xmid[loop]+xmid[loop+1])/2.
            if cubic:
                yout = csfit(xcalc)
            else:
                yout = numpy.interp(xcalc, xvalues, yvalues)
        else:
            xwork = numpy.copy(xvalues[inds1[0][-1]:inds2[0][0]])
            ywork = numpy.copy(yvalues[inds1[0][-1]:inds2[0][0]])
            if xwork[0] != xmid[loop]:
                xwork = numpy.insert(xwork, 0, xmid[loop])
                ywork = numpy.insert(ywork, 0, ymid[loop])
            if xwork[-1] != xmid[loop+1]:
                xwork = numpy.append(xwork, xmid[loop+1])
                ywork = numpy.append(ywork, ymid[loop+1])
            xrange = xmid[loop+1]-xmid[loop]
            if cubic:
                newxwork = numpy.zeros(len(xwork)*(noversample+1),
                                       dtype=numpy.float32)
                for n1 in range(len(xwork)):
                    newxwork[(noversample+1)*n1] = xwork[n1]
                    if n1+1 < len(xwork):
                        delx = (xwork[n1+1] - xwork[n1])/(noversample+1)
                    for l1 in range(noversample):
                        newxwork[(noversample+1)*n1+l1+1] = \
                            newxwork[(noversample+1)*n1]+(l1+1)*delx
                newxwork[-1] = xwork[-1]
                for l1 in range(noversample):
                    newxwork[(noversample+1)*(len(xwork)-1)+l1+1] = \
                        xwork[-2]+delx*l1
                newywork = csfit(newxwork)
                yout = numpy.trapz(newywork, newxwork)/xrange
            else:
                yout = numpy.trapz(ywork, xwork)/xrange
        ynew[loop] = yout
    return ynew
