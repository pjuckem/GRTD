'''
Script to process template and glacial category & coarse frac files for PEST.

This script is to be used within the batch file used to run the model.
Note that matplotlib & arcpy and associated functions have been removed to allow
this script to run on the cluster.  It is expected that this file may be copied into
the Modflow/PEST working directory.

Input and output to this script is entered at the bottom, below
the "if __name__=='__main__':" line.
'''

from __future__ import division

__authors__= 'Brian Clark, Daniel Feinstein, and Paul Juckem'
__maintainer__= 'NAWQA'

import os, csv, sys
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
#import arcpy

class fitpts():
    '''
    fitpts Class.
    Fits an exponential curve through three user-defined points of coarse
    fraction vs property value.
    Curently expects <prop>FunctionInput.csv file present in script directory
    produced by kmaker.calcStats function.

    Parameters
    ----------
    p : string, hk, vk, vani, or por
    inputfile : string <optional>, path and filename of function input values
    outdir : string <optional>, path for saving output figures

    Attributes
    ----------

    Methods
    -------

    See also
    --------

    Notes
    -----

    Examples
    --------
    >>> fit = fitpts('vani')
    >>> coeff = fit.getCoeff()
    >>> fit.graphit()

    '''

    def __init__(self, p, inputfile='', outdir=os.getcwd()):
        self.prop = p
        self.outdir = outdir
        if inputfile == '':
            f = open('{}FunctionInput.csv'.format(self.prop),'r')
        else:
            f = open(inputfile,'r')

        c = csv.reader(f)
        self.lay = []
        self.cat = []
        self.x = []
        self.y = []
        self.popt = []
        self.pcov = []
        self.eqtype = []
        self.mult = []
        self.catidx = dict()
        self.expectedk = []
        c.next() #skip header
        for i,l in enumerate(c):
            self.lay.append(l[0])
            self.cat.append(l[1])
            # self.catidx[int(l[1])]=i
            self.catidx[float(l[0].strip() + '.' + l[1].strip())] = i
            self.x.append([np.float64(z) + 0.000001 for z in l[2:5]])
            self.x[-1] = np.array(self.x[-1])  # generate a list of arrays
            self.y.append([np.float64(z) for z in l[5:]])
            self.y[-1] = np.array(self.y[-1])
            self.expectedk.append(float(l[6]))

            self.mult.append(1.)

        f.close()

    def funcTesting(self,xi, *args):
        # return a *(xi**b)
        # return a * np.log(xi)+b
        # return a * np.exp(b*xi)
        # return a * np.power(xi,b)
        # return a*np.exp(b*xi) + c #exponential variation
        # return args[0] * np.exp(args[1] * xi) + args[2] #exponential variation
        return args[0] * xi + args[2] #linear


# These are for solutionType = 'interp1d' as created by Paul Juckem
    def scipy_interp(self, xi, yi, xnew, kind):  # perform the interpolation to the new values
        if kind == 'cubic':
            if len(xi) <4:
                print ('Use of "cubic" interpolation (by specifying kind=cubic for solverType=interp1d) requires \n'
                'at least 4 pairs of coarse fraction:property values.  We limit this to 3 pairs. \nPlease '
                'select an alternative interpolation method (eg: kind = quadratic provides a curved interpolation.)')
                sys.exit()
        xnew[np.where(xnew < xi.min())] = xi.min()  # No extrapolation!!
        xnew[np.where(xnew > xi.max())] = xi.max()  # No extrapolation!!
        intrp = interp1d(xi, yi, kind)
        return intrp(xnew)

    def curve1(self,xi,yi):  # Alternative function using scipy_interp
        self.eqtype.append(self.scipy_interp)
        poptm = xi
        pcovm = yi  # popt is xvals; pcov is yvals (hijacking variables)
        return (poptm, pcovm)

# These are for solutionType = 'curve_fit' as created by Brian Clark
    def funcPower(self, xi, *args):
        return(args[0] * np.power(xi, args[1]))

    def funcExp(self,xi,*args):
        return args[0] * np.exp(args[1] * xi) + args[2] #exponential variation
        # return args[0] * xi + args[1] #linear

    def funcLinear(self,xi,*args):
        return args[0] * xi + args[1] #linear

    def funcLogLinear(self,xi,*args):
        xi = np.log(xi)
        return args[0] * xi + args[1] #linear

    def curve0(self,xi,yi):  #####  Start processing  #####

        if self.prop == 'por':
            cutoff = 0.00001  # force use of Linear interpolation; other eqs can't handle curvature of input
            #cutoff = 0.0001
            strt = -10
            stp = 10
            nstep = 1
        else:
            cutoff = 1.
            strt = 500
            stp = -500
            nstep = -9

        fres = cutoff + 0.1
        if xi[0] == 0.:
            xi[0] = 0.0001
        xi = np.array(xi)

        if fres > cutoff:
            for pi in range(strt, stp, nstep):
            # for pi in range(1, 2):
                # actually found that seed value matters - and also the increment,
                # so looping over seed values for the first coefficent seems to provide
                # the best fit with an exponential
                if self.prop == 'por':
                    p0 = [1,pi,1]
                else:
                    p0 = [pi,1,1]
                try:
                    poptm, pcov = curve_fit(self.funcExp, xi, yi, maxfev=5000, p0=p0)
                    resids = yi - self.funcExp(xi, *poptm)
                    fres = sum(resids**2)
                    if fres <= cutoff:
                        self.eqtype.append(self.funcExp)
                        break
                except:
                    pass

        if fres > cutoff:
            p0 = [1,1]
            poptm, pcov = curve_fit(self.funcPower, xi, yi, maxfev=3000, p0=p0)
            resids = yi - self.funcPower(xi, *poptm)
            fres = sum(resids**2)
            if fres <= cutoff:
                self.eqtype.append(self.funcPower)        
        
        if fres > cutoff:
            p0 = [1,1]
            poptm, pcov = curve_fit(self.funcLogLinear, xi, yi, maxfev=3000, p0=p0)
            resids = yi - self.funcLogLinear(xi, *poptm)
            fres = sum(resids**2)
            if fres <= cutoff:
                self.eqtype.append(self.funcLogLinear)

        if fres > cutoff:
            # if all else fails, default to linear fit.
            p0 = [1,1]
            poptm, pcov = curve_fit(self.funcLinear, xi, yi, maxfev=3000, p0=p0)
            self.eqtype.append(self.funcLinear)

        return poptm, pcov

    def getCoeff(self):
        for i in range(len(self.x)):
            if solutionType.lower() == 'curvefit':
                if self.y[i][0] == self.y[i][1] == self.y[i][2]:
                    self.popt.append([0.0,self.y[i][0]])
                    self.eqtype.append(self.funcLinear)
                else:
                    self.popt.append(self.curve0(self.x[i],self.y[i])[0])  # Brian's original curve_fit functions
                                                                           # popt is coeff

            elif solutionType.lower() == 'interp1d':
                self.popt.append(self.curve1(self.x[i], self.y[i])[0])  # PFJ's scipy_interp1d function
                self.pcov.append(self.curve1(self.x[i], self.y[i])[1])  # popt is xvals; pcov is yvals (hijacking variables)
            else:
                print('error in the specified solutionType')
                sys.exit()
        return (self.popt, self.pcov, self.mult, self.catidx, self.eqtype)  # everything is now tied to the catindex

    def getExpectedK(self):
        c = []
        c.append([int(z) for z in self.cat[:]])
        c = c[0]
        return dict(zip(c, self.expectedk))
        
    def graphit(self):
        # mpl.rcParams['font.family'] = 'Univers 57 Condensed'
        for i in range(len(self.x)):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.plot(self.x[i],self.y[i], 'ko', label='points', markersize=5)
            x0 = np.linspace(min(self.x[i]), max(self.x[i]),20)
            if self.y[i][0] == self.y[i][1] == self.y[i][2]:
                y0 = [self.y[i][0]] * 20
                plt.plot(x0, y0, 'r-', label="variant exp")
            else:
                if solutionType.lower() == 'curvefit':
                    plt.plot(x0, self.eqtype[i](x0, *self.popt[i]), 'r-', label="variant exp")
                elif solutionType.lower() == 'interp1d':
                    plt.plot(x0, self.eqtype[i](self.x[i], self.y[i], x0, kind), 'r-', label="variant exp")
                else:
                    print('error in the specified solutionType')
                    sys.exit()

            for tick in ax.yaxis.get_major_ticks():
              tick.label1.set_fontsize(10)
            for tick in ax.xaxis.get_major_ticks():
              tick.label1.set_fontsize(10)
            #plt.xlabel('Coarse fraction (multiplier={})'.format(self.mult[i]),size=10)
            plt.xlabel('Coarse fraction)', size=10)
            if self.prop.lower() == 'hk':
                ylab = 'Hydraulic conductivity'
            elif self.prop.lower() == 'vk':
                ylab = 'Vertical Hydraulic conductivity'
            elif self.prop.lower() == 'vani':
                ylab = 'Vertical anisotropy'
            elif self.prop.lower() == 'por':
                ylab = 'Porosity'
            plt.ylabel(ylab, size=10)
            eqsign = ''
            if self.prop == 'vani':
                eqsign = ''
            if self.eqtype[i] == self.funcExp:
                plt.suptitle('{}\nLayer {}\nCategory {}\n$y = {}{:6.4} \\times e^{{{:6.4} \\times x}}+ {:6.4}$'.format(ylab,self.lay[i],
                                                                self.cat[i],eqsign,self.popt[i][0],
                                                                self.popt[i][1],self.popt[i][2]),x=0.30, y=0.87)

            elif self.eqtype[i] == self.funcLinear:
                plt.suptitle('{}\nLayer {}\nCategory {}\n$y = {:6.4}x + {:6.4}$'.format(ylab,self.lay[i],
                                                                self.cat[i],self.popt[i][0],
                                                                self.popt[i][1]),x=0.30, y=0.87)

            elif self.eqtype[i] == self.funcLogLinear:
                plt.suptitle('{}\nLayer {}\nCategory {}\n$y = {:6.4} * log(x) + {:6.4}$'.format(ylab,self.lay[i],
                                                                self.cat[i],self.popt[i][0],
                                                                self.popt[i][1]),x=0.30, y=0.87)
            elif self.eqtype[i] == self.funcPower:
                plt.suptitle('{}\nLayer {}\nCategory {}\n$y = {:6.4} * x^{{:6.4}}$'.format(ylab,self.lay[i],
                                                                self.cat[i],self.popt[i][0],
                                                                self.popt[i][1]), x=0.30, y=0.87)
            elif self.eqtype[i] == self.scipy_interp:
                plt.suptitle('{}\nLayer {}\nCategory {}\n{} Numpy 1D interpolation'.format(ylab, self.lay[i],
                                                                self.cat[i], kind), x=0.30,y=0.87)

            plt.savefig(os.path.join(self.outdir, 'lay{}cat{}prop{}.png'.format(self.lay[i],self.cat[i],self.prop)))
            plt.close()
        print('curve plots saved at ' + self.outdir)

    def plotProps(self,arr,p,ilay):
        plt.imshow(arr)
        plt.colorbar()
        plt.savefig(os.path.join(self.outdir, p+'lay'+str(ilay)+'.png'))
        plt.close()
        print(p + 'lay' + str(ilay) + '.png saved at ' + self.outdir)
    

class pestTools():
    '''
    pestTools Class

    Parameters
    ----------
    nlayg : int, the number of glacial layers.  Obviously, this can't
    exceed the total number of layers.
    modelfc: string, full path and name of the feature class of model grid points with glacial
    and coarse fraction attributes.
    outpath: string, full path of location to create text files of glacial and coarse fraction arrays.
    nrow: integer, number of rows in model.
    ncol: integer, number of columns in model.
    skiplay: not yet supported.
    glacDictDescriptors: dictionary, key/value combinations of glacial unit number and description,
    optional.
    plotarray: boolean, if true, writes png's of curve fit program and pushes other plots
    to the screen.
    savePropArray: boolean, if true, saves output array to text file.

    Attributes
    ----------

    Methods
    -------

    See also
    --------

    Notes
    -----

    Examples
    --------

    '''

    def __init__(self, nlayg, outpath, nrow, ncol, skiplay=[],
                 glacDictDescriptors=None, plotarray=False, savePropArray=True):
        self.nlayg = nlayg
        self.outpath = outpath
        self.nrow = nrow
        self.ncol = ncol
        self.skiplay = skiplay
        self.glacDictDescriptors = glacDictDescriptors
        self.plotarray = plotarray
        self.savePropArray = savePropArray


    def applyKfunctionNumpy(self, propType, prop_file, arrayOut, userProp='', ):
        t = time.time()

        for p in propType:
            print 'Applying functions for {}'.format(p)
            eq = fitpts(p, prop_file)
            if solutionType == 'curvefit':
                coeff, junk, mult, catidx, eqType = eq.getCoeff()
            elif solutionType == 'interp1d':
                xvals, yvals, mult, catidx, eqType = eq.getCoeff()

            if self.plotarray:
                eq.graphit()
            for i in range(self.nlayg):
                coarseArray = np.loadtxt(os.path.join(self.outpath, 'coarse' + str(i+1) + '.dat'))
                glacialArray = np.loadtxt(os.path.join(self.outpath,'glacialNum' + str(i+1) + '.dat'), dtype=np.int)
                glacCats = np.unique(glacialArray) #distinct list of glacial categories
                if all(x != i+1 for x in self.skiplay):  # if not a layer to be skipped, proceed
                    for g in glacCats:
                        if solutionType.lower() == 'curvefit':
                            coarseArray[np.where(glacialArray == g)] = eqType[catidx[float(str(i+1)+'.'+str(g))]](coarseArray[np.where(glacialArray == g)]/100*mult[catidx[float(str(i+1)+'.'+str(g))]],*coeff[catidx[float(str(i+1)+'.'+str(g))]])
                        if solutionType.lower() == 'interp1d':
                            coarseArray[np.where(glacialArray == g)] = eqType[catidx[float(str(i + 1) + '.' + str(g))]](
                                xvals[catidx[float(str(i+1)+'.'+str(g))]], yvals[catidx[float(str(i+1)+'.'+str(g))]],
                                coarseArray[np.where(glacialArray == g)] / 100 * mult[catidx[float(str(i + 1) + '.' + str(g))]],
                                kind)

                else:  # only perform this loop if some layers are specified in the "skip_layer" list.  May want to merge PFJ's script into this else statement.
                    catlist = []
                    proplist = []
                    cklist = 'no'
                    try:
                        prop_input = open(userProp,'r')
                    except:
                        print 'No file with user defined values of min, max, expected properties by ZONE was supplied.'
                        print 'For an example, see ../glackmaker/pykmaker/glacKmaker/propValInput.csv'
                        sys.exit()
                    lines = csv.reader(prop_input)
                    lines.next()
                    for j in lines:
                        if i + 1 == int(j[0]):  # must specify K per category for EACH layer.
                            catlist.append(int(j[1]))
                            if p == 'hk':
                                proplist.append(float(j[3]))
                            elif (p == 'vk' or p == 'vani'):
                                proplist.append(float(j[7]))
                            elif (p == 'por'):
                                proplist.append(float(j[10]))
                            cklist = 'yes'
                    prop_input.close()
                    if cklist == 'yes':
                         catkdict = dict(zip(catlist, proplist))
                    else:
                        print '\nMin K, expected K, and max K per glacial category must specified for EACH layer of the model. \n' \
                              'Please update {}, and re-run the \"applyKfunction\" function.'.format(prop_file)
                        sys.exit()

                    for g in glacCats:
                        coarseArray[np.where(glacialArray == g)] = catkdict[g]  # unpack expected property values based on keys of glac. category


                if self.plotarray:
                    eq.plotProps(coarseArray,p,i+1)
                if self.savePropArray:
                    np.savetxt(os.path.join(arrayOut, p + '_Layer_' + str(i+1) + '.ref'),
                               coarseArray, delimiter=' ')

        print('completed in ' + str(time.time() - t) + ' seconds')


if __name__=='__main__':

    ''' Note that this script is designed to work as part of a PEST run.  It is to be called
    by a batch file, which is run as part of a PEST run.'''

    # path to save text version of glacial and coarse fraction arrays from feature class
    #outpath = '../../2_input/A_modflow/general'
    #outpath = 'D:/PFJData2/Projects/Lake_Team/ANVIL/MODFLOW/working/PEST'
    outpath = '../Models/FWP5L_hK/'

    # pest directory
    #tplDir = '../../4_pest/3G_hetK_HF/pstfiles'
    #tplDir = '../../pstfiles'  # relative to the copied 1_inputscripts/pest directory copied into the worker# directory.
    #tplDir = '../../2_input/A_modflow/KmakerOutput'  # temporary, just to get this working with a single forward run.
    #tplDir = 'D:/PFJData2/Projects/Lake_Team/ANVIL/MODFLOW/working/PEST/'
    tplDir = outpath

    # output directory for external array referenced by the UPW
    # the output file will have include the property abbreviation (hk)
    # and '_Layer_n' where n is the layer and end in .ref - according to flopy structure
    #arrayOut = '../../2_input/A_modflow/KmakerOutput'
    #arrayOut = 'D:/PFJData2/Projects/Lake_Team/ANVIL/MODFLOW/working/PEST/'
    arrayOut = outpath


    # property file(s) generated by kmaker - will be overwritten with tpl file and PEST
    #Kprop_file = os.path.join(tplDir, 'hkFunctionInput.csv')
    #Vprop_file = os.path.join(tplDir, 'vaniFunctionInput.csv')
    #porprop_file = os.path.join(tplDir, 'porFunctionInput.csv')
    porprop_file = os.path.join(tplDir, 'porFunctionInput_SM96_V7.csv')

    skiplay = []
    nlayg = 3
    nrow = 930
    ncol = 650
    #propType = ['hk', 'vani']
    propType = ['por']
    solutionType = 'interp1d'   # Flag for which interpolation method to use.
                                # Options:  1. 'interp1d' = use scipy.interpolate.interp1d as implemented by Paul Juckem.
                                #           2. 'curvefit' = use equation optimization method implemented by Brian Clark.
    # if 'interp1d', can optionally specify the kind of interpolation (linear, quadratic, etc., see the
    # docs for interp1d here: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d  )
    # default is linear.  Can not use "Cubic" with this code, as Cubic requires >3 points, but we limit the CF to property
    # relationship to just 3 pairs.
    kind = 'quadratic'


    kcpest = pestTools(nlayg, outpath, nrow, ncol,
                      skiplay=skiplay, glacDictDescriptors=None, plotarray=True, savePropArray=True)

    # primary method - hopefully to be used with PEST
    #kcpest.applyKfunctionNumpy(['hk'], Kprop_file, arrayOut, userProp='')
    #kcpest.applyKfunctionNumpy(['vani'], Vprop_file, arrayOut, userProp='')
    kcpest.applyKfunctionNumpy(['por'], porprop_file, arrayOut, userProp='')
