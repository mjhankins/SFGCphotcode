#!/usr/bin/env python

#import all required packages
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import subprocess

import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits,ascii
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.table import join
from astropy.coordinates import SkyCoord

from photutils.aperture import SkyCircularAperture,SkyCircularAnnulus,aperture_photometry 
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
from photutils.background import Background2D, MedianBackground, SExtractorBackground, MMMBackground
from photutils.utils import calc_total_error

from regions import read_ds9
from regions import read_ds9, write_ds9, CircleSkyRegion

from photutils.psf import DAOGroup, IntegratedGaussianPRF
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.psf import BasicPSFPhotometry
from photutils.datasets import make_gaussian_sources_image
from astropy.nddata import Cutout2D

import scipy.optimize as opt

#import configuration for selected file
from config import wavelength, segdetsig, finddetsig, bkgbox #import additional common paramters
from config import dpath, dpathalt, ds9path #import additional common paramters

from config import *

interactive=False


#cutout size to use
size=11



def Gaussian2D(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x,y)=xdata_tuple
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2))) + offset
    return g.ravel()

def Gaussian2Dfixedpos(xdata_tuple, amplitude, sigma_x, sigma_y, theta, offset):
    (x,y)=xdata_tuple  
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = amplitude*np.exp( - (a*((x-size/2.0)**2) + 2*b*(x-size/2.0)*(y-size/2.0) + c*((y-size/2.0)**2))) + offset
    return g.ravel()
    
def doGaussian2DModel(Tab,size,plot=False):
    cutouts=[]
    residuals=[]
    varcuts=[]
    params=[]
    for i in range(len(Tab)):

        psloc=(int(np.round(Tab['xcentroid'][i])),int(np.round(Tab['ycentroid'][i])))
        cutout = Cutout2D(data_bkgsub, psloc, (size, size))
        varcut = Cutout2D(varmap, psloc, (size, size))
        # add some noise to the data and try to fit the data generated beforehand
        initial_guess = (Tab["aperture_sum_3.5as"][i]/62.,int(size/2.0),int(size/2.0),4,4,0,0)
        
        #replace any nans with the median value
        cutout.data[np.isnan(cutout.data)]=np.nanmedian(cutout.data)

        # Create x and y indices
        x = np.linspace(0, size-1, size)
        y = np.linspace(0, size-1, size)
        x, y = np.meshgrid(x, y)
        
        #try fit with all parameters allowed to vary
        try:
            popt, pcov = opt.curve_fit(Gaussian2D, (x, y), cutout.data.ravel(), p0=initial_guess)
            qflag=0
        except RuntimeError:
            popt=[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
            
        data_fitted = Gaussian2D((x, y), *popt)
        
        
        #if fit fails then try running with fixed position     
        if popt[0] is np.nan:
            ig2 = (Tab["aperture_sum_3.5as"][i]/62.,4,4,0,0.15)
            try:
                popt, pcov = opt.curve_fit(Gaussian2Dfixedpos, (x, y), cutout.data.ravel(), p0=ig2)
                data_fitted = Gaussian2Dfixedpos((x, y), *popt)
                qflag=1
            except RuntimeError:
                popt=[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
                qflag=3
                
        '''
        xfit=int(size/2.0)-popt[1]
        yfit=int(size/2.0)-popt[2]

        #if it seems like we're at the wrong source try refittingn at center  
        if np.abs(xfit)>6. or np.abs(yfit)>6.:
            ig2 = (Tab["aperture_sum_3.5as"][i]/62.,4,4,0,0.15)
            try:
                popt, pcov = opt.curve_fit(Gaussian2Dfixedpos, (x, y), cutout.data.ravel(), p0=ig2)
                data_fitted = Gaussian2Dfixedpos((x, y), *popt)
                qflag=2
            except RuntimeError:
                #popt=[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
                donothing=1
                qflag=4
        '''
            
        residual = cutout.data - data_fitted.reshape(size, size)
        fiterr=np.sqrt(np.diagonal(pcov)*np.sum(residual**2)/(np.shape(cutout.data.ravel())[0]-5))
        
        if qflag==0:
            
            volume = 2*np.pi*np.sqrt(np.abs(popt[3]))*np.sqrt(np.abs(popt[4]))
            volerr = 2*np.pi*np.sqrt((fiterr[3]/popt[3])+(fiterr[4]/popt[4]))
            xfwhm=np.abs(popt[3])*2.355
            xfwhmerr=np.abs(fiterr[3])*2.355
            yfwhm=np.abs(popt[4])*2.355
            yfwhmerr=np.abs(fiterr[4])*2.355
            fitAmp=popt[0]
            fitAmperr=fiterr[0]
            xfit=int(size/2.0)-popt[1]
            xfiterr=fiterr[1]
            yfit=int(size/2.0)-popt[2]
            yfiterr=fiterr[2]
            theta=popt[5]
            theta_err=fiterr[5]
            dcOffset=popt[6] 
            dco_err=fiterr[6]
        else:

            volume = 2*np.pi*np.sqrt(np.abs(popt[1]))*np.sqrt(np.abs(popt[2]))
            volerr = 2*np.pi*np.sqrt((fiterr[1]/popt[1])+(fiterr[2]/popt[2]))
            xfwhm=np.abs(popt[1])*2.355
            xfwhmerr=np.abs(fiterr[1])*2.355
            yfwhm=np.abs(popt[2])*2.355
            yfwhmerr=np.abs(fiterr[2])*2.355
            fitAmp=popt[0]
            fitAmperr=fiterr[0]
            xfit=0
            xfiterr=0
            yfit=0
            yfiterr=0
            theta=popt[3]
            theta_err=fiterr[3]
            dcOffset=popt[4]
            dco_err=fiterr[4]
        
        if plot and popt[0] is not np.nan:
            fig, ax = plt.subplots(1, 1)
            #ax.hold(True)
            ax.imshow(cutout.data, cmap=plt.cm.jet, origin='lower',
                extent=(x.min(), x.max(), y.min(), y.max()))
            ax.contour(x, y, data_fitted.reshape(size,size), 4, colors='w')
            plt.show()


        model=[volume,volerr,xfwhm,xfwhmerr,yfwhm,yfwhmerr,fitAmp,fitAmperr,xfit,xfiterr,yfit,yfiterr,theta,theta_err,
               dcOffset,dco_err,qflag]
        
        cutouts.append(cutout.data)
        residuals.append(residual)
        varcuts.append(varcut.data)
        params.append(model)
        
        
    return cutouts,residuals,varcuts,np.array(params)



#this is main
for info in field._registry:
    filename=info.filename
    name=info.name
    m1cut=info.m1cut
    m2lims=info.m2lims
    m3lims=info.m3lims
    
    print('\nScript running on field: ', name)
    
    
    try:
        os.chdir(dpath)
    except:
        os.chdir(dpathalt)
    
    #import data - unpack fits file into header, image (data), varance (varmap), and exposure time (tmap)
    hdu=fits.open(filename)
    header=hdu[0].header
    
    #pull the first image plane
    ims=hdu[0].data
    
    #use the first image plane shape to determine how to unpack the rest of the data
    if len(np.shape(ims))==2:
        data=ims
        varmap=hdu[1].data
        tmap=hdu[2].data
    elif len(np.shape(ims))==3:
        data=ims[0]
        varmap=ims[1]
        tmap=ims[2]
        hdu[0].header['NAXIS']=2 #hack to make non-standard WCS work with astropy
    hdu.close()
    
    #define wcs object for header
    wcs=WCS(hdu[0].header)
    
    #create pixel error map by taking sqrt of variance map
    #errormap=np.sqrt(varmap)
    
    #load in the catalog files if they exist
    if os.path.isfile(name+'_'+str(wavelength)+'um_segCat.dat'):
        segTab=ascii.read(name+'_'+str(wavelength)+'um_segCat.dat')
    else:
        segTab=None
    
    if os.path.isfile(name+'_'+str(wavelength)+'um_daoCat.dat'):   
        daoTab=ascii.read(name+'_'+str(wavelength)+'um_daoCat.dat')
    else:
        daoTab=None
        
    #Source coordinates need to be in the form of skycoord objects to create apertures. Ascii tables save as strings so this code puts them back in the right form. 
    scseg=[]
    if segTab is not None:
        sourcecoords=segTab['sky_centroid']
    
        for coord in sourcecoords:
            pos=coord.find(",")
            ra=coord[:pos]
            dec=coord[pos+1:]
            scobj=SkyCoord(ra,dec,unit=u.deg)
            scseg.append(scobj)
    
        segTab['skycoords']=scseg
        sourcesseg=segTab['skycoords']
    

    scdao=[]
    if daoTab is not None:
        sourcecoords=daoTab['sky_centroid']
        for coord in sourcecoords:
            pos=coord.find(",")
            ra=coord[:pos]
            dec=coord[pos+1:]
            scobj=SkyCoord(ra,dec,unit=u.deg)
            scdao.append(scobj)
        
        daoTab['skycoords']=scdao
        sourcesdao=daoTab['skycoords']
    
    
    #check if user defined ds9 file exists
    if os.path.isfile(name+'_'+str(wavelength)+'um_usrCat.dat'):
        usrTab=ascii.read(name+'_'+str(wavelength)+'um_usrCat.dat')
    else:
        usrTab=None
    
    scusr=[]
    if usrTab is not None:
        sourcecoords=usrTab['sky_centroid']
    
        for coord in sourcecoords:
            pos=coord.find(",")
            ra=coord[:pos]
            dec=coord[pos+1:]
            scobj=SkyCoord(ra,dec,unit=u.deg)
            scusr.append(scobj)
    
        usrTab['skycoords']=scusr
        sourcesusr=usrTab['skycoords']
                
        usersources=True
        print('Number of user defined DS9 sources found: ', len(clist))  
    else:
        usersources=False
        print('No user defined DS9 sources found')
        

    #create background model for image using median method
    bkg_estimator = MedianBackground() #MMMBackground() #SExtractorBackground() #MedianBackground()
    bkg_data = Background2D(data,(bkgbox, bkgbox), filter_size=(3, 3),bkg_estimator=bkg_estimator,edge_method='pad')
    bkg_rms=bkg_data.background_rms
    bkg=bkg_data.background 

    #create background subtracted image
    data_bkgsub = data - bkg
    
    
    if segTab is not None:
        #perform psf photometry
        cuts,resid,vcuts,modpar = doGaussian2DModel(segTab,size,plot=interactive)
        
        #add results to table
        segTab['Cutouts']=cuts
        segTab['VarCutout']=vcuts
        segTab['Residuals']=resid
        segTab['PSF_FLUX']=modpar[:,0]
        segTab['PSF_FLUX_ERR']=modpar[:,1]
        segTab['PSF_XFWHM']=modpar[:,2]
        segTab['PSF_XFWHM_ERR']=modpar[:,3]
        segTab['PSF_YFWHM']=modpar[:,4]
        segTab['PSF_yFWHM_ERR']=modpar[:,5]
        segTab['PSF_FitAmp']=modpar[:,6]
        segTab['PSF_FitAmpERR']=modpar[:,7]
        segTab['PSF_XFIT']=modpar[:,8]
        segTab['PSF_XFIT_ERR']=modpar[:,9]
        segTab['PSF_YFIT']=modpar[:,10]
        segTab['PSF_YFIT_ERR']=modpar[:,11]
        segTab['PSF_THETA']=modpar[:,12]
        segTab['PSF_THETA_ERR']=modpar[:,13]
        segTab['PSF_DCOFF']=modpar[:,14]
        segTab['PSF_DCO_ERR']=modpar[:,15]
        segTab['PSF_QFLAG']=modpar[:,16]
        
        #write out the resulting tables to file
        segTab.write(name+'_'+str(wavelength)+'um_segCat_psf.fits',overwrite=True)
        ascii.write(segTab,name+'_'+str(wavelength)+'um_segCat_psf.dat',overwrite=True)

    if daoTab is not None:
        #perform psf photometry
        cuts,resid,vcuts,modpar = doGaussian2DModel(daoTab,size,plot=interactive)
        
        #add results to table
        daoTab['Cutouts']=cuts
        daoTab['Residuals']=resid
        daoTab['VarCutout']=vcuts
        daoTab['PSF_FLUX']=modpar[:,0]
        daoTab['PSF_FLUX_ERR']=modpar[:,1]
        daoTab['PSF_XFWHM']=modpar[:,2]
        daoTab['PSF_XFWHM_ERR']=modpar[:,3]
        daoTab['PSF_YFWHM']=modpar[:,4]
        daoTab['PSF_yFWHM_ERR']=modpar[:,5]
        daoTab['PSF_FitAmp']=modpar[:,6]
        daoTab['PSF_FitAmpERR']=modpar[:,7]
        daoTab['PSF_XFIT']=modpar[:,8]
        daoTab['PSF_XFIT_ERR']=modpar[:,9]
        daoTab['PSF_YFIT']=modpar[:,10]
        daoTab['PSF_YFIT_ERR']=modpar[:,11]
        daoTab['PSF_THETA']=modpar[:,12]
        daoTab['PSF_THETA_ERR']=modpar[:,13]
        daoTab['PSF_DCOFF']=modpar[:,14]
        daoTab['PSF_DCO_ERR']=modpar[:,15]
        daoTab['PSF_QFLAG']=modpar[:,16]
            
        #write out the resulting tables to file
        daoTab.write(name+'_'+str(wavelength)+'um_daoCat_psf.fits',overwrite=True)
        ascii.write(daoTab,name+'_'+str(wavelength)+'um_daoCat_psf.dat',overwrite=True)
    
    if usersources:
        #perform psf photometry
        cuts,resid,vcuts,modpar = doGaussian2DModel(daoTab,size,plot=interactive)

        #add results to table
        usrTab['Cutouts']=cuts
        usrTab['Residuals']=resid
        usrTab['VarCutout']=vcuts
        usrTab['PSF_FLUX']=modpar[:,0]
        usrTab['PSF_FLUX_ERR']=modpar[:,1]
        usrTab['PSF_XFWHM']=modpar[:,2]
        usrTab['PSF_XFWHM_ERR']=modpar[:,3]
        usrTab['PSF_YFWHM']=modpar[:,4]
        usrTab['PSF_yFWHM_ERR']=modpar[:,5]
        usrTab['PSF_FitAmp']=modpar[:,6]
        usrTab['PSF_FitAmpERR']=modpar[:,7]
        usrTab['PSF_XFIT']=modpar[:,8]
        usrTab['PSF_XFIT_ERR']=modpar[:,9]
        usrTab['PSF_YFIT']=modpar[:,10]
        usrTab['PSF_YFIT_ERR']=modpar[:,11]
        usrTab['PSF_THETA']=modpar[:,12]
        usrTab['PSF_THETA_ERR']=modpar[:,13]
        usrTab['PSF_DCOFF']=modpar[:,14]
        usrTab['PSF_DCO_ERR']=modpar[:,15]
        usrTab['PSF_QFLAG']=modpar[:,16]
        
        usrTab.write(name+'_'+str(wavelength)+'um_usrCat_psf.fits',overwrite=True)
        ascii.write(usrTab,name+'_'+str(wavelength)+'um_usrCat_psf.dat',overwrite=True)
    