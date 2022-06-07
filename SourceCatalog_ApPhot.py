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
from astropy.table import Table, join, vstack
from astropy.wcs import WCS
from astropy import units as u
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D

from photutils.aperture import SkyCircularAperture,SkyCircularAnnulus,aperture_photometry 
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
from photutils.background import Background2D, MedianBackground, SExtractorBackground, MMMBackground
from photutils.utils import calc_total_error
from photutils.morphology import data_properties
from photutils import make_source_mask

from regions import read_ds9
from regions import read_ds9, write_ds9, CircleSkyRegion

#import configuration for selected file
from config import wavelength, bkgbox, cutsize, radii, r_in, r_out #import additional common paramters
from config import dpath, dpathalt, ds9path #import additional common paramters

from config import *

from FORCASTphot import performApPhoto, fitshapes


#Specify options for script

interactive=False


####------------------------Start of main-------------------------------------
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
    wcsmap=WCS(hdu[0].header)
    
    #create pixel error map by taking sqrt of variance map
    errormap=np.sqrt(varmap)
    
    
    #------------------Updated Background estimation method---------------------------
    #create initial background model for building source mask
    bkg_estimator = MMMBackground()  #Alternates -  SExtractorBackground() or MedianBackground()
    bkg_data = Background2D(data,(bkgbox,bkgbox),bkg_estimator=bkg_estimator,edge_method='pad')
    #bkg_rms=bkg_data.background_rms
    bkg=bkg_data.background
    
    tmapnorm=tmap/np.max(tmap) #creating a normalized exposure time map for the mask
    maskTPS=np.where(tmapnorm<0.05,tmapnorm,0).astype('bool')

    #create masked array for the background subtracted data
    data_ma = np.ma.masked_array(data, mask=maskTPS)
    
    mask_3sigma = make_source_mask(data_ma-bkg, nsigma=3, npixels=3, dilate_size=3, filter_fwhm=3)

    data_ma2 = np.ma.masked_array(data, mask=mask_3sigma)
    
    #create updated background model detected sources masked
    bkg_data = Background2D(data_ma2,(bkgbox,bkgbox),bkg_estimator=bkg_estimator,edge_method='pad')
    bkg_rms=bkg_data.background_rms
    bkg=bkg_data.background

    #create background subtracted image
    data_bkgsub = data - bkg
    
    
    #specify radii to use with source measurements
    #radii = [4, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 8, 10, 10.5, 12] #aperture radii to use in photoemtry - units are pixels
    #r_in = 12  #inner radius for background annulus - units are pixels #12 or 15
    #r_out = 20  #outer radius for background annulus - units are pixels #20 or 25
    
    #Start by doing photometry on the combined source list
    #load in source lists if they exist
    if os.path.isfile(name+'_'+str(wavelength)+'um_CombinedSources.fits'):
        combTab=Table.read(name+'_'+str(wavelength)+'um_CombinedSources.fits')
    else:
        combTab=None

    #Get Source coordinates from table
    if combTab is not None:
        sourcesAll=combTab['sky_centroid'] 

    if combTab is not None:
        CombPhotTable=performApPhoto(data_bkgsub,errormap,tmap,header,sourcesAll,radii,r_in,r_out,plot=interactive)
        
        #add additonal information to table
        CombPhotTable['Field']=name
        CombPhotTable['wv']=wavelength

        #display the photometry table
        #CombPhotTable
    else:
        print('No sources found in Combined Source List')
    
    
    #merge Tables
    mtComb = join(combTab, CombPhotTable, keys='id')
    
    #add shape parameters to table
    mtComb=fitshapes(data_bkgsub,mtComb,cutouts=True,cutsize=cutsize) #optional plot=True for diagnostic plots
    
    #write out catalog 
    mtComb.write(name+'_'+str(wavelength)+'um_CombCat.fits', overwrite=True)

    
    
    '''
    
    #------------- Do photometry and processing on other files (optional)---------
    

    #get seg table from detection step
    if os.path.isfile(name+'_'+str(wavelength)+'um_seg.fits'):
        segTab=Table.read(name+'_'+str(wavelength)+'um_seg.fits')
    else:
        segTab=None

    #Get Source coordinates from table
    if segTab is not None:
        sourcesseg=segTab['sky_centroid']

    if segTab is not None:
        SegPhotTable=performApPhoto(data_bkgsub,tmap,wcsmap,sourcesseg,radii,r_in,r_out,plot=interactive)
    
        #display the photometry table
        if interactive:
            print('\nPhotometry table for segmentation map sources.')
            print(SegPhotTable)
            
        #fix id keywords in tables so theu can be merged
        try:
            segTab.rename_column('label', 'id')
        except:
            print('Do nothing because keyword is already changed.')
            
        #merge Tables
        mtSeg = join(segTab, SegPhotTable, keys='id')
        
        #add shape parameters to table
        mtSeg=fitshapes(data_bkgsub,mtSeg) #optional plot=True for diagnostic plots
        
        #write out the resulting tables to file
        mtSeg.write(name+'_'+str(wavelength)+'um_segCat.fits',overwrite=True)

        
        
        
    if os.path.isfile(name+'_'+str(wavelength)+'um_dao.fits'):   
        daoTab=Table.read(name+'_'+str(wavelength)+'um_dao.fits')
    else:
        daoTab=None

    if daoTab is not None:
        sourcesdao=daoTab['sky_centroid']
        
    if daoTab is not None:
        DaoPhotTable=performApPhoto(data_bkgsub,tmap,wcsmap,sourcesdao,radii,r_in,r_out,plot=interactive)
    
        #display the table
        if interactive:
            print('\nPhotometry table for DAOfind sources.')
            print(DaoPhotTable)
            
        #merge Tables
        mtDao = join(daoTab, DaoPhotTable, keys='id')
        
        #add shape parameters to table
        mtDao=fitshapes(data_bkgsub,mtDao) #optional plot=True for diagnostic plots
            
        #write out the resulting tables to file
        mtDao.write(name+'_'+str(wavelength)+'um_daoCat.fits',overwrite=True)

    '''