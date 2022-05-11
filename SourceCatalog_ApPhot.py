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
from astropy.table import join, Table
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

#import configuration for selected file
from config import wavelength, segdetsig, finddetsig, bkgbox #import additional common paramters
from config import dpath, dpathalt, ds9path #import additional common paramters

from config import *


#Specify options for script

interactive=False

####-------------------------functions-------------------------------------

def performApPhoto(data,tmap,wcs,sourceCoords,radii,rin,rout,plot=True):

	#create aperture objects for all specified radii 
	apertures =[SkyCircularAperture(sourceCoords, r=r*0.786*u.arcsec) for r in radii]

	#do aperture photometry on data using defined apertures 
	phot_table = aperture_photometry(data, apertures,wcs=wcs,error=errormap,method='exact')

	#now estimate local backgrounds using background annulus
	annulus_aperture = SkyCircularAnnulus(sourceCoords, r_in=rin*0.786*u.arcsec, r_out=rout*0.786*u.arcsec) #define annulus

	#convert to pixel coords for calcs and plotting
	pix_aperture = apertures[1].to_pixel(wcs) #only use one of the apertures for this. The default is the 2nd in the list
	pix_annulus_aperture = annulus_aperture.to_pixel(wcs)

	#store area value of annulus in case we want it later
	phot_table['pixAnnArea']=pix_annulus_aperture.area

    
    
    
	#now do robust statistics on the background annuli
	#create lists to store information for later
	bkg_median=[]
	bkg_mean=[]
	bkg_std=[]
	appmasks=[]

	#create mask array for the annuli
	annulus_masks = pix_annulus_aperture.to_mask(method='exact')

	#for each of the annuli loop through and calculate stats using sigma cliped stats
	for mask in annulus_masks:
		annulus_data = mask.multiply(data)
		maskdata=mask.data

		#do statistics
		annulus_data_1d = annulus_data[maskdata > 0]
		meansc, median_sigclip, stdsc = sigma_clipped_stats(annulus_data_1d)
		bkg_median.append(median_sigclip)
		bkg_mean.append(meansc)
		bkg_std.append(stdsc)
		appmasks.append(mask.data)

	#store values in numpy arrays
	bkg_median = np.array(bkg_median)
	bkg_mean = np.array(bkg_mean)
	bkg_std = np.array(bkg_std)

	#add columns for background information and also background subtracted apertures
	phot_table['ann_bkg_med'] = bkg_median
	phot_table['ann_bkg_mean'] = bkg_mean 
	phot_table['ann_bkg_std'] = bkg_std 


    
	#information from exposure time maps
	#create lists to store information for later
	texp_mean=[]
	texp_med=[]
	texpmasks=[]

	#create mask array for the exp time map apertures
	ap_masks = pix_aperture.to_mask(method='exact')

	#for each of the annuli loop through and calculate stats using sigma cliped stats
	for mask in ap_masks:
		ap_texp = mask.multiply(tmap)
        
		#do statistics
		ap_texp_1d = ap_texp[mask.data > 0]
		meansc, median_sigclip, stdsc = sigma_clipped_stats(ap_texp_1d)
		texp_med.append(median_sigclip)
		texp_mean.append(meansc)
		texpmasks.append(mask.data)

	phot_table['texp_med'] = texp_med
	phot_table['texp_mean'] = texp_mean
    


	#caclulate background subtracted photometry and snr for each source
	for i in range(0,len(radii)):
		#coloumn names and from the phot table and new columns we'll add
		cname1='aperture_sum_'+str(i)
		cname2='aperture_sum_err_'+str(i)
		newcol1='aper_sum_bkgsub_'+str(radii[i])+'pix'
		newcol2='aper_snr_'+str(radii[i])+'pix'
		newcol3='aper_area_'+str(radii[i])+'pix'
        
		#get pixel aperture areas for calculations
		pixel_ap = apertures[i].to_pixel(wcs)
		pixarea=pixel_ap.area
		#get background subracted photo by subtracting median annulus value from each pixel in aperture
		phot_table[newcol1]=(phot_table[cname1]/pixarea-phot_table['ann_bkg_med'])*pixarea
		#calculate SNR following equation from forcast photometry cookbook -https://sofia-data-analysis-cookbooks.readthedocs.io/en/latest/FORCAST-photometry_detailed.html
		phot_table[newcol2]=phot_table[newcol1]/np.sqrt(2*np.pi*(radii[i]*phot_table['ann_bkg_std']/phot_table[newcol1])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025) 
		phot_table[newcol3]=pixarea
		
		#rename aperture columns in table to be more descriptive
		rename1='aperture_sum_'+str(radii[i])+'pix'
		rename2='aperture_sum_err_'+str(radii[i])+'pix'
		phot_table.rename_column(cname1, rename1)
		phot_table.rename_column(cname2, rename2)    
    
	#add additonal information for wavelength and which field 
	phot_table['Field']=name
	phot_table['wv']=wavelength
    
    
	if plot:
		#show figure with apertures overlayed
		plt.figure(figsize=(8,8))
		norm = simple_norm(data, 'sqrt', percent=99)
		plt.imshow(data, norm=norm, interpolation='nearest',origin='lower')
		plt.colorbar()

		ap_patches = pix_aperture.plot(color='white', lw=2,
	                          label='Photometry aperture')
		ann_patches = pix_annulus_aperture.plot(color='red', lw=2,
		                                    label='Background annulus')
		handles = (ap_patches[0], ann_patches[0])
		plt.legend(loc='best', facecolor='#458989', labelcolor='white',
		           handles=handles, prop={'weight': 'bold', 'size': 11})
		plt.show()
	return phot_table



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
    
    
    #create background model for image using median method
    bkg_estimator = MedianBackground() #MMMBackground() #SExtractorBackground() #MedianBackground()
    bkg_data = Background2D(data,(bkgbox, bkgbox), filter_size=(5, 5),bkg_estimator=bkg_estimator,edge_method='pad') #smaller box?, 20x20, 25x25?
    bkg_rms=bkg_data.background_rms
    bkg=bkg_data.background 

    #create background subtracted image for photometry
    data_bkgsub = data - bkg
    
    #specify radii to use with source measurements
    radii = [4,7,12] #aperture radii to use in photoemtry - units are pixels
    r_in = 12  #inner radius for background annulus - units are pixels
    r_out = 20  #outer radius for background annulus - units are pixels
    
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
        CombPhotTable=performApPhoto(data_bkgsub,tmap,wcsmap,sourcesAll,radii,r_in,r_out,plot=interactive)

        #display the photometry table
        #CombPhotTable
    else:
        print('No sources found in Combined Source List')
    
    
    #merge Tables
    mtComb = join(combTab, CombPhotTable, keys='id')
    
    #mtComb=doPSFphoto(data_bkgsub,bkg,mtComb,2.0)
    
    #write out catalog 
    mtComb.write(name+'_'+str(wavelength)+'um_CombCat.fits', overwrite=True)

    
    
    
    
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
        
        #mtSeg=doPSFphoto(data_bkgsub,bkg,mtSeg,2.0)
        
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
        
        #mtDao=doPSFphoto(data_bkgsub,bkg,mtDao,2.0)
            
        #write out the resulting tables to file
        mtDao.write(name+'_'+str(wavelength)+'um_daoCat.fits',overwrite=True)

    