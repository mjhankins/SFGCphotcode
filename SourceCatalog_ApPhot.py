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

interactive=False


def performApPhoto(data,wcs,sourceCoords,radii,plot=True):
	#first let's do some simple annulus extractions...   
	apertures =[SkyCircularAperture(sourceCoords, r=r*u.arcsec) for r in radii] # for pixels: r=r*u.arcsec*0.786

	#do aperture photometry on data using defined apertures 
	phot_table = aperture_photometry(data, apertures,wcs=wcs,error=errormap,method='exact')

	#display phot table
	#phot_table

	#now try photometry with local background subtraction
	aperture2 =SkyCircularAperture(sourceCoords, r=6*u.arcsec) #define aperture
	annulus_aperture = SkyCircularAnnulus(sourceCoords, r_in=10*u.arcsec, r_out=16*u.arcsec) #define annulus

	#convert to pixel coords for plotting
	pix_aperture = aperture2.to_pixel(wcs)
	pix_annulus_aperture = annulus_aperture.to_pixel(wcs)
    
	#print(pix_aperture)

	if plot:
		#show figure with apertures overlayed
		plt.figure(figsize=(8,8))
		norm = simple_norm(data, 'sqrt', percent=99)
		plt.imshow(data, norm=norm, interpolation='nearest',origin='lower')
		plt.colorbar()
		#plt.xlim(40, 140)
		#plt.ylim(50, 125)

		ap_patches = pix_aperture.plot(color='white', lw=2,
 	                          label='Photometry aperture')
		ann_patches = pix_annulus_aperture.plot(color='red', lw=2,
		                                    label='Background annulus')
		handles = (ap_patches[0], ann_patches[0])
		plt.legend(loc='best', facecolor='#458989', labelcolor='white',
		           handles=handles, prop={'weight': 'bold', 'size': 11})
		plt.show()

	#now lets do robust statistics on the background annuli


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
		#print(np.shape(mask))
        
		#this is a bit of debugging to handle if the mask array is the wrong shape
		if np.shape(mask.data)[0]==41:
			maskdata=mask.data[:-1,:]
		else:
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

	#do aperture photometry
	phot_table2 = aperture_photometry(data, aperture2,wcs=wcs,error=errormap,method='exact') #

	#add columns for background information and also background subtracted apertures
	phot_table2['ann_bkg_med'] = bkg_median
	phot_table2['ann_bkg_mean'] = bkg_mean 
	phot_table2['ann_bkg_std'] = bkg_std 
	phot_table2['aper_sum_bkgsub_6as'] = (phot_table2['aperture_sum']/pix_aperture.area - phot_table2['ann_bkg_med'])* pix_aperture.area 
	phot_table2['aper_sum_bkgsub_err_6as'] = phot_table2['aperture_sum_err'] # should this be modified by bkgsub in some way?
	#not sure if the above is right for the error array...



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
		#print(np.shape(mask))
        
		##this is a bit of debugging to handle if the mask array is the wrong shape
		#if np.shape(mask.data)[0]==41:
		#	maskdata=mask.data[:-1,:]
		#else:
		#	maskdata=mask.data
        
		#do statistics
		ap_texp_1d = ap_texp[mask.data > 0]
		meansc, median_sigclip, stdsc = sigma_clipped_stats(ap_texp_1d)
		texp_med.append(median_sigclip)
		texp_mean.append(meansc)
		texpmasks.append(mask.data)

    
	phot_table2['texp_med'] = texp_med
	phot_table2['texp_mean'] = texp_mean
    
    
	'''
	#calculate sky noise for 6 pixel aperture
	skynoise=np.sqrt(phot_table2['ann_bkg_med']*pix_aperture.area)
	sna=np.array(skynoise)
	masknan=np.isnan(sna)
	sna[masknan]=0.0
	skynoise=sna**2

	#store per pix sky noise for later
	phot_table2['skynoise_pix']=skynoise/pix_aperture.area

	#shot noise from the source
	sourcenoise=phot_table2['aper_sum_bkgsub_6as']

	#thermal noise from camera (from error map)
	thermalnoise=phot_table2['aperture_sum_err']
    


	#compute total noise 
	#totalnoise=np.sqrt(sourcenoise+thermalnoise+skynoise) #all noise sources
	#totalnoise=np.sqrt(thermalnoise+skynoise) # no shot noise -> For some reason this seems to give more 'reasonable' values. Need to think about why this is a bit more...
	totalnoise=np.sqrt((thermalnoise+skynoise)*(1+pix_aperture.area/pix_annulus_aperture.area)) #modified to account for pixel stats
	'''  
    
	r_ap=6.0/0.768 #aperture in arcsecs converted to pixels
    
	#from https://sofia-data-analysis-cookbooks.readthedocs.io/en/latest/FORCAST-photometry_detailed.html
	#totalnoise=np.sqrt(2*np.pi*(np.nanmax(tmap)/phot_table2['texp_mean'])*(r_ap*imgbkg/phot_table2['aper_sum_bkgsub_6as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
    
	totalnoise=np.sqrt(2*np.pi*(r_ap*phot_table2['ann_bkg_std']/phot_table2['aper_sum_bkgsub_6as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)

    
	#SNR calc for 6 pixel aperture
	phot_table2['aper_snr_6as']=phot_table2['aper_sum_bkgsub_6as']/totalnoise

    
	#save these for later use
	phot_table2['pixApArea']=pix_aperture.area
	phot_table2['pixAnnArea']=pix_annulus_aperture.area
    
    
	merged_table = join(phot_table, phot_table2, keys='id')

	return merged_table

	
def modTabCol2(merged_table):
	merged_table.remove_columns(['xcenter_1','ycenter_1','xcenter_2','ycenter_2','sky_center_1','sky_center_2','aperture_sum','aperture_sum_err'])
	
	#rename some columns to avoid possible confusion
	merged_table.rename_column('aperture_sum_0', 'aperture_sum_3.5as')
	merged_table.rename_column('aperture_sum_err_0', 'aperture_sum_err_3.5as')
	merged_table.rename_column('aperture_sum_1', 'aperture_sum_3.75as')
	merged_table.rename_column('aperture_sum_err_1', 'aperture_sum_err_3.75as')
	merged_table.rename_column('aperture_sum_2', 'aperture_sum_4.0as')
	merged_table.rename_column('aperture_sum_err_2', 'aperture_sum_err_4.0as')
	merged_table.rename_column('aperture_sum_3', 'aperture_sum_4.25as')
	merged_table.rename_column('aperture_sum_err_3', 'aperture_sum_err_4.25as')
	merged_table.rename_column('aperture_sum_4', 'aperture_sum_4.5as')
	merged_table.rename_column('aperture_sum_err_4', 'aperture_sum_err_4.5as')
	merged_table.rename_column('aperture_sum_5', 'aperture_sum_4.75as')
	merged_table.rename_column('aperture_sum_err_5', 'aperture_sum_err_4.75as')
	merged_table.rename_column('aperture_sum_6', 'aperture_sum_5.0as')
	merged_table.rename_column('aperture_sum_err_6', 'aperture_sum_err_5.0as')
	merged_table.rename_column('aperture_sum_7', 'aperture_sum_5.25as')
	merged_table.rename_column('aperture_sum_err_7', 'aperture_sum_err_5.25as')
	merged_table.rename_column('aperture_sum_8', 'aperture_sum_5.5as')
	merged_table.rename_column('aperture_sum_err_8', 'aperture_sum_err_5.5as')
	merged_table.rename_column('aperture_sum_9', 'aperture_sum_5.75as')
	merged_table.rename_column('aperture_sum_err_9', 'aperture_sum_err_5.75as')
	merged_table.rename_column('aperture_sum_10', 'aperture_sum_6.0as')
	merged_table.rename_column('aperture_sum_err_10', 'aperture_sum_err_6.0as') 
	merged_table.rename_column('aperture_sum_11', 'aperture_sum_10as')
	merged_table.rename_column('aperture_sum_err_11', 'aperture_sum_err_10as')
	
	
	#compute area for the different size apertures 
	ap350area=merged_table['pixApArea']*(3.5/6.)**2
	ap375area=merged_table['pixApArea']*(3.75/6.)**2
	ap400area=merged_table['pixApArea']*(4.0/6.)**2
	ap425area=merged_table['pixApArea']*(4.25/6.)**2
	ap450area=merged_table['pixApArea']*(4.5/6.)**2
	ap475area=merged_table['pixApArea']*(4.75/6.)**2
	ap500area=merged_table['pixApArea']*(5.0/6.)**2
	ap525area=merged_table['pixApArea']*(5.25/6.)**2
	ap550area=merged_table['pixApArea']*(5.5/6.)**2
	ap575area=merged_table['pixApArea']*(5.75/6.)**2
	ap600area=merged_table['pixApArea']*(6.0/6.)**2
	ap10area=merged_table['pixApArea']*(10./6.)**2
	
	
	#calculate local bkg subtracted photometry for the other apertures 
	merged_table['aper_sum_bkgsub_3.5as']=(merged_table['aperture_sum_3.5as']/ap350area-merged_table['ann_bkg_med'])*ap350area
	merged_table['aper_sum_bkgsub_3.75as']=(merged_table['aperture_sum_3.75as']/ap375area-merged_table['ann_bkg_med'])*ap375area
	merged_table['aper_sum_bkgsub_4.0as']=(merged_table['aperture_sum_4.0as']/ap400area-merged_table['ann_bkg_med'])*ap400area
	merged_table['aper_sum_bkgsub_4.25as']=(merged_table['aperture_sum_4.25as']/ap425area-merged_table['ann_bkg_med'])*ap425area
	merged_table['aper_sum_bkgsub_4.5as']=(merged_table['aperture_sum_4.5as']/ap450area-merged_table['ann_bkg_med'])*ap450area
	merged_table['aper_sum_bkgsub_4.75as']=(merged_table['aperture_sum_4.75as']/ap475area-merged_table['ann_bkg_med'])*ap475area
	merged_table['aper_sum_bkgsub_5.0as']=(merged_table['aperture_sum_5.0as']/ap500area-merged_table['ann_bkg_med'])*ap500area
	merged_table['aper_sum_bkgsub_5.25as']=(merged_table['aperture_sum_5.25as']/ap525area-merged_table['ann_bkg_med'])*ap525area
	merged_table['aper_sum_bkgsub_5.5as']=(merged_table['aperture_sum_5.5as']/ap550area-merged_table['ann_bkg_med'])*ap550area
	merged_table['aper_sum_bkgsub_5.75as']=(merged_table['aperture_sum_5.75as']/ap575area-merged_table['ann_bkg_med'])*ap575area
	merged_table['aper_sum_bkgsub_6.0as']=(merged_table['aperture_sum_6.0as']/ap600area-merged_table['ann_bkg_med'])*ap600area
	merged_table['aper_sum_bkgsub_10as']=(merged_table['aperture_sum_10as']/ap10area-merged_table['ann_bkg_med'])*ap10area
	'''
	#calculate snr for each aperture
	merged_table['aper_snr_3.5as']=merged_table['aper_sum_bkgsub_3.5as']/np.sqrt((merged_table['aperture_sum_err_3.5as']+merged_table['skynoise_pix']*ap350area)*(1+ap350area/merged_table['pixAnnArea']))
	merged_table['aper_snr_3.75as']=merged_table['aper_sum_bkgsub_3.75as']/np.sqrt((merged_table['aperture_sum_err_3.75as']+merged_table['skynoise_pix']*ap375area)*(1+ap375area/merged_table['pixAnnArea']))
	merged_table['aper_snr_4.0as']=merged_table['aper_sum_bkgsub_4.0as']/np.sqrt((merged_table['aperture_sum_err_4.0as']+merged_table['skynoise_pix']*ap400area)*(1+ap400area/merged_table['pixAnnArea']))
	merged_table['aper_snr_4.25as']=merged_table['aper_sum_bkgsub_4.25as']/np.sqrt((merged_table['aperture_sum_err_4.25as']+merged_table['skynoise_pix']*ap425area)*(1+ap425area/merged_table['pixAnnArea']))
	merged_table['aper_snr_4.5as']=merged_table['aper_sum_bkgsub_4.5as']/np.sqrt((merged_table['aperture_sum_err_4.5as']+merged_table['skynoise_pix']*ap450area)*(1+ap450area/merged_table['pixAnnArea']))
	merged_table['aper_snr_4.75as']=merged_table['aper_sum_bkgsub_4.75as']/np.sqrt((merged_table['aperture_sum_err_4.75as']+merged_table['skynoise_pix']*ap475area)*(1+ap475area/merged_table['pixAnnArea']))
	merged_table['aper_snr_5.0as']=merged_table['aper_sum_bkgsub_5.0as']/np.sqrt((merged_table['aperture_sum_err_5.0as']+merged_table['skynoise_pix']*ap500area)*(1+ap500area/merged_table['pixAnnArea']))
	merged_table['aper_snr_5.25as']=merged_table['aper_sum_bkgsub_5.25as']/np.sqrt((merged_table['aperture_sum_err_5.25as']+merged_table['skynoise_pix']*ap525area)*(1+ap525area/merged_table['pixAnnArea']))
	merged_table['aper_snr_5.5as']=merged_table['aper_sum_bkgsub_5.5as']/np.sqrt((merged_table['aperture_sum_err_5.5as']+merged_table['skynoise_pix']*ap550area)*(1+ap550area/merged_table['pixAnnArea']))
	merged_table['aper_snr_5.75as']=merged_table['aper_sum_bkgsub_5.75as']/np.sqrt((merged_table['aperture_sum_err_5.75as']+merged_table['skynoise_pix']*ap575area)*(1+ap575area/merged_table['pixAnnArea']))
	merged_table['aper_snr_6.0as']=merged_table['aper_sum_bkgsub_6.0as']/np.sqrt((merged_table['aperture_sum_err_6.0as']+merged_table['skynoise_pix']*ap600area)*(1+ap600area/merged_table['pixAnnArea']))
	merged_table['aper_snr_10as']=merged_table['aper_sum_bkgsub_10as']/np.sqrt((merged_table['aperture_sum_err_10as']+merged_table['skynoise_pix']*ap10area)*(1+ap10area/merged_table['pixAnnArea']))
	'''
    
    
	#calculate snr for each aperture - updated using https://sofia-data-analysis-cookbooks.readthedocs.io/en/latest/FORCAST-photometry_detailed.html
	merged_table['aper_snr_3.5as']=merged_table['aper_sum_bkgsub_3.5as']/np.sqrt(2*np.pi*((3.5/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_3.5as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025) 
	merged_table['aper_snr_3.75as']=merged_table['aper_sum_bkgsub_3.75as']/np.sqrt(2*np.pi*((3.75/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_3.75as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_4.0as']=merged_table['aper_sum_bkgsub_4.0as']/np.sqrt(2*np.pi*((4.0/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_4.0as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_4.25as']=merged_table['aper_sum_bkgsub_4.25as']/np.sqrt(2*np.pi*((4.25/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_4.25as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_4.5as']=merged_table['aper_sum_bkgsub_4.5as']/np.sqrt(2*np.pi*((4.5/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_4.5as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_4.75as']=merged_table['aper_sum_bkgsub_4.75as']/np.sqrt(2*np.pi*((4.75/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_4.75as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_5.0as']=merged_table['aper_sum_bkgsub_5.0as']/np.sqrt(2*np.pi*((5.0/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_5.0as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_5.25as']=merged_table['aper_sum_bkgsub_5.25as']/np.sqrt(2*np.pi*((5.25/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_5.25as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_5.5as']=merged_table['aper_sum_bkgsub_5.5as']/np.sqrt(2*np.pi*((5.5/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_5.5as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_5.75as']=merged_table['aper_sum_bkgsub_5.75as']/np.sqrt(2*np.pi*((5.75/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_5.75as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_6.0as']=merged_table['aper_sum_bkgsub_6.0as']/np.sqrt(2*np.pi*((6.0/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_6.0as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	merged_table['aper_snr_10as']=merged_table['aper_sum_bkgsub_10as']/np.sqrt(2*np.pi*((10.0/0.768)*merged_table['ann_bkg_std']/merged_table['aper_sum_bkgsub_10as'])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025)
	
    
	#calculate max snr for all apertures
	snr_values=np.array(merged_table['aper_snr_3.5as','aper_snr_3.75as','aper_snr_4.0as','aper_snr_4.25as','aper_snr_4.5as','aper_snr_4.75as','aper_snr_5.0as','aper_snr_5.25as','aper_snr_5.5as','aper_snr_5.75as','aper_snr_6.0as','aper_snr_10as'])
	snr_values.dtype=np.float
	snr_values=np.reshape(snr_values, (-1,12))
	maxsnr=np.nanmax(snr_values,axis=1)
	merged_table['aper_snr_max']=maxsnr
	
	#add additonal information for wavelength and which field 
	merged_table['Field']='C7'+name
	merged_table['wv']=wavelength
	
	#display table
	return merged_table

#new addition to include FWHM measurements for DAO
def doPSFphoto(image,bkgmodel,sourceTable,sigma_init,plotting=False): 
	 #create initial guess positions for fitting routine
	
	initTab=Table()
	initTab['x_0']=sourceTable['xcentroid']
	initTab['y_0']=sourceTable['ycentroid']
	initTab['flux']=sourceTable['aper_sum_bkgsub_5.0as']
	
	try:
		initTab['sigma_0']=sourceTable['fwhm']
	except:
		initTab['sigma_0']=2.5
	
	daogroup = DAOGroup(crit_separation=8)
	mmm_bkg = MedianBackground() #MedianBackground()#MMMBackground()
	fitter = LevMarLSQFitter()
	gaussian_prf = IntegratedGaussianPRF(sigma=sigma_init)
	gaussian_prf.sigma.fixed = False
	mmm_bkg.sigma_clip.fixed = False
	
	basic_phot_obj=BasicPSFPhotometry(group_maker=daogroup, bkg_estimator=mmm_bkg, psf_model=gaussian_prf, fitter=fitter, fitshape=(11, 11))
	
	results = basic_phot_obj(image,init_guesses=initTab) #must provide initial guesses as an astropy table with columns x_0 and y_0 in pixel coords
	
	#results=photB_results[photB_results['flux_fit']>0]
	#results=results[(results['sigma_fit']>0) & (results['sigma_fit']<14) ]
	
	sourceTable['PSF_fwhm']=results['sigma_fit']*2.355
	#sourceTable['PSF_fwhm_unc']=results['sigma_unc']*2.355 #note sigma_unc isn't consistently returned for some reason...
	sourceTable['PSF_Flux_1D']=results['flux_fit']
	#sourceTable['PSF_Flux_unc_1D']=results['flux_unc']
	#sourceTable['psfSNR']=results['flux_fit']/results['flux_unc']
	
	'''
	#construct model residuals
	sources = Table()
	sources['flux'] = results['flux_fit']
	sources['x_mean'] = results['x_fit']
	sources['y_mean'] = results['y_fit']
	sources['x_stddev'] = results['sigma_fit']
	sources['y_stddev'] = sources['x_stddev']
	#sources['theta'] = [0] * 2
	
	modelimage = make_gaussian_sources_image(np.shape(data), sources)
	
	#residual=data-modelimage
	residual=data-(modelimage+bkgmodel) #include bkg model   
	    
	chivals=[]
	
	for source in results:
		spos=(np.int(source['x_fit']),np.int(source['y_fit']))
		cutout=Cutout2D(residual,spos,15)
		#print(np.sum(cutout.data))
		fudge=100.
		chivals.append((np.nansum(cutout.data)/21.)**2*100)
		
	sourceTable['psfFitChi']=chivals
	#print(chivals)
	
	#rename flux initial guess in table 
	results.rename_column('flux', 'flux_init')
	
	if plotting==True:
		#a few diagnostic plots to examine background contributions to aperture flux
		fig, (ax1,ax2) = plt.subplots(1, 2,figsize=(15,5.5))
		
		ax1.set_title('Data')
		p1=ax1.imshow(data,origin='lower')
		fig.colorbar(p1, ax=ax1)
		
		
		ax2.set_title('Residual Map from DAOfind')
		p2=ax2.imshow(DAOresidual,origin='lower',vmin=-0.2,vmax=0.2)
		fig.colorbar(p2, ax=ax2)
	    
	'''
	#return photB_results, residual#, chivals
	return sourceTable

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
    wcsmap=WCS(hdu[0].header)
    
    #create pixel error map by taking sqrt of variance map
    errormap=np.sqrt(varmap)
    
    #load in the catalog files if they exist
    if os.path.isfile(name+'_'+str(wavelength)+'um_seg.fits'):
        #segTab=ascii.read(name+'_'+str(wavelength)+'um_seg.dat')
        segTab=Table.read(name+'_'+str(wavelength)+'um_seg.fits')
    else:
        segTab=None
    
    if os.path.isfile(name+'_'+str(wavelength)+'um_dao.fits'):   
        #daoTab=ascii.read(name+'_'+str(wavelength)+'um_dao.dat')
        daoTab=Table.read(name+'_'+str(wavelength)+'um_dao.fits')
    else:
        daoTab=None
        
    '''    
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
    '''
    #get source coordinates
    if segTab is not None:
        sourcesseg=segTab['sky_centroid']
        
    if daoTab is not None:
        sourcesdao=daoTab['sky_centroid']
    
    #check if user defined ds9 file exists
    if os.path.isfile(name+'_ds9.reg'):
        sourcesDS9=read_ds9(name+'_ds9.reg')
    
        clist=[]
        
        for source in sourcesDS9:
            if source.visual['color']=='green':
                sc=source.center
                ra=sc.ra.value
                dec=sc.dec.value
                a=(ra,dec)
                clist.append(a)
                
        ds9sc=SkyCoord(clist,unit=u.deg)
        usersources=True
        print('Number of user defined DS9 sources found: ', len(clist))  
    else:
        usersources=False
        print('No user defined DS9 sources found')
        
    #radii = [2,4,6,8,10] #define aperture radii & construct apertures (line below) - old
    radii = [3.5,3.75,4.0,4.25,4.5,4.75,5,5.25,5.5,5.75,6.0,10] #define aperture radii & construct apertures (line below) - new

    #create background model for image using median method
    bkg_estimator = MedianBackground() #MMMBackground() #SExtractorBackground() #MedianBackground()
    bkg_data = Background2D(data,(bkgbox, bkgbox), filter_size=(3, 3),bkg_estimator=bkg_estimator,edge_method='pad')
    bkg_rms=bkg_data.background_rms
    bkg=bkg_data.background 

    #create background subtracted image
    data_bkgsub = data - bkg
    
    
    if segTab is not None:
        SegPhotTable=performApPhoto(data,wcsmap,sourcesseg,radii,plot=interactive)
    
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
        merged_table_seg = join(segTab, SegPhotTable, keys='id')
        
        mtSeg=modTabCol2(merged_table_seg)
        
        #mtSeg=doPSFphoto(data_bkgsub,bkg,mtSeg,2.0)
        
        #write out the resulting tables to file
        #ascii.write(mtSeg, name+'_'+str(wavelength)+'um_segCat.dat', overwrite=True)
        mtSeg.write(name+'_'+str(wavelength)+'um_segCat.fits',overwrite=True)


    if daoTab is not None:
        DaoPhotTable=performApPhoto(data,wcsmap,sourcesdao,radii,plot=interactive)
    
        #display the table
        if interactive:
            print('\nPhotometry table for DAOfind sources.')
            print(DaoPhotTable)
            
        #merge Tables
        merged_table_dao = join(daoTab, DaoPhotTable, keys='id')
        
        mtDao=modTabCol2(merged_table_dao)
        
        #mtDao=doPSFphoto(data_bkgsub,bkg,mtDao,2.0)
            
        #write out the resulting tables to file
        #ascii.write(mtDao, name+'_'+str(wavelength)+'um_daoCat.dat', overwrite=True)
        mtDao.write(name+'_'+str(wavelength)+'um_daoCat.fits',overwrite=True)

    
    if usersources:
        UserPhotTable=performApPhoto(data,ds9sc,radii,plot=interactive)
    
        #display the table
        if interactive:
            print('\nPhotometry table for user defined sources in ds9.')
            print(UserPhotTable)

        mtds9=modTabCol2(UserPhotTable)
        
        #mtds9=doPSFphoto(data_bkgsub,bkg,mtds9,2.0)
 
        #ascii.write(mtds9, name+'_'+str(wavelength)+'um_usrCat.dat', overwrite=True)
        usrCat.write(name+'_'+str(wavelength)+'um_usrCat.fits', overwrite=True)
    