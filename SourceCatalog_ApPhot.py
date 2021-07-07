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
	phot_table2['aper_sum_bkgsub_6as'] = (phot_table2['aperture_sum']/pix_aperture.area - phot_table2['ann_bkg_med'])* pix_aperture.area 
	phot_table2['aper_sum_bkgsub_err_6as'] = phot_table2['aperture_sum_err'] # should this be modified by bkgsub in some way?
	#not sure if the above is right for the error array...

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
    
	#save these for later use
	phot_table2['pixApArea']=pix_aperture.area
	phot_table2['pixAnnArea']=pix_annulus_aperture.area

	#compute total noise 
	#totalnoise=np.sqrt(sourcenoise+thermalnoise+skynoise) #all noise sources
	#totalnoise=np.sqrt(thermalnoise+skynoise) # no shot noise -> For some reason this seems to give more 'reasonable' values. Need to think about why this is a bit more...
	totalnoise=np.sqrt((thermalnoise+skynoise)*(1+pix_aperture.area/pix_annulus_aperture.area)) #modified to account for pixel stats

	#SNR calc for 6 pixel aperture
	phot_table2['aper_snr_6as']=phot_table2['aper_sum_bkgsub_6as']/totalnoise
    
	merged_table = join(phot_table, phot_table2, keys='id')

	return merged_table

	
def modTabCol(merged_table):
	merged_table.remove_columns(['xcenter_1','ycenter_1','xcenter_2','ycenter_2','sky_center_1','sky_center_2','aperture_sum','aperture_sum_err'])
	
	#rename some columns to avoid possible confusion
	merged_table.rename_column('aperture_sum_0', 'aperture_sum_2as')
	merged_table.rename_column('aperture_sum_err_0', 'aperture_sum_err_2as')
	merged_table.rename_column('aperture_sum_1', 'aperture_sum_4as')
	merged_table.rename_column('aperture_sum_err_1', 'aperture_sum_err_4as')
	merged_table.rename_column('aperture_sum_2', 'aperture_sum_6as')
	merged_table.rename_column('aperture_sum_err_2', 'aperture_sum_err_6as')
	merged_table.rename_column('aperture_sum_3', 'aperture_sum_8as')
	merged_table.rename_column('aperture_sum_err_3', 'aperture_sum_err_8as')
	merged_table.rename_column('aperture_sum_4', 'aperture_sum_10as')
	merged_table.rename_column('aperture_sum_err_4', 'aperture_sum_err_10as')
	
	
	#compute area for the different size apertures 
	ap2area=merged_table['pixApArea']*(2./6.)**2
	ap4area=merged_table['pixApArea']*(4./6.)**2
	#ap6area=merged_table['pixApArea']*(6./6.)**2
	ap8area=merged_table['pixApArea']*(8./6.)**2
	ap10area=merged_table['pixApArea']*(10./6.)**2
	
	
	#calculate local bkg subtracted photometry for the other apertures 
	merged_table['aper_sum_bkgsub_2as']=(merged_table['aperture_sum_2as']/ap2area-merged_table['ann_bkg_med'])*ap2area
	merged_table['aper_sum_bkgsub_4as']=(merged_table['aperture_sum_4as']/ap4area-merged_table['ann_bkg_med'])*ap4area
	#merged_table['aper_sum_bkgsub_6as']=(merged_table['aperture_sum_6as']/ap6area-merged_table['ann_bkg_med'])*ap6area
	merged_table['aper_sum_bkgsub_8as']=(merged_table['aperture_sum_8as']/ap8area-merged_table['ann_bkg_med'])*ap8area
	merged_table['aper_sum_bkgsub_10as']=(merged_table['aperture_sum_10as']/ap10area-merged_table['ann_bkg_med'])*ap10area
	
	#calculate snr for each aperture
	merged_table['aper_snr_2as']=merged_table['aper_sum_bkgsub_2as']/np.sqrt((merged_table['aperture_sum_err_2as']+merged_table['skynoise_pix']*ap2area)*(1+ap2area/merged_table['pixAnnArea']))
	merged_table['aper_snr_4as']=merged_table['aper_sum_bkgsub_4as']/np.sqrt((merged_table['aperture_sum_err_4as']+merged_table['skynoise_pix']*ap4area)*(1+ap4area/merged_table['pixAnnArea']))
	#merged_table['aper_snr_6as']=merged_table['aper_sum_bkgsub_6as']/np.sqrt((merged_table['aperture_sum_err_6as']+merged_table['skynoise_pix']*ap6area)*(1+ap6area/merged_table['pixAnnArea']))
	merged_table['aper_snr_8as']=merged_table['aper_sum_bkgsub_8as']/np.sqrt((merged_table['aperture_sum_err_8as']+merged_table['skynoise_pix']*ap8area)*(1+ap8area/merged_table['pixAnnArea']))
	merged_table['aper_snr_10as']=merged_table['aper_sum_bkgsub_10as']/np.sqrt((merged_table['aperture_sum_err_10as']+merged_table['skynoise_pix']*ap10area)*(1+ap10area/merged_table['pixAnnArea']))
	
	#calculate max snr for all apertures
	snr_values=np.array(merged_table['aper_snr_2as','aper_snr_4as','aper_snr_6as','aper_snr_8as','aper_snr_10as'])
	snr_values.dtype=np.float
	snr_values=np.reshape(snr_values, (-1,5))
	maxsnr=np.nanmax(snr_values,axis=1)
	merged_table['aper_snr_max']=maxsnr
	
	#add additonal information for wavelength and which field 
	merged_table['Field']='C7'+name
	merged_table['wv']=wavelength
	
	#return the table
	return merged_table



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
    if os.path.isfile(name+'_'+str(wavelength)+'um_seg.dat'):
        segTab=ascii.read(name+'_'+str(wavelength)+'um_seg.dat')
    else:
        segTab=None
    
    if os.path.isfile(name+'_'+str(wavelength)+'um_dao.dat'):   
        daoTab=ascii.read(name+'_'+str(wavelength)+'um_dao.dat')
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
        
    radii = [2,4,6,8,10] #define aperture radii & construct apertures (line below)
    
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
        
        mtSeg=modTabCol(merged_table_seg)
        
        #write out the resulting tables to file
        ascii.write(mtSeg, name+'_'+str(wavelength)+'um_segCat.dat', overwrite=True)

    if daoTab is not None:
        DaoPhotTable=performApPhoto(data,wcsmap,sourcesdao,radii,plot=interactive)
    
        #display the table
        if interactive:
            print('\nPhotometry table for DAOfind sources.')
            print(DaoPhotTable)
            
        #merge Tables
        merged_table_dao = join(daoTab, DaoPhotTable, keys='id')
        
        mtDao=modTabCol(merged_table_dao)
            
        #write out the resulting tables to file
        ascii.write(mtDao, name+'_'+str(wavelength)+'um_daoCat.dat', overwrite=True)
    
    if usersources:
        UserPhotTable=performApPhoto(data,ds9sc,radii,plot=interactive)
    
        #display the table
        if interactive:
            print('\nPhotometry table for user defined sources in ds9.')
            print(UserPhotTable)

        mtds9=modTabCol(UserPhotTable)
 
        ascii.write(mtds9, name+'_'+str(wavelength)+'um_usrCat.dat', overwrite=True)
    