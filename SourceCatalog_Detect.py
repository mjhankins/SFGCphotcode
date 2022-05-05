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
from astropy.coordinates import SkyCoord,Angle 
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord,skycoord_to_pixel
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.table import join, Table

from photutils.aperture import SkyCircularAperture,SkyCircularAnnulus,aperture_photometry, CircularAperture
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
from photutils.background import Background2D, MedianBackground, SExtractorBackground, MMMBackground
from photutils.utils import calc_total_error
from photutils import DAOStarFinder,IRAFStarFinder


from regions import read_ds9, write_ds9, CircleSkyRegion

#import configuration for selected file
from config import wavelength, segdetsig, finddetsig, bkgbox #import additional common paramters
from config import dpath, dpathalt, ds9path #import additional common paramters

from config import *

interactive=False

for info in field._registry:
    name=info.name
    filename=info.filename
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
    
    
    '''
    #take a quick look at the maps that were loaded in 
    plt.figure()
    plt.title('Data')
    plt.imshow(data,origin='lower',interpolation='none')
    plt.colorbar()
    #plt.clim(0.0,0.1)
    plt.show()
    
    plt.figure()
    plt.title('Variance Map')
    plt.imshow(varmap,origin='lower',interpolation='none')
    plt.colorbar()
    plt.show()
    
    plt.figure()
    plt.title('Exposure Time Map')
    plt.imshow(tmap,origin='lower',interpolation='none')
    plt.colorbar()
    plt.show()
    '''
    
    #create pixel error map by taking sqrt of variance map
    errormap=np.sqrt(varmap)
    
    #create mask for edges of field where less integration time was collected
    tmapnorm=tmap/np.max(tmap) #normalize the exposure time map
    mask=np.where(tmapnorm<m1cut,tmapnorm,0).astype('bool') #create mask for any locations with less than 50% of max exposure time -Can be modified as needed
    #the above can be adjusted if there are obvious sources near the edges of the map
    
    '''
    #plot any of the mask to verify
    plt.figure()
    plt.imshow(mask,origin='lower',interpolation='none')
    plt.colorbar()
    plt.show()
    '''
    
    #create background model for image using median method
    bkg_estimator = MedianBackground() #MMMBackground() #SExtractorBackground() #MedianBackground()
    bkg_data = Background2D(data,(bkgbox, bkgbox), filter_size=(3, 3),bkg_estimator=bkg_estimator,edge_method='pad') #smaller box?, 20x20, 25x25?
    bkg_rms=bkg_data.background_rms
    bkg=bkg_data.background 
    
    #create background subtracted image
    data_bkgsub = data - bkg
    
    #set detection threshold for source finding based on modeled background rms
    threshold = segdetsig*bkg_rms
    
    '''
    #plot the data and background model
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3,figsize=(14,10))
    ax1.set_title('Data')
    
    #set the image limits for the plots
    minval=-0.1
    maxval=0.1
    
    ax1.imshow(data,origin='lower',vmin=minval,vmax=maxval)
    ax2.set_title('Background Model')
    ax2.imshow(bkg,origin='lower',vmin=minval,vmax=maxval)
    ax3.set_title('Data - Bkg. Model')
    ax3.imshow(data_bkgsub,origin='lower',vmin=minval,vmax=maxval)
    plt.show()
    '''

    
    
    #do source detection. A 3x3 FWHM gaussian is used to smooth image prior to detection
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data_bkgsub, threshold, mask=mask, npixels=5, kernel=kernel)
    
    #removed labels that exist in masked region
    if m2lims is not None:
        mask2=np.zeros(np.shape(mask))
        for lim in m2lims:
            mask2[lim[0]:lim[1],lim[2]:lim[3]]=1
        segm.remove_masked_labels(mask2.astype('bool'))

    #lets take a look at deblending sources
    if segm is not None:
        segm_deblend = deblend_sources(data_bkgsub, segm, npixels=5,kernel=kernel, nlevels=64,contrast=0.001)
    else:
        segm_deblend=None #if segm is empty pass it on to avoid errors

    #remove any sources that should be masked (mask3) 
    if m3lims is not None:
        mask3=np.zeros(np.shape(mask))
        for lim in m3lims:
            mask3[lim[0]:lim[1],lim[2]:lim[3]]=1
        segm_deblend.remove_masked_labels(mask3.astype('bool'))
    
    '''
    #make plot of segmentation image to show detected sources side by side with data
    norm = ImageNormalize(stretch=SqrtStretch())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
    ax1.imshow(data_bkgsub, origin='lower', cmap='Greys_r', norm=norm)
    ax1.set_title('Data')
    cmap = segm.make_cmap(seed=123)
    ax2.imshow(segm, origin='lower', cmap=cmap, interpolation='nearest')
    ax2.set_title('Deblend Seg. Image')
    plt.show()
    '''

    
    #now lets look at building a catalog from the deblended segmentation map
    if segm is not None:
        catwerr = SourceCatalog(data_bkgsub, segm_deblend,background=bkg,error=errormap,wcs=wcsmap)
        columns=['label','xcentroid','ycentroid','sky_centroid','background_centroid', 'background_mean','background_sum','background','area',
            'semimajor_sigma','semiminor_sigma','orientation', 'cxx','cxy','cyy','eccentricity','ellipticity','elongation',
            'fwhm','kron_flux','kron_fluxerr','kron_radius','segment_flux','segment_fluxerr']#,'kron_aperture']
        tbl2 = catwerr.to_table(columns)
    
        #calculate statistics for background cutouts in table
        segbkg_median=[]
        segbkg_mean=[]
        segbkg_std=[]
    
        #loop through each cutout and use sigma_cliped_stats to get mean, median, and std
        for i in range (0,len(tbl2['background'])):
            bkgdata=tbl2['background'][i]
            meansc, median_sigclip, stdsc = sigma_clipped_stats(bkgdata)
            segbkg_median.append(median_sigclip)
            segbkg_mean.append(meansc)
            segbkg_std.append(stdsc)
    
        #add the above calculated information to our table
        tbl2['segbkg_mean_sc']=segbkg_mean
        tbl2['segbkg_median_sc']=segbkg_median
        tbl2['segbkg_std_sc']=segbkg_std  
    
        #remove the 2d background array to tidy up the table
        tbl2.remove_column('background')

    
        #calculate noise stats
        #sky noise from background
        skynoise=np.sqrt(tbl2['segbkg_median_sc']*tbl2['area']/u.pix**2)
        #replace any nan values with zero
        sna=np.array(skynoise)
        masknan=np.isnan(sna)
        sna[masknan]=0.0
        skynoise=sna**2
    
        #shot noise from the source
        sourcenoise=tbl2['segment_flux']
    
        #thermal noise from camera (from error map)
        thermalnoise=tbl2['segment_fluxerr']
    
        #total noise
        #totalnoise=np.sqrt(sourcenoise+thermalnoise+skynoise) #includes all noise sources
        totalnoise=np.sqrt(thermalnoise+skynoise) #no shot noise -> For some reason this seems to work much better for the apertures. Need to think about why this is a bit more...
    
        #calculate SNR for the segments
        tbl2['segmentSNR']=tbl2['segment_flux']/ totalnoise
        
        #change format of columns to save fewer decimal places
        for col in tbl2.colnames:
            if col!='sky_centroid': #skip sky centroid since its problematic in this context
                tbl2[col].info.format = '%.4G'
                
        #print the number of sources found
        print('Number of segmentation map sources found: ', len(tbl2))
        
        #convert Qtable to table so we can write as fits
        ntbl2=Table(tbl2)

        #write out the resulting table to file
        #ascii.write(tbl2, name+'_'+str(wavelength)+'um_seg.dat', overwrite=True)
        ntbl2.write(name+'_'+str(wavelength)+'um_seg.fits',overwrite=True)
    else:
        tbl2=None
        #write out the resulting table to file
        #ascii.write([], name+'_'+str(wavelength)+'um_seg.dat', overwrite=True)
      
    
    #lets produce a new background map that's more optimized for point sources to use with DAOfind...
    bkg_estimator_ps = MMMBackground() #MMMBackground() #SExtractorBackground() #MedianBackground()
    bkg_data_ps = Background2D(data,(8, 8), filter_size=(3, 3),bkg_estimator=bkg_estimator_ps,edge_method='pad')
    bkgps_rms=bkg_data_ps.background_rms
    bkgps=bkg_data_ps.background 

    #create background subtracted image
    data_bkgsub_ps = data - bkgps
    
    #do statistics on image - daofind requires a single value for threshold rather than a 2d map
    #print('Seg. Map mean/median Threshold/3.0 for comparison: ', (np.mean(threshold)/3.,np.median(threshold)/3.))

    #mean, median, std = sigma_clipped_stats(data_bkgsub, sigma=3.0)  
    #print('Seg Map Bkgsub data mean, med, std for comparison: ',(mean, median, std))  

    mean, median, std = sigma_clipped_stats(data_bkgsub_ps, sigma=3.0)  
    #print('DAOFind Data mean, med, std: ',(mean, median, std))
   

    #Before star finding lets mask areas in the image that were also masked in the seg. map
    
    #make high cut on edges of fields for DAOphot - note that this is quite a bit larger than the seg map
    mask4=np.where(tmapnorm<0.75,tmapnorm,0).astype('bool')
    
    #create a single mask from the combination of the 3 earlier mask
    if m3lims is not None:
        maskcombine=(mask4 == 1) | (mask2 == 1) | (mask3 == 1)
    elif m2lims is not None:
        maskcombine=(mask4 == 1) | (mask2 == 1)
    else:
        maskcombine=(mask4 == 1)
    
    #create masked array for the background subtracted data
    data_bkgsub_ma = np.ma.masked_array(data_bkgsub_ps, mask=maskcombine)
    
    
    #now run starfinder to find sources 
    daofind = DAOStarFinder(fwhm=4.2, threshold=finddetsig*std)
    DAOsources = daofind(data_bkgsub_ma,mask=maskcombine)
    
    StarFinder = IRAFStarFinder(fwhm=4.2, threshold=finddetsig*std)
    IRAFsources = StarFinder(data_bkgsub_ma,mask=maskcombine)

    #plot data with apertures on detected sources
    if DAOsources is not None:
        Dpositions = np.transpose((DAOsources['xcentroid'], DAOsources['ycentroid']))
        Dapertures = CircularAperture(Dpositions, r=4.)
    if IRAFsources is not None:
        Ipositions = np.transpose((IRAFsources['xcentroid'], IRAFsources['ycentroid']))
        Iapertures = CircularAperture(Ipositions, r=4.)
    if tbl2 is not None:
        Spositions = np.transpose((tbl2['xcentroid'], tbl2['ycentroid']))
        Sapertures = CircularAperture(Spositions, r=4.)
    
    '''
    norm = ImageNormalize(stretch=SqrtStretch())
    
    fig, (ax1, ax2,ax3) = plt.subplots(1, 3, figsize=(16, 10))
    ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
    ax1.set_title('Seg. map Sources')
    Sapertures.plot(color='red', lw=1.5, alpha=0.5,axes=ax1)
    
    ax2.imshow(data, origin='lower', cmap='Greys_r',norm=norm)
    ax2.set_title('DAOfind Sources')
    Dapertures.plot(color='magenta', lw=1.5, alpha=0.5,axes=ax2)
    
    ax3.imshow(data, origin='lower', cmap='Greys_r',norm=norm)
    ax3.set_title('IRAF StarFind Sources')
    Iapertures.plot(color='cyan', lw=1.5, alpha=0.5,axes=ax3)
    
    #plt.show()
    '''
    
    
    #lets add ra, dec coordinates of sources to DAOsources table
    if DAOsources is None:
        Nsources=0
    else:
        Nsources=len(DAOsources['id'])
        
    scs=[]
    
    for i in range(0,Nsources):
        xcoord=DAOsources['xcentroid'][i]
        ycoord=DAOsources['ycentroid'][i]
        sc=pixel_to_skycoord(xcoord,ycoord,wcsmap)
        scs.append(sc)
        
    if DAOsources is not None:
        DAOsources['sky_centroid']=scs
    
        #change format of columns to save fewer decimal places
        for col in DAOsources.colnames:
            if col!='sky_centroid': #skip sky centroid since its problematic in this context
                DAOsources[col].info.format = '%.4G'
            
    #print the number of sources found
    if DAOsources is not None:
        print('Number of DAOfind sources found: ', len(DAOsources))
        #write out the resulting table to file
        #ascii.write(DAOsources, name+'_'+str(wavelength)+'um_dao.dat', overwrite=True)
        
        #convert Qtable to table so we can write as fits
        nDAOsources=Table(DAOsources)
        nDAOsources.write(name+'_'+str(wavelength)+'um_dao.fits',overwrite=True)
    else:
        print('Number of DAOfindSources found: 0')
    
    
    #set size of regions in ds9
    radius = Angle(0.00083333, u.deg) #must be in degrees - current value is r=3"
       
    #create ds9 regions file for segmentation map sources
    #start by getting lists of coords from table
    if tbl2 is not None:
        sourcecoords=tbl2['sky_centroid']
    
    
        #loop through and create region instances for each source
        regions=[]
        for i in range(0,len(sourcecoords)):
            region = CircleSkyRegion(sourcecoords[i], radius)
            regions.append(region)
        
        #write out region file
        write_ds9(regions, name+'_seg.reg')
    
    
        #change the color of the regions to red - no built in way to do this in regions package :-/
        with open(name+'_seg.reg', 'r+') as f:
            text = f.read()
            text = re.sub(r'\)', r') # color=red', text)
            f.seek(0)
            f.write(text)
            f.truncate()

    
    regions=[]
    #create ds9 region file for DAO sources
    if DAOsources is not None:
        sourcecoords=DAOsources['sky_centroid']
    
        
        for i in range(0,len(sourcecoords)):
            region = CircleSkyRegion(sourcecoords[i], radius)
            regions.append(region)
        
    #write out region file
    write_ds9(regions, name+'_DF.reg')
    

    #change the color of the regions to cyan - no built in way to do this in regions package :-/
    with open(name+'_DF.reg', 'r+') as f:
        text2 = f.read()
        text2 = re.sub(r'\)', r') # color=cyan', text2)
        f.seek(0)
        f.write(text2)
        f.truncate()
    

    #create new region file including both sets of sources
    newregtext=text+text2[45:]
    
    with open(name+'.reg', 'w+') as f:
        f.seek(0)
        f.write(newregtext)
        f.truncate()
    
    
    #open ds9 and load in the region file that was just saved
    if interactive:
        if ds9path is not None:
            p1=subprocess.Popen(ds9path+' '+filename+ ' -mode region -regions load '+name+'.reg', shell=True) 
            p1.wait()

    
    
    
    
    
