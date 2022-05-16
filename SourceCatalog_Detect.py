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
from astropy.table import join, Table, vstack

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

####-------------------------functions-------------------------------------

#Now add RA, DEC coordinates of sources to table
def addSkyCentroid(tab,wcsmap):
    if 'sky_centroid' in tab.columns:
        print('this table already has sky centroid column!')
    else:
        #Now add RA, DEC coordinates of sources to table
        if len(tab)>0:

            Nsources=len(tab)
            scs=[]

            for i in range(0,Nsources):
                xcoord=tab['xcentroid'][i]
                ycoord=tab['ycentroid'][i]
                sc=pixel_to_skycoord(xcoord,ycoord,wcsmap)
                scs.append(sc)

            tab['sky_centroid']=scs
        
    return tab

# functions for creating ds9 region files
def makeDS9reg(tab,radius,outname,color=None):
    #get source coordinates from table
    sourcecoords=tab['sky_centroid']

    #set size of regions 
    radius = Angle(radius, u.deg) 

    #loop through and create region instances for each source
    regions=[]
    for i in range(0,len(sourcecoords)):
        region = CircleSkyRegion(sourcecoords[i], radius)
        regions.append(region)
        
    #write out region file
    write_ds9(regions, outname)
    
    if color is not None:
        #change the color of the regions to red - no built in way to do this in regions package :-/
        with open(outname, 'r+') as f:
            text = f.read()
            text = re.sub(r'\)', r') # color='+str(color), text)
            f.seek(0)
            f.write(text)
            f.truncate()
    return text

#create ds9 file by combining two region files
def makeCombDS9file(text1,text2,outname):
    newregtext=text1+text2[45:]

    with open(outname, 'w+') as f:
        f.seek(0)
        f.write(newregtext)
        f.truncate()
        
#function to find the 'not' indcies for sources in table
def findNOTindex(tab,index):
    allpos=np.linspace(0,len(tab)-1,len(tab),dtype=np.int64) 
    notindex=list(set(allpos)-set(index))
    return notindex

####----------------------- start of code-----------------------------------------
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
    
    #create pixel error map by taking sqrt of variance map
    errormap=np.sqrt(varmap)
    
    if interactive:
        #plot the data and background model
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3,figsize=(14,10))

        ax1.set_title('Data')
        im1=ax1.imshow(data,origin='lower')
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im1, cax=cax, orientation='vertical')
        ax1.set_axis_off()

        ax2.set_title('Variance Map')
        im2=ax2.imshow(varmap,origin='lower')
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im2, cax=cax, orientation='vertical')
        ax2.set_axis_off()

        ax3.set_title('Exposure Time Map')
        im3=ax3.imshow(tmap,origin='lower')
        divider = make_axes_locatable(ax3)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im3, cax=cax, orientation='vertical')
        ax3.set_axis_off()
        plt.show()
    
    #create background model for image using median method
    bkg_estimator = MedianBackground() #MMMBackground() #SExtractorBackground() #MedianBackground()
    bkg_data = Background2D(data,(bkgbox, bkgbox), filter_size=(5, 5),bkg_estimator=bkg_estimator,edge_method='pad')
    bkg_rms=bkg_data.background_rms
    bkg=bkg_data.background 
    
    #create background subtracted image
    data_bkgsub = data - bkg
    
    #set detection threshold for source finding based on modeled background rms
    threshold = segdetsig*bkg_rms
    
    if interactive:
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
        
    
    ###----------------------- DAO/IRAF Star Finders --------------------------------
    
    #First, create a mask for the data
    #Make a cut based on exposure time. Note this is more agressive than SEG because we can't use 2D background with DAO or IRAF
    tmapnorm=tmap/np.max(tmap) #creating a normalized exposure time map for the mask
    maskTPS=np.where(tmapnorm<0.75,tmapnorm,0).astype('bool')
    
    mask2=np.zeros(np.shape(maskTPS))
    mask3=np.zeros(np.shape(maskTPS))

    #create a single mask from the combination of mask parameters specified in the config file
    if m3lims is not None:
        maskPS=(maskTPS == 1) | (mask2 == 1) | (mask3 == 1)
    elif m2lims is not None:
        maskPS=(maskTPS == 1) | (mask2 == 1)
    else:
        maskPS=(maskTPS == 1)

    #create masked array for the background subtracted data
    data_bkgsub_ma = np.ma.masked_array(data_bkgsub, mask=maskPS)

    #get standard devaition of image to specify detection threshold
    #std = np.median(bkg_rms) # use rms map value from background model - must be single value becuase methods don't allow for passing an array

    #Alternatively - use sigma_clipped_stats on the image to estimate background for source finding
    mean, median, std = sigma_clipped_stats(data_bkgsub, sigma=3.0)      
        
    #now run starfinder routines to find sources 
    daofind = DAOStarFinder(fwhm=4, threshold=finddetsig*std)
    DAOsources = daofind(data_bkgsub_ma,mask=maskPS)

    StarFinder = IRAFStarFinder(fwhm=4, threshold=finddetsig*std)
    IRAFsources = StarFinder(data_bkgsub_ma,mask=maskPS)    
        
    if interactive:
        #plot data with apertures on detected sources
        if DAOsources is not None:
            Dpositions = np.transpose((DAOsources['xcentroid'], DAOsources['ycentroid']))
            Dapertures = CircularAperture(Dpositions, r=4.)

        if IRAFsources is not None:
            Ipositions = np.transpose((IRAFsources['xcentroid'], IRAFsources['ycentroid']))
            Iapertures = CircularAperture(Ipositions, r=4.)


        norm = ImageNormalize(stretch=SqrtStretch())

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))

        ax1.imshow(data_bkgsub, origin='lower', cmap='Greys_r',norm=norm)
        ax1.set_title('DAOfind Sources')
        if DAOsources is not None:
            Dapertures.plot(color='magenta', lw=1.5, alpha=0.5,axes=ax1)

        ax2.imshow(data_bkgsub, origin='lower', cmap='Greys_r',norm=norm)
        ax2.set_title('IRAF StarFind Sources')
        if IRAFsources is not None:
            Iapertures.plot(color='cyan', lw=1.5, alpha=0.5,axes=ax2)

        plt.show()
        
    #add skycoords to source tables
    DAOsources=addSkyCentroid(DAOsources,wcsmap)
    IRAFsources=addSkyCentroid(IRAFsources,wcsmap)
    
    #print the number of sources found
    if DAOsources is not None:
        print('Number of DAOfind sources found: ', len(DAOsources))
      
        #convert Qtable to table so we can write as fits
        #nDAOsources=Table(DAOsources)
        #nDAOsources.write(name+'_'+str(wavelength)+'um_dao.fits',overwrite=True)

        #save a ds9 regions file for the sources
        txt2=makeDS9reg(DAOsources,0.00083333,name+'_daoFind.reg','magenta') #r=3"
    else:
        print('Number of DAOfind Sources found: 0')
        txt2=[]
        
        #print the number of sources found
    if IRAFsources is not None:
        print('Number of IRAFfind sources found: ', len(IRAFsources))
      
        #convert Qtable to table so we can write as fits
        #nDAOsources=Table(DAOsources)
        #nDAOsources.write(name+'_'+str(wavelength)+'um_iraf.fits',overwrite=True)

        #save a ds9 regions file for the sources
        txt3=makeDS9reg(IRAFsources,0.00083333,name+'_irafFind.reg','cyan') #r=3"
    else:
        print('Number of IRAFfind Sources found: 0')
        txt2=[]
        
        
    ###------------------------- SEG MAP ----------------------------------    
        
    
    #create mask for 'bad' parts of the image (e.g., edges of field where dithers don't overlap well)
    mask=np.where(tmapnorm<m1cut,tmapnorm,0).astype('bool') #create mask for any locations less than specified fraction of total integration time
    #note this can be adjusted in the config file parameters

    #Use a gaussian to smooth image prior to running the segmentation map - (FWHM 3x3, but can be adjusted in the line below)
    sigma = 3.0 * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()


    #Create Segmentation Map and detect sources
    segm = detect_sources(data_bkgsub, threshold, mask=mask, npixels=5, kernel=kernel)
    
    #remove any that exist in masked region as defined in the config file
    if m2lims is not None:
        mask2=np.zeros(np.shape(mask))
        for lim in m2lims:
            mask2[lim[0]:lim[1],lim[2]:lim[3]]=1
        segm.remove_masked_labels(mask2.astype('bool'))
    if m3lims is not None:
        mask3=np.zeros(np.shape(mask))
        for lim in m3lims:
            mask3[lim[0]:lim[1],lim[2]:lim[3]]=1
        segm.remove_masked_labels(mask3.astype('bool'))

    #lets take a look at deblending sources
    if segm is not None:
        segm_deblend = deblend_sources(data_bkgsub, segm, npixels=5,kernel=kernel, nlevels=64,contrast=0.001)
    else:
        segm_deblend=None #if segm is empty pass it on to avoid errors

    if interactive:
        #make plot of segmentation image to show detected sources side by side with data
        norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
        ax1.imshow(data_bkgsub, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        cmap = segm.make_cmap(seed=123)
        ax2.imshow(segm, origin='lower', cmap=cmap, interpolation='nearest')
        ax2.set_title('Deblend Seg. Image')
        plt.show()

    #Build a table from the deblended segmentation map
    if segm is not None:    
        catwerr = SourceCatalog(data_bkgsub, segm_deblend,background=bkg,error=errormap,wcs=wcsmap)
        #columns=['label','xcentroid','ycentroid','sky_centroid','background_centroid', 'background_mean','background_sum','background', 'area','semimajor_sigma','semiminor_sigma','orientation', 'cxx','cxy','cyy','eccentricity','ellipticity','elongation','fwhm','kron_flux','kron_fluxerr','kron_radius','segment_flux','segment_fluxerr']#Specify table columns
        columns=['label','xcentroid','ycentroid','sky_centroid']
        SEGsources= catwerr.to_table(columns)
                 
        #print the number of sources found
        print('Number of segmentation map sources found: ', len(SEGsources))
        
        #convert Qtable to Table so it can be written as fits
        #SEGsources=Table(SEGsources)

        #write out the resulting table to file
        #SEGsources.write(name+'_'+str(wavelength)+'um_seg.fits',overwrite=True)

        #save a ds9 regions file for the sources
        txt1=makeDS9reg(SEGsources,0.00083333,name+'_segFind.reg','red') #r=3"
    else:
        SEGsources=None
        txt1=[]
      
    

    #plot data with apertures on detected sources
    if DAOsources is not None:
        Dpositions = np.transpose((DAOsources['xcentroid'], DAOsources['ycentroid']))
        Dapertures = CircularAperture(Dpositions, r=4.)
    if IRAFsources is not None:
        Ipositions = np.transpose((IRAFsources['xcentroid'], IRAFsources['ycentroid']))
        Iapertures = CircularAperture(Ipositions, r=4.)
    if SEGsources is not None:
        Spositions = np.transpose((SEGsources['xcentroid'], SEGsources['ycentroid']))
        Sapertures = CircularAperture(Spositions, r=4.)
    
    if interactive:
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

        plt.show()
    
   
    #open ds9 and load in the region file that was just saved
    if interactive:
        if ds9path is not None:
            
            savename=name+'.reg'
            makeCombDS9file(t1,t2,t3,savename)  
    
            p1=subprocess.Popen(ds9path+' '+filename+ ' -mode region -regions load '+name+'.reg', shell=True) 
            p1.wait()

            
    #read in user defined sources if they exist
    UserFile=name+'_ds9.reg'
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

        ds9sc=SkyCoord(clist,unit=u.deg,frame='fk5')
        usersources=True
        print('Number of user defined DS9 sources found: ', len(clist))  
    else:
        usersources=False
        print('No user defined DS9 sources found')       
     
    
        #plot data with sources found in each method used
    if interactive:
        #plot data with sources found in each method used
        if DAOsources is not None:
            Dpositions = np.transpose((DAOsources['xcentroid'], DAOsources['ycentroid']))
            Dapertures = CircularAperture(Dpositions, r=4.)

        if IRAFsources is not None:
            Ipositions = np.transpose((IRAFsources['xcentroid'], IRAFsources['ycentroid']))
            Iapertures = CircularAperture(Dpositions, r=4.)

        if SEGsources is not None:
            Spositions = np.transpose((SEGsources['xcentroid'], SEGsources['ycentroid']))
            Sapertures = CircularAperture(Spositions, r=4.)

        if usersources:
            Upositions=np.transpose(skycoord_to_pixel(ds9sc,wcsmap))
            Uapertures = CircularAperture(Upositions, r=4.)

        norm = ImageNormalize(stretch=SqrtStretch())

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 12))

        ax1.imshow(data, origin='lower', cmap='Greys_r',norm=norm)
        ax1.set_title('DAOfind Sources')
        if DAOsources is not None:
            Dapertures.plot(color='magenta', lw=1.5, alpha=0.5,axes=ax1)

        ax2.imshow(data, origin='lower', cmap='Greys_r',norm=norm)
        ax2.set_title('IRAFfind Sources')
        if DAOsources is not None:
            Dapertures.plot(color='cyan', lw=1.5, alpha=0.5,axes=ax2)

        ax3.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
        ax3.set_title('Seg. map Sources')
        if SEGsources is not None:
            Sapertures.plot(color='red', lw=1.5, alpha=0.5,axes=ax3)

        ax4.imshow(data, origin='lower', cmap='Greys_r',norm=norm)
        ax4.set_title('User Defined Sources')
        if usersources:
            Uapertures.plot(color='yellow', lw=1.5, alpha=0.5,axes=ax4)

        plt.show()
    
    ### ------------ label source discovery method -----------------
    
    #now create a combined list of all sources found by each method
    if (DAOsources is not None)&(IRAFsources is not None):   
    
        #get source coordinates from DAO and SEG
        sourcesDao=DAOsources['sky_centroid']
        sourcesIRAF=IRAFsources['sky_centroid']

        #crossmatch source lists to look for duplication
        idx,rdx, d2d, d3d = sourcesDao.search_around_sky(sourcesIRAF, 3*u.arcsec)
        print('Number of crossmatched DAO/IRAF sources found: ', len(idx))

        #create tables for sources - start with crossmatched sources
        if len(rdx)>0:
            tab1=Table()
            tab1['sky_centroid']=sourcesDao[rdx]
            tab1['type']=100 #crossmatch IRAF & DAO
            t1=True
        else:
            t1=False

        #Find unique sources in list
        notrdx=findNOTindex(sourcesDao,rdx)
        if len(notrdx)>0:
            tab2=Table()
            tab2['sky_centroid']=sourcesDao[notrdx]
            tab2['type']=200 #DAO
            t2=True
        else:
            t2=False

        #Find unique sources in list
        notidx=findNOTindex(sourcesIRAF,idx)
        if len(notidx)>0:
            tab3=Table()
            tab3['sky_centroid']=sourcesIRAF[notidx]
            tab3['type']=300 #IRAF
            t3=True
        else:
            t3=False

        if t1&t2&t3:
            #create combined list
            mergeTab=vstack([tab1,tab2,tab3])
        elif t1 is False:
            mergeTab=vstack([tab2,tab3])
        elif t2 is False:
            mergeTab=vstack([tab1,tab3])
        elif t3 is False:
            mergeTab=vstack([tab1,tab2])
        else:
            print('Throw error! position 1')   
            
    elif IRAFsources is None:
        mergeTab=Table()
        mergeTab['sky_centroid']=DAOsources['sky_centroid']
        mergeTab['type']=200
        
    elif DAOsources is None:
        mergeTab=Table()
        mergeTab['sky_centroid']=IRAFsources['sky_centroid']
        mergeTab['type']=300
        
    else:
        print('Danger, unexpected behavior...')
    
    
    
    #now create a combined list of all sources found by each method
    if (mergeTab is not None)&(SEGsources is not None):
        #get source coordinates from DAO and SEG
        Findsources=mergeTab['sky_centroid']
        sourcesSeg=SEGsources['sky_centroid']

        #crossmatch source lists to look for duplication
        idx,rdx, d2d, d3d = Findsources.search_around_sky(sourcesSeg, 3*u.arcsec)
        print('Number of crossmatched Finder sources and Seg sources: ', len(idx))

        #create tables for sources - start with crossmatched sources
        if len(rdx)>0:
            tab4=Table()
            tab4['sky_centroid']=mergeTab['sky_centroid'][rdx]
            tab4['type']=mergeTab['type'][rdx]+10
            t4=True
        else:
            t4=False

        #Find unique sources in list
        notrdx=findNOTindex(mergeTab,rdx)
        if len(notrdx)>0:
            tab5=Table()
            tab5['sky_centroid']=mergeTab['sky_centroid'][notrdx]
            tab5['type']=mergeTab['type'][notrdx]+20
            t5=True
        else:
            t5=False

        #Find unique sources in list
        notidx=findNOTindex(sourcesSeg,idx)
        if len(notidx)>0:
            tab6=Table()
            tab6['sky_centroid']=sourcesSeg[notidx]
            tab6['type']=10
            t6=True
        else:
            t6=False

        if t4&t5&t6:
            #create combined table
            CombTab=vstack([tab4,tab5,tab6])
        elif t6 is False:
            CombTab=vstack([tab4,tab5])
        elif t5 is False:
            CombTab=vstack([tab4,tab6])
        elif t4 is False:
            CombTab=vstack([tab5,tab6])
        else:
            print('Throw error! position 2') 

            
    elif SEGsources is None:
        CombTab=mergeTab
        
    elif Findsources is None:
        CombTab=Table()
        CombTab['sky_centroid']=SEGsources['sky_centroid']
        CombTab['type']=10
        
    else:
        print('Danger, unexpected behavior...')

        
    if usersources:
        #now check against user sources to ensure no duplication there. 
        mergeList1=CombTab['sky_centroid']

        idx, rdx, d2d, d3d = mergeList1.search_around_sky(ds9sc, 3*u.arcsec)
        print('Number of crossmatched sources found: ', len(idx))    
        
        if len(rdx)>0:
            tab7=Table()
            tab7['sky_centroid']=CombTab['sky_centroid'][rdx]
            tab7['type']=CombTab['type'][rdx]+1
            t7=True
        else:
            t7=False

        notrdx=findNOTindex(mergeList1,rdx)
        if len(notrdx)>0:
            tab8=Table()
            tab8['sky_centroid']=CombTab['sky_centroid'][notrdx]
            tab8['type']=CombTab['type'][notrdx]+0 #not a user source
            t8=True
        else:
            t8=False


        notidx=findNOTindex(ds9sc,idx)
        if len(notidx)>0:
            tab9=Table()
            tab9['sky_centroid']=ds9sc[notidx]
            tab8['type']=1
            t9=True
        else:
            t9=False



        if t7&t8&t9:
            #create combined table
            CombTabNew=vstack([tab7,tab8,tab9])
        elif t9 is False:
            CombTabNew=vstack([tab7,tab8])
        elif t8 is False:
            CombTabNew=vstack([tab7,tab9])
        elif t7 is False:
            CombTabNew=vstack([tab8,tab9])
        else:
            print('Throw error! position 3')

    else:
        CombTabNew=CombTab
        print('No user Sources Defined')


    if CombTab is not None:
        #add id column to table before saving
        CombTabNew['id']= np.arange(1, len(CombTabNew) + 1)
    
        #write out the resulting table to file
        CombTabNew.write(name+'_'+str(wavelength)+'um_CombinedSources.fits',overwrite=True)
    

