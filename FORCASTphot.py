#import all required packages
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import subprocess

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.coordinates import SkyCoord, Angle 
from astropy.coordinates.matching import _get_cartesian_kdtree
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord,skycoord_to_pixel
from astropy import units as u
from astropy.stats import SigmaClip,sigma_clipped_stats
from astropy.table import join, Table, vstack
from astropy.nddata import Cutout2D

from photutils.aperture import SkyCircularAperture,SkyCircularAnnulus,aperture_photometry, CircularAperture
from photutils.background import Background2D, MedianBackground, SExtractorBackground, MMMBackground
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
from photutils import DAOStarFinder,IRAFStarFinder, make_source_mask
from photutils.morphology import data_properties
from photutils import make_source_mask

from regions import read_ds9, write_ds9, CircleSkyRegion

from config import *

# ----------------- newest stuff --------------------
from astropy.modeling import models, fitting
import warnings
from astropy.utils.exceptions import AstropyUserWarning




def extractModel(psfmodel, datacutout, unccutout,wv,correctErrors=True):
    #Normalize psfmodel
    psfmodel = psfmodel / np.sum(psfmodel)

    #Compute best-fit flux for model
    model_flux = np.sum(psfmodel * datacutout) / np.sum(np.square(psfmodel))
    
    unccutout=unccutout/2.0 #chi2 too small - adjust error array estimate to compensate
    
    #Compute chi2 of model fit
    if wv==25:
        deg_freedom = np.size(datacutout) - 14 
    else:
        deg_freedom = np.size(datacutout) - 21
    chi2 = np.sum(np.square(datacutout - psfmodel * model_flux) / np.square(unccutout)) / deg_freedom 
    
    if (correctErrors and chi2>2.0):
        re_unccutout=unccutout*np.sqrt(chi2)
        model_flux_err = np.sqrt(np.sum(np.square(psfmodel) * np.square(re_unccutout))) / np.sum(np.square(psfmodel))
        #chi2 = np.sum(np.square(datacutout - psfmodel * model_flux) / np.square(re_unccutout)) / deg_freedom    #don't update values
    else:
        model_flux_err = np.sqrt(np.sum(np.square(psfmodel) * np.square(unccutout))) / np.sum(np.square(psfmodel))
        
    return model_flux, model_flux_err, chi2




def modelSources(data,errorimg,tab,header,cutouts=False,cutsize=25):
    
    wavelen=tab['wv'][0]
    
    if wavelen==25:
        #print('Wavelength 25')
        csize=17
        psfwhm=4.5
        maxfwhm=10
        vexlim=12
    elif wavelen==37:
        csize=21
        psfwhm=5.5
        maxfwhm=12
        vexlim=14
        
    
    csize=17
    xfitm=[]
    yfitm=[]
    mamp=[]
    alpha=[]
    alphaerr=[]
    gamma=[]
    gammaerr=[]
    fwhm=[]
    fwhmerr=[]
    xfitg=[]
    yfitg=[]
    xfiterr=[]
    yfiterr=[]
    gamp=[]
    xfwhm=[]
    xfwhmerr=[]
    yfwhm=[]
    yfwhmerr=[]
    elong=[]
    orient=[]
    modelflux=[]
    modelfluxerr=[]
    modelSNRs=[]
    modelchi2=[]
    gfitwarn=[]
    mfitwarn=[]
    gmodelflux=[]
    gmodelfluxerr=[]
    gmodelSNRs=[]
    gmodelchi2=[]
    modeltype=[]
    #fitwarn=[]
    fdist=[]
    vex=[]

    cimages=np.zeros((len(tab),csize,csize))
    gresidual=np.zeros((len(tab),csize,csize))
    mresidual=np.zeros((len(tab),csize,csize))
    cutimgs=np.zeros((len(tab),cutsize,cutsize)) #these are the display cutouts
    pcovmatrix=[]
    
    
    i=0
    for source in tab:
        
        #not x/y centroid came from fitsources which is depreciated now...
        #try:
        #    spos=(source['xcentroid'].value,source['ycentroid'].value)
        #except:
        #    print('centroids are missing from table. Using x/y center instead')
        #    spos=(source['xcenter'].value,source['ycenter'].value)
        
        spos=(source['xcenter'].value,source['ycenter'].value)
        
        cut_img=Cutout2D(data,spos,csize,mode='partial',fill_value=0.0,copy=True)
        unc_img=Cutout2D(errorimg,spos,csize,mode='partial',fill_value=0.0,copy=True)
        
        cimg=cut_img.data
        uimg=unc_img.data
        
        #Background subtraction based on annulus background
        cimg=cimg-tab['ann_bkg_med'][i] #could do this better with a 2D bkg model... but larger cutout needed.
        
        y, x = np.mgrid[:csize, :csize]
        
        # Fit the data using astropy.modeling
        p1_init = models.Moffat2D(x_0=np.int32(csize/2),y_0=np.int32(csize/2))
        fit_p1 = fitting.LevMarLSQFitter(calc_uncertainties=True)
        
        
        with warnings.catch_warnings(record=True) as w:      
            
            p1 = fit_p1(p1_init, x, y, cimg)
            #print(w[0].message) #debugging to check error messages
            #fitwarn.append(len(w))
            

        if p1.cov_matrix is not None:
            p1cov=np.diag(p1.cov_matrix.cov_matrix)
            #xfiterr.append(p1cov[1])  #in pixels
            #yfiterr.append(p1cov[2])  #in pixels
            gammaerr.append(p1cov[3])  #in pixels
            alphaerr.append(p1cov[4])  #alpha is unit-less
            grammaer=p1cov[3] #in pixels
            mfit=0
            xfitm.append(p1.x_0.value) #in pixels
            yfitm.append(p1.y_0.value) #in pixels
            mamp.append(p1.amplitude.value)
            alpha.append(p1.alpha.value) #unitless
            gamma.append(p1.gamma.value) #in pixels
        else:
            #xfiterr.append(np.nan)
            #yfiterr.append(np.nan)
            gammaerr.append(np.nan)
            grammaer=np.nan
            alphaerr.append(np.nan)
            mfit=1
            xfitm.append(np.nan) #in pixels
            yfitm.append(np.nan) #in pixels
            mamp.append(np.nan)
            alpha.append(np.nan) #unitless
            gamma.append(np.nan) #in pixels
        
        mmodel=p1(x, y)
        mresid=cimg.data - mmodel
        
        mflux,mfluxerr,chi2m=extractModel(mmodel,cimg,uimg,wavelen)
        
        
        r_effective=(p1.fwhm)/2.355*1.645 #use 1.645 sigma? could be adjusted...
        #etamoffat=np.sqrt(2*np.pi*(r_effective*tab['ann_bkg_std'][i]/mflux)**2+(header['ERRCALF']/header['CALFCTR'])**2+0.005625+(mfluxerr/mflux)**2)
        #mfluxunc=etamoffat*mflux
        #modelSNR=mflux/(np.sqrt(2*np.pi)*r_effective*tab['ann_bkg_std'][i])
        
        mfluxunc=np.sqrt(2*np.pi)*r_effective*tab['ann_bkg_std'][i]
        modelSNR=mflux/mfluxunc
        
        
        p2_init = models.Gaussian2D(x_mean=np.int32(csize/2),y_mean=np.int32(csize/2))
        fit_p2 = fitting.LevMarLSQFitter(calc_uncertainties=True)
        
        with warnings.catch_warnings(record=True) as w2:
            p2 = fit_p2(p2_init, x, y, cimg)
            #gfitwarn.append(len(w2))
            
        if p2.cov_matrix is not None:
            xfwhm.append(p2.x_fwhm*0.768) # in arcseconds
            yfwhm.append(p2.y_fwhm*0.768) # in arcseconds
            gamp.append(p2.amplitude.value)
            xfitg.append(p2.x_mean.value)
            yfitg.append(p2.y_mean.value)
            p2cov=np.diag(p2.cov_matrix.cov_matrix)
            xfiterr.append(p2cov[1])
            yfiterr.append(p2cov[2])
            xfwhmerr.append(p2cov[3]*0.768) #in arcsec
            yfwhmerr.append(p2cov[4]*0.768) #in arcsec
            gfwhmerr=np.sqrt((p2cov[3])**2+(p2cov[4])**2) #in pixels
            gfit=0
            gfwhm=np.sqrt(p2.x_fwhm*p2.y_fwhm) # calculate an effective 1D FWHM for the gaussian model   
        else:
            xfwhm.append(np.nan)
            yfwhm.append(np.nan)
            gamp.append(np.nan)
            xfitg.append(np.nan)
            yfitg.append(np.nan)
            xfiterr.append(np.nan)
            yfiterr.append(np.nan)
            xfwhmerr.append(np.nan) #was -1
            yfwhmerr.append(np.nan) #was -1
            gfit=1
            gfwhm=np.nan
            gfwhmerr=np.nan #was -1
        
        
        gmodel=p2(x, y)
        
        if not np.sum(gmodel)>0: #make dummy psf if model psf is array of zeros 
            gmodel[7:11,7:11]=1
            gfit=1

        
        gresid=cimg - gmodel   
        
        gflux,gfluxerr,chi2g=extractModel(gmodel,cimg,uimg,wavelen)
        
        
        
        r_effective=gfwhm/2.355*1.645 #use 1.645 sigma? could be adjusted...
        #etagauss=np.sqrt(2*np.pi*(r_effective*tab['ann_bkg_std'][i]/gflux)**2+(header['ERRCALF']/header['CALFCTR'])**2+0.005625+(gfluxerr/gflux)**2)
        #gfluxunc=etagauss*gflux
        #gmodelSNR=gflux/(np.sqrt(2*np.pi)*r_effective*tab['ann_bkg_std'][i])
        
        gfluxunc=np.sqrt(2*np.pi)*r_effective*tab['ann_bkg_std'][i]
        gmodelSNR=gflux/gfluxunc
        
        if cutouts:
            c2=Cutout2D(data,spos,cutsize,mode='partial',fill_value=0.0,copy=True)
            cutimgs[i,:,:]=c2.data
            gresidual[i,:,:]=gresid
            mresidual[i,:,:]=mresid
            
        
        
                
        if p2.x_fwhm>=p2.y_fwhm:
            elong.append(p2.x_fwhm/p2.y_fwhm)
        elif p2.x_fwhm<p2.y_fwhm:
            elong.append(p2.y_fwhm/p2.x_fwhm)
        else:
            elong.append(np.nan)
        
        
        #initial classification based on model fit performance
        if grammaer/p1.gamma.value<0.8:
            if chi2m<(1.1*chi2g) and elong[-1]<1.25 and (p1.fwhm*0.768)<psfwhm:
                modelflag='m'          
            else:
                if (gfwhmerr/gfwhm)<1.0:
                    modelflag='g'
                else:
                    modelflag='a'      
        elif (gfwhmerr/gfwhm)<1.0:
            modelflag='g'
        else:
            modelflag='a'
        
        #adjust model assessment based on chi squared values
        if modelflag=='g':
            if chi2m<(1.1*chi2g) and elong[-1]<=1.25 and (p1.fwhm*0.768)<psfwhm and grammaer/p1.gamma.value<0.8:
                modelflag='m' 
                
        if modelflag=='m' and chi2m>4:
            if chi2g<=4 and (gfwhmerr/gfwhm)<1.0:
                modelflag='g'
            else:
                modelflag='a'
                
        if modelflag=='g' and chi2g>4:
            modelflag='a'
            
        
        #calculate source dist from center and use for flagging
        msdist=np.sqrt((p2.x_mean.value-np.int32(csize/2))**2+(p2.y_mean.value-np.int32(csize/2))**2)
        fdist.append(msdist)
        
        
        if modelflag!='a' and msdist>4:
            modelflag='a'
            
        
        if modelflag!='a' and gfwhm>maxfwhm:
                modelflag='a'
                

        
        #retroactively update model flag if fit wasn't good
        if modelflag=='g' and gfit==1:
            modelflag='a'
        if modelflag=='m' and mfit==1:
            modelflag='a'
            
            
        
        #use moffat FWHM on sources marked moffat. 
        #if gfwhm<=vexlim:
        
        if modelflag=='m':
            fx=p1.fwhm*0.768
            fxerr=(grammaer/p1.gamma.value)*p1.fwhm*0.768
        elif modelflag=='g':
            fx=gfwhm*0.768
            fxerr=gfwhmerr*0.768
        elif modelflag=='a':
            if gfit==0:
                fx=gfwhm*0.768
            else:
                fx=np.nan
            fxerr=np.nan
        
        
        #append the FWHM values
        if fx == np.nan:
            fwhm.append(fx)
        elif fx<(cutsize-5):
            fwhm.append(fx)
        else:
            fwhm.append(cutsize-5)
            
        fwhmerr.append(fxerr)
        
        
        #flag if very extended
        if fx>vexlim:
            vex.append(1)
        else:
            vex.append(0)
                    
        
        
        '''# use gaussian FWHM for all fwhm measures
        if gfit==0:
            fwhm.append(gfwhm*0.768)
            
            if (gfwhmerr/gfwhm)<1.0:
                fwhmerr.append(gfwhmerr*0.768)
            else:
                fwhmerr.append(np.nan)
                #gfit=1
        else:
            fwhm.append(np.nan)
            fwhmerr.append(np.nan)
        '''
            
        
        mfitwarn.append(mfit)
        gfitwarn.append(gfit)
        
        modeltype.append(modelflag)
        
        if gfit==0:
            orient.append(p2.theta.value)
            gmodelflux.append(gflux)
            gmodelfluxerr.append(gfluxunc)
            gmodelSNRs.append(gmodelSNR)
            gmodelchi2.append(chi2g)
        else:
            orient.append(np.nan)
            gmodelflux.append(np.nan)
            gmodelfluxerr.append(np.nan)
            gmodelSNRs.append(np.nan)
            gmodelchi2.append(np.nan)            
            
        if mfit==0:    
            modelflux.append(mflux)
            modelfluxerr.append(mfluxunc)
            modelSNRs.append(modelSNR)
            modelchi2.append(chi2m)
        else:
            modelflux.append(np.nan)
            modelfluxerr.append(np.nan)
            modelSNRs.append(np.nan)
            modelchi2.append(np.nan)

        #sourcepixels.append(sourcepix)
        i=i+1
        
            
    #elongation
    tab['fwhm']=fwhm   #was 'fit_fwhm'
    tab['fwhm_err']=fwhmerr #was 'fit_fwhm_err'
    tab['elong']=elong #was 'fit_elong'
    tab['mfit_x0']=xfitm
    tab['mfit_y0']=yfitm
    tab['mfit_amp']=mamp
    tab['mfit_gamma']=gamma
    tab['mfit_gamma_err']=gammaerr
    tab['mfit_alpha']=alpha
    tab['mfit_alpha_err']=alphaerr
    tab['gfit_x0']=xfitg
    tab['gfit_y0']=yfitg
    tab['gfit_x0_err']=xfiterr
    tab['gfit_y0_err']=yfiterr
    tab['gfit_amp']=gamp
    tab['gfit_xfwhm']=xfwhm
    tab['gfit_xfwhm_err']=xfwhmerr
    tab['gfit_yfwhm']=yfwhm
    tab['gfit_yfwhm_err']=yfwhmerr
    tab['gfit_orient']=orient
    tab['FluxMoffat2D']=modelflux
    tab['FluxMoffatErr']=modelfluxerr
    tab['mModelSNR']=modelSNRs
    tab['Moffat2DChi2']=modelchi2
    tab['FluxGauss2D']=gmodelflux
    tab['FluxGauss2DErr']=gmodelfluxerr
    tab['gModelSNR']=gmodelSNRs
    tab['Gauss2DChi2']=gmodelchi2
    tab['MFitFlag']=mfitwarn
    tab['GFitFlag']=gfitwarn
    #tab['SourcePix']=sourcepixels
    tab['cutouts']=cutimgs
    tab['gresid']=gresidual
    tab['mresid']=mresidual
    tab['fit_dist']=fdist
    tab['vexFlag']=vex
        
    tab['BestModel']=modeltype
    
    return tab





'''
def extractModel(psfmodel, datacutout, unccutout):
	#Normalize psfmodel
	psfmodel = psfmodel / np.sum(psfmodel)
	
	#Compute best-fit flux for model
	model_flux = np.sum(psfmodel * datacutout) / np.sum(np.square(psfmodel))
	#Compute error for flux model
	model_flux_err = np.sqrt(np.sum(np.square(psfmodel) * np.square(unccutout)))/ np.sum(np.square(psfmodel))/10. #fudge factor...
	
	#Compute chi2 of model fit
	deg_freedom = np.size(datacutout) - 1 
	chi2 = np.sum(np.square(datacutout - psfmodel * model_flux) / np.square(unccutout)) / deg_freedom
	
	return model_flux, model_flux_err, chi2



def modelSources(data,errorimg,tab,header,cutouts=False,cutsize=25):
    csize=17
    xfit=[]
    yfit=[]
    xfiterr=[]
    yfiterr=[]
    amp=[]
    alpha=[]
    alphaerr=[]
    gamma=[]
    gammaerr=[]
    fwhm=[]
    xfwhm=[]
    yfwhm=[]
    elong=[]
    orient=[]
    modelflux=[]
    modelfluxerr=[]
    modelSNRs=[]
    modelchi2=[]
    fitwarn=[]
    gmodelflux=[]
    gmodelfluxerr=[]
    gmodelSNRs=[]
    gmodelchi2=[]
    gfitwarn=[]
    modeltype=[]
    sourcepixels=[]
    cimages=np.zeros((len(tab),csize,csize))
    residuals=np.zeros((len(tab),csize,csize))
    cutimgs=np.zeros((len(tab),cutsize,cutsize)) #these are the display cutouts
    #pcovmatrix=[]
    
    
    i=0
    for source in tab:
        spos=(source['xcentroid'].value,source['ycentroid'].value)
        
        cut_img=Cutout2D(data,spos,csize,mode='partial',fill_value=0.0,copy=True)
        unc_img=Cutout2D(errorimg,spos,csize,mode='partial',fill_value=0.0,copy=True)
        
        cimg=cut_img.data
        uimg=unc_img.data
        
        #Background subtraction based on annulus background
        cimg=cimg-tab['ann_bkg_med'][i] #could do this better with a 2D bkg model... but larger cutout needed.
        
        y, x = np.mgrid[:csize, :csize]
        
        # Fit the data using astropy.modeling
        p1_init = models.Moffat2D(x_0=np.int32(csize/2),y_0=np.int32(csize/2))
        fit_p1 = fitting.LevMarLSQFitter(calc_uncertainties=True)
        
        
        with warnings.catch_warnings(record=True) as w:      
            
            p1 = fit_p1(p1_init, x, y, cimg)
            #print(w[0].message) #debugging to check error messages
            fitwarn.append(len(w))
        
        if p1.cov_matrix is not None:
            p1cov=np.diag(p1.cov_matrix.cov_matrix)
            xfiterr.append(p1cov[1])
            yfiterr.append(p1cov[2])
            gammaerr.append(p1cov[3])
            alphaerr.append(p1cov[4])
            #pcovmatrix.append(p1cov)
        else:
            xfiterr.append(-1)
            yfiterr.append(-1)
            gammaerr.append(-1)
            alphaerr.append(-1)
            #pcovmatrix.append([-1, -1, -1, -1, -1])
        
        mmodel=p1(x, y)
        resid=cimg.data - mmodel
        
        mflux,mfluxerr,chi2m=extractModel(mmodel,cimg,uimg)
        
        r_effective=(p1.fwhm)/2.0
        noisecalc=np.sqrt(2*np.pi*(r_effective*tab['ann_bkg_std'][i]/mflux)**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025+(mfluxerr/mflux)**2)
        modelSNR=1.0*mflux/noisecalc
        
        
        
        p2_init = models.Gaussian2D(x_mean=np.int32(csize/2),y_mean=np.int32(csize/2))
        fit_p2 = fitting.LevMarLSQFitter(calc_uncertainties=True)
        
        with warnings.catch_warnings(record=True) as w2:
            p2 = fit_p2(p2_init, x, y, cimg)
            
            gfitwarn.append(len(w2))
        
        gmodel=p2(x, y)
        gresid=cimg - gmodel
        
        gflux,gfluxerr,chi2g=extractModel(gmodel,cimg,uimg)
        
        r_effective=np.sqrt(p2.x_fwhm*p2.y_fwhm)/2.0
        noisecalc=np.sqrt(2*np.pi*(r_effective*tab['ann_bkg_std'][i]/gflux)**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025+(gfluxerr/gflux)**2)
        gmodelSNR=1.0*gflux/noisecalc
        
        if cutouts:
            c2=Cutout2D(data,spos,cutsize,mode='partial',fill_value=0.0,copy=True)
            cutimgs[i,:,:]=c2
            
        try:    
            pix_aperture=CircularAperture((p1.x_0.value,p1.y_0.value), r=r_effective*1.2)
        except:
            pix_aperture=CircularAperture((17/.2,17.2), r=r_effective*1.2)
            
        source_mask = pix_aperture.to_mask(method='exact')
        source_data = source_mask.multiply(data)
        sourcedata=source_mask.data
        sourcepix=sum(sum(sourcedata>5*tab['ann_bkg_std'][i]))
        

        xfit.append(p1.x_0.value)
        yfit.append(p1.y_0.value)
        amp.append(p1.amplitude.value)
        alpha.append(p1.alpha.value)
        gamma.append(p1.gamma.value)
        fwhm.append(p1.fwhm*0.768) # in arcseconds
        xfwhm.append(p2.x_fwhm*0.768) # in arcseconds
        yfwhm.append(p2.y_fwhm*0.768) # in arcseconds
        
        if p2.x_fwhm<p2.y_fwhm:
            elong.append(p2.x_fwhm/p2.y_fwhm)
        else:
            elong.append(p2.y_fwhm/p2.x_fwhm)
        
        orient.append(p2.theta.value)
        modelflux.append(mflux)
        modelfluxerr.append(mfluxerr)
        modelSNRs.append(modelSNR)
        modelchi2.append(chi2m)
        
        gmodelflux.append(gflux)
        gmodelfluxerr.append(gfluxerr)
        gmodelSNRs.append(gmodelSNR)
        gmodelchi2.append(chi2g)
        sourcepixels.append(sourcepix)
        
        
        if chi2m<=2.3:
            mtype='M'
            residuals[i,:,:]=resid
        elif chi2g<=2.3:
            mtype='G'
            residuals[i,:,:]=gresid
        else:
            mtype='A'     
        
        #cimages[i,:,:]=cimg
        #if chi2m/chi2g < 1.05:
        #    if chi2m<=2.3:
        #        mtype = 'M'
        #        residuals[i,:,:]=resid
        #    else:
        #        mtype= 'A'
        #else:
        #    if chi2g<=2.3:
        #        mtype = 'G'
        #        residuals[i,:,:]=gresid
        #    else:
        #        mtype = 'A'
                
        modeltype.append(mtype)
        i=i+1
        
        
    
    #elongation
    tab['fit_x0']=xfit
    tab['fit_y0']=yfit
    tab['fit_x0_err']=xfiterr
    tab['fit_y0_err']=yfiterr
    tab['fit_amp']=amp
    tab['fit_gamma']=gamma
    tab['fit_gamma_err']=gammaerr
    tab['fit_alpha']=alpha
    tab['fit_alpha_err']=alphaerr
    tab['fit_fwhm']=fwhm
    tab['fit_xfwhm']=xfwhm
    tab['fit_yfwhm']=yfwhm
    tab['fit_elong']=elong
    tab['fit_orient']=orient
    tab['FluxMoffat2D']=modelflux
    tab['FluxMoffatErr']=modelfluxerr
    tab['ModelSNR']=modelSNRs
    tab['Moffat2DChi2']=modelchi2
    tab['FluxGauss2D']=gmodelflux
    tab['FluxGauss2DErr']=gmodelfluxerr
    tab['gModelSNR']=gmodelSNRs
    tab['Gauss2DChi2']=gmodelchi2
    tab['FitWarningFlag']=fitwarn
    tab['BestModel']=modeltype
    tab['SourcePix']=sourcepixels
    #tab['DataCutout']=cimages
    tab['ModelResidual']=residuals
    tab['cutouts']=cutimgs
    #tab['Model_pcov']=pcovmatrix
    
    return tab


'''


def performApPhoto(data,errormap,tmap,header,sourceCoords,radii,rin,rout,plot=True):
    
    #define wcs object for header
    wcs=WCS(header)
    
    wv=header['WAVELNTH']

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
    bkgflag=[]
    
    if wv==25.2:
        stdcheck=0.015 #25um
    elif wv==37.1:
        stdcheck=0.09 #37um
    else:
        print('Filter for '+str(wv)+' is not set up in photometry function...')
        print('Defaulting value to 0, which will result in all data backround measurements being sigma clipped - this may not be desirable so please be cautious!')
        stdcheck=0.0

    #create mask array for the annuli
    annulus_masks = pix_annulus_aperture.to_mask(method='exact')

    #for each of the annuli loop through and calculate stats using sigma cliped stats
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        maskdata=mask.data

        try:
            #do statistics with sig clip based on noise in annulus
            annulus_data_1d = annulus_data[maskdata > 0]
            stdval=np.nanstd(annulus_data_1d)
            
            if stdval>stdcheck: 
                meansc, median_sigclip, stdsc = sigma_clipped_stats(annulus_data_1d,maxiters=3)
                bkg_median.append(median_sigclip)
                bkg_mean.append(meansc)
                bkg_std.append(stdsc)
                bkgflag.append(1)
            else:
                meansc, median_sigclip, stdsc = sigma_clipped_stats(annulus_data_1d,sigma=5,maxiters=1)
                bkg_median.append(median_sigclip)
                bkg_mean.append(meansc)
                bkg_std.append(stdsc)
                bkgflag.append(0)
                
            appmasks.append(mask.data)
        except TypeError:
            bkg_median.append(np.nanmedian(annulus_data_1d))
            bkg_mean.append(np.nanmean(annulus_data_1d))
            bkg_std.append(np.nanstd(annulus_data_1d))
            appmasks.append(np.nan)
            
            if np.nanstd(annulus_data_1d)>stdcheck:
                bkgflag.append(1)
            else:
                bkgflag.append(0)
            
            
        #bkg_std.append(np.nanstd(annulus_data_1d)) #only for testing...
                
    #store values in numpy arrays
    bkg_median = np.array(bkg_median)
    bkg_mean = np.array(bkg_mean)
    bkg_std = np.array(bkg_std)
    bkg_flag=np.array(bkgflag)

    #add columns for background information and also background subtracted apertures
    phot_table['ann_bkg_med'] = bkg_median
    phot_table['ann_bkg_mean'] = bkg_mean 
    phot_table['ann_bkg_std'] = bkg_std 
    phot_table['bkg_flag'] = bkg_flag

    
    #information from exposure time maps
    #create lists to store information for later
    texp_mean=[]
    texp_med=[]
    texpmasks=[]
    edflag=[]
    
    #get maximum of exposure time map
    texpmax=np.nanmax(tmap)

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
        
        if meansc < 0.8*texpmax:
            edflag.append(1)
        else:
            edflag.append(0)

    #phot_table['texp_med'] = texp_med
    phot_table['texp_mean'] = texp_mean
    phot_table['edge_flag'] = edflag
    

    #caclulate background subtracted photometry and snr for each source
    for i in range(0,len(radii)):
        #coloumn names and from the phot table and new columns we'll add
        cname1='aperture_sum_'+str(i)
        cname2='aperture_sum_err_'+str(i)
        newcol1='aper_sum_bkgsub_'+str(radii[i])+'pix'
        newcol2='aper_'+str(radii[i])+'pix_noise'
        newcol3='aper_snr_'+str(radii[i])+'pix'
        newcol4='aper_area_'+str(radii[i])+'pix'
        
        #get pixel aperture areas for calculations
        pixel_ap = apertures[i].to_pixel(wcs)
        pixarea=pixel_ap.area
        #get background subracted photo by subtracting median annulus value from each pixel in aperture
        phot_table[newcol1]=(phot_table[cname1]/pixarea-phot_table['ann_bkg_med'])*pixarea
        #calculate noise following equation from forcast photometry cookbook -https://sofia-data-analysis-cookbooks.readthedocs.io/en/latest/FORCAST-photometry_detailed.html - Note there is an added error from the aperture extraction included (cname2)
        
        
        
        #treat SNR as % error
        #etavals=np.sqrt(2*np.pi*(radii[i]*phot_table['ann_bkg_std']/phot_table[newcol1])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.005625+(phot_table[cname2]/phot_table[cname1])**2) #upped from 5% to 7.5%
        #phot_table[newcol2]=etavals*phot_table[newcol1] #note eta is relative error
        
        phot_table[newcol2]=np.sqrt(2*np.pi)*radii[i]*phot_table['ann_bkg_std']
        
        
        #calculate SNR - old      
        #phot_table[newcol3]=phot_table[newcol1]/phot_table[newcol2]
        #phot_table[newcol3]=phot_table[newcol1]/(np.sqrt(2*np.pi*(radii[i]*phot_table['ann_bkg_std']/phot_table[newcol1])**2+(header['ERRCALF']/header['CALFCTR'])**2+(phot_table[cname2]*2.0/phot_table[cname1])**2))
        
        phot_table[newcol3]=phot_table[newcol1]/(np.sqrt(2*np.pi)*radii[i]*phot_table['ann_bkg_std'])
        
        phot_table[newcol4]=pixarea
        
        phot_table['ECF']=header['ERRCALF']/header['CALFCTR']

        #rename aperture columns in table to be more descriptive
        rename1='aperture_sum_'+str(radii[i])+'pix'
        rename2='aperture_sum_err_'+str(radii[i])+'pix'
        phot_table.rename_column(cname1, rename1)
        phot_table.rename_column(cname2, rename2)    
    
    
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
    #plt.savefig("detectimg.png",dpi=600)    
    plt.show()

    return phot_table



'''
def fitshapes(image,tab,plot=False,cutouts=False,cutsize=25):
    columns = ['xcentroid', 'ycentroid','fwhm' ,'semimajor_sigma','semiminor_sigma', 'orientation','elongation','covar_sigx2','covar_sigy2']
    
    #initialize table with correct column formatting by creating dummy table
    from photutils.datasets import make_4gaussians_image
    data = make_4gaussians_image()[43:79, 76:104]
    cat = data_properties(data)
    tbl = cat.to_table(columns=columns)
    tbl.remove_row(0) #remove entry from dummy table
      
    if cutouts:
        cimgs=np.zeros((len(tab),cutsize,cutsize))
        i=0
    
    #fix keywords in table if they don't match what is expected
    if 'xcentroid' not in tab.columns:
        tab.rename_column('xcenter', 'xcentroid')
        tab.rename_column('ycenter', 'ycentroid')
    
    #loop through sources and fit shapes
    for source in tab:
        #create data cutout around source centroid position
        spos=(np.int64(source['xcentroid']),np.int64(source['ycentroid']))
        cutout=Cutout2D(image,spos,13,mode='partial',fill_value=0.0,copy=True)
        
        #do addtional bkg subtraciton here? - think not for now...
        mean, median, std = sigma_clipped_stats(cutout.data, sigma=3.0)
        c1 = cutout.data - median
        
        #get shape fitting results and store to table
        cat = data_properties(c1,background=std)
        temp = cat.to_table(columns=columns)
        tbl=vstack([tbl,temp])
        
        if cutouts:
            c2=Cutout2D(image,spos,cutsize,mode='partial',fill_value=0.0,copy=True)
            
            #shapeinfo=np.shape(c2.data)
            #if shapeinfo==(cutsize,cutsize): #only store cutout if its correct size
            #try:
            #    cimgs[i,:,:]=c2.data
            #except:
            #    print('cutout generator failed...')

            cimgs[i,:,:]=c2.data
            i=i+1 #iterate
            
        #optional plots for visual diagnostics
        if plot==True:
            position = (cat.xcentroid, cat.ycentroid)
            r = 1.0  # approximate isophotal extent
            a = cat.semimajor_sigma.value * r
            b = cat.semiminor_sigma.value * r
            theta = cat.orientation.to(u.rad).value
            apertures = EllipticalAperture(position, a, b, theta=theta)
            plt.imshow(c1, origin='lower', cmap='viridis',interpolation=None)
            apertures.plot(color='r')
            plt.show()
    
    #if needed, rename duplicate columns from seg table
    if 'fwhm' in tab.columns:
        tab.rename_column('fwhm', 'fwhm_seg')
        #print('this already has FWHM column!')
    if 'semimajor_sigma' in tab.columns:
        tab.rename_column('semimajor_sigma', 'semimajor_sigma_seg')
    if 'semiminor_sigma' in tab.columns:
        tab.rename_column('semiminor_sigma', 'semiminor_sigma_seg')        
    if 'orientation' in tab.columns:
        tab.rename_column('orientation', 'orientation_seg')        
        
    #add fitted shape parameters to table
    tab['fit_x_cent']=tbl['xcentroid']
    tab['fit_y_cent']=tbl['ycentroid']        
    tab['fwhm']=tbl['fwhm']    
    tab['semimajor_sigma']=tbl['semimajor_sigma']
    tab['semiminor_sigma']=tbl['semiminor_sigma']
    tab['orientation']=tbl['orientation']
    tab['elongation']=tbl['elongation']
    tab['covar_sigx2']=tbl['covar_sigx2']
    tab['covar_sigy2']=tbl['covar_sigy2']
    
    tab['fit_dist']=np.sqrt((tbl['xcentroid']-6.5)**2+(tbl['ycentroid']-6.5)**2)
    
    if cutouts:
        tab['cutouts']=cimgs

    #return input table with new columns from shape fitting
    return tab
'''


def CombineFieldResults(CatName,wavelength): # #Name options- 'CombCat', 'SegCat', 'DaoCat'
    #pick a field to start with so we can get the table structure. 
    #The files for the field cannot be empty or an error message will display
    startfield=Field01.name

    if os.path.exists(startfield+'_'+str(wavelength)+'um_'+CatName+'.fits'):
        cat=Table.read(startfield+'_'+str(wavelength)+'um_'+CatName+'.fits')

        if 'sky_centroid' in cat.columns:

            #fix for issue with combining tables with objects as columns
            cat['RA(J2000)']=cat['sky_centroid'].ra.to_string(u.hour,precision=2)
            cat['DEC(J2000)']=cat['sky_centroid'].dec.to_string(u.degree,precision=2)
            cat['glon']=cat['sky_centroid'].galactic.l.deg
            cat['glat']=cat['sky_centroid'].galactic.b.deg
            cat.remove_column('sky_centroid')
            
        elif 'sky_center' in cat.columns:
            
            #fix for issue with combining tables with objects as columns
            cat['RA(J2000)']=cat['sky_center'].ra.to_string(u.hour,precision=2)
            cat['DEC(J2000)']=cat['sky_center'].dec.to_string(u.degree,precision=2)
            cat['glon']=cat['sky_center'].galactic.l.deg
            cat['glat']=cat['sky_center'].galactic.b.deg
            cat.remove_column('sky_center')

    else:
        print('Error... Must pick another file to start with')
        
        
    #loop through all the saved photometry tables for individual fields and append them together
    for info in field._registry:
        name=info.name
        
        if wavelength==25:
            filename=info.file25
        elif wavelength==37:
            filename=info.file37
        else:
            print('Error! Wavlength must be 25 or 37...')
            sys.exit(1)

        #print('\nLoading in photometry data from field: ', name)

        if name is not startfield:
            
            #print(name)

            if os.path.exists(name+'_'+str(wavelength)+'um_'+CatName+'.fits'):
                newtab1=Table.read(name+'_'+str(wavelength)+'um_'+CatName+'.fits')

                if 'sky_centroid' in newtab1.columns:
                    #fix for issue with combining tables with objects as columns
                    newtab1['RA(J2000)']=newtab1['sky_centroid'].ra.to_string(u.hour,precision=2)
                    newtab1['DEC(J2000)']=newtab1['sky_centroid'].dec.to_string(u.degree,precision=2)
                    newtab1['glon']=newtab1['sky_centroid'].galactic.l.deg
                    newtab1['glat']=newtab1['sky_centroid'].galactic.b.deg
                    newtab1.remove_column('sky_centroid')
                elif 'sky_center' in newtab1.columns:
                    #fix for issue with combining tables with objects as columns
                    newtab1['RA(J2000)']=newtab1['sky_center'].ra.to_string(u.hour,precision=2)
                    newtab1['DEC(J2000)']=newtab1['sky_center'].dec.to_string(u.degree,precision=2)
                    newtab1['glon']=newtab1['sky_center'].galactic.l.deg
                    newtab1['glat']=newtab1['sky_center'].galactic.b.deg
                    newtab1.remove_column('sky_center')

                #combine tables
                cat=vstack([cat,newtab1])


    #print the table sizes to get source counts
    #print('Raw number of combinded sources: ', len(cat))
    
    #sort table entries by gal. lat. 
    cat.sort('glon')
    
    #fix to wrap galactic coordinates
    alttabidx2=cat['glon']>180
    alttabidx1=cat['glon']<180
    c1=cat[alttabidx2]
    c2=cat[alttabidx1]
    catnew=vstack([c1,c2])
    
    #rename field based id to anothter name to avoid confusion
    #catnew.rename_column('id', 'old_id')
    
    #remove original source 'id' values because they are field based and not useful in combined table
    catnew.remove_column('id')

    #add new ids to the master catalog
    catnew['id']=np.linspace(1,len(catnew),len(catnew),dtype=np.int64)
    
    #re-add "sky_centroid" column 
    catnew['sky_centroid']=SkyCoord(catnew['RA(J2000)'],catnew['DEC(J2000)'],unit=(u.degree,u.hour))

    #return combined catalog
    return catnew
    
#source coords should be skycoord objects, sep should be an angular quantity e.g., 4*u.arcsec
def remove_duplicates(cat, sep):
    
    #check cat type to specify which SNR column to use
    if "aper_snr_4pix" in cat.columns:
        colname='aper_snr_4pix'
    elif "aper_snr_5.5pix" in cat.columns:
        colname='aper_snr_5.5pix'
    
    #get sky centroid
    sourcecoords=cat['sky_centroid']

    #create KD Tree
    kdt=_get_cartesian_kdtree(sourcecoords)
    
    #search radius
    r = (2 * np.sin(Angle(sep) / 2.0)).value #search radius
    
    counter=0
    #create array to hold matches
    matchstore=np.zeros(len(cat),dtype=np.int64)
    
    #loop through and find all matching sources within crossmatch radius
    for i, matches in enumerate(kdt.query_ball_tree(kdt,r)):
        if len(matches)>1:
            counter+=1
            
            for match in matches:
                #print(match)
                
                if match == matches[0]:
                    other=matches[1]
                else:
                    other=matches[0]
                
                matchstore[match]=cat['id'][other]
                
    #store possible crossmatches
    cat['selfXmatch']=matchstore
    print('number of likely duplicates: ', counter)
    
    

    #create keep and remove lists
    keep=[]
    remove=[]
    
    #loop through and mark sources to keep and remove based on which has the greater aperture photometry snr
    for row in cat:
        if row['selfXmatch']>0:
            #print(row['Master_id'])
            sidx1 = row['id']
            sidx2 = row['selfXmatch']

            row1=cat[cat['id']==sidx1]
            row2=cat[cat['id']==sidx2]
            
            snr1=row1[colname].data[0]
            if len(row2[colname].data)>0:
                snr2=row2[colname].data[0]
            else:
                snr2=0
                
            if np.ma.is_masked(snr1):
                snr1=-1
                
            if np.ma.is_masked(snr2):
                snr2=-1
                
            if snr1>snr2:
                keep.append(sidx1)
            else:
                remove.append(sidx1)
    rid=[]
    
    for sid in remove:
        arr=cat['id']==sid
        tval=np.where(arr)[0]
        rid.append(tval.data[0])
        
    if len(rid)>0:
        cat.remove_rows(rid)
    
    return cat, counter    
    
    
    
    
#create catalog names for sources based on galactic coords
def createNames(tab):
    #get lat and long values
    latvals=tab['glat'].value
    lonvals=tab['glon'].value
    #convert to strings
    strings1 = ["%.3f" % lonvals for lonvals in lonvals] #updated to 4 decimals because 3 isn't sufficient to distinguish some sources from one another
    strings2 = ["%.4f" % latvals for latvals in latvals] #updated to 4 decimals because 3 isn't sufficient to distinguish some sources from one another
    #update string 2 to include '+' signs where needed
    s2update=[]
    for string in strings2:
        if not string.startswith('-'):
            string='+'+string
        s2update.append(string)
    #create names array
    names=['SFGC'+strings1+s2update for strings1,s2update in zip(strings1,s2update)]
    #add names to table
    tab['SourceID']=names
    #return table with names added
    return tab
    


#find the 'not' indcies for sources in table
def findNOTindex(tab,index):
    allpos=np.linspace(0,len(tab)-1,len(tab),dtype=np.int64) 
    notindex=list(set(allpos)-set(index))
    return notindex



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

'''
# A few useful functions for creating ds9 files
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
'''


def makeDS9reg(savename, table, radius, color='green',label=None):
    
    
    if 'sky_centroid' in table.columns:    
        scs=table['sky_centroid']
    elif 'skycoords' in table.columns:
        scs=table['skycoords']
    else:
        print('Table is missing colum for sky coordinates. Please add column \'sky_centroid\' ')

    decimals = 1 # number of decimal places to display with labels
    
    #loop through and create region instances for each source
    regions=[]
    for i in range(0,len(scs)):
        region = CircleSkyRegion(scs[i], radius)
        regions.append(region)
            #write out region file
    write_ds9(regions, savename)
    
    
    hld=""
    with open(savename, 'r+') as f:
        lineno=0
        for line in f:
            if lineno<=(len(table)+1):
                if label is not None:
                    text=line.replace(r')', r') # color='+color+' text={'+str(np.round(table[label][lineno-2],decimals))+'}')
                else:
                    text=line.replace(r')', r') # color='+color)
                hld=hld+text
                lineno=lineno+1
        f.seek(0)
        f.write(hld)
        f.truncate()
        
    return hld
        
        
        

def CombDS9files(text1,text2,outname):
    newregtext=text1+text2[45:]

    with open(outname, 'w+') as f:
        f.seek(0)
        f.write(newregtext)
        f.truncate()


    
    
    
    
        
        

