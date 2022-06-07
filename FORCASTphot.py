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



def performApPhoto(data,errormap,tmap,header,sourceCoords,radii,rin,rout,plot=True):
    
    #define wcs object for header
    wcs=WCS(header)

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

        try:
            #do statistics
            annulus_data_1d = annulus_data[maskdata > 0]
            meansc, median_sigclip, stdsc = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
            bkg_mean.append(meansc)
            bkg_std.append(stdsc)
            appmasks.append(mask.data)
        except TypeError:
            bkg_median.append(np.nanmedian(annulus_data_1d))
            bkg_mean.append(np.nanmean(annulus_data_1d))
            bkg_std.append(np.nanstd(annulus_data_1d))
            appmasks.append(np.nan)
                
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
        newcol2='aper_'+str(radii[i])+'pix_noise'
        newcol3='aper_snr_'+str(radii[i])+'pix'
        newcol4='aper_area_'+str(radii[i])+'pix'
        
        #get pixel aperture areas for calculations
        pixel_ap = apertures[i].to_pixel(wcs)
        pixarea=pixel_ap.area
        #get background subracted photo by subtracting median annulus value from each pixel in aperture
        phot_table[newcol1]=(phot_table[cname1]/pixarea-phot_table['ann_bkg_med'])*pixarea
        #calculate noise following equation from forcast photometry cookbook -https://sofia-data-analysis-cookbooks.readthedocs.io/en/latest/FORCAST-photometry_detailed.html - Note there is an added error from the aperture extraction included (cname2)  
        phot_table[newcol2]=np.sqrt(2*np.pi*(radii[i]*phot_table['ann_bkg_std']/phot_table[newcol1])**2+(header['ERRCALF']/header['CALFCTR'])**2+0.0025+phot_table[cname2]**2)
        #calculate SNR 
        phot_table[newcol3]=phot_table[newcol1]/phot_table[newcol2]
        phot_table[newcol4]=pixarea

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
    plt.show()

    return phot_table




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



def CombineFieldResults(CatName,wavelength): # #Name options- 'CombCat', 'SegCat', 'DaoCat'
    #pick a field to start with so we can get the table structure. 
    #The files for the field cannot be empty or an error message will display
    startfield=Field01.name

    if os.path.exists(startfield+'_'+str(wavelength)+'um_'+CatName+'.fits'):
        cat=Table.read(startfield+'_'+str(wavelength)+'um_'+CatName+'.fits')

        if 'sky_centroid' in cat.columns:

            #fix for issue with combining tables with objects as columns
            cat['RA(J2000)']=cat['sky_centroid'].ra
            cat['DEC(J2000)']=cat['sky_centroid'].dec
            cat['glon']=cat['sky_centroid'].galactic.l.deg
            cat['glat']=cat['sky_centroid'].galactic.b.deg
            cat.remove_column('sky_centroid')
            
        elif 'sky_center' in cat.columns:
            
            #fix for issue with combining tables with objects as columns
            cat['RA(J2000)']=cat['sky_center'].ra
            cat['DEC(J2000)']=cat['sky_center'].dec
            cat['glon']=cat['sky_center'].galactic.l.deg
            cat['glat']=cat['sky_center'].galactic.b.deg
            cat.remove_column('sky_center')

    else:
        print('Error... Must pick another file to start with')
        
        
    #loop through all the saved photometry tables for individual fields and append them together
    for info in field._registry:
        filename=info.filename
        name=info.name

        #print('\nLoading in photometry data from field: ', name)

        if name is not startfield:

            if os.path.exists(name+'_'+str(wavelength)+'um_'+CatName+'.fits'):
                newtab1=Table.read(name+'_'+str(wavelength)+'um_'+CatName+'.fits')

                if 'sky_centroid' in newtab1.columns:
                    #fix for issue with combining tables with objects as columns
                    newtab1['RA(J2000)']=newtab1['sky_centroid'].ra
                    newtab1['DEC(J2000)']=newtab1['sky_centroid'].dec
                    newtab1['glon']=newtab1['sky_centroid'].galactic.l.deg
                    newtab1['glat']=newtab1['sky_centroid'].galactic.b.deg
                    newtab1.remove_column('sky_centroid')
                elif 'sky_center' in newtab1.columns:
                    #fix for issue with combining tables with objects as columns
                    newtab1['RA(J2000)']=newtab1['sky_center'].ra
                    newtab1['DEC(J2000)']=newtab1['sky_center'].dec
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
    catnew['sky_centroid']=SkyCoord(catnew['RA(J2000)'],catnew['DEC(J2000)'],unit=u.deg)

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
    
    #initialize counter
    counter=0
    
    #create array to hold matches
    matchstore=np.zeros(len(cat),dtype=np.int64)
    
    #loop through and find all matching sources within crossmatch radius
    for i, matches in enumerate(kdt.query_ball_tree(kdt,r)):
        if len(matches)>1:
            counter+=1

            for match in matches:
                if match!=matches[(match+1)%len(matches)]:
                    matchstore[match]=matches[(match+1)%len(matches)]+1
                else:
                    matchstore[match]=matches[(match+2)%len(matches)]+1

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
            idx1 = row['id']
            idx2 = row['selfXmatch']

            row1=cat[cat['id']==idx1]
            row2=cat[cat['id']==idx2]
            
            snr1=row1[colname].data[0]
            snr2=row2[colname].data[0]
           
            if snr1>snr2:
                keep.append(idx1)
                remove.append(idx2)
            else:
                keep.append(idx2)
                remove.append(idx1)

    #remove duplicates in lists
    keep=list(set(keep))
    remove=list(set(remove))
    
    #remove duplicates in table
    removeIdx=np.array(remove)-1 #fix zero/one initialization issue
    cat.remove_rows(removeIdx)
    
    return cat
    
    
    
#create catalog names for sources based on galactic coords
def createNames(tab):
    #get lat and long values
    latvals=tab['glat'].value
    lonvals=tab['glon'].value
    #convert to strings
    strings1 = ["%.4f" % lonvals for lonvals in lonvals]
    strings2 = ["%.4f" % latvals for latvals in latvals]
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
        
        
        

def CombDS9files(text1,text2,text3,outname):
    newregtext=text1+text2[45:]+text3[45:]

    with open(outname, 'w+') as f:
        f.seek(0)
        f.write(newregtext)
        f.truncate()



