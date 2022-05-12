#!/usr/bin/env python

#import all required packages
import os
import numpy as np
import re
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.table import join, vstack, Table
from astropy.coordinates import SkyCoord, search_around_sky, Angle

from photutils.aperture import SkyCircularAperture,SkyCircularAnnulus,aperture_photometry 
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
from photutils.background import Background2D, MedianBackground, SExtractorBackground, MMMBackground
from photutils.utils import calc_total_error

from astropy.coordinates.matching import _get_cartesian_kdtree
from regions import read_ds9, write_ds9, CircleSkyRegion

#import configuration for selected file
#from config import wavelength, segdetsig, finddetsig, bkgbox #import additional common paramters
from config import dpath, dpathalt, ds9path #import additional common paramters
from config import *

import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


interactive=False

SegMap=True



####-------------------------functions-------------------------------------

def CombineFieldResults(CatName,wavelength): # #Name options- 'CombCat', 'SegCat', 'DaoCat'
    #pick a field to start with so we can get the table structure. 
    #The files for the field cannot be empty or an error message will display
    startfield=FieldA.name

    if os.path.exists(startfield+'_'+str(wavelength)+'um_'+CatName+'.fits'):
        cat=Table.read(startfield+'_'+str(wavelength)+'um_'+CatName+'.fits')

        #fix for issue with combining tables with objects as columns
        cat['RA(J2000)']=cat['sky_centroid'].ra
        cat['DEC(J2000)']=cat['sky_centroid'].dec
        cat.remove_column('sky_centroid')

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

                #fix for issue with combining tables with objects as columns
                newtab1['RA(J2000)']=newtab1['sky_centroid'].ra
                newtab1['DEC(J2000)']=newtab1['sky_centroid'].dec
                newtab1.remove_column('sky_centroid')

                #combine tables
                cat=vstack([cat,newtab1])

    #re-add "sky_centroid" column 
    cat['sky_centroid']=SkyCoord(cat['RA(J2000)'],cat['DEC(J2000)'],unit=u.deg)

    #print the table sizes to get source counts
    #print('Raw number of combinded sources: ', len(cat))
    
    #rename field based id to anothter name to avoid confusion
    cat.rename_column('id', 'old_id')

    #add new ids to the master catalog
    cat['id']=np.linspace(1,len(cat),len(cat),dtype=np.int64)

    #return combined catalog
    return cat
    

#source coords should be skycoord objects, sep should be an angular quantity e.g., 4*u.arcsec
def remove_duplicates(cat, sep):
	sourcecoords=cat['sky_centroid']

	#create KD Tree
	kdt=_get_cartesian_kdtree(sourcecoords)
    
	#search radius
	r = (2 * np.sin(Angle(sep) / 2.0)).value #search radius
    
	#initialize counter
	counter=0
    
	#create array to hold matches
	matchstore=np.zeros(len(cat),dtype=np.int)
    
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
			
			snr1=row1['aper_snr_4pix'].data[0]
			snr2=row2['aper_snr_4pix'].data[0]
			           
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


def makeDS9file(savename, table, radius, color='green', labelon=False,label='ColumnName'):
    scs=table['sky_centroid']
    
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
            if lineno<len(table):
                if labelon:
                    text=line.replace(r')', r') # color='+color+' text={'+str(np.round(table[label][lineno-2],decimals))+'}')
                else:
                    text=line.replace(r')', r') # color='+color)
                hld=hld+text
                lineno=lineno+1
        f.seek(0)
        f.write(hld)
        f.truncate()

#find the 'not' indcies for sources in table
def findNOTindex(tab,index):
    allpos=np.linspace(0,len(tab)-1,len(tab),dtype=np.int64) 
    notindex=list(set(allpos)-set(index))
    return notindex    
    
####----------------------Start of main----------------------------------    


#change directory to where data is
try:
    os.chdir(dpath)
except:
    os.chdir(dpathalt)


#create initial master catalog at 25 um
mastercat25=CombineFieldResults('CombCat',25)


#Now lets look at possible duplications that may exist because of overlapping data coverage
print('Number of 25 um sources with duplicates included', len(mastercat25))
mastercat25=remove_duplicates(mastercat25, 3.0*u.arcsec)
print('Number of 25 um sources with duplicates removed', len(mastercat25))


if interactive:
    #look at SNR distrobution
    binlist=np.linspace(1,40,20)
    plt.figure(figsize=(12,7))
    plt.title('Histogram of 4pix Ap SNR for 25 um Sources')
    plt.hist(mastercat25['aper_snr_4pix'],bins=binlist)
    plt.xlabel('SNR of Source')
    plt.ylabel('Number of Sources')
    #plt.xlim(0,25)
    plt.show()

#lets examine SNR cuts...
snrcut4pix=mastercat25['aper_snr_4pix']>=5.0 #max snr in all computed apertures must be gtreq to 5
mcat25snrcut=mastercat25[snrcut4pix] #apply snr cut

print("Number of Catalog Sources after SNR cut: ", len(mcat25snrcut))

#set size of regions 
r = Angle(0.00083333, u.deg) #must be in degrees - current value is r=3"
#write out ds9 files

#write files
makeDS9file('mastercatComb_4pixSNRselect_25m_NoLabel.reg', mcat25snrcut, r, color='green')
mcat25snrcut.write('masterCat_25um_Step3.fits',overwrite=True)    

###--------------------------------------------------------------------------------------
#create initial master catalog at 37 um
mastercat37=CombineFieldResults('CombCat',37)

#Now lets look at possible duplications that may exist because of overlapping data coverage
print('Number of 37 um sources with duplicates included', len(mastercat37))
mastercat37=remove_duplicates(mastercat37, 3.0*u.arcsec)
print('Number of 37 um sources with duplicates removed', len(mastercat37))


if interactive:
    #look at SNR distrobution
    binlist=np.linspace(1,40,20)
    plt.figure(figsize=(12,7))
    plt.title('Histogram of 4pix Ap SNR for 37 um Sources')
    plt.hist(mastercat37['aper_snr_4pix'],bins=binlist)
    plt.xlabel('SNR of Source')
    plt.ylabel('Number of Sources')
    #plt.xlim(0,37)
    plt.show()

#lets examine SNR cuts...
snrcut4pix=mastercat37['aper_snr_4pix']>=5.0 #max snr in all computed apertures must be gtreq to 5
mcat37snrcut=mastercat37[snrcut4pix] #apply snr cut

print("Number of Catalog Sources after SNR cut: ", len(mcat37snrcut))    

#write files
makeDS9file('mastercatComb_4pixSNRselect_37m_NoLabel.reg', mcat37snrcut, r, color='green')
mcat37snrcut.write('masterCat_37um_Step3.fits',overwrite=True)  
    
###--------------------------------------------------------------------------------------  
    
#get source coordinates from both tables
sources25=mcat25snrcut['sky_centroid']
sources37=mcat37snrcut['sky_centroid']

#crossmatch source lists to look for duplication
idx,rdx, d2d, d3d = sources25.search_around_sky(sources37, 3*u.arcsec) #use larger search radius?
print('Number of crossmatched 25/37 sources found: ', len(idx))    
    
matched25=mcat25snrcut[rdx]
matched37=mcat37snrcut[idx]    
    
notrdx=findNOTindex(mcat25snrcut,rdx)
only25=mcat25snrcut[notrdx]

notidx=findNOTindex(mcat37snrcut,idx)
only37=mcat37snrcut[notidx]    
    
    
t1=Table()
t1['RA(J2000)']=matched25['RA(J2000)']
t1['DEC(J2000)']=matched25['DEC(J2000)']
t1['25um Flux (Jy)']=matched25['aper_sum_bkgsub_4pix']
t1['fwhm25']=matched25['fwhm']
t1['25um SNR']=matched25['aper_snr_4pix']
t1['37um Flux (Jy)']=matched37['aper_sum_bkgsub_4pix']
t1['37um SNR']=matched37['aper_snr_4pix']
t1['fwhm37']=matched37['fwhm']


t2=Table()
t2['RA(J2000)']=only25['RA(J2000)']
t2['DEC(J2000)']=only25['DEC(J2000)']
t2['25um Flux (Jy)']=only25['aper_sum_bkgsub_4pix']
t2['25um SNR']=only25['aper_snr_4pix']
t2['fwhm25']=only25['fwhm']
t2['37um Flux (Jy)']=-1
t2['37um SNR']=-1
t2['fwhm37']=-1


t3=Table()
t3['RA(J2000)']=only37['RA(J2000)']
t3['DEC(J2000)']=only37['DEC(J2000)']
t3['25um Flux (Jy)']=-1
t3['25um SNR']=-1
t3['fwhm25']=-1
t3['37um Flux (Jy)']=only37['aper_sum_bkgsub_4pix']
t3['37um SNR']=only37['aper_snr_4pix']
t3['fwhm37']=only37['fwhm']    
    
mastercat=vstack((t1,t2,t3))    

#change format of columns to save fewer decimal places
for col in mastercat.colnames:
    if col!='RA(J2000)' and col!='DEC(J2000)' and col!='type': #skip columns that aren't relevant
        mastercat[col].info.format = '%.4G'

#write out final catalog
mastercat.write('masterCat_step3_final.fits',overwrite=True)



    
###--------------------------------------------------------------------------------------  

if SegMap==True:
    #create combined field segment map source catalog at 25 um
    Segcat25=CombineFieldResults('segCat',25)

    #Remove duplications from overlapping fields
    SegCat25=remove_duplicates(Segcat25, 3.0*u.arcsec)

    #lets examine SNR cuts...
    snrcut4pix=SegCat25['aper_snr_4pix']>=5.0 #max snr in all computed apertures must be gtreq to 5
    Segcat25snr=SegCat25[snrcut4pix] #apply snr cut
    print("Number of 25 um Seg Catalog Sources after SNR cut: ", len(Segcat25snr))

    #write ds9 file
    makeDS9file('Segcat_4pixSNRselect_25m_NoLabel.reg', Segcat25snr, r, color='green')

    #write fits file
    Segcat25snr.write('SegCat_CombinedFields_snrcut.fits',overwrite=True)

    #create combined field segment map source catalog at 37 um
    Segcat37=CombineFieldResults('segCat',37)

    #Remove duplications from overlapping fields
    SegCat37=remove_duplicates(Segcat37, 3.0*u.arcsec)

    #lets examine SNR cuts...
    snrcut4pix=SegCat37['aper_snr_4pix']>=5.0 #max snr in all computed apertures must be gtreq to 5
    Segcat37snr=SegCat37[snrcut4pix] #apply snr cut
    print("Number of 37 um Seg Catalog Sources after SNR cut: ", len(Segcat37snr))

    #write ds9 file
    makeDS9file('Segcat_4pixSNRselect_37m_NoLabel.reg', Segcat37snr, r, color='green')

    #write fits file
    Segcat37snr.write('SegCat_CombinedFields_snrcut.fits',overwrite=True)
    
    
    
    
    
    
    
    
    