#config file for running different fields Contains all relevant parameters to be loaded in by the jupyter notebook

import numpy

#data path
dpath='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/SgrBandNearby'
dpathalt='E:\Documents\SofiaLegacyProgram\CAL_files\CAL_R'

#file name
fname='F0588_FO_IMA_07018936_FORF253_CAL_0008-0041_Field_13.fits'

#Name of field
field='Field13'

#wavelength of data
wavelength=25.2

#normalized exposure time map cutoff
tmapnormcut=0.2

#use mask2
usemask2=True

#limits for mask 2
#m2x1=0
#m2x2=330
#m2y1=250
#m2y2=300

mask2=np.zeros(np.shape(mask))  #mask 2 is applied to the first segmentation map instance for source detection
#mask2[0:300,0:100]=1




#use mask3
usemask3=True

#limits for mask 3
#m3x1=0
#m3x2=300
#m3y1=0
#m3y2=70

mask3=np.zeros(np.shape(mask))  #mask 3 is applied to the deblended segmentation map for source detection - this may or may not be needed if the first 2 masks work ok
#mask3[0:330,0:70]=1

#bkg esitmator box size
bkgbox=20