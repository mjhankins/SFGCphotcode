#config file for running different fields Contains all relevant parameters to be loaded in by the jupyter notebook

#data path
dpath='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/SgrBandNearby'
dpathalt='E:\Documents\SofiaLegacyProgram\CAL_files\CAL_R'

#file name
fname='F0591_FO_IMA_07018935_FORF253_CAL_0039-0061_Field_12.fits'

#Name of field
field='Field12'

#wavelength of data
wavelength=25.2

#normalized exposure time map cutoff
tmapnormcut=0.1

#use mask2 - applied to first segmentation map
usemask2=True
# list of regions to define mask2 in form [y1:y2,x1,x2]
m2limits=[[0,50,0,350],[0,350,0,40],[300,350,40,80]]

#use mask3
usemask3=False
# list of regions to define mask3 in form [y1:y2,x1,x2]
#m3limits=[[290,370,0,370]

#bkg esitmator box size
bkgbox=20