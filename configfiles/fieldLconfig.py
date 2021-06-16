#config file for running different fields Contains all relevant parameters to be loaded in by the jupyter notebook

#data path
dpath='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/SgrAandNearby'
dpathalt='E:\Documents\SofiaLegacyProgram\CAL_files\CAL_C'

#file name
fname='F0589_FO_IMA_07018923_FORF253_CAL_0050-0079_Field_L.fits'

#Name of field
field='FieldL'

#wavelength of data
wavelength=25.2

#normalized exposure time map cutoff
tmapnormcut=0.5

#use mask2
usemask2=False
# list of regions to define mask2 in form [y1:y2,x1,x2]
#m2limits=[[290,370,0,370]

#use mask3
usemask3=False
# list of regions to define mask3 in form [y1:y2,x1,x2]
#m3limits=[[290,370,0,370]

#bkg esitmator box size
bkgbox=20