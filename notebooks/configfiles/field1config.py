#config file for running different fields Contains all relevant parameters to be loaded in by the jupyter notebook

#data path
dpath='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/SgrCandNearby'
dpathalt='E:\Documents\SofiaLegacyProgram\CAL_files\CAL_R'

#file name
fname='F0592_FO_IMA_0701891_FORF253_CAL_0095-0118_FIELD_1.fits'

#Name of field
field='Field1'

#wavelength of data
wavelength=25.2

#normalized exposure time map cutoff
tmapnormcut=0.2

#use mask2 - applied to first segmentation map
usemask2=True
# list of regions to define mask2 in form [y1:y2,x1,x2]
m2limits=[[280,360,0,370],[0,350,0,50],[200,350,300,370]]

#use mask3
usemask3=False
# list of regions to define mask3 in form [y1:y2,x1,x2]
#m3limits=[[290,370,0,370]

#bkg esitmator box size
bkgbox=20
