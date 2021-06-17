#config file for running different fields Contains all relevant parameters to be loaded in by the jupyter notebook

#data path
dpath='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/SgrAandNearby'
dpathalt='E:\Documents\SofiaLegacyProgram\CAL_files\CAL_C'

#file name
fname='F0590_FO_IMA_07018918_FORF253_CAL_0117-0141_Field_G.fits'

#Name of field
field='FieldG'

#wavelength of data
wavelength=25.2

#normalized exposure time map cutoff
tmapnormcut=0.01

#use mask2
usemask2=True
# list of regions to define mask2 in form [y1:y2,x1,x2]
m2limits=[[0,100,0,250],[0,90,250,350],[350,380,0,200],[150,230,300,380]] #possibly also [150,250,0,30]

#use mask3
usemask3=False
# list of regions to define mask3 in form [y1:y2,x1,x2]
#m3limits=[[290,370,0,370]

#bkg esitmator box size
bkgbox=20