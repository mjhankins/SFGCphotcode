#configuration file for the photometry code

#some common parameters for all files
dpath='/home/matt/Data/SLP_GC_data/CAL_files/all'  #path to data
dpathalt='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/all'   #alternate path to data

ds9path='/home/matt/anaconda3/bin/ds9' #path to envoke ds9. Use None if you don't want to open ds9 at the end. If you are unsure use the command 'which ds9' inside of your python environment to find the correct ds9 install. 

wavelength=25  #wavelength of the observations - 25 or 37
segdetsig=3.0  #detection threshold sigma for segmentation map
finddetsig=3.5 #detection threshold sigma for daofinder and starfinder
bkgbox=13 #size of box used for background model

class mosaic():
    _registry = []
    
    # init method or constructor 
    def __init__(self, name, filename,m1cut,m2lims,m3lims):
        self.name=name
        self.filename=filename
        self.m1cut=m1cut
        self.m2lims=m2lims
        self.m3lims=m3lims
        self._registry.append(self)

Mosaic25=mosaic(
name='Mosaic25',  #name of field
filename='F0217_FO_IMA_70030015_FORF253_MOS_0001-0348_final_MATT_Corrected.fits', #file name
m1cut=0.01, # mask cutoff for normalized exposure time map
m2lims=None, #list of manually defined masked regions in the form [y1,y2,x1,x2]
m3lims=None) #additional mask if needed in the form [y1,y2,x1,x2]

Mosaic37=mosaic(
name='Mosaic37',  #name of field
filename='F0217_FO_IMA_70030016_FORF371_MOS_0001-0348_final_MATT_Corrected.fits', #file name
m1cut=0.01, # mask cutoff for normalized exposure time map
m2lims=None, #list of manually defined masked regions in the form [y1,y2,x1,x2]
m3lims=None) #additional mask if needed in the form [y1,y2,x1,x2]

class field():
    _registry = []
    
    # init method or constructor 
    def __init__(self, name, filename,m1cut,m2lims,m3lims):
        self.name=name
        self.filename=filename
        self.m1cut=m1cut
        self.m2lims=m2lims
        self.m3lims=m3lims
        self._registry.append(self)

#file names for field1 
file25='F0592_FO_IMA_0701891_FORF253_CAL_0095-0118.fits'
file37='F0592_FO_IMA_0701891_FORF371_CAL_0095-0118.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field1
Field1=field(
name='Field1',  #name of field
filename=fname, #file name
m1cut=0.2, # mask cutoff for normalized exposure time map
m2lims=[[280,360,0,370],[0,350,0,50],[200,350,300,370]], #list of manually defined masked regions in the form [y1,y2,x1,x2]
m3lims=None) #additional mask if needed in the form [y1,y2,x1,x2]

#file names for field2 
file25='F0592_FO_IMA_0701892_FORF253_CAL_0120-0141.fits'
file37='F0592_FO_IMA_0701892_FORF371_CAL_0120-0141.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field2
Field2=field(
name='Field2',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field3 
file25='F0595_FO_IMA_0701893_FORF253_CAL_0283-0303.fits'
file37='F0595_FO_IMA_0701893_FORF371_CAL_0283-0303.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Feild3
Field3=field(
name='Field3',
filename=fname,
m1cut=0.5,
m2lims=[[0,150,250,300]],
m3lims=None)

#file names for field5 
file25='F0595_FO_IMA_0701895_FORF253_CAL_0259-0282.fits'
file37='F0595_FO_IMA_0701895_FORF371_CAL_0259-0282.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field5
Field5=field(
name='Field5',
filename=fname,
m1cut=0.5,
m2lims=[[0,330,250,300]],
m3lims=[[0,330,0,70]])

#file names for field6
file25='F0593_FO_IMA_0701896_FORF253_CAL_0099-0124.fits'
file37='F0593_FO_IMA_0701896_FORF371_CAL_0099-0124.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field6
Field6=field(
name='Field6',
filename=fname,
m1cut=0.1,
m2lims=[[290,370,0,370],[0,370,290,380],[0,30,0,370],[0,370,0,50]],
m3lims=None)

#file names for field7 
file25='F0593_FO_IMA_0701897_FORF253_CAL_0125-0136.fits'
file37='F0593_FO_IMA_0701897_FORF371_CAL_0125-0136.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field7
Field7=field(
name='Field7',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field8 
file25='F0591_FO_IMA_0701898_FORF253_CAL_0006-0038.fits'
file37='F0591_FO_IMA_0701898_FORF371_CAL_0006-0038.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field8
Field8=field(
name='Field8',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field9 
file25='F0588_FO_IMA_0701899_FORF253_CAL_0042-0081.fits'
file37='F0588_FO_IMA_0701899_FORF371_CAL_0042-0081.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field9
Field9=field(
name='Field9',
filename=fname,
m1cut=0.2,
m2lims=[[0,230,0,40],[0,50,0,50],[0,300,280,320]],
m3lims=None)

#file names for field10 
file25='F0588_FO_IMA_07018910_FORF253_CAL_0082-0111.fits'
file37='F0588_FO_IMA_07018910_FORF371_CAL_0082-0111.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field10
Field10=field(
name='Field10',
filename=fname,
m1cut=0.2,
m2lims=None,
m3lims=None)

#file names for field11
file25='F0588_FO_IMA_07018911_FORF253_CAL_0114-0124.fits'
file37='F0588_FO_IMA_07018911_FORF371_CAL_0114-0124.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field11
Field11=field(
name='Field11',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field12
file25='F0591_FO_IMA_07018935_FORF253_CAL_0039-0061.fits'
file37='F0591_FO_IMA_07018935_FORF371_CAL_0039-0061.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field12
Field12=field(
name='Field12',
filename=fname,
m1cut=0.1,
m2lims=[[0,50,0,350],[0,350,0,40],[300,350,40,80]],
m3lims=None)

#file names for field13
file25='F0588_FO_IMA_07018936_FORF253_CAL_0008-0041.fits'
file37='F0588_FO_IMA_07018936_FORF371_CAL_0008-0041.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field13
Field13=field(
name='Field13',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field A
file25='F0588_FO_IMA_07018912_FORF253_CAL_0301-0325.fits'
file37='F0588_FO_IMA_07018912_FORF371_CAL_0301-0325.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field A
FieldA=field(
name='FieldA',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field B
file25='F0591_FO_IMA_07018913_FORF253_CAL_0088-0114.fits'
file37='F0591_FO_IMA_07018913_FORF371_CAL_0088-0114.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field B
FieldB=field(
name='FieldB',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field C
file25='F0588_FO_IMA_07018914_FORF253_CAL_0326-0348.fits'
file37='F0588_FO_IMA_07018914_FORF371_CAL_0326-0348.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field C
FieldC=field(
name='FieldC',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field D
file25='F0592_FO_IMA_07018915_FORF253_CAL_0276-0303.fits'
file37='F0592_FO_IMA_07018915_FORF371_CAL_0276-0303.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field D
FieldD=field(
name='FieldD',
filename=fname,
m1cut=0.2,
m2lims=[[200,250,0,50]],
m3lims=None)

#file names for field E
file25='F0590_FO_IMA_07018916_FORF253_CAL_0142-0161.fits'
file37='F0590_FO_IMA_07018916_FORF371_CAL_0142-0161.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field E
FieldE=field(
name='FieldE',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field F
file25='F0592_FO_IMA_07018917_FORF253_CAL_0068-0094.fits'
file37='F0592_FO_IMA_07018917_FORF371_CAL_0068-0094.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field F
FieldF=field(
name='FieldF',
filename=fname,
m1cut=0.3,
m2lims=None,
m3lims=None)

#file names for field G
file25='F0590_FO_IMA_07018918_FORF253_CAL_0117-0141.fits'
file37='F0590_FO_IMA_07018918_FORF371_CAL_0117-0141.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field G
FieldG=field(
name='FieldG',
filename=fname,
m1cut=0.01,
m2lims=[[0,100,0,250],[0,90,250,350],[350,380,0,200],[150,230,300,380]],
m3lims=None)

#file names for field H
file25='F0589_FO_IMA_07018919_FORF253_CAL_0017-0049.fits'
file37='F0589_FO_IMA_07018919_FORF371_CAL_0017-0049.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field H
FieldH=field(
name='FieldH',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field I
file25='F0592_FO_IMA_07018920_FORF253_CAL_0005-0036.fits'
file37='F0592_FO_IMA_07018920_FORF371_CAL_0005-0036.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field I
FieldI=field(
name='FieldI',
filename=fname,
m1cut=0.3,
m2lims=None,
m3lims=None)

#file names for field K
file25='F0590_FO_IMA_07018922_FORF253_CAL_0087-0116.fits'
file37='F0590_FO_IMA_07018922_FORF371_CAL_0087-0116.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field K
FieldK=field(
name='FieldK',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field L
file25='F0589_FO_IMA_07018923_FORF253_CAL_0050-0079.fits'
file37='F0589_FO_IMA_07018923_FORF371_CAL_0050-0079.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field L
FieldL=field(
name='FieldL',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field M
file25='F0593_FO_IMA_07018924_FORF253_CAL_0071-0098.fits'
file37='F0593_FO_IMA_07018924_FORF371_CAL_0071-0098.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field M
FieldM=field(
name='FieldM',
filename=fname,
m1cut=0.2,
m2lims=[[230,330,0,80]],
m3lims=None)

#file names for field O
file25='F0590_FO_IMA_07018926_FORF253_CAL_0057-0086.fits'
file37='F0590_FO_IMA_07018926_FORF371_CAL_0057-0086.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field O
FieldO=field(
name='FieldO',
filename=fname,
m1cut=0.4,
m2lims=None,
m3lims=None)

#file names for field P
file25='F0589_FO_IMA_07018927_FORF253_CAL_0138-0169.fits'
file37='F0589_FO_IMA_07018927_FORF371_CAL_0138-0169.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field P
FieldP=field(
name='FieldP',
filename=fname,
m1cut=0.3,
m2lims=[[350,380,0,370]],
m3lims=None)

#file names for field Q
file25='F0590_FO_IMA_07018928_FORF253_CAL_0024-0055.fits'
file37='F0590_FO_IMA_07018928_FORF371_CAL_0024-0055.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field Q
FieldQ=field(
name='FieldQ',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field R
file25='F0589_FO_IMA_07018929_FORF253_CAL_0080-0106.fits'
file37='F0589_FO_IMA_07018929_FORF371_CAL_0080-0106.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field R
FieldR=field(
name='FieldR',
filename=fname,
m1cut=0.2,
m2lims=None,
m3lims=None)

#file names for field S
file25='F0592_FO_IMA_07018930_FORF253_CAL_0037-0067.fits'
file37='F0592_FO_IMA_07018930_FORF371_CAL_0037-0067.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field S
FieldS=field(
name='FieldS',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field T
file25='F0589_FO_IMA_07018931_FORF253_CAL_0108-0137.fits'
file37='F0589_FO_IMA_07018931_FORF371_CAL_0108-0137.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field T
FieldT=field(
name='FieldT',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field U
file25='F0591_FO_IMA_07018932_FORF253_CAL_0062-0087.fits'
file37='F0591_FO_IMA_07018932_FORF371_CAL_0062-0087.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field U
FieldU=field(
name='FieldU',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field V
file25='F0593_FO_IMA_07018933_FORF253_CAL_0039-0070.fits'
file37='F0593_FO_IMA_07018933_FORF371_CAL_0039-0070.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field V
FieldV=field(
name='FieldV',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field W
file25='F0593_FO_IMA_07018934_FORF253_CAL_0005-0038.fits'
file37='F0593_FO_IMA_07018934_FORF371_CAL_0005-0038.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field W
FieldW=field(
name='FieldW',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#file names for field X
file25='F0594_FO_IMA_07018937_FORF253_CAL_0018-0047.fits'
file37='F0594_FO_IMA_07018937_FORF371_CAL_0018-0047.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field X
FieldX=field(
name='FieldX',
filename=fname,
m1cut=0.2,
m2lims=[[0,40,0,340],[0,100,300,340]],
m3lims=None)

#file names for field Y
file25='F0594_FO_IMA_07018938_FORF253_CAL_0048-0075.fits'
file37='F0594_FO_IMA_07018938_FORF371_CAL_0048-0075.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for Field Y
FieldY=field(
name='FieldY',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#other fields from earlier cycles
#File names for H north
file25='F0217_FO_IMA_70030015_FORF253_CAL_0032-0057.fits'
file37='F0217_FO_IMA_70030016_FORF371_CAL_0058-0071.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for H_north
Hnorth=field(
name='H_North',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#File Names for H South
file25='F0217_FO_IMA_70030012_FORF253_CAL_0084-0093.fits'
file37='F0217_FO_IMA_70030013_FORF371_CAL_0094-0100.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for H_south
Hsouth=field(
name='H_South',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#File Names for ArchE
file25='F0225_FO_IMA_70030026_FORF253_CAL_0027-0044.fits'
file37='F0225_FO_IMA_70030028_FORF371_CAL_0062-0077.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for ArchE
ArchE=field(
name='ArchE',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#File Names for ArchNE
file25='F0225_FO_IMA_70030023_FORF253_CAL_0078-0095.fits'
file37='F0324_FO_IMA_70040067_FORF371_CAL_0030-0108.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for ArchNE
ArchNE=field(
name='ArchNE',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#File Names for ArchNW
file25='F0224_FO_IMA_70030020_FORF253_CAL_0033-0056.fits'
file37='F0224_FO_IMA_70030021_FORF371_CAL_0058-0080.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for ArchNW
ArchNW=field(
name='ArchNW',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#File Names for ArchSE
file25='F0227_FO_IMA_70030029_FORF253_CAL_0059-0077.fits'
file37='F0227_FO_IMA_70030030_FORF371_CAL_0078-0097.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for ArchSE
#ArchSE parameters
ArchSE=field(
name='ArchSE',
filename=fname,
m1cut=0.5,
m2lims=None,
m3lims=None)

#File Names for ArchW
file25='F0225_FO_IMA_70030017_FORF371_CAL_0004-0026.fits'
file37='F0225_FO_IMA_70030017_FORF253_CAL_0004-0026.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for ArchW
ArchW=field(
name='ArchW',
filename=fname,
m1cut=0.5,
m2lims=[[250,350,150,250]],
m3lims=None)



#File Names for cycle9 Field A
file25='F0754_FO_IMA_0902161_FORF253_CAL_0176-0198.fits'
file37='F0754_FO_IMA_0902161_FORF371_CAL_0176-0198.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for A9
FieldA9=field(
name='FieldA9',
filename=fname,
m1cut=0.5,
m2lims=None, #[[250,350,150,250]],
m3lims=None)

#File Names for cycle9 Field B
file25='F0754_FO_IMA_0902162_FORF253_CAL_0127-0175.fits'
file37='F0754_FO_IMA_0902162_FORF371_CAL_0127-0175.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for B9
FieldB9=field(
name='FieldB9',
filename=fname,
m1cut=0.5,
m2lims=None, #[[250,350,150,250]],
m3lims=None)


#File Names for cycle9 Field C
file25='F0754_FO_IMA_0902163_FORF253_CAL_0077-0126.fits'
file37='F0754_FO_IMA_0902163_FORF371_CAL_0077-0126.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for C9
FieldC9=field(
name='FieldC9',
filename=fname,
m1cut=0.5,
m2lims=None, #[[250,350,150,250]],
m3lims=None)


#File Names for cycle9 Field F
file25='F0753_FO_IMA_0902166_FORF253_CAL_0192-0227.fits'
file37='F0753_FO_IMA_0902166_FORF371_CAL_0192-0227.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for F9
FieldF9=field(
name='FieldF9',
filename=fname,
m1cut=0.5,
m2lims=None, #[[250,350,150,250]],
m3lims=None)


#File Names for cycle9 Field I
file25='F0755_FO_IMA_0902169_FORF253_CAL_0133-0182.fits'
file37='F0755_FO_IMA_0902169_FORF253_CAL_0133-0182.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for I9
FieldI9=field(
name='FieldI9',
filename=fname,
m1cut=0.5,
m2lims=None, #[[250,350,150,250]],
m3lims=None)


#File Names for cycle9 Field L
file25='F0753_FO_IMA_09021612_FORF253_CAL_0163-0187.fits'
file37='F0753_FO_IMA_09021612_FORF371_CAL_0163-0187.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for L9
FieldL9=field(
name='FieldL9',
filename=fname,
m1cut=0.5,
m2lims=None, #[[250,350,150,250]],
m3lims=None)


#File Names for cycle9 Field N
file25='F0755_FO_IMA_09021614_FORF253_CAL_0002-0132.fits'
file37='F0755_FO_IMA_09021614_FORF371_CAL_0002-0132.fits'
if wavelength==25:
    fname=file25
else:
    fname=file37
#parameters defining the field object for N9
FieldN9=field(
name='FieldN9',
filename=fname,
m1cut=0.5,
m2lims=None, #[[250,350,150,250]],
m3lims=None)


