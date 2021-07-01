#configuration file for the photometry code

#some common parameters for all files
dpath='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/all'  #path to data
dpathalt='E:\Documents\SofiaLegacyProgram\CAL_files\CAL_C'   #alternate path to data

ds9path='C:\\Users\\mhankins1\\ds9.exe' #path to envoke ds9. Use None if you don't want to open ds9 at the end

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
file25='F0592_FO_IMA_0701891_FORF253_CAL_0095-0118_FIELD_1.fits'
file37='F0592_FO_IMA_0701891_FORF371_CAL_0095-0118_FIELD_1.fits'
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
file25='F0592_FO_IMA_0701892_FORF253_CAL_0120-0141_FIELD_2.fits'
file37='F0592_FO_IMA_0701892_FORF371_CAL_0120-0141_FIELD_2.fits'
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
file25='F0595_FO_IMA_0701893_FORF253_CAL_0283-0303_FIELD_3.fits'
file37='F0595_FO_IMA_0701893_FORF371_CAL_0283-0303_FIELD_3.fits'
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
file25='F0595_FO_IMA_0701895_FORF253_CAL_0259-0282_FIELD_5.fits'
file37='F0595_FO_IMA_0701895_FORF371_CAL_0259-0282_FIELD_5.fits'
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
file25='F0593_FO_IMA_0701896_FORF253_CAL_0099-0124_Field_6.fits'
file37='F0593_FO_IMA_0701896_FORF371_CAL_0099-0124_Field_6.fits'
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
file25='F0593_FO_IMA_0701897_FORF253_CAL_0125-0136_Field_7.fits'
file37='F0593_FO_IMA_0701897_FORF371_CAL_0125-0136_Field_7.fits'
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
file25='F0591_FO_IMA_0701898_FORF253_CAL_0006-0038_Field_8.fits'
file37='F0591_FO_IMA_0701898_FORF371_CAL_0006-0038_Field_8.fits'
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
file25='F0588_FO_IMA_0701899_FORF253_CAL_0042-0081_Field_9.fits'
file37='F0588_FO_IMA_0701899_FORF371_CAL_0042-0081_Field_9.fits'
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
file25='F0588_FO_IMA_07018910_FORF253_CAL_0082-0111_Field_10.fits'
file37='F0588_FO_IMA_07018910_FORF371_CAL_0082-0111_Field_10.fits'
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
file25='F0588_FO_IMA_07018911_FORF253_CAL_0114-0124_Field_11.fits'
file37='F0588_FO_IMA_07018911_FORF371_CAL_0114-0124_Field_11.fits'
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
file25='F0591_FO_IMA_07018935_FORF253_CAL_0039-0061_Field_12.fits'
file37='F0591_FO_IMA_07018935_FORF371_CAL_0039-0061_Field_12.fits'
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
file25='F0588_FO_IMA_07018936_FORF253_CAL_0008-0041_Field_13.fits'
file37='F0588_FO_IMA_07018936_FORF371_CAL_0008-0041_Field_13.fits'
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
file25='F0588_FO_IMA_07018912_FORF253_CAL_0301-0325_Field_A.fits'
file37='F0588_FO_IMA_07018912_FORF371_CAL_0301-0325_Field_A.fits'
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
file25='F0591_FO_IMA_07018913_FORF253_CAL_0088-0114_Field_B.fits'
file37='F0591_FO_IMA_07018913_FORF371_CAL_0088-0114_Field_B.fits'
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
file25='F0588_FO_IMA_07018914_FORF253_CAL_0326-0348_Field_C.fits'
file37='F0588_FO_IMA_07018914_FORF371_CAL_0326-0348_Field_C.fits'
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
file25='F0592_FO_IMA_07018915_FORF253_CAL_0276-0303_Field_D.fits'
file37='F0592_FO_IMA_07018915_FORF371_CAL_0276-0303_Field_D.fits'
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
file25='F0590_FO_IMA_07018916_FORF253_CAL_0142-0161_Field_E.fits'
file37='F0590_FO_IMA_07018916_FORF371_CAL_0142-0161_Field_E.fits'
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
file25='F0592_FO_IMA_07018917_FORF253_CAL_0068-0094_Field_F.fits'
file37='F0592_FO_IMA_07018917_FORF371_CAL_0068-0094_Field_F.fits'
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
file25='F0590_FO_IMA_07018918_FORF253_CAL_0117-0141_Field_G.fits'
file37='F0590_FO_IMA_07018918_FORF371_CAL_0117-0141_Field_G.fits'
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
file25='F0589_FO_IMA_07018919_FORF253_CAL_0017-0049_Field_H.fits'
file37='F0589_FO_IMA_07018919_FORF371_CAL_0017-0049_Field_H.fits'
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
file25='F0592_FO_IMA_07018920_FORF253_CAL_0005-0036_Field_I.fits'
file37='F0592_FO_IMA_07018920_FORF371_CAL_0005-0036_Field_I.fits'
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
file25='F0590_FO_IMA_07018922_FORF253_CAL_0087-0116_Field_K.fits'
file37='F0590_FO_IMA_07018922_FORF371_CAL_0087-0116_Field_K.fits'
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
file25='F0589_FO_IMA_07018923_FORF253_CAL_0050-0079_Field_L.fits'
file37='F0589_FO_IMA_07018923_FORF371_CAL_0050-0079_Field_L.fits'
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
file25='F0593_FO_IMA_07018924_FORF253_CAL_0071-0098_Field_M.fits'
file37='F0593_FO_IMA_07018924_FORF371_CAL_0071-0098_Field_M.fits'
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
file25='F0590_FO_IMA_07018926_FORF253_CAL_0057-0086_Field_O.fits'
file37='F0590_FO_IMA_07018926_FORF371_CAL_0057-0086_Field_O.fits'
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
file25='F0589_FO_IMA_07018927_FORF253_CAL_0138-0169_Field_P.fits'
file37='F0589_FO_IMA_07018927_FORF371_CAL_0138-0169_Field_P.fits'
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
file25='F0590_FO_IMA_07018928_FORF253_CAL_0024-0055_Field_Q.fits'
file37='F0590_FO_IMA_07018928_FORF371_CAL_0024-0055_Field_Q.fits'
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
file25='F0589_FO_IMA_07018929_FORF253_CAL_0080-0106_Field_R.fits'
file37='F0589_FO_IMA_07018929_FORF371_CAL_0080-0106_Field_R.fits'
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
file25='F0592_FO_IMA_07018930_FORF253_CAL_0037-0067_Field_S.fits'
file37='F0592_FO_IMA_07018930_FORF371_CAL_0037-0067_Field_S.fits'
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
file25='F0589_FO_IMA_07018931_FORF253_CAL_0108-0137_Field_T.fits'
file37='F0589_FO_IMA_07018931_FORF371_CAL_0108-0137_Field_T.fits'
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
file25='F0591_FO_IMA_07018932_FORF253_CAL_0062-0087_Field_U.fits'
file37='F0591_FO_IMA_07018932_FORF371_CAL_0062-0087_Field_U.fits'
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
file25='F0593_FO_IMA_07018933_FORF253_CAL_0039-0070_Field_V.fits'
file37='F0593_FO_IMA_07018933_FORF371_CAL_0039-0070_Field_V.fits'
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
file25='F0593_FO_IMA_07018934_FORF253_CAL_0005-0038_Field_W.fits'
file37='F0593_FO_IMA_07018934_FORF371_CAL_0005-0038_Field_W.fits'
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
file25='F0594_FO_IMA_07018937_FORF253_CAL_0018-0047_Field_X.fits'
file37='F0594_FO_IMA_07018937_FORF371_CAL_0018-0047_Field_X.fits'
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
file25='F0594_FO_IMA_07018938_FORF253_CAL_0048-0075_Field_Y.fits'
file37='F0594_FO_IMA_07018938_FORF371_CAL_0048-0075_Field_Y.fits'
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
file25='F0217_FO_IMA_70030015_FORF253_CAL_0032-0057_Hnorth.fits'
file37='F0217_FO_IMA_70030016_FORF371_CAL_0058-0071_Hnorth.fits'
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
file25='F0217_FO_IMA_70030012_FORF253_CAL_0084-0093_Hsouth.fits'
file37='F0217_FO_IMA_70030013_FORF371_CAL_0094-0100_Hsouth.fits'
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
file25='F0225_FO_IMA_70030026_FORF253_CAL_0027-0044_ArchE.fits'
file37='F0225_FO_IMA_70030028_FORF371_CAL_0062-0077_ArchE.fits'
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
file25='F0225_FO_IMA_70030023_FORF253_CAL_0078-0095_ArchNE.fits'
file37='F0324_FO_IMA_70040067_FORF371_CAL_0030-0108_ArchNE.fits'
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
file25='F0224_FO_IMA_70030020_FORF253_CAL_0033-0056_ArchNW.fits'
file37='F0224_FO_IMA_70030021_FORF371_CAL_0058-0080_ArchNW.fits'
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
file25='F0227_FO_IMA_70030029_FORF253_CAL_0059-0077_ArchSE.fits'
file37='F0227_FO_IMA_70030030_FORF371_CAL_0078-0097_ArchSE.fits'
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
file25='F0225_FO_IMA_70030017_FORF371_CAL_0004-0026_ArchW.fits'
file37='F0225_FO_IMA_70030017_FORF253_CAL_0004-0026_ArchW.fits'
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




