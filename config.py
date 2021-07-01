#listing of the defined files
def __dir__():
    #return ['Field1', 'Field2', 'Field3', 'Field5', 'Field6', 'Field7', 'Field8', 'Field9', 'Field10', 'Field11', 'Field12', 'Field13', ]
    return ['Mosaic25','Mosaic37','Field1', 'Field10', 'Field11', 'Field12', 'Field13', 'Field2', 'Field3', 'Field5', 'Field6', 'Field7', 'Field8', 'Field9', 'FieldA', 'FieldArchE', 'FieldArchNE', 'FieldArchNW', 'FieldArchSE', 'FieldArchW', 'FieldB', 'FieldC', 'FieldD', 'FieldE', 'FieldF', 'FieldG', 'FieldH', 'FieldHnorth', 'FieldHsouth', 'FieldI', 'FieldK', 'FieldL', 'FieldM', 'FieldO', 'FieldP', 'FieldQ', 'FieldR', 'FieldS', 'FieldT', 'FieldU', 'FieldV', 'FieldW', 'FieldX', 'FieldY']

#some common parameters for all files
dpath='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/all'  #path to data
dpathalt='E:\Documents\SofiaLegacyProgram\CAL_files\CAL_C'   #alternate path to data
detectsigma=3.0  #detection threshold sigma
bkgbox=13 #size of box used for background model
ds9path='C:\\Users\\mhankins1\\ds9.exe' #path to envoke ds9. Use None if you don't want to open ds9 at the end

class field():
    # init method or constructor 
    def __init__(self, name, filename,m1cut,m2lims,m3lims):
        self.name=name
        self.filename=filename
        self.m1cut=m1cut
        self.m2lims=m2lims
        self.m3lims=m3lims
        

Mosaic25=field(
name='Mosaic25',  #name of field
filename='F0217_FO_IMA_70030015_FORF253_MOS_0001-0348_final_MATT_Corrected.fits', #file name
m1cut=0.01, # mask cutoff for normalized exposure time map
m2lims=None, #list of manually defined masked regions in the form [y1,y2,x1,x2]
m3lims=None) #additional mask if needed in the form [y1,y2,x1,x2]

Mosaic37=field(
name='Mosaic37',  #name of field
filename='F0217_FO_IMA_70030016_FORF371_MOS_0001-0348_final_MATT_Corrected.fits', #file name
m1cut=0.01, # mask cutoff for normalized exposure time map
m2lims=None, #list of manually defined masked regions in the form [y1,y2,x1,x2]
m3lims=None) #additional mask if needed in the form [y1,y2,x1,x2]

#field1 parameters
Field1=field(
name='Field1',  #name of field
filename='F0592_FO_IMA_0701891_FORF253_CAL_0095-0118_FIELD_1.fits', #file name
m1cut=0.2, # mask cutoff for normalized exposure time map
m2lims=[[280,360,0,370],[0,350,0,50],[200,350,300,370]], #list of manually defined masked regions in the form [y1,y2,x1,x2]
m3lims=None) #additional mask if needed in the form [y1,y2,x1,x2]

#field2 parameters
Field2=field(
name='Field2',
filename='F0592_FO_IMA_0701892_FORF253_CAL_0120-0141_FIELD_2.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#field3 parameters
Field3=field(
name='Field3',
filename='F0595_FO_IMA_0701893_FORF253_CAL_0283-0303_FIELD_3.fits',
m1cut=0.5,
m2lims=[[0,150,250,300]],
m3lims=None)

#field5 parameters
Field5=field(
name='Field5',
filename='F0595_FO_IMA_0701895_FORF253_CAL_0259-0282_FIELD_5.fits',
m1cut=0.5,
m2lims=[[0,330,250,300]],
m3lims=[[0,330,0,70]])

#field6 parameters
Field6=field(
name='Field6',
filename='F0593_FO_IMA_0701896_FORF253_CAL_0099-0124_Field_6.fits',
m1cut=0.1,
m2lims=[[290,370,0,370],[0,370,290,380],[0,30,0,370],[0,370,0,50]],
m3lims=None)

#field7 parameters
Field7=field(
name='Field7',
filename='F0593_FO_IMA_0701897_FORF253_CAL_0125-0136_Field_7.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#field8 parameters
Field8=field(
name='Field8',
filename='F0591_FO_IMA_0701898_FORF253_CAL_0006-0038_Field_8.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#field9 parameters
Field9=field(
name='Field9',
filename='F0588_FO_IMA_0701899_FORF253_CAL_0042-0081_Field_9.fits',
m1cut=0.2,
m2lims=[[0,230,0,40],[0,50,0,50],[0,300,280,320]],
m3lims=None)


#field10 parameters
Field10=field(
name='Field10',
filename='F0588_FO_IMA_07018910_FORF253_CAL_0082-0111_Field_10.fits',
m1cut=0.2,
m2lims=None,
m3lims=None)


#field11 parameters
Field11=field(
name='Field11',
filename='F0588_FO_IMA_07018911_FORF253_CAL_0114-0124_Field_11.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#field12 parameters
Field12=field(
name='Field12',
filename='F0591_FO_IMA_07018935_FORF253_CAL_0039-0061_Field_12.fits',
m1cut=0.1,
m2lims=[[0,50,0,350],[0,350,0,40],[300,350,40,80]],
m3lims=None)

#field13 parameters
Field13=field(
name='Field13',
filename='F0588_FO_IMA_07018936_FORF253_CAL_0008-0041_Field_13.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldA parameters
FieldA=field(
name='FieldA',
filename='F0588_FO_IMA_07018912_FORF253_CAL_0301-0325_Field_A.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldB parameters
FieldB=field(
name='FieldB',
filename='F0591_FO_IMA_07018913_FORF253_CAL_0088-0114_Field_B.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldC parameters
FieldC=field(
name='FieldC',
filename='F0588_FO_IMA_07018914_FORF253_CAL_0326-0348_Field_C.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldD parameters
FieldD=field(
name='FieldD',
filename='F0592_FO_IMA_07018915_FORF253_CAL_0276-0303_Field_D.fits',
m1cut=0.2,
m2lims=[[200,250,0,50]],
m3lims=None)

#fieldE parameters
FieldE=field(
name='FieldE',
filename='F0590_FO_IMA_07018916_FORF253_CAL_0142-0161_Field_E.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldF parameters
FieldF=field(
name='FieldF',
filename='F0592_FO_IMA_07018917_FORF253_CAL_0068-0094_Field_F.fits',
m1cut=0.3,
m2lims=None,
m3lims=None)

#fieldG parameters
FieldG=field(
name='FieldG',
filename='F0590_FO_IMA_07018918_FORF253_CAL_0117-0141_Field_G.fits',
m1cut=0.01,
m2lims=[[0,100,0,250],[0,90,250,350],[350,380,0,200],[150,230,300,380]],
m3lims=None)

#fieldH parameters
FieldH=field(
name='FieldH',
filename='F0589_FO_IMA_07018919_FORF253_CAL_0017-0049_Field_H.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldI parameters
FieldI=field(
name='FieldI',
filename='F0592_FO_IMA_07018920_FORF253_CAL_0005-0036_Field_I.fits',
m1cut=0.3,
m2lims=None,
m3lims=None)

#fieldK parameters
FieldK=field(
name='FieldK',
filename='F0590_FO_IMA_07018922_FORF253_CAL_0087-0116_Field_K.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldL parameters
FieldL=field(
name='FieldL',
filename='F0589_FO_IMA_07018923_FORF253_CAL_0050-0079_Field_L.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldM parameters
FieldM=field(
name='FieldM',
filename='F0593_FO_IMA_07018924_FORF253_CAL_0071-0098_Field_M.fits',
m1cut=0.2,
m2lims=[[230,330,0,80]],
m3lims=None)

#fieldO parameters
FieldO=field(
name='FieldO',
filename='F0590_FO_IMA_07018926_FORF253_CAL_0057-0086_Field_O.fits',
m1cut=0.4,
m2lims=None,
m3lims=None)

#fieldP parameters
FieldP=field(
name='FieldP',
filename='F0589_FO_IMA_07018927_FORF253_CAL_0138-0169_Field_P.fits',
m1cut=0.3,
m2lims=[[350,380,0,370]],
m3lims=None)

#fieldQ parameters
FieldQ=field(
name='FieldQ',
filename='F0590_FO_IMA_07018928_FORF253_CAL_0024-0055_Field_Q.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldR parameters
FieldR=field(
name='FieldR',
filename='F0589_FO_IMA_07018929_FORF253_CAL_0080-0106_Field_R.fits',
m1cut=0.2,
m2lims=None,
m3lims=None)

#fieldS parameters
FieldS=field(
name='FieldS',
filename='F0592_FO_IMA_07018930_FORF253_CAL_0037-0067_Field_S.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldT parameters
FieldT=field(
name='FieldT',
filename='F0589_FO_IMA_07018931_FORF253_CAL_0108-0137_Field_T.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldT parameters
FieldT=field(
name='FieldT',
filename='F0589_FO_IMA_07018931_FORF253_CAL_0108-0137_Field_T.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldU parameters
FieldU=field(
name='FieldU',
filename='F0591_FO_IMA_07018932_FORF253_CAL_0062-0087_Field_U.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldV parameters
FieldV=field(
name='FieldV',
filename='F0593_FO_IMA_07018933_FORF253_CAL_0039-0070_Field_V.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldW parameters
FieldW=field(
name='FieldW',
filename='F0593_FO_IMA_07018934_FORF253_CAL_0005-0038_Field_W.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldW parameters
FieldW=field(
name='FieldW',
filename='F0593_FO_IMA_07018934_FORF253_CAL_0005-0038_Field_W.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#fieldX parameters
FieldX=field(
name='FieldX',
filename='F0594_FO_IMA_07018937_FORF253_CAL_0018-0047_Field_X.fits',
m1cut=0.2,
m2lims=[[0,40,0,340],[0,100,300,340]],
m3lims=None)

#fieldY parameters
FieldY=field(
name='FieldY',
filename='F0594_FO_IMA_07018938_FORF253_CAL_0048-0075_Field_Y.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#other fields from earlier cycles
#H North parameters
FieldHnorth=field(
name='H_North',
filename='F0217_FO_IMA_70030015_FORF253_CAL_0032-0057_Hnorth.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#H South parameters
FieldHsouth=field(
name='H_South',
filename='F0217_FO_IMA_70030012_FORF253_CAL_0084-0093_Hsouth.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#ArchE parameters
FieldArchE=field(
name='ArchE',
filename='F0225_FO_IMA_70030026_FORF253_CAL_0027-0044_ArchE.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#ArchNE parameters
FieldArchNE=field(
name='ArchNE',
filename='F0225_FO_IMA_70030023_FORF253_CAL_0078-0095_ArchNE.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#ArchNW parameters
FieldArchNW=field(
name='ArchNW',
filename='F0224_FO_IMA_70030020_FORF253_CAL_0033-0056_ArchNW.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#ArchSE parameters
FieldArchSE=field(
name='ArchSE',
filename='F0227_FO_IMA_70030029_FORF253_CAL_0059-0077_ArchSE.fits',
m1cut=0.5,
m2lims=None,
m3lims=None)

#ArchW parameters
FieldArchW=field(
name='ArchW',
filename='F0225_FO_IMA_70030017_FORF253_CAL_0004-0026_ArchW.fits',
m1cut=0.5,
m2lims=[[250,350,150,250]],
m3lims=None)