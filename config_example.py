#configuration file for the SOFIA/FORCAST Galactic Center photometry code

#some common parameters for all files
dpath='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/update'  #path to data
dpathalt='C:/Users/mhankins1/Documents/SofiaLegacyProgram/SOFIA_Cycle7_Data/all'    #alternate path to data

ds9path='C:/Users/mhankins1/ds9.exe' #path to envoke ds9. Use None if you don't want to open ds9 at the end

#--------------------------------------------------------------------------------------------------------

class wv():
    _registry = []
    
    # init method or constructor 
    def __init__(self, name, wavelength, bkgbox, finddetsig, segdetsig, radii, r_in,r_out, cutsize):
        self.name=name
        self.wavelength=wavelength
        self.bkgbox=bkgbox
        self.finddetsig=finddetsig
        self.segdetsig=segdetsig
        self.radii=radii
        self.r_in=r_in
        self.r_out=r_out
        self.cutsize=cutsize
        self._registry.append(self)
        
        
#wavelength specific parameters
F252=wv(
name='F252',
wavelength=25,
bkgbox=9, #size of box used for background model
finddetsig=5.0, #detection threshold sigma for daofinder
segdetsig=3.0, #detection threshold sigma for segmentation map
radii = [4, 8, 12], #aperture radii to use in photoemtry - units are pixels
r_in  = 12, #inner radius for background annulus
r_out = 20, #outer radius for background annulus
cutsize=25) #size of source cutout

#wavelength specific parameters
F371=wv(
name='F371',
wavelength=37,
bkgbox=11, #size of box used for background model
finddetsig=5.0, #detection threshold sigma for daofinder
segdetsig=3.0, #detection threshold sigma for segmentation map
radii = [5.5, 10.5, 14], #aperture radii to use in photoemtry - units are pixels
r_in  = 14, #inner radius for background annulus
r_out = 22, #outer radius for background annulus
cutsize=31) #size of source cutout



class mosaic():
    _registry = []
    
    # init method or constructor 
    def __init__(self, name, filename, m1cut,m2lims):
        self.name=name
        self.filename=filename
        self.m1cut=m1cut
        self.m2lims=m2lims
        self._registry.append(self)

Mosaic25=mosaic(
name='Mosaic25',  #name of field
filename='F0217_FO_IMA_70030015_FORF253_MOS_0001-0348_final_MATT_Corrected.fits', #file name
m1cut=0.01, # mask cutoff for normalized exposure time map
m2lims=None) #list of manually defined masked regions in the form [y1,y2,x1,x2]

Mosaic37=mosaic(
name='Mosaic37',  #name of field
filename='F0217_FO_IMA_70030016_FORF371_MOS_0001-0348_final_MATT_Corrected.fits', #file name
m1cut=0.01, # mask cutoff for normalized exposure time map
m2lims=None) #list of manually defined masked regions in the form [y1,y2,x1,x2]

class field():
    _registry = []
    
    # init method or constructor 
    def __init__(self, name, file25, file37, m1cut, m2lims):
        self.name=name
        self.file25=file25
        self.file37=file37
        self.m1cut=m1cut
        self.m2lims=m2lims
        self._registry.append(self)
        
        
#file names for new field01 - old cycle 7, field 11
Field01=field(
name='Field01',
file25='F0588_FO_IMA_07018911_FORF253_CAL_0114-0124.fits',
file37='F0588_FO_IMA_07018911_FORF371_CAL_0114-0124.fits',
m1cut=0.5,
m2lims=None)


#file names for new field02 - old cycle 7, field 10
Field02=field(
name='Field02',
file25='F0588_FO_IMA_07018910_FORF253_CAL_0082-0111.fits',
file37='F0588_FO_IMA_07018910_FORF371_CAL_0082-0111.fits',
m1cut=0.2,
m2lims=None)


#file names for new field03 - old cycle 7, field 9
Field03=field(
name='Field03',
file25='F0588_FO_IMA_0701899_FORF253_CAL_0042-0081.fits',
file37='F0588_FO_IMA_0701899_FORF371_CAL_0042-0081.fits',
m1cut=0.2,
m2lims=[[0,230,0,40],[0,50,0,50],[0,300,280,320]])


#file names for new field04 - old cycle 7 field 8
Field04=field(
name='Field04',
file25='F0591_FO_IMA_0701898_FORF253_CAL_0006-0038.fits',
file37='F0591_FO_IMA_0701898_FORF371_CAL_0006-0038.fits',
m1cut=0.5,
m2lims=None)


#file names for new field05 - old cycle 7, field 12
Field05=field(
name='Field05',
file25='F0591_FO_IMA_07018935_FORF253_CAL_0039-0061.fits',
file37='F0591_FO_IMA_07018935_FORF371_CAL_0039-0061.fits',
m1cut=0.1,
m2lims=[[0,50,0,350],[0,350,0,40],[300,350,40,80]])


#file names for new field06 - old cycle 7, field 13
Field06=field(
name='Field06',
file25='F0588_FO_IMA_07018936_FORF253_CAL_0008-0041.fits',
file37='F0588_FO_IMA_07018936_FORF371_CAL_0008-0041.fits',
m1cut=0.5,
m2lims=None)


#File Names for new field07 - old cycle 9, field A
Field07=field(
name='Field07',
file25='F0754_FO_IMA_0902161_FORF253_CAL_0176-0198.fits',
file37='F0754_FO_IMA_0902161_FORF371_CAL_0176-0198.fits',
m1cut=0.5,
m2lims=None)

#File Names for new field08  - old cycle 9, field B
Field08=field(
name='Field08',
file25='F0754_FO_IMA_0902162_FORF253_CAL_0127-0175.fits',
file37='F0754_FO_IMA_0902162_FORF371_CAL_0127-0175.fits',
m1cut=0.5,
m2lims=None)


#File Names for new field09 -  - old cycle 9, field C
Field09=field(
name='Field09',
file25='F0754_FO_IMA_0902163_FORF253_CAL_0077-0126.fits',
file37='F0754_FO_IMA_0902163_FORF371_CAL_0077-0126.fits',
m1cut=0.5,
m2lims=None)


#file names for new field10 - old cycle 7 field 7
Field10=field(
name='Field10',
file25='F0593_FO_IMA_0701897_FORF253_CAL_0125-0136.fits',
file37='F0593_FO_IMA_0701897_FORF371_CAL_0125-0136.fits',
m1cut=0.5,
m2lims=None)


#file names for new field11 - old cycle 7, field 6
Field11=field(
name='Field11',
file25='F0593_FO_IMA_0701896_FORF253_CAL_0099-0124.fits',
file37='F0593_FO_IMA_0701896_FORF371_CAL_0099-0124.fits',
m1cut=0.1,
m2lims=[[290,370,0,370],[0,370,290,380],[0,30,0,370],[0,370,0,50]])


#file names for new field 12 - old cycle 7  field U
Field12=field(
name='Field12',
file25='F0591_FO_IMA_07018932_FORF253_CAL_0062-0087.fits',
file37='F0591_FO_IMA_07018932_FORF371_CAL_0062-0087.fits',
m1cut=0.5,
m2lims=None)


#file names for new field13 - old cycle 7 field T
Field13=field(
name='Field13',
file25='F0589_FO_IMA_07018931_FORF253_CAL_0108-0137.fits',
file37='F0589_FO_IMA_07018931_FORF371_CAL_0108-0137.fits',
m1cut=0.5,
m2lims=None)


#file names for new field14 - old cycle 7 field V
Field14=field(
name='Field14',
file25='F0593_FO_IMA_07018933_FORF253_CAL_0039-0070.fits',
file37='F0593_FO_IMA_07018933_FORF371_CAL_0039-0070.fits',
m1cut=0.5,
m2lims=None)


#file names for new field15 - old cycle 7 field R
Field15=field(
name='Field15',
file25='F0589_FO_IMA_07018929_FORF253_CAL_0080-0106.fits',
file37='F0589_FO_IMA_07018929_FORF371_CAL_0080-0106.fits',
m1cut=0.2,
m2lims=None)


#file names for new field16 - old cycle 7 field W
Field16=field(
name='Field16',
file25='F0593_FO_IMA_07018934_FORF253_CAL_0005-0038.fits',
file37='F0593_FO_IMA_07018934_FORF371_CAL_0005-0038.fits',
m1cut=0.5,
m2lims=None)


#file names for new field17 - old cycle 7 field S
Field17=field(
name='Field17',
file25='F0592_FO_IMA_07018930_FORF253_CAL_0037-0067.fits',
file37='F0592_FO_IMA_07018930_FORF371_CAL_0037-0067.fits',
m1cut=0.5,
m2lims=None)


#file names for new field18 - old cycle 7, field Y
Field18=field(
name='Field18',
file25='F0594_FO_IMA_07018938_FORF253_CAL_0048-0075.fits',
file37='F0594_FO_IMA_07018938_FORF371_CAL_0048-0075.fits',
m1cut=0.5,
m2lims=None)


#File Names for new field19 - old ArchNE
Field19=field(
name='Field19',
file25='F0225_FO_IMA_70030023_FORF253_CAL_0078-0095.fits',
file37='F0324_FO_IMA_70040067_FORF371_CAL_0030-0108.fits',
m1cut=0.5,
m2lims=None)


#file names for new field 20 - old cycle 7 field Q
Field20=field(
name='Field20',
file25='F0590_FO_IMA_07018928_FORF253_CAL_0024-0055.fits',
file37='F0590_FO_IMA_07018928_FORF371_CAL_0024-0055.fits',
m1cut=0.5,
m2lims=None)


#file names for new field21 - old cycle 7 field X
Field21=field(
name='Field21',
file25='F0594_FO_IMA_07018937_FORF253_CAL_0018-0047.fits',
file37='F0594_FO_IMA_07018937_FORF371_CAL_0018-0047.fits',
m1cut=0.2,
m2lims=[[0,40,0,340],[0,100,300,340]])


#File Names for new field22 - old ArcheNW
Field22=field(
name='Field22',
file25='F0224_FO_IMA_70030020_FORF253_CAL_0033-0056.fits',
file37='F0224_FO_IMA_70030021_FORF371_CAL_0058-0080.fits',
m1cut=0.5,
m2lims=None)


#File Names for new field 23 - old ArchE
Field23=field(
name='Field23',
file25='F0225_FO_IMA_70030026_FORF253_CAL_0027-0044.fits',
file37='F0225_FO_IMA_70030028_FORF371_CAL_0062-0077.fits',
m1cut=0.5,
m2lims=None)


#file names for new field24 - old cycle 7 field P
Field24=field(
name='Field24',
file25='F0589_FO_IMA_07018927_FORF253_CAL_0138-0169.fits',
file37='F0589_FO_IMA_07018927_FORF371_CAL_0138-0169.fits',
m1cut=0.3,
m2lims=[[350,380,0,370]])


#file names for new field25 - old cycle 7 field O
Field25=field(
name='Field25',
file25='F0590_FO_IMA_07018926_FORF253_CAL_0057-0086.fits',
file37='F0590_FO_IMA_07018926_FORF371_CAL_0057-0086.fits',
m1cut=0.4,
m2lims=None)


#File Names for new field26 - old ArchSE
Field26=field(
name='Field26',
file25='F0227_FO_IMA_70030029_FORF253_CAL_0059-0077.fits',
file37='F0227_FO_IMA_70030030_FORF371_CAL_0078-0097.fits',
m1cut=0.5,
m2lims=None)


#File Names for new field27 - old ArchW
Field27=field(
name='Field27',
file25='F0225_FO_IMA_70030017_FORF253_CAL_0004-0026.fits',
file37='F0225_FO_IMA_70030017_FORF371_CAL_0004-0026.fits',
m1cut=0.5,
m2lims=[[250,350,150,250]])


#file names for new field28 - old cycle 7 field M
Field28=field(
name='Field28',
file25='F0593_FO_IMA_07018924_FORF253_CAL_0071-0098.fits',
file37='F0593_FO_IMA_07018924_FORF371_CAL_0071-0098.fits',
m1cut=0.2,
m2lims=[[230,330,0,80]])


#file names for new field29- old cycle 7 field L
Field29=field(
name='Field29',
file25='F0589_FO_IMA_07018923_FORF253_CAL_0050-0079.fits',
file37='F0589_FO_IMA_07018923_FORF371_CAL_0050-0079.fits',
m1cut=0.5,
m2lims=None)


#File names for new field30 - old H North
Field30=field(
name='Field30',
file25='F0217_FO_IMA_70030015_FORF253_CAL_0032-0057.fits',
file37='F0217_FO_IMA_70030016_FORF371_CAL_0058-0071.fits',
m1cut=0.5,
m2lims=None)


#file names for new field31 - old cycle 7 field K
Field31=field(
name='Field31',
file25='F0590_FO_IMA_07018922_FORF253_CAL_0087-0116.fits',
file37='F0590_FO_IMA_07018922_FORF371_CAL_0087-0116.fits',
m1cut=0.5,
m2lims=None)


#file names for new field32 - old cycle 7 field I
Field32=field(
name='Field32',
file25='F0592_FO_IMA_07018920_FORF253_CAL_0005-0036.fits',
file37='F0592_FO_IMA_07018920_FORF371_CAL_0005-0036.fits',
m1cut=0.3,
m2lims=None)


#file names for new field33 - old cycle 7 field G
Field33=field(
name='Field33',
file25='F0590_FO_IMA_07018918_FORF253_CAL_0117-0141.fits',
file37='F0590_FO_IMA_07018918_FORF371_CAL_0117-0141.fits',
m1cut=0.01,
m2lims=[[0,100,0,250],[0,90,250,350],[350,380,0,200],[150,230,300,380]])


#file names for new field34 - old cycle 7 field H
Field34=field(
name='Field34',
file25='F0589_FO_IMA_07018919_FORF253_CAL_0017-0049.fits',
file37='F0589_FO_IMA_07018919_FORF371_CAL_0017-0049.fits',
m1cut=0.5,
m2lims=None)


#File Names for new field35 -  old cycle 9 field F
Field35=field(
name='Field35',
file25='F0753_FO_IMA_0902166_FORF253_CAL_0192-0227.fits',
file37='F0753_FO_IMA_0902166_FORF371_CAL_0192-0227.fits',
m1cut=0.5,
m2lims=None)


#File Names for new field 36 - old H South
Field36=field(
name='Field36',
file25='F0217_FO_IMA_70030012_FORF253_CAL_0084-0093.fits',
file37='F0217_FO_IMA_70030013_FORF371_CAL_0094-0100.fits',
m1cut=0.5,
m2lims=None)


#file names for new field37 - old cycle 7 field F
Field37=field(
name='Field37',
file25='F0592_FO_IMA_07018917_FORF253_CAL_0068-0094.fits',
file37='F0592_FO_IMA_07018917_FORF371_CAL_0068-0094.fits',
m1cut=0.3,
m2lims=None)


#file names for new field38 - old cycle 7 field E
Field38=field(
name='Field38',
file25='F0590_FO_IMA_07018916_FORF253_CAL_0142-0161.fits',
file37='F0590_FO_IMA_07018916_FORF371_CAL_0142-0161.fits',
m1cut=0.5,
m2lims=None)


#file names for new field39 - old cycle 7 field B
Field39=field(
name='Field39',
file25='F0591_FO_IMA_07018913_FORF253_CAL_0088-0114.fits',
file37='F0591_FO_IMA_07018913_FORF371_CAL_0088-0114.fits',
m1cut=0.5,
m2lims=None)


#file names for new field40 - old cycle 7 field D
Field40=field(
name='Field40',
file25='F0592_FO_IMA_07018915_FORF253_CAL_0276-0303.fits',
file37='F0592_FO_IMA_07018915_FORF371_CAL_0276-0303.fits',
m1cut=0.2,
m2lims=[[200,250,0,50]])


#file names for new field41 -  old cycle 7 field C
Field41=field(
name='Field41',
file25='F0588_FO_IMA_07018914_FORF253_CAL_0326-0348.fits',
file37='F0588_FO_IMA_07018914_FORF371_CAL_0326-0348.fits',
m1cut=0.5,
m2lims=None)

#file names for new field42 - old cycle 7 field A
Field42=field(
name='Field42',
file25='F0588_FO_IMA_07018912_FORF253_CAL_0301-0325.fits',
file37='F0588_FO_IMA_07018912_FORF371_CAL_0301-0325.fits',
m1cut=0.5,
m2lims=None)


#file names for new field43 - old cycle 7 field S
Field43=field(
name='Field43',  
file25='F0595_FO_IMA_0701895_FORF253_CAL_0259-0282.fits',
file37='F0595_FO_IMA_0701895_FORF371_CAL_0259-0282.fits',
m1cut=0.5,
m2lims=[[0,330,250,300],[0,330,0,70]])


#File Names for new field44 - old cycle 9 field I
Field44=field(
name='Field44',
file25='F0755_FO_IMA_0902169_FORF253_CAL_0133-0182.fits',
file37='F0755_FO_IMA_0902169_FORF253_CAL_0133-0182.fits',
m1cut=0.5,
m2lims=None)


#file names for new field45 - old cycle 7 field 3
Field45=field(
name='Field45', 
file25='F0595_FO_IMA_0701893_FORF253_CAL_0283-0303.fits',
file37='F0595_FO_IMA_0701893_FORF371_CAL_0283-0303.fits',
m1cut=0.5,
m2lims=[[0,150,250,300]])


#file names for new field46 - old cycle 9 field J
Field46=field(
name='Field46', 
file25='F0877_FO_IMA_09021610_FORF253_CAL_0190-0246.fits',
file37='F0877_FO_IMA_09021610_FORF371_CAL_0190-0246.fits',
m1cut=0.5,
m2lims=None)


#file names for new field47 - old cycle 9 field K
Field47=field(
name='Field47', 
file25='F0877_FO_IMA_09021611_FORF253_CAL_0247-0296.fits',
file37='F0877_FO_IMA_09021611_FORF371_CAL_0247-0296.fits',
m1cut=0.5,
m2lims=None)


#File Names for new field48 - old cycle 9 field L 
Field48=field(
name='Field48',
file25='F0873_FO_IMA_09021612_FORF253_CAL_0211-0255.fits', #latest files
file37='F0873_FO_IMA_09021612_FORF371_CAL_0211-0255.fits', #latest files
m1cut=0.5,
m2lims=None)


#File Names for new field49 - old cycle 9 field N
Field49=field(
name='Field49',
file25='F0755_FO_IMA_09021614_FORF253_CAL_0002-0132.fits',
file37='F0755_FO_IMA_09021614_FORF371_CAL_0002-0132.fits',
m1cut=0.5,
m2lims=None)


#file names for new field50 - old cycle 7 field 2
Field50=field(
name='Field50', 
file25='F0592_FO_IMA_0701892_FORF253_CAL_0120-0141.fits',
file37='F0592_FO_IMA_0701892_FORF371_CAL_0120-0141.fits',
m1cut=0.5,
m2lims=None)


#file names for new field51 - old cycle 9 field O
Field51=field(
name='Field51', 
file25='F0876_FO_IMA_09021615_FORF253_CAL_0606-0638.fits',
file37='F0876_FO_IMA_09021615_FORF371_CAL_0606-0638.fits',
m1cut=0.5,
m2lims=None)


#file names for new field52 - old cycle 7 field 1
Field52=field(
name='Field52', 
file25='F0592_FO_IMA_0701891_FORF253_CAL_0095-0118.fits',
file37='F0592_FO_IMA_0701891_FORF371_CAL_0095-0118.fits',
m1cut=0.2,
m2lims=[[280,360,0,370],[0,350,0,50],[200,350,300,370]])




