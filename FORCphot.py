



#Now add RA, DEC coordinates of sources to table
def addSkyCentroid(tab,wcsmap):
    if 'sky_centroid' in tab.columns:
        print('this table already has sky centroid column!')
    else:
        #Now add RA, DEC coordinates of sources to table
        if len(tab)>0:

            Nsources=len(tab)
            scs=[]

            for i in range(0,Nsources):
                xcoord=tab['xcentroid'][i]
                ycoord=tab['ycentroid'][i]
                sc=pixel_to_skycoord(xcoord,ycoord,wcsmap)
                scs.append(sc)

            tab['sky_centroid']=scs
        
    return tab



# A few useful functions for creating ds9 files
def makeDS9reg(tab,radius,outname,color=None):
    #get source coordinates from table
    sourcecoords=tab['sky_centroid']

    #set size of regions 
    radius = Angle(radius, u.deg) 

    #loop through and create region instances for each source
    regions=[]
    for i in range(0,len(sourcecoords)):
        region = CircleSkyRegion(sourcecoords[i], radius)
        regions.append(region)
        
    #write out region file
    write_ds9(regions, outname)
    
    if color is not None:
        #change the color of the regions to red - no built in way to do this in regions package :-/
        with open(outname, 'r+') as f:
            text = f.read()
            text = re.sub(r'\)', r') # color='+str(color), text)
            f.seek(0)
            f.write(text)
            f.truncate()
    return text


def makeCombDS9file(text1,text2,text3,outname):
    newregtext=text1+text2[45:]+text3[45:]

    with open(outname, 'w+') as f:
        f.seek(0)
        f.write(newregtext)
        f.truncate()













