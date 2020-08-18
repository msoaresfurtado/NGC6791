import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import style
import os.path
import matplotlib.patheffects as path_effects
import seaborn as sns
from astropy.io import fits
from astropy.table import Table

#astropy statistics
from astropy.stats import median_absolute_deviation
#astrobase
import astrobase
from astrobase import periodbase, checkplot


fitsfilesin=sorted(glob.glob('/Users/melinda/Dropbox/Research_NGC6791/Colman_stitched2/*.fits'))

import time
start_time = time.time()
print(start_time)
for i in range(0,16): 
    hdulist = fits.open(fitsfilesin[i])
    hdu = hdulist[0]
    print('File ',i)
    print(hdu.header['KEPLERID'])
    outfile1='/Users/melinda/Dropbox/Research_NGC6791/Colman_PKL_files/'+str(hdu.header['KEPLERID'])+'.pkl'
    if os.path.exists(outfile1):
        print("Skipping existing file: ",outfile1)
        continue
    else:
        binarytable = Table(hdulist[1].data)
        fluxes=1+binarytable['CORRECTED FLUX']
        times=binarytable['TIME']
        errs=0.01*fluxes
        print(round((time.time() - start_time),1),'seconds | ',round((time.time() - start_time)/60.0,1),'minutes')
        ls = periodbase.pgen_lsp(times,fluxes,errs,magsarefluxes=True)
        print(round((time.time() - start_time),1),'seconds | ',round((time.time() - start_time)/60.0,1),'minutes')
        print('LS done')
        bls=periodbase.bls_parallel_pfind(times,fluxes,errs,magsarefluxes=True)
        print(round((time.time() - start_time),1),'seconds | ',round((time.time() - start_time)/60.0,1),'minutes')
        print('bls done')
        pdm = periodbase.stellingwerf_pdm(times,fluxes,errs,magsarefluxes=True)
        print(round((time.time() - start_time),1),'seconds | ',round((time.time() - start_time)/60.0,1),'minutes')
        print('pdm done')
        cpf = checkplot.checkplot_pickle([ls,pdm,bls],times, fluxes,errs,magsarefluxes=True,\
                outfile=outfile1,minbinelems=1, plotxlim=(-1,1),\
                objectinfo={'objectid': hdulist[0].header['HLSPTARG'],'ra': hdulist[0].header['RA_OBJ'],'decl': hdulist[0].header['DEC_OBJ'],'ndet': hdulist[0].header['TELAPSE']})
        print(round((time.time() - start_time),1),'seconds | ',round((time.time() - start_time)/60.0,1),'minutes')
        print('____________________________________')

print('Now we will quickly generate PNG files from the PKL files')

#Creating png files from pkl files
fitsfilesin=sorted(glob.glob('/Users/melinda/Dropbox/Research_NGC6791/Colman_PKL_files/*.pkl'))

for i in range(0,len(fitsfilesin)):
    ID=fitsfilesin[i].split(".pkl",1)[0]
    ID=ID.split("files/",1)[1]
    outfile1='/Users/melinda/Dropbox/Research_NGC6791/Colman_PKL_files/'+str(ID)+'.pkl'
    outfile2='/Users/melinda/Dropbox/Research_NGC6791/Colman_PNG_files/'+str(ID)+'.png'
    if os.path.exists(outfile2):
        print("Skipping existing file: ",outfile2)
        continue
    else:
        cpfpng = checkplot.checkplot_pickle_to_png(outfile1,outfile2)
        print(outfile1,'\n',outfile2)
        print(round((time.time() - start_time),1),'seconds | ',round((time.time() - start_time)/60.0,1),'minutes')
        print('____________________________________')

print('Now we will make the Master PNG')

df=pd.read_table('ColmanTable_NGC6791_FullTable.dat',delimiter=' ')
df_BSS=df[(df.Gmag<17)&(df.bp_rp<1.3)&df.bp_rp>0.5]
PNGfilesin=sorted(glob.glob('/Users/melinda/Dropbox/Research_NGC6791/Colman_PNG_files/*.png'))

sns.set(font_scale=2.0)
sns.set_style({'font.family': 'Times New Roman'})

for i in range(0,len(PNGfilesin)):
    ID=PNGfilesin[i].split(".png",1)[0]
    ID=ID.split("files/",1)[1]
    print(ID)
    outfile3='/Users/melinda/Dropbox/Research_NGC6791/Colman_PNG_Full_files/'+str(ID)+'.png'
    if os.path.exists(outfile3):
        print("Skipping existing file: ",outfile3)
        continue
    else:
        idx = (df['KEPID'] == float(ID)).idxmax()
        fig = plt.figure(figsize=(30,40),frameon=False)
        ax0 = plt.subplot2grid((5, 3), (0, 0), colspan=3,rowspan=3)
        ax1 = plt.subplot2grid((5, 3), (3, 0), colspan=1,rowspan=2)
        ax2 = plt.subplot2grid((5, 3), (3,1))
        ax3 = plt.subplot2grid((5, 3), (3,2))
        ax4 = plt.subplot2grid((5, 3), (4, 1))

        img1 = mpimg.imread(PNGfilesin[i])
        ax0.imshow(img1)
        ax0.axis('off')

        ax1.scatter(df.bp_rp,df.Gmag,c='grey',s=2.0,alpha=0.5,label='All LCs --- Coleman')
        ax1.scatter(df.bp_rp[(df.Gmag<17)&(df.bp_rp<1.3)&df.bp_rp>0.5],\
                    df.Gmag[(df.Gmag<17)&(df.bp_rp<1.3)&df.bp_rp>0.5],\
                    c='blue',s=40,label='BSS --- Coleman')
        if ID in df_BSS['KEPID']:
            ax1.annotate('Identified as a BSS', xy=(-0.8, 21.75), size=28,color='green')
            ax1.scatter(df.bp_rp[df.KEPID==float(ID)],df.Gmag[df.KEPID==float(ID)],marker='*',color='green',s=400,label=str(ID))
        else:
            ax1.annotate('Not identified as a BSS', xy=(-0.8, 21.75), size=28,color='red')
            ax1.scatter(df.bp_rp[df.KEPID==float(ID)],df.Gmag[df.KEPID==float(ID)],marker='*',color='red',s=400,label=str(ID))
        ax1.set_xlim(-1,3)
        ax1.set_ylim(22,12)
        ax1.set_ylabel('Gmag')
        ax1.set_xlabel('B-R')
        ax1.legend()

        ax2.scatter(df.pmra,df.pmdec,c='grey',s=2.0,alpha=0.5)
        dist=round(1/(df.parallax[idx]*1.0e-3),1)
        tmp='d = '+str(dist)+' pc'
        if ID in df_BSS['KEPID']:
            ax2.scatter(df.pmra[df.KEPID==float(ID)],df.pmdec[df.KEPID==float(ID)],marker='*',color='green',s=400)
            tmp='d ='+str(dist)+'pc'
            ax2.annotate(tmp, xy=(-1.75, -4), size=28,color='green')
        else:
            ax2.scatter(df.pmra[df.KEPID==float(ID)],df.pmdec[df.KEPID==float(ID)],marker='*',color='red',s=400)
            ax2.annotate(tmp, xy=(-1.75, -4), size=28,color='red')
        ax2.set_ylabel('pmdec')
        ax2.set_xlabel('pmra')

        ax3.scatter(df.Mass,df.astrometric_excess_noise,c='grey',s=2.0,alpha=0.5)
        mass=df.Mass[idx]
        aen=round(df.astrometric_excess_noise[idx],4)
        tmp='Mass = '+str(mass)+' Msun'
        tmp2='AEN = '+str(aen)
        if ID in df_BSS['KEPID']:
            ax3.scatter(df.Mass[df.KEPID==float(ID)],df.astrometric_excess_noise[df.KEPID==float(ID)],\
                        marker='*',color='green',s=400)
            ax3.annotate(tmp, xy=(0.75, 1.4), size=28,color='green')
            ax3.annotate(tmp2, xy=(0.75, 1.2), size=28,color='green')
        else:
            ax3.scatter(df.Mass[df.KEPID==float(ID)],df.astrometric_excess_noise[df.KEPID==float(ID)],\
                    marker='*',color='red',s=400)
            ax3.annotate(tmp, xy=(0.75, 1.4), size=28,color='red')
            ax3.annotate(tmp2, xy=(0.75, 1.2), size=28,color='red')
        ax3.set_ylabel('astrometric_excess_noise')
        ax3.set_xlabel('Mass')

        ax4.scatter(df.RA,df.DEC,c='grey',s=2.0,alpha=0.5)
        if ID in df_BSS['KEPID']:
            ax4.scatter(df.RA[df.KEPID==float(ID)],df.DEC[df.KEPID==float(ID)],marker='*',color='green',s=400)
        else:
            ax4.scatter(df.RA[df.KEPID==float(ID)],df.DEC[df.KEPID==float(ID)],marker='*',color='red',s=400)
        ax4.set_ylabel('DEC')
        ax4.set_xlabel('RA')


        fig.savefig(outfile3, dpi=150)
        print(outfile3)
        print(round((time.time() - start_time),1),'seconds | ',round((time.time() - start_time)/60.0,1),'minutes')
        print('____________________________________')

print('Everything is done!')
