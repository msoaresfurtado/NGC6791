import glob
import time
import pandas as pd
import numpy as np
#Live dangerously
import warnings
warnings.filterwarnings("ignore")
from astropy.table import Table
from astropy.io import fits
from astropy import coordinates
from astropy import units as u
from astroquery.vizier import Vizier
from astroquery.mast import Catalogs
v = Vizier()

t = Table.read('ColmanTable_NGC6791.txt', format='ascii')

df = t.to_pandas()
df=df.drop(columns=['Quarter', 'EXPOSURE','CLUSTER'])
df=df.rename(columns={"RA_OBJ": "RA", "DEC_OBJ": "DEC", "KEPLERID":"KEPID"})



df=df.assign(GaiaID=0,KIC=0,TIC=0,\
            Vmag=0.0,Jmag=0.0,Kmag=0.0,kepmag=0.0,Imag=0.0,Rmag=0.0,Bmag=0.0,Gmag=0.0,bp_rp=0.0,\
            pmra=0.0,pmdec=0.0,pmra_error=0.0,pmdec_error=0.0,radial_velocity=0.0,
            rv_template_logg=0.0,rv_template_fe_h=0.0,\
            Rad=0.0,Mass=0.0,Lum=0.0,logg=0.0,LClass=np.nan,astrometric_excess_noise=0.0,\
            Teff=0.0,e_Dist=0,a_g_val=0,\
            parallax=0.0,parallax_error=0.0,ProbCluster=0.0,phot_variable_flag=0)

a1=0
a2=len(df)

start_time = time.time()
for i in range(a1,a2):
    coords=str(df.RA[i])+' '+str(df.DEC[i])
    c = coordinates.SkyCoord(coords,unit=('deg','deg'),frame='icrs')
    #first fill in any values found in GAIA
    gaia_cat = Catalogs.query_region(coords, radius=0.0005,catalog="Gaia", version=2)
    if len(gaia_cat)>0:
        df.GaiaID[i]=gaia_cat['source_id'][0]
        df.parallax[i]=gaia_cat['parallax'][0]
        df.parallax_error[i]=gaia_cat['parallax_error'][0]
        df.pmra[i]=gaia_cat['pmra'][0]
        df.pmra_error[i]=gaia_cat['pmra_error'][0]
        df.pmdec[i]=gaia_cat['pmdec'][0]
        df.pmdec_error[i]=gaia_cat['pmdec_error'][0]
        df.astrometric_excess_noise[i]=gaia_cat['astrometric_excess_noise'][0]
        df.Rad[i]=gaia_cat['radius_val'][0]
        df.Lum[i]=gaia_cat['lum_val'][0]
        df.Teff[i]=gaia_cat['teff_val'][0]
        df.Gmag[i]=gaia_cat['phot_g_mean_mag'][0]
        df.bp_rp[i]=gaia_cat['bp_rp'][0]
        df.radial_velocity[i]=gaia_cat['radial_velocity'][0]
        df.a_g_val[i]=gaia_cat['a_g_val'][0]
        df.rv_template_logg[i]=gaia_cat['rv_template_logg'][0]
        df.rv_template_fe_h[i]=gaia_cat['rv_template_fe_h'][0]
    #now let us search other catalogs for other magnitudes and info
    result = v.query_region(c, radius=0.0005*u.deg)
    print(len(result),'total catalogs for source',i)
    #collect all the values and choose the one with the best precision
    Vmags, Jmags, Kmags, Imags, Rmags, Bmags,Rads,Masses,loggs,Lums\
    =[],[],[],[],[],[],[],[],[],[]
    for ii in result:
        if ('KIC' in ii.columns):
            if ('-' not in str(ii['KIC'][0])):
                df.KIC[i]=ii['KIC'][0]
        if (df.TIC[i]==0) & ('TIC' in (ii.columns)):
            df.TIC[i]=ii['TIC'][0]
        if ('Vmag' in (ii.columns)):
            Vmags.append(ii['Vmag'][0])
        if ('Jmag' in (ii.columns)):
            Jmags.append(ii['Jmag'][0])
        if ('jmag' in (ii.columns)):
            Jmags.append(ii['jmag'][0])
        if ('Kmag' in (ii.columns)):
            Kmags.append(ii['Kmag'][0])
        if ('kmag' in (ii.columns)):
            Kmags.append(ii['kmag'][0])
        if ('kepmag' in (ii.columns)):
            df.kepmag[i]=str(ii['kepmag'][0])
        if ('Imag' in (ii.columns)):
            Imags.append(ii['Imag'][0])
        if ('Rmag' in (ii.columns)):
            Rmags.append(ii['Rmag'][0])
        if ('Bmag' in (ii.columns)):
            Bmags.append(ii['Bmag'][0])
        if ('Proba' in (ii.columns)):
            df.ProbCluster[i]=ii['Proba'][0]
        ###
        if ('Mass' in (ii.columns)) and (df.Mass[i]==0):
            Masses.append(ii['Mass'][0])
        if ('Lum' in (ii.columns)) and (df.Lum[i]==0):
            Lums.append(ii['Lum'][0])
        if ('logg' in (ii.columns)):
            loggs.append(ii['logg'][0])
        if ('LClass' in (ii.columns)):
            df.LClass[i]=ii['LClass'][0]
        if ('Rad' in (ii.columns)) and (df.Rad[i]==0):
            Rads.append(ii['Rad'][0])
        if ('Teff' in (ii.columns)) and (df.Teff[i]==0):
            Teffs.append(ii['Teff'][0])
    df.Vmag[i]=str(round(np.nanmedian(Vmags),2))
    df.Jmag[i]=str(round(np.nanmedian(Jmags),2))
    df.Kmag[i]=str(round(np.nanmedian(Kmags),2))
    df.Kmag[i]=str(round(np.nanmedian(Kmags),2))
    df.Imag[i]=str(round(np.nanmedian(Imags),2))
    df.Rmag[i]=str(round(np.nanmedian(Rmags),2))
    df.Bmag[i]=str(round(np.nanmedian(Bmags),2))
    if len(Masses)>1:
        df.Mass[i]=str(round(np.nanmedian(Masses),2))
    df.logg[i]=str(round(np.nanmedian(loggs),2))
    if len(Rads)>1:
        df.Rad[i]=str(round(np.nanmedian(Rad),2))
    if len(Teffs)>1:
        df.Teffs[i]=str(round(np.nanmedian(Teffs),2))
    print('| source ', i, '| ', round((time.time() - start_time),1),'seconds | ',round((time.time() - start_time)/60.0,1),'minutes | ',\
         round((time.time() - start_time)/3600.0,1),'hours')
    print('____________________________________________________________________________________________________________')

t2 = Table.from_pandas(df)
output='ColmanTable_NGC6791_FullTable.dat'
t2.write(output,format='ascii',overwrite=True)
