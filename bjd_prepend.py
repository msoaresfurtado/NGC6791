import pyfits, glob
from astropy.table import Table
from astropy.io import fits

tpf_files= sorted(glob.glob("/Volumes/rmathieu/NGC6791/kepler/*lpd-targ.fits"))
print('total tpf files:',len(tpf_files))
print(tpf_files[0])

t = Table(names=('ID', 'RA', 'DEC'), dtype=('S18','f4','f4'))

for z in range(0,len(tpf_files)):
    hdulist = fits.open(tpf_files[z])
    hdulist_HEADER=pyfits.open(tpf_files[z])
    prihdr = hdulist[0].header
    prihdr_HEAD=hdulist_HEADER[1].data
    cols = hdulist[1].columns
    print(tpf_files[z])
    ID=tpf_files[z].split('-')[0]
    ID=ID.split('kplr00')[1]
    scidata = hdulist[1].data
    K2CADENCENO=scidata['CADENCENO']
    K2TIME=scidata['TIME'] #units of BJD - 2454833
    RA=prihdr['RA_OBJ']
    DEC=prihdr['DEC_OBJ']
    print(ID,RA,DEC)
    t.add_row((ID,RA,DEC))
    output='BJD_extracted/CADENCE_TIME_FILES.out'
    t.write(output,format='ascii',overwrite=True)
    fileout='BJD_extracted/CADENCE_TIME_'+ID+'_'+str(RA)+'_'+str(DEC)+'.out'
    print(fileout)
    f = open(fileout, 'w')
    for i in range(0,len(K2TIME)):
        stringout=str(i)+' '+str(K2TIME[i]) + ' ' +str(K2CADENCENO[i])+'\n'
        print (stringout)
        f.write(stringout)
    f.close()


    
