import glob
from astropy.io import fits
import numpy as np



files = sorted(glob.glob("/Volumes/rmathieu/NGC6791/kepler/*lpd-targ.fits"))
                         
x_positions = [] # These lists are to later figure out relative positions of frames
y_positions = []
ras = []
decs = []
names = []

for i in range(len(files)):
    # First, get appropriate part of the name
    name = files[i].split('-')[0]
    name = name[-9:]
    names.append(name)
    print(name) # The unique name/number associated with this particular stamp
    f = fits.open(files[i]) #open fits file
    #print f.info()
    header = f[1].header
    x_positions.append(header['1CRV4P'])
    y_positions.append(header['2CRV4P'])
    ras.append(header['RA_OBJ'])
    decs.append(header['DEC_OBJ'])
    table = f[1].data #get appropriate table
    f.close()

    BJD = [] #array to catch BJD
    cadence = []
    for j in range(len(table)): #loop and get BJD's
        BJD.append(table[j][0])
        cadence.append(table[j][2])

    # Write these numbers to a file
    print(BJD[5])
    with open("BJD_extracted/BJD_" + name + ".txt","w") as fwrite:
        fwrite.write("# Cadence_number zero_index_cadnum    BJD\n")
        for i in range(len(BJD)):
            fwrite.write(str(cadence[i]) + "             " +\
                             str(i) + "                   " +\
                             str(BJD[i])+"\n")


print(x_positions)
print(y_positions)

# Plotting just to check that that the lowest x_positions is also the lowest RA.
#   It is.
#import matplotlib.pyplot as plt
plt.scatter(x_positions,y_positions)
for i in range(len(ras)):
    plt.text(x_positions[i],y_positions[i],ras[i],horizontalalignment='left')
    plt.text(x_positions[i],y_positions[i],names[i],horizontalalignment='left')
plt.savefig("temp.pdf")

x_subtract = 512
y_subtract = 20

with open("BJD_extracted/xy_ranges_by_TPF.txt","w") as fwrite:
    fwrite.write("#name   xmin     xmax     ymin     ymax\n")
    for i in range(len(names)):
        xmin = x_positions[i] - x_subtract
        ymin = y_positions[i] - y_subtract
        fwrite.write(names[i] + "%8d %8d %8d %8d\n" % (xmin, xmin+50, ymin, ymin+50))
        #fwrite.write(names[i] + "      " + str(xmin) + "      " +\
        #             str(xmin+50) + "     " + str(ymin) + "    " +\
        #             str(ymin+50) + "\n")
