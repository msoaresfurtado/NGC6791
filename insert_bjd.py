# This file will insert the cadence number and BJD into all the output files
# in subtraction_photometry_output/ so that, after running grcollect,
# we will have an accurate measure of the time of each observation for use
# for making the light curves.  grcollect isn't sorting things in proper
# time order for whatever reason, and we're going to need the BJD eventually
# anyway.

import numpy as np
import glob as glob

## First, get the BJD's for all the target pixel files
BJD_files = glob.glob("BJD_extracted/BJD_*.txt")

BJD_dicts_dict = {}
for i in range(len(BJD_files)):
    print(BJD_files[i])
    key = int(BJD_files[i].split("/")[-1].split("_")[1].split(".")[0])
    cadence_no, zeroindex_cno, BJD = np.genfromtxt(BJD_files[i],
                                                   unpack=True)
    cadence_no = [int(val) for val in cadence_no]
    zeroindex_cno = [int(val) for val in zeroindex_cno]
    BJD_dict_justthisfile = {}
    for j in range(len(zeroindex_cno)):
        BJD_dict_justthisfile[zeroindex_cno[j]] = BJD[j]
    BJD_dicts_dict[key] = BJD_dict_justthisfile


## Now, get the ranges of each of the target pixel files to use
ranges = np.genfromtxt("BJD_extracted/xy_ranges_by_TPF.txt",
                       dtype=int)

# Assemble a list to be used as a quick match for stars to the respective 
# TPF number
quick_TPF_match = [[None, None, None], [None, None, None], [None, None, None],
                   [None, None, None], [None, None, None], [None, None, None]]

for i in range(len(ranges)):
    quick_TPF_match[ranges[i][1]/50][ranges[i][3]/50] = ranges[i][0]

###########################################
##################
#  All right, now to get the subtraction photometry files, open then, and
#  modify
#
#####################################

photometry_output_files = glob.glob("subtraction_photometry_output/" + 
                                    "photometry_output_*.out")
photometry_output_files = [f for f in photometry_output_files if "BJDprepended" not in f] # to cut out the 
#previously prepended files in the same directory

print(len(photometry_output_files), " photometry output files to process")

for p_file in photometry_output_files:
    one_indexed_cad_num = p_file.split(".")[0].split("_")[-1]
    zero_index = int(one_indexed_cad_num) - 1
    with open(p_file,"r") as f:
        lines = f.readlines()
    #print lines

    # modify the comment line
    if lines[3][0] != '#' or 'ID' not in lines[3] or 'X Coord.' not in lines[3]:
        raise RuntimeError("Incorrect comment line found!!!")

    lines[3] = lines[3][:22] + "BJD          cad. no. " + lines[3][22:]

    # Now to actually edit the stuff
    for j in range(4,len(lines)):
        first_space = lines[j].index(' ')
        splits = lines[j].split()
        BJD = BJD_dicts_dict[quick_TPF_match[int(float(splits[1]))/50][int(float(splits[2]))/50]][zero_index]
        #print BJD
        lines[j] = lines[j][:first_space] + "  " + str(BJD).ljust(13) + "  " + one_indexed_cad_num.rjust(4) + " " +\
            lines[j][first_space:]

    output_file = "subtraction_photometry_output/photometry_output_BJDprepended_" +\
        str(one_indexed_cad_num) + ".out"

    with open(output_file,"w") as f:
        for line in lines:
            f.write(line)
print('Huzzah, you have incorporated the BJD dates! Make sure these look sane!')
