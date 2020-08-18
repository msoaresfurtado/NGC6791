import glob

photometry_output_files = glob.glob("subtraction_photometry_output/" +
                                    "photometry_output_*.out")
photometry_output_files = [f for f in photometry_output_files if "BJDprepended" not in f] # to cut out the           
#previously prepended files in the same directory                                                                    
print(len(photometry_output_files), " photometry output files to process")
