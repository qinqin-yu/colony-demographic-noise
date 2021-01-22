msd_avg: msd averaged across tracks for each window length
msd_ste: standard error of the mean of the msd, calculated from all tracks a given window length
ensemble_msd: values for each track, rows are for each track, columns are window lengths
ensemble_ste: standard error of the mean for each track, rows are for each track, columns are window lengths
all_L: the array of window sizes in microns 
filename='t'+zero-indexed position on 96-well plate (across rows then columns, i.e. 0=A1, 1=A2, ..., 12=B1, etc)
foldername=plate index (ex. DE1 or mreB) + replicate (a, b, c, d)