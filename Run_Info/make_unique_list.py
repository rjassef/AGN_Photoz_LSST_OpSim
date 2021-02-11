#!/usr/bin/env python

import numpy as np
import run_infos as ri

#Create the families object.
families = ri.FamilyInfo()

#Read all the families info.
families.read_summary_csv(csv_file='all_summaries_2021_02_09.csv')

#Now, run through all the families getting the OpSim names.
OpSims = list()
for f in families.list_of_families():
    for findex in families.family_info(f).index:
        OpSims.append(findex)

#Filter out the repeated names.
OpSims = np.unique(OpSims)

#Write the results.
np.savetxt("Current_list_of_OpSims.txt", OpSims, '%s')
