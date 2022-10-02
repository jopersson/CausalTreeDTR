import numpy as np
import pandas as pd
import h5py
import csv

path = 'Data/all_hourly_data.h5'
hf = h5py.File(path, 'r')

# dynamic data frames
keys = [key for key in hf.keys()]

# codes are patient IDs
codes = hf['codes']
codes = np.array(codes['axis1_label0'])

patients = hf['patients']["table"]
patientkeys = [key for key in patients.keys()]


interventions = hf['interventions']
intervention_name = np.array(interventions['block0_items'])
interventions = np.array(interventions['block0_values']) 

vitals_labs_mean = hf['vitals_labs_mean']
vitals_labs_mean_name =  np.array(vitals_labs_mean['block0_items_level0'])
vitals_labs_mean = np.array(vitals_labs_mean['block0_values']) 

interventions = np.concatenate([codes.reshape(len(codes),1), interventions], axis=1)
vitals_labs_mean = np.concatenate([codes.reshape(len(codes),1), vitals_labs_mean], axis=1)




if __name__ == "__main__":
    print("start")
    np.savetxt("patients.csv", patients, delimiter=",", fmt = '%s')
    np.savetxt("interventions.csv", interventions, delimiter=",")
    print("interventions done")
    np.savetxt("vitals_labs_mean.csv", vitals_labs_mean, delimiter=",")
    print("vitals done")
    np.savetxt("codes.csv", codes, delimiter=",")
    print("codes done")
    np.savetxt("intervention_name.csv", intervention_name, delimiter=",", fmt='%s')
    np.savetxt("vitals_labs_mean_name.csv", vitals_labs_mean_name, delimiter=",", fmt='%s')
    print("labels done")
