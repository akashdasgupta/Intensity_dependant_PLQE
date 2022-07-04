import numpy as np
import os
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d as inter
import scipy.constants as sci
from copy import deepcopy
from scipy.signal import find_peaks
from scipy.integrate import simpson
from scipy import stats

### Loads calibration (hard coded for now): 
cal_wavel, cal_values = ([],[])

with open(f".\calibration_data\cal_QEpro_steel_200um_2022-05-05.txt",'r') as file:
    reader = csv.reader(file, delimiter = ' ')
    for row in reader:
        cal_wavel.append(float(row[0]))
        cal_values.append(float(row[1]))

calibration = inter(cal_wavel,cal_values, fill_value='extrapolate')
del cal_wavel
del cal_values


def process_csv(path):
    powers = []
    exposures = []
    for _,_,filenames in os.walk(path):
        for filename in filenames:
            power, exposure = filename.split('.csv')[0].split('_')
            powers.append(float(power))
            exposures.append(float(exposure))
            
        break

    all_data = []
    for filename, exposure in zip(filenames, exposures):
        wavel = []
        counts = []
        with open(f"{path}/{filename}",'r') as file:
            reader = csv.reader(file)
            for row in reader:
                wavel.append(float(row[0]))
                counts.append((float(row[1]) / exposure)*calibration(float(row[0]))*(1e-6) / (sci.h * sci.c / float(row[0]))) 
        wavel = np.array(wavel)
        counts = np.array(counts)
        
        data = np.vstack((wavel, counts)).T
        all_data.append(data)
#     print(exposures)
    return np.array(powers), all_data


def load_data(datapath, feature_cutoff=450, correct_empty=True):
    
    # creates database to hold data:
    measurement_db = {}
    # Hard coded: All three needed to calculate:
    measurement_types = ['empty', 'in', 'out']
    
    
    for measurement_type in measurement_types:
        # Sorted short measurement:
        short_power, short_data = process_csv(f"{datapath}/{measurement_type}/short")
        short_data = [i for _,i in sorted(zip(short_power, short_data))]
        short_power = np.array(sorted(short_power))

        # sorted long measurement:
        long_power, long_data = process_csv(f"{datapath}/{measurement_type}/long")
        long_data = [i for _,i in sorted(zip(long_power, long_data))]
        long_power = np.array(sorted(long_power))

        # sorted short reference:
        dark_s_power, dark_s_data = process_csv(f"{datapath}/{measurement_type}/dark_s")
        dark_s_data = [i for _,i in sorted(zip(dark_s_power, dark_s_data))]
        dark_s_power = np.array(sorted(dark_s_power))

        # sorted long reference: 
        dark_l_power, dark_l_data = process_csv(f"{datapath}/{measurement_type}/dark_l")
        dark_l_data = [i for _,i in sorted(zip(dark_l_power, dark_l_data))]
        dark_l_power = np.array(sorted(dark_l_power))

        
        # Dark correction: 
        for i in range(len(short_power)):
            short_counts = short_data[i][:,1]
            short_dark_counts = dark_s_data[i][:,1]

            short_data[i][:,1] = short_counts - short_dark_counts

            long_counts = long_data[i][:,1]
            long_dark_counts = dark_l_data[i][:,1]

            long_data[i][:,1] = long_counts - long_dark_counts
        
        # Average power from 4 measurements:
        power = (short_power + long_power + dark_s_power + dark_l_power) / 4

        measurement_db[measurement_type] = [power, short_data, long_data]
    
    if correct_empty:
        #####
        # empty correction: For removing stray peaks present because of the laser
        # Scales an empty measurement so features are same size
        #####

        normalisation_db = {} # Holds empty to measurement scaling factor (asumes stray features scale linearly wrt each other)

        # Fits peaks to find feature: 
        power, short_datas, long_datas = measurement_db['empty']
        peaks, properties = find_peaks(np.log(short_datas[len(power)-2][:,1]),  
                                           prominence=1, width=10)

        # selects 1st feature past specified wavelength:
        for peak in peaks:
            peak_lambda = short_datas[len(power)-2][:,0][peak]
            if peak_lambda >=feature_cutoff:
                correction_index = peak
                break
        # Stores scaling factor: 
        for measurement in measurement_types:
            power, short_datas, long_datas = measurement_db[measurement]

            short_lowest_peak = np.mean(short_datas[len(power)-2][correction_index-1:correction_index+1 ,1])
            long_lowest_peak = np.mean(long_datas[len(power)-2][correction_index-1:correction_index+1 ,1])

            normalisation_db[measurement] = (short_lowest_peak, long_lowest_peak)

        # IN substraction: 

        measurement_db['in_empty_subtracted'] = deepcopy(measurement_db['in'])
        # short:
        for i in range(len( measurement_db['in_empty_subtracted'][0])):

            empty_correction_factor = normalisation_db['in'][0] / normalisation_db['empty'][0]
            empty_scaled = measurement_db['empty'][1][i][:,1] * empty_correction_factor 
            measurement_db['in_empty_subtracted'][1][i][:,1] -= empty_scaled
        # long   
        for i in range(len( measurement_db['in_empty_subtracted'][0])):

            empty_correction_factor = normalisation_db['in'][1] / normalisation_db['empty'][1]
            empty_scaled = measurement_db['empty'][2][i][:,1] * empty_correction_factor 
            measurement_db['in_empty_subtracted'][2][i][:,1] -= empty_scaled


        # OUT substraction:
        measurement_db['out_empty_subtracted'] = deepcopy(measurement_db['out'])

        # short
        for i in range(len( measurement_db['out_empty_subtracted'][0])):

            empty_correction_factor = normalisation_db['out'][0] / normalisation_db['empty'][0]
            empty_scaled = measurement_db['empty'][1][i][:,1] * empty_correction_factor 
            measurement_db['out_empty_subtracted'][1][i][:,1] -= empty_scaled
        # long
        for i in range(len( measurement_db['out_empty_subtracted'][0])):

            empty_correction_factor = normalisation_db['out'][1] / normalisation_db['empty'][1]
            empty_scaled = measurement_db['empty'][2][i][:,1] * empty_correction_factor 
            measurement_db['out_empty_subtracted'][2][i][:,1] -= empty_scaled


    return measurement_db



def calculate_PLQE(db, laser_range, PL_range):

    PLQE = []
    PLQE_no_out = []
    
    A = []
    LA =[]
    L_in = []
    PL_in = []
    PL_out = []
    
    for i in range(len(db['in'][0])):
        laser_min, laser_max = laser_range
        PL_min, PL_max = PL_range

        A.append(1 - simpson(db['in'][1][i][:,1][np.argmin(abs(db['in'][1][i][:,0]-laser_min)): 
                                                 np.argmin(abs(db['in'][1][i][:,0]-laser_max))], 
                            db['in'][1][i][:,0][np.argmin(abs(db['in'][1][i][:,0]-laser_min)): 
                                                np.argmin(abs(db['in'][1][i][:,0]-laser_max))])/simpson(db['out'][1][i][:,1][np.argmin(abs(db['out'][1][i][:,0]-laser_min)): 
                                                                                                                             np.argmin(abs(db['out'][1][i][:,0]-laser_max))], 
                                                                                                     db['out'][1][i][:,0][np.argmin(abs(db['out'][1][i][:,0]-laser_min)): 
                                                                                                                          np.argmin(abs(db['out'][1][i][:,0]-laser_max))]))
        LA.append(simpson(db['empty'][1][i][:,1][np.argmin(abs(db['empty'][1][i][:,0]-laser_min)): 
                                                 np.argmin(abs(db['empty'][1][i][:,0]-laser_max))], 
                      db['empty'][1][i][:,0][np.argmin(abs(db['empty'][1][i][:,0]-laser_min)): 
                                             np.argmin(abs(db['empty'][1][i][:,0]-laser_max))]))
        
        L_in.append(simpson(db['in'][1][i][:,1][np.argmin(abs(db['in'][1][i][:,0]-laser_min)): 
                                                np.argmin(abs(db['in'][1][i][:,0]-laser_max))], 
                      db['in'][1][i][:,0][np.argmin(abs(db['in'][1][i][:,0]-laser_min)): 
                                          np.argmin(abs(db['in'][1][i][:,0]-laser_max))]))
        
        PL_in.append(simpson(db['in_empty_subtracted'][2][i][:,1][np.argmin(abs(db['in_empty_subtracted'][2][i][:,0]-PL_min)): 
                                                                  np.argmin(abs(db['in_empty_subtracted'][2][i][:,0]-PL_max))], 
                         db['in_empty_subtracted'][2][i][:,0][np.argmin(abs(db['in_empty_subtracted'][2][i][:,0]-PL_min)): 
                                                              np.argmin(abs(db['in_empty_subtracted'][2][i][:,0]-PL_max))]))

        
        PL_out.append(simpson(db['out_empty_subtracted'][2][i][:,1][np.argmin(abs(db['out_empty_subtracted'][2][i][:,0]-PL_min)): 
                                                                    np.argmin(abs(db['out_empty_subtracted'][2][i][:,0]-PL_max))], 
                         db['out_empty_subtracted'][2][i][:,0][np.argmin(abs(db['out_empty_subtracted'][2][i][:,0]-PL_min)): 
                                                               np.argmin(abs(db['out_empty_subtracted'][2][i][:,0]-PL_max))]))
        
    A = np.array(A)
    LA =np.array(LA)
    L_in = np.array(L_in)
    PL_in = np.array(PL_in)
    PL_out = np.array(PL_out)
    
    
    PLQE = (PL_in - (1-A)*PL_out)/(LA*A) 
    PLQE_no_out = PL_in/LA

    return db['in'][0], np.array(PLQE),  np.array(PLQE_no_out)



# Voc, radiative limit:

bandgap_v = []
voc_rad = []
with open(r"./calibration_data/voc_rad_eg.csv", "r") as file:
    reader = csv.reader(file)
    for row in reader:
        bandgap_v.append(float(row[0]))
        voc_rad.append(float(row[1]))
vocradf = inter(bandgap_v, voc_rad)
del bandgap_v
del voc_rad

def interpolate_gaps(values, limit=None):
    """
    Fill gaps using linear interpolation, optionally only fill gaps up to a
    size of `limit`.
    """
    values = np.asarray(values)
    i = np.arange(values.size)
    valid = np.isfinite(values)
    filled = np.interp(i, i[valid], values[valid])

    if limit is not None:
        invalid = ~valid
        for n in range(1, limit+1):
            invalid[:-n] &= invalid[n:]
        filled[invalid] = np.nan

    return filled



def ideality_fit(QFLS, num_suns, sunmin, sunmax):
    lnmin = np.log(sunmin)
    lnmax = np.log(sunmax)

    ln_curr =  np.log(num_suns)
    
    minindex = np.argmin(abs(ln_curr-lnmin))
    maxindex = np.argmin(abs(ln_curr-lnmax))

    m,c = stats.linregress(ln_curr[minindex:maxindex], QFLS[minindex:maxindex])[0:2]
    return sci.e*m/(sci.k*293), (m,c)
            