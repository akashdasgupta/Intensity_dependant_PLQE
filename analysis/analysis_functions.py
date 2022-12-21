import numpy as np
import os
import csv
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d as inter
import scipy.constants as sci
from copy import deepcopy
from scipy.signal import find_peaks
from scipy.integrate import simps as simpson
from scipy import stats
from scipy.optimize import curve_fit

# A nice color pallet:
cpal = ["#b30000", "#7c1158", "#4421af", "#1a53ff", "#0d88e6", "#00b7c7", "#5ad45a", "#8be04e", "#ebdc78"]

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

# Manule noise fitting for our oceanoptics spectrometr
# Measured 10/2022

def readnoise(int_time_ms):
    y0=2.07658
    A1=6.17822
    t1=7673.19756
    A2=0.24656
    t2=377.04338

    return y0 + A1*(1 - np.exp(-int_time_ms/t1)) + A2*(1 - np.exp(-int_time_ms/t2))
###


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
    all_raw_data= []
    all_noise = []

    for filename, exposure in zip(filenames, exposures):
        wavel = []
        counts = []
        raw_counts = []
        noise = []

        with open(f"{path}/{filename}",'r') as file:
            reader = csv.reader(file)
            for row in reader:
                wavel.append(float(row[0]))
                counts.append((float(row[1]) / exposure)*calibration(float(row[0]))*(1e-6) / (sci.h * sci.c / float(row[0]))) 


                raw_counts.append(float(row[1]))
                noise.append(readnoise(exposure))

        wavel = np.array(wavel)
        counts = np.array(counts)
        raw_counts = np.array(raw_counts)
        noise = np.array(noise)
        
        data = np.vstack((wavel, counts)).T
        data_raw = np.vstack((wavel, raw_counts)).T
        data_noise = np.vstack((wavel, noise)).T
        
        all_data.append(data)
        all_raw_data.append(data_raw)
        all_noise.append(data_noise)
    
    
    all_data = [i for _,i in sorted(zip(powers, all_data))]
    all_raw_data = [i for _,i in sorted(zip(powers, all_raw_data))]
    all_noise = [i for _,i in sorted(zip(powers, all_noise))]
    powers = np.array(sorted(powers))

    return powers, all_data,all_raw_data,all_noise


def load_data(datapath, feature_cutoff=450, correct_empty=True, custom_feature_location=None):
    
    # creates database to hold data:
    measurement_db = {}
    noise_db = {}
    # Hard coded: All three needed to calculate:
    measurement_types = ['empty', 'in', 'out']
    
    
    for measurement_type in measurement_types:
        # Load Data:
        short_power, short_data, short_raw,short_noise  = process_csv(f"{datapath}/{measurement_type}/short")
        long_power, long_data,long_raw,long_noise = process_csv(f"{datapath}/{measurement_type}/long")
        dark_s_power, dark_s_data, dark_s_raw, dark_s_noise= process_csv(f"{datapath}/{measurement_type}/dark_s")
        dark_l_power, dark_l_data,dark_l_raw, dark_l_noise = process_csv(f"{datapath}/{measurement_type}/dark_l")

        # S2n arrays
        short_s2r = deepcopy(short_noise)
        long_s2r =deepcopy(long_noise)
        
        # Dark correction: 
        for i in range(len(short_power)):
            short_counts = short_data[i][:,1]
            short_dark_counts = dark_s_data[i][:,1]
            short_data[i][:,1] = short_counts - short_dark_counts

            long_counts = long_data[i][:,1]
            long_dark_counts = dark_l_data[i][:,1]
            long_data[i][:,1] = long_counts - long_dark_counts


            # Short S2n:
            short_signal = short_raw[i][:,1]-dark_s_raw[i][:,1]
            short_readnoise = np.sqrt(short_noise[i][:,1]**2 + dark_s_noise[i][:,1]**2)
            short_s2r[i][:,1] = short_signal/short_readnoise

            long_signal = long_raw[i][:,1]-dark_l_raw[i][:,1]
            long_readnoise = np.sqrt(long_noise[i][:,1]**2 + dark_l_noise[i][:,1]**2)
            long_s2r[i][:,1] = long_signal/long_readnoise
        
        # Average power from 4 measurements:
        power = (short_power + long_power + dark_s_power + dark_l_power) / 4

        measurement_db[measurement_type] = [power, short_data, long_data]
        noise_db[measurement_type] = [power, short_s2r, long_s2r]
    
    if correct_empty:
        #####
        # empty correction: For removing stray peaks present because of the laser
        # Scales an empty measurement so features are same size
        #####

        normalisation_db = {} # Holds empty to measurement scaling factor (asumes stray features scale linearly wrt each other)

        # Fits peaks to find feature: 
        power, short_datas, long_datas = measurement_db['empty']
        
        if not custom_feature_location:
            peaks, properties = find_peaks(np.log(short_datas[len(power)-2][:,1]),  
                                               prominence=1, width=10)

            # selects 1st feature past specified wavelength:
            for peak in peaks:
                peak_lambda = short_datas[len(power)-2][:,0][peak]
                if peak_lambda >=feature_cutoff:
                    correction_index = peak
                    break
        else: 
            correction_index = np.argmin([abs(i-custom_feature_location) for i in short_datas[0][:,0]])
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


    return measurement_db, noise_db


def remove_noisy_data(measurement_db, noise_db, PL_range, min_s2rn=1, prnt_name=','):
    wavel_min, wavel_max = PL_range
    wavels = measurement_db['in'][1][0][:,0]

    min_index = np.argmin([abs(i-wavel_min) for i in wavels])
    max_index = np.argmin([abs(i-wavel_max) for i in wavels])
    powers, short_noisees, long_noisees = noise_db['in']

    db = {}
    for key in measurement_db.keys():
        db[key] = []
    
    cutoff = 0
    for i, (power, short_noise, long_noise) in enumerate(zip(powers, short_noisees, long_noisees)):
        s2rn = np.mean(long_noise[min_index:max_index, 1])
        if s2rn < min_s2rn:
            print (f'Sample {prnt_name}: Signal-to-read-noise of {round(s2rn,2)} at power {power} was under the min value of: {min_s2rn}. Data discarded.')
            cutoff = i
    
    for key in measurement_db.keys():
        for entry in measurement_db[key]:
            db[key].append(deepcopy(entry)[cutoff:])
    return db




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


def round_sig(x, sig=2):
    return round(x, sig-int(np.floor(np.log10(abs(x))))-1)

def ideality_fit(QFLS, num_suns, sunmin, sunmax):
    lnmin = np.log(sunmin)
    lnmax = np.log(sunmax)

    ln_curr =  np.log(num_suns)
    
    minindex = np.argmin(abs(ln_curr-lnmin))
    maxindex = np.argmin(abs(ln_curr-lnmax))

    mask = ~np.isnan(ln_curr[minindex:maxindex]) & ~np.isnan(QFLS[minindex:maxindex])
    m,c = stats.linregress(ln_curr[minindex:maxindex][mask], QFLS[minindex:maxindex][mask])[0:2]
    return sci.e*m/(sci.k*293), (m,c)
            

def JV_forfit(V,J0, n):
    # assuming Jsc of 1, from numsuns scaling
    return 1 - J0*(np.exp(sci.e*V/(n*sci.k*298)))