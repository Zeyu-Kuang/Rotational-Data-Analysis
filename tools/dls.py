
# coding: utf-8

# In[ ]:

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from progressbar import ProgressBar

def add_var(value1, sigma1, value2, sigma2):
    """
    This function adds two uncorrelated varibles with their values and standard error, sigmas
    It returns the total value and sigma
    """
    sigma_tot = 1 / (np.sqrt(1 / sigma1**2 + 1 / sigma2**2))
    value_tot = (value1 / sigma1**2 + value2 / sigma2**2) * sigma_tot**2
    return value_tot, sigma_tot

def sub_var(value1, sigma1, value2, sigma2):
    """
    This function subtracts two uncorrelated varibles with their values and standard error, sigmas
    'variable 1 - variable 2'.
    It returns the total value and sigma
    """
    return add_var(value1, sigma1, -value2, sigma2)

def quick_read(runNumber):
    """
    Read data without any specification other than the run number. Return a DataFrame object contains all the data.
    runNumber could be an integer or a list.
    If runNumber is an integer, it specifies the number of one run to read. 
    If runNumber is a list, it specifies the numbers of runs to read and concatenate them together.
    For instance, if you want to extract data from run 11: quick_read(11). Extract data from run 8 and 9: quick_read([8,9]).
    """
    
    # if runNumber is an integer, read data directly
    if isinstance(runNumber, int):
        # read parameters
        if runNumber == 8:
            run_location = './Data/File_location_run_8.txt' 
            lst = 4055.8764  # local sidereal time for run 8
            rf = 360  # rotational frequency for run 8
        elif runNumber == 9:
            run_location = './Data/File_location_run_9.txt'
            lst = 1131905 + 4055.8764  # local sidereal time for run 9
            rf = 360  # rotational frequency for run 9
        elif runNumber == 10:
            run_location = './Data/File_location_run_10.txt' 
            lst = 5903354 + 4055.8764  # local sidereal time for run 10
            rf = 420  # rotational frequency for run 10
        elif runNumber == 11:
            run_location = './Data/File_location_run_11.txt' 
            lst = 7986836 + 4055.8764  # local sidereal time for run 11
            rf = 320  # rotational frequency for run 11
        else:
            raise ValueError('Run number could be 8, 9, 10, or 11')

        data_df = read_data(run_location, lst, rf)  # read data stored in the file
        data_df = data_df.reset_index(drop='True')  # reset index

        return data_df
    
    # if runNumber is a list, break it into integer. 
    if isinstance(runNumber, list):   
        # reading one by one 
        data_df = pd.DataFrame()  # storing all data
        for each in runNumber:
            data_df = pd.concat([data_df,quick_read(each)])
        
        data_df = data_df.reset_index(drop='True')  # reset index
        
        return data_df
            
    
    

def read_data(location_file, local_time, rotation_frequency):
    """
    This function reads the raw data file, whose pathname are stored in the location_file.
    Return: DataFrame obejct that contains all the raw data.
    local_time (s) should be added to the time tag so that it starts from the local sidereal frame. 
    """
    # read locations of data file
    with open(location_file) as f:
        files = f.readlines()
    files = [x.strip() for x in files]  # remove the '\n' in each line.
    
    # read and concatenate each data file
    data_df = pd.DataFrame()  # for fast reading and storing data
    for location in files:
        data = pd.read_csv(location, sep='\t', header=None, index_col=False)
        data_df = pd.concat([data_df, data])  # concatenate each data
    
    # change time to local sidereal time
    data_df[1] = data_df[1] + local_time  
    
    # change voltage to fractional frequency [unitless]
    factor_v2f = 6.65726 * 1e-8  # calibration factor when rotation frequency is 360
    data_df[3] = data_df[3] * factor_v2f * rotation_frequency / 360
    
    return data_df

def demodulate(subset):
    """
    Given a subset of signal [index, timetag, angle, delta_V], fitting it by equation (stated in OLS section). 
    Return parameters include [theta, standard deviation, averaged time]
    """
    # initialization
    N = subset[:,0]
    t = subset[:,1]
    phi_R = subset[:,2]
    delta_V = subset[:,3]
    
    # creating matrix X
    X = np.zeros([subset.shape[0],5])
    X[:,0] = 1
    X[:,1] = np.sin(phi_R)
    X[:,2] = np.cos(phi_R)
    X[:,3] = np.sin(2*phi_R)
    X[:,4] = np.cos(2*phi_R)
    
    # perform OLS
    varcov = np.linalg.inv(np.dot(X.T, X))  # inverted inner product of X
    beta = np.dot(varcov, np.dot(X.T,delta_V))  # fitting parameters
    
    res = np.inner(beta.T, X) - delta_V  # residual
    dres = np.std(res)  # standard deviation of the residual
    
    V_beta = dres * dres * varcov  # variance matrix of fitting parameters
    sigma_beta = np.sqrt(np.diag(V_beta))  # standard deviation for fitting parameters
    
    # calculate the mean time for all these parameters
    t_mean = np.mean(t)  # average time for the subset
    
    # calculate the number of samples inside the subset
    N_sample = np.size(N)
    
    # packing all parameters in one vector
    params = pd.DataFrame([[beta, sigma_beta, t_mean, N_sample]],columns=['beta','sigma','time','Nsamples'])
    
    return params

def select_param(params_df, param):
    """
    Given a params dataframe, select a certain param. Output the result as 2D narray
    """
    param_df = params_df[param]  # choose a parameter
    param_na = param_df.values  # change it to narray
    param_na = np.array(param_na.tolist())  # change it to 2D narray (tricky part)
    
    return param_na
    
def DLS1(data_na, subset_rot):
    """
    Perform the first stage of DLS. 
    data_na: raw data (numpy array).
    subset_rot: number of rotation per subset has.
    return: parameters output from DLS1 as DataFrame object
    """
    # initialization of looping parameters
    [data_rows, data_cols] = np.shape(data_na)  
    #subset_rot = 2500  # number of subset_rot to demodulate
    rotCounter = 0  # initialize the rotation counter
    index_start = 0  # initialize the starting index for slicing a subset

    # initialize variables for storing demodulated data
    params_df = pd.DataFrame()
    pbar = ProgressBar()  # showing progress

    # looping through row of the data 
    for i in pbar(range(1,data_rows)):
        delta_angle = data_na[i][2] - data_na[i-1][2]
        if delta_angle < -4: 
            # there is a jump of the angle, this means a rotation 
            rotCounter += 1
        if rotCounter == subset_rot:
            # calculate the subset till rotCounter reaches subset_rot
            index_end = i  # the ending index for slicing a subset ('index_end' is the first candidate of the second subset)
            subset = data_na[index_start : index_end]  # creating a subset for the demodulation (slicing from 'index_start' to 'index_end - 1')
            rotCounter = 0  # reset the counter to zero
            index_start = i  # reset the starting index 

            # demodulate the subset
            params = demodulate(subset)

            # storing demodulated parameters
            params_df = pd.concat([params_df, params])

    params_df = params_df.reset_index(drop='True')  # reset index
    return params_df

def DLS2(params_df, w_s, par=0):
    """
    Perform 2nd stage of DLS. Demodulate over frequency w_s (2*w_s).
    It recieves the output of DLS1 and return SME coefficients
    A list is returned. Its first element is the value of SME coefficients as a dictionary. Its second element is also a dictionary, stores the standard deviation of the SME coeffieicnet.
    
    Note: if par==1, then the parameters of the DLS2 will be returned as well
    """
    ## initalization
    A_cols = 4
    A_rows = params_df.shape[0]
    A = np.zeros((A_rows, A_cols))
    A[:,0] = params_df.index.values
    A[:,1] = select_param(params_df, 'time')

    A[:,2] = A[:,1] * w_s

    A_sin2w = A.copy()
    A_sin2w[:,3] = select_param(params_df, 'beta')[:,3]  # coefficeint of sin(2w) [unitless]

    A_cos2w = A.copy()
    A_cos2w[:,3] = select_param(params_df, 'beta')[:,4]  # coefficient of cos(2w) [unitless]

    ## Calculating

    para_S = demodulate(A_sin2w)
    para_C = demodulate(A_cos2w)

    ## Storing result

    DLS2_co_values = {}

    DLS2_co_values['S_c_0'] = para_S['beta'][0][0]
    DLS2_co_values['S_s_ws'] = para_S['beta'][0][1]
    DLS2_co_values['S_c_ws'] = para_S['beta'][0][2]
    DLS2_co_values['S_s_2ws'] = para_S['beta'][0][3]
    DLS2_co_values['S_c_2ws'] = para_S['beta'][0][4]

    DLS2_co_values['C_c_0'] = para_C['beta'][0][0]
    DLS2_co_values['C_s_ws'] = para_C['beta'][0][1]
    DLS2_co_values['C_c_ws'] = para_C['beta'][0][2]
    DLS2_co_values['C_s_2ws'] = para_C['beta'][0][3]
    DLS2_co_values['C_c_2ws'] = para_C['beta'][0][4]

    DLS2_co_sigma = {}

    DLS2_co_sigma['S_c_0'] = para_S['sigma'][0][0]
    DLS2_co_sigma['S_s_ws'] = para_S['sigma'][0][1]
    DLS2_co_sigma['S_c_ws'] = para_S['sigma'][0][2]
    DLS2_co_sigma['S_s_2ws'] = para_S['sigma'][0][3]
    DLS2_co_sigma['S_c_2ws'] = para_S['sigma'][0][4]

    DLS2_co_sigma['C_c_0'] = para_C['sigma'][0][0]
    DLS2_co_sigma['C_s_ws'] = para_C['sigma'][0][1]
    DLS2_co_sigma['C_c_ws'] = para_C['sigma'][0][2]
    DLS2_co_sigma['C_s_2ws'] = para_C['sigma'][0][3]
    DLS2_co_sigma['C_c_2ws'] = para_C['sigma'][0][4]

    ## Calculate EMS coefficients

    X = 58.014 * 2 * np.pi / 360  # colatitude of the laborotory (radius)
    sin_X = np.sin(X)
    cos_X = np.cos(X)

    EMS_co_values = {}
    EMS_co_sigma = {}

    EMS_co_values['c_T_Q'] = - DLS2_co_values['C_c_0'] / (4 * sin_X**2)
    EMS_co_sigma['c_T_Q'] = - DLS2_co_sigma['C_c_0'] / (4 * sin_X**2)

    EMS_co_values['c_T_Y_1'] = - DLS2_co_values['C_c_ws'] / (8 * sin_X)
    EMS_co_sigma['c_T_Y_1'] = - DLS2_co_sigma['C_c_ws'] / (8 * sin_X)
    EMS_co_values['c_T_Y_2'] = DLS2_co_values['S_s_ws'] / (8 * sin_X * cos_X)
    EMS_co_sigma['c_T_Y_2'] = DLS2_co_sigma['S_s_ws'] / (8 * sin_X * cos_X)
    EMS_co_values['c_T_Y'], EMS_co_sigma['c_T_Y'] = add_var(EMS_co_values['c_T_Y_1'], EMS_co_sigma['c_T_Y_1'],
                                                            EMS_co_values['c_T_Y_2'], EMS_co_sigma['c_T_Y_2'])

    EMS_co_values['c_T_X_1'] = - DLS2_co_values['C_s_ws'] / (8 * sin_X * cos_X)
    EMS_co_sigma['c_T_X_1'] = - DLS2_co_sigma['C_s_ws'] / (8 * sin_X * cos_X)
    EMS_co_values['c_T_X_2'] = - DLS2_co_values['S_c_ws'] / (8 * sin_X)
    EMS_co_sigma['c_T_X_2'] = - DLS2_co_sigma['S_c_ws'] / (8 * sin_X)
    EMS_co_values['c_T_X'], EMS_co_sigma['c_T_X'] = add_var(EMS_co_values['c_T_X_1'], EMS_co_sigma['c_T_X_1'],
                                                            EMS_co_values['c_T_X_2'], EMS_co_sigma['c_T_X_2'])

    EMS_co_values['c_T_Z_1'] = DLS2_co_values['C_s_2ws'] / (2 * (1 + cos_X)**2)
    EMS_co_sigma['c_T_Z_1'] = DLS2_co_sigma['C_s_2ws'] / (2 * (1 + cos_X)**2)
    EMS_co_values['c_T_Z_2'] = DLS2_co_values['S_c_2ws'] / (2 * (1 + cos_X)**2)
    EMS_co_sigma['c_T_Z_2'] = DLS2_co_sigma['S_c_2ws'] / (2 * (1 + cos_X)**2)
    EMS_co_values['c_T_Z'], EMS_co_sigma['c_T_Z'] = add_var(EMS_co_values['c_T_Z_1'], EMS_co_sigma['c_T_Z_1'],
                                                            EMS_co_values['c_T_Z_2'], EMS_co_sigma['c_T_Z_2'])

    EMS_co_values['c_T_M'], EMS_co_sigma['c_T_M'] = add_var(DLS2_co_values['C_c_2ws'], DLS2_co_sigma['C_c_2ws'], 
                                                            DLS2_co_values['S_s_2ws'], DLS2_co_sigma['S_s_2ws']) / (4 * (cos_X - 1)**2)
    EMS_co_values['c_T_u'], EMS_co_sigma['c_T_u'] = sub_var(DLS2_co_values['C_c_2ws'], DLS2_co_sigma['C_c_2ws'], 
                                                            DLS2_co_values['S_s_2ws'], DLS2_co_sigma['S_s_2ws']) / (4 * (cos_X + 1)**2)

    # change the EMS sigma to absolute 
    for key, values in EMS_co_sigma.items():
        EMS_co_sigma[key] = abs(values)
    
    # return parameters as well if par is set to 1
    if par == 1:
        return EMS_co_values, EMS_co_sigma, DLS2_co_values, DLS2_co_sigma
    
    return EMS_co_values, EMS_co_sigma