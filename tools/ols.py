import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


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


def OLS(dataset, update=0):
    """
    Given a dataset of signal which should be a numpy array with the form [index, sidereal_time, angle, delta_f], fitting it by equation (stated in OLS section). 
    Return parameters include [theta, standard deviation, averaged time]
    """
    # initialization
    N = dataset[:,0]
    t = dataset[:,1]
    phi_R = dataset[:,2]
    delta_f = dataset[:,3]
    
    # pre-calculation
    w_s = 2 * np.pi / (23*3600 + 56*60 + 4) # sidereal frequency
    phi_S = w_s * t  # sidereal phase
    
    # creating matrix X
    X = np.zeros([dataset.shape[0],11])
    X[:,0] = 1
    X[:,1] = np.sin(2*phi_R)
    X[:,2] = np.cos(2*phi_R)
    
    X[:,3] = np.sin(2*phi_R + phi_S)
    X[:,4] = np.cos(2*phi_R + phi_S)
    
    X[:,5] = np.sin(2*phi_R - phi_S)
    X[:,6] = np.cos(2*phi_R - phi_S)
    
    X[:,7] = np.sin(2*phi_R + 2*phi_S)
    X[:,8] = np.cos(2*phi_R + 2*phi_S)
    
    X[:,9] = np.sin(2*phi_R - 2*phi_S)
    X[:,10] = np.cos(2*phi_R - 2*phi_S)
    
    
    # perform OLS
    varcov = np.linalg.inv(np.dot(X.T, X))  # inverted inner product of X
    beta = np.dot(varcov, np.dot(X.T,delta_f))  # fitting parameters
    
    res = np.inner(beta.T, X) - delta_f  # residual
    dres = np.std(res)  # standard deviation of the residual
    
    V_beta = dres * dres * varcov  # variance matrix of fitting parameters
    sigma_beta = np.sqrt(np.diag(V_beta))  # standard deviation for fitting parameters
    
    # calculate the mean time for all these parameters
    t_mean = np.mean(t)  # average time for the dataset
    
    # calculate the number of samples inside the dataset
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

def ols_slicing(data_df,N_slice):
    """
    input data_df contains the raw data.
    N_slice is the number of subdata you want to have.
    
    output params that contains parameters for each subdata.
    """
    
    # slice 
    N_data = data_df.shape[0]
    N_piece = int(N_data / N_slice)
    params = pd.DataFrame()  # for storing params in each sub data

    for i in range(N_slice):
        # slicing each sub_data
        start = i * N_piece
        stop = (i+1) * N_piece
        next_stop = (i+2) * N_piece
        if next_stop > N_data:
            sub_data_df = data_df[start:]
        else:
            sub_data_df = data_df[start:stop]

        # do OLS on each sub_data
        sub_data = sub_data_df.values  # change to numpy array
        sub_param = OLS(sub_data)
        params = pd.concat((params, sub_param))

    params.reset_index(drop='True')
    return params
   
class var():
    """
    This class is for calculating the addition between variables with standard deviation, sigma.
    """
    def __init__(self, value, sigma):
        self.value = value
        self.sigma = sigma
        
    def __add__(self, other):
        if isinstance(other, var):
            sum_sigma = 1 / (np.sqrt(1 / self.sigma**2 + 1 / other.sigma**2))
            sum_value = (self.value / self.sigma**2 + other.value / other.sigma**2) * sum_sigma**2
            return var(sum_value, sum_sigma)
        else:
            print('Operation not supported!')
            return None
        
    def __sub__(self, other):
        sub_other = var(-other.value, other.sigma)
        
        return self.__add__(sub_other)
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return var(self.value * other, self.sigma * abs(other))
        else:
            print('Operation not supported!')
            return None
        
    def __turediv__(self, other):
        return __mul__(self, 1 / other)
    
    __rmul__ = __mul__
            
def add_params(params):
    """
    This function recieves the output parameter dataframe from function slice.
    It returns a var object that contains all summed up parameters.
    """
    # initilization
    beta = select_param(params, 'beta')
    sigma = select_param(params, 'sigma')
    
    # create var object
    N_params = 11
    param_sum = [None] * N_params
    for i in range(N_params):
        # for wach parameter
        beta_vi = beta[:, i]
        sigma_vi = sigma[:, i]
        # sum different pieces of one parameter up (with sigma)
        param_sum[i] = var(beta_vi[0], sigma_vi[0])
        for j in range(1, beta_vi.shape[0]):
            param_sum[i] = param_sum[i] + var(beta_vi[j], sigma_vi[j])
            
    return param_sum

def cal_co(param_sum):
    """
    This function calculates the EMS coefficient from the param_sum variable.
    """
    X = 58.014 * 2 * np.pi / 360  # colatitude of the laborotory (radius)
    sin_X = np.sin(X)
    cos_X = np.cos(X)

    # calculate the coefficients of EMS
    c_T_Q =  1 / (-4 * sin_X**2) * param_sum[2]
    c_T_X_1 = 1 / (-4 * (1 + cos_X) * sin_X) * param_sum[3]
    c_T_X_2 = 1 / (4 * (cos_X - 1) * sin_X) * param_sum[5]
    c_T_Y_1 = 1 / (-4 * (1 + cos_X) * sin_X) * param_sum[4]
    c_T_Y_2 = 1 / (4 * (cos_X - 1) * sin_X) * param_sum[6]
    c_T_Z = 1 / (2 * (1 + cos_X)**2) * param_sum[7]
    c_T_M = 1 / (2 * (cos_X - 1)**2) * param_sum[10]
    c_T_u = 1 / (2 * (1 + cos_X)**2) * param_sum[8]

    c_T_X = c_T_X_1 + c_T_X_2
    c_T_Y = c_T_Y_1 + c_T_Y_2
    
    return [c_T_Q, c_T_X, c_T_Y, c_T_Z, c_T_M, c_T_u]
