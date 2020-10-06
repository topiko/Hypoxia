# Scripts to calculate the TCP (cell survival fraction)
# after given dose has been given:
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from input_utils import res_dir, K_names, get_cmd_args_parser, terminated_flag, dose_data
from run_all import active_vessels_arr, K0s, rids
from plot_utils import plot_D99, figure_folder_name, figw
OERm = 3

def OMF(K, alpha=.01, beta=.001, Krm=1):

    return OER(K, alpha, beta, Krm)/OERm

def OER(K, alpha=.01, beta=.001, Krm=1):
    # [Powathil 2012] Jakob has different!
    return (OERm*K + Krm)/(K + Krm)

def S(K, dose, alpha=.01, beta=.001, Krm=1):
    # Cell survival farction after dose:
    return np.exp(-alpha*OMF(K, alpha, beta, Krm)*dose - beta*(OMF(K, alpha, beta, Krm)*dose)**2)

def integrate_S(K_arr, dose, alpha, beta, Krm):

    total_area = K_arr[:, 1].sum()
    dS_frac = np.zeros(len(K_arr))
    for j in range(len(K_arr)):
        dS_frac[j] = S(K_arr[j, 0], dose, alpha=alpha, beta=beta, Krm=Krm)*K_arr[j, 1]/total_area

    print('Dose: {} --> Surv frac: {}'.format(dose, dS_frac.sum()))
    return dS_frac.sum()

def get_D99_dose(data, alpha, beta, Krm):
    # Get dose [Gy] required for 99 cell destruction:

    K_arr= np.zeros((len(data)//3, 2))
    for i in range(len(data)//3):
        cell_data = data[3*i:3*(i+1)]
        K_arr[i,0] = cell_data[:,4].mean()
        K_arr[i,1] = cell_data[0,5]

    #total_area = K_arr[:, 1].sum()
    #dS_frac = np.zeros(len(K_arr))

    def get_S_diff(dose):

        intS = integrate_S(K_arr, dose, alpha, beta, Krm)
        #for j in range(len(K_arr)):
        #    dS_frac[j] = S(K_arr[j, 0], dose, alpha=alpha, beta=beta, Km=Km)*K_arr[j, 1]/total_area

        #print('Dose: {} --> Surv frac: {}'.format(dose, dS_frac.sum()))
        return (intS-0.01)**2

    dose = fmin(get_S_diff, 10, ftol=.0001, xtol=.01)[0] #, args=(), x

    return dose

if __name__ == '__main__':


    if len(sys.argv) == 1:
        doses = np.linspace(0, 200, 100)
        K = 20
        plt.plot(doses, S(K, doses))
        plt.xlabel('Dose [Gy]')
        plt.ylabel('Survival fraction')
        plt.semilogy()
        plt.show()
    else:

        parser = get_cmd_args_parser()
        args = parser.parse_args()

        params = vars(args)
        res_dir_f = res_dir.format(**params)
        Krm = params['K_rm']
        alpha = params['alpha']  # Powathil 2012 /see ../Refs/Powathil2012
        beta = params['beta'] #.03 # Powathil 2012
        oereq = params['oereq']

        for K0 in K0s:

            params['K_0'] = K0
            path = dose_data.format(**params) #, alpha=alpha, beta=beta, oereq='powathil')
            try:
                np.load(path)
                print('Found data in: {}'.format(path))
            except FileNotFoundError:

                doses_arr = np.zeros((len(active_vessels_arr), len(rids)))
                for i, active_vessels in enumerate(active_vessels_arr):
                    params['active_vessels'] = active_vessels

                    doses = np.zeros(len(rids))

                    rids_run = [0] if active_vessels==100 else rids
                    for j, rid in enumerate(rids_run):
                        params['run_idx'] = rid
                        print(params)
                        res_dir_f = res_dir.format(**params)
                        data = np.load(res_dir_f + K_names(terminated_flag)[1])

                        doses[j:] = get_D99_dose(data, alpha, beta, Krm)
                    doses_arr[i, :] = doses

                np.save(path, doses_arr)


        plot_D99(figure_folder_name.format(**params), params, figw, figh=3)

