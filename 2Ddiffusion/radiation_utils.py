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

def OMF(K, alpha=.01, beta=.001, Km=1):

    return OER(K, alpha, beta, Km)/OERm

def OER(K, alpha=.01, beta=.001, Km=1):
    # [Powathil 2012] Jakob has different!
    return (OERm*K + Km)/(K + Km)

def S(K, dose, alpha=.01, beta=.001, Km=1):
    # Cell survival farction after dose:
    return np.exp(-alpha*OMF(K, alpha, beta, Km)*dose - beta*(OMF(K, alpha, beta, Km)*dose)**2)


def get_D99_dose(data, alpha, beta, Km):
    # Get dose [Gy] required for 99 cell destruction:

    K_arr= np.zeros((len(data)//3, 2))
    for i in range(len(data)//3):
        cell_data = data[3*i:3*(i+1)]
        K_arr[i,0] = cell_data[:,4].mean()
        K_arr[i,1] = cell_data[0,5]

    total_area = K_arr[:, 1].sum()

    dS_frac = np.zeros(len(K_arr))

    def get_S_diff(dose):
        for j in range(len(K_arr)):
            dS_frac[j] = S(K_arr[j, 0], dose, alpha=alpha, beta=beta, Km=Km)*K_arr[j, 1]/total_area

        print('Dose: {} --> Surv frac: {}'.format(dose, dS_frac.sum()))
        return (dS_frac.sum()-0.01)**2

    dose = fmin(get_S_diff, 40, ftol=.001, xtol=.1)[0] #, args=(), x
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
        Km = params['K_m']
        alpha = .3 # Powathil 2012 /see ../Refs/Powathil2012
        beta = .03 # Powathil 2012

        for K0 in K0s:

            params['K_0'] = K0
            try:
                path = dose_data.format(**params)
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

                        doses[j:] = get_D99_dose(data, alpha, beta, Km)
                    doses_arr[i, :] = doses


                np.save(dose_data.format(**params), doses_arr)


        plot_D99(figure_folder_name.format(**params), params, figw, figh=2)

