import sys
import argparse

res_dir = 'Data/diffusion_iz={layer_idx:02d}_np={npixels}_activeV={active_vessels:02d}_K0={K_0:g}_D={D:g}_Km={K_m:g}_C0={C_0:g}_rid={run_idx:03d}/'

dose_data = 'Data/D99_doses_layer_idx={layer_idx}_K0={K_0}_D={D}_Km={K_m}_C0={C_0}_alpha={alpha}_beta={beta}_eq={oereq}.npy'

terminated_flag = 'fin'

def K_names(t):
    # Takes care of the bit crptic K array name genration...
    if isinstance(t, str):
        if t == 'fin':
            K_name = 'K_arr_fin.'
        else:
            t = int(t)

    if isinstance(t, float):
        t = int(1000*t)
    if isinstance(t, int):
        K_name = 'K_arr_t={}ms.'.format(t)

    return K_name + 'pdf', K_name + 'npy'

def get_cmd_args_parser():

    parser = argparse.ArgumentParser('Hypoxia')

    parser.add_argument('-layer_idx', type=int, default=1)
    parser.add_argument('-npixels', type=int, default=100)
    parser.add_argument('-active_vessels', metavar='Active vessel percentage', type=int, default=50)
    parser.add_argument('-K_0', type=float, default=40)
    parser.add_argument('-D', type=float, default=2000)
    parser.add_argument('-K_m', type=float, default=1)
    parser.add_argument('-C_0', type=float, default=15)
    parser.add_argument('-run_idx', type=int, default=0)
    parser.add_argument('-alpha', type=float, default=0.3)
    parser.add_argument('-beta', type=float, default=0.03)
    parser.add_argument('-K_rm', type=float, default=3)
    parser.add_argument('-oereq', type=str, default='powathil')

    return parser
