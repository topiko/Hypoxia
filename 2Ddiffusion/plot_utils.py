import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from input_utils import res_dir, K_names, get_cmd_args_parser


test = False
show = False

error_bar_dict = {'elinewidth':1, 'capsize':2}
figw = 3
figh = 5.5
figure_folder_name = 'Pictures_layer_idx={layer_idx}/' #.format(**params)
def plot_oxygen_hmap(data, fig, ax):

    data[:,1:4] *= 1e6 # convert to mum
    patches = []
    colors = []
    N = 100 if test else 1000000

    for i in range(len(data)//3)[:N]:
        cell_data = data[3*i:3*(i+1)]
        tri = Polygon(cell_data[:,1:3])

        colors.append(cell_data[:,4].mean())
        patches.append(tri)

    p = PatchCollection(patches)
    p.set_clim([0, 50])

    p.set_array(np.array(colors))

    ax.add_collection(p)

    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("top", size="5%", pad=0.1)

    # Colorbar.
    cbar = plt.colorbar(p, cax=color_axis, orientation="horizontal")
    cbar.set_label('K [$mmHg$]', labelpad=1)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')

    #fig.colorbar(p, ax=ax)
    ax.set_xlim([data[:, 1].min(), data[:,1].max()])
    ax.set_ylim([data[:, 2].min(), data[:,2].max()])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\mu m$')
    ax.set_ylabel(r'$\mu m$')

    return ax

def plot_oxygen_hist(datas, ax, nbins=50):

    # Buiöd histogram of oxygen level in tissue.
    # Each cell in tissue is weighted by its area.


    bins = np.linspace(0, 50, nbins+1)
    oxygens = np.zeros((nbins, len(datas)))
    for n, data in enumerate(datas):
        N = int(len(data)//3)
        area = np.zeros(N)
        oxygen = np.zeros(N)
        for i in range(len(data)//3):
            cell_data = data[3*i:3*(i+1)]

            oxygen[i] = cell_data[:,4].mean()
            area[i] = cell_data[0, 5]*1e6 # form m^2 --> mm^2

        oxygens[:, n], _ = np.histogram(oxygen, bins=bins, weights=area)

    if (oxygens.sum(axis=0).std() > .0001):
        raise ValueError('Areas should be equal... std={}'.format(oxygens.sum(axis=0).std()))

    width = bins[1] - bins[0]
    ax.bar(bins[:-1], oxygens.mean(axis=1),
           width=width,
           yerr=oxygens.std(axis=1),
           align='edge',
           alpha=.7,
           error_kw=error_bar_dict)

    ax.set_xlabel(r'$K$, [$mmHg$]')
    ax.set_ylabel(r'Area, [$mm^2$]')

    return ax

def plot_K_arr(data, t, params, res_dir, figw=3, figh=4):

    title = r'Active Vessels={active_vessels}%, $K_0={K_0:g}$'.format(**params) + '\n' \
            + '$D={D:g}$, $C_0={C_0:g}$, $K_m={K_m:g}$'.format(**params)

    fig_name = res_dir + title.replace(' ', '').replace(',', '_').replace('$', '').replace('%', '').replace('\n', '_') + '.png'
    if os.path.isfile(fig_name):
        print('Figure already exists. Exiting')
        return

    fig = plt.figure(figsize=(figw, figh)) #, constrained_layout=True))

    gs = fig.add_gridspec(3, 1)
    ax1 = fig.add_subplot(gs[:2, 0])
    ax2 = fig.add_subplot(gs[2,0])

    if isinstance(data, list):
        ax1 = plot_oxygen_hmap(data[0], fig, ax1)
    elif isinstance(data, np.ndarray):
        ax1 = plot_oxygen_hmap(data, fig, ax1)
        data = [data]

    ax2 = plot_oxygen_hist(data, ax2)

    plt.suptitle(title, fontsize=12)

    fig.subplots_adjust(hspace=.6)
    plt.gcf().subplots_adjust(right=0.95, left=.22, top=.85)
    plt.savefig(fig_name, dpi=250)

    if show: plt.show()

def plot_D99(save_folder, params, figw=4, figh=3):

    from run_all import active_vessels_arr, K0s
    from input_utils import dose_data

    fig, ax = plt.subplots(figsize=(figw, figh))
    for K0 in K0s:
        params['K_0'] = K0
        doses_arr = np.load(dose_data.format(**params))
        ax.errorbar(active_vessels_arr,
                    doses_arr.mean(axis=1),
                    yerr=doses_arr.std(axis=1),
                    label=r'K_0={}'.format(K0),
                    **error_bar_dict)

    ax.set_xlabel('Active vessels [%]')
    ax.set_ylabel('Dose [Gy] for D99')

    plt.gcf().subplots_adjust(right=0.95, left=.22, bottom=.20)
    plt.legend(frameon=False)
    plt.savefig(save_folder + 'D99.png', dpi=250)

if __name__ == '__main__':


    t = 'fin'

    parser = get_cmd_args_parser()
    args = parser.parse_args()

    params = vars(args)
    res_dir_f = res_dir.format(**params)

    print(vars(args))
    rids = []
    for entity in os.listdir('Data/'):
        if entity.startswith(res_dir_f.replace('Data/', '').split('_rid')[0]):
            rids.append(int(entity.split('rid=')[-1]))

    if len(rids) == 0:
        raise FileNotFoundError('Not a single run found.')

    datas = []
    for run_id in rids:
        params['run_idx'] = run_id
        res_dir_f = res_dir.format(**params)
        datas.append(np.load(res_dir_f + K_names(t)[1]))

    folder_name = figure_folder_name.format(**params)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)



    plot_K_arr(datas, t, params, folder_name, figw, figh)
