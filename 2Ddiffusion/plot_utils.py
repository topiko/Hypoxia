import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from input_utils import res_dir, K_names, get_cmd_args_parser, terminated_flag


test = False #True #False
show = test #False

error_bar_dict = {'elinewidth':1, 'capsize':2}
show_step = 'first' # 'full' False
figw = 3
figh = 5.5
figure_folder_name = 'Pictures_layer_idx={layer_idx}/' #.format(**params)

def plot_oxygen_hmap(data, fig, ax):

    data_scaled = data.copy()
    data_scaled[:,1:4] *= 1e6 # convert to mum
    patches = []
    colors = []
    N = 10000 if test else 1000000

    for i in range(len(data)//3)[:N]:
        cell_data = data_scaled[3*i:3*(i+1)]
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
    ax.set_xlim([data_scaled[:, 1].min(), data_scaled[:,1].max()])
    ax.set_ylim([data_scaled[:, 2].min(), data_scaled[:,2].max()])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\mu m$')
    ax.set_ylabel(r'$\mu m$')

    return ax

def plot_K_on_diag(datas, ax3):

    def get_K_on_diag_(data):
        rK = []
        for i in range(len(data)//3):
            cell_data = data[3*i:3*(i+1)]
            v1 = cell_data[0,1:4]
            v2 = cell_data[1,1:4]
            v3 = cell_data[2,1:4]

            x = cell_data[:, 1]
            y = cell_data[:, 2]

            signs = np.sign(x-y)
            if (1 in signs) and (-1 in signs):
                # on x==y line
                x = x.mean()
                y = y.mean()
                r = np.sqrt(x**2 + y**2)
                rK.append([r, cell_data[:,4].mean()])

        rK = np.array(rK)
        return rK[rK[:,0].argsort(), :]


    for data in datas:
        rK = get_K_on_diag_(data)
        rK = rK[len(rK)//2:]
        ax3.plot(rK[:, 0], rK[:,1])
        idx_max = rK[:, 1].argmax() #[-1]
        val_max = rK[:, 1].max()
        ax3.vlines(rK[idx_max, 0], 0, val_max)
        idx_half = ((rK[:, 1] - rK[idx_max, 1]/2)**2).argmin()
        ax3.vlines(rK[idx_half, 0], 0, rK[idx_half, 1],
                   label=r'half dist = {:.3f}\mu m'.format((-rK[idx_max, 0] + rK[idx_half, 0])*1e6))
        ax3.legend(frameon=False)

def plot_oxygen_hist(datas, params, ax, nbins=25):

    # Build histogram of oxygen level in tissue.
    # Each cell in tissue is weighted by its area.

    if params['layer_idx'] == -1:
        bins = np.linspace(0, 50, nbins+1)
    else:
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
            #v1 = cell_data[0,1:4]
            #v2 = cell_data[1,1:4]
            #v3 = cell_data[2,1:4]
            #print(v1-v2)
            #print(v1-v3)
            #print(area[i] - np.linalg.norm(np.cross((v1-v2), (v1-v3)))/2*1e6)
        oxygens[:, n], _ = np.histogram(oxygen, bins=bins, weights=area)

    print('Areas:')
    print(oxygens.sum(axis=0))

    if (oxygens.sum(axis=0).std() > .0001):
        raise ValueError('Areas should be equal... std={}'.format(oxygens.sum(axis=0).std()))

    width = bins[1] - bins[0]
    ax.bar(bins[:-1],
           oxygens.mean(axis=1),
           width=width,
           yerr=oxygens.std(axis=1),
           align='edge',
           alpha=.7,
           error_kw=error_bar_dict)
    if show_step:
        for i in range(len(datas)):
            color = 'red' if i == 0 else 'black'
            if (show_step == 'full') or (i == 0):
                ax.step(bins,
                        np.insert(oxygens[:, i], 0, oxygens[0, i]),
                        where='pre',
                        color=color,
                        lw=1,
                        alpha=.5)


    ax.set_xlabel(r'$K$, [$mmHg$]')
    ax.set_ylabel(r'Area, [$mm^2$]')

    return ax

def plot_K_arr(data, params, res_dir, figw=3, figh=4):

    title = r'Active Vessels={active_vessels}%, $K_0={K_0:g}$, $K_m={K_m:g}$'.format(**params)
    #+ '\n' \
    #        + '$D={D:g}$, $C_0={C_0:g}$, $K_m={K_m:g}$'.format(**params)

    fig_name = res_dir + title.replace(' ', '').replace(',', '_').replace('$', '').replace('%', '').replace('\n', '_') + '.png'
    if os.path.isfile(fig_name) and False: # and (not test):
        print('Figure already exists. Exiting')
        return

    fig = plt.figure(figsize=(figw, figh)) #, constrained_layout=True))

    if params['layer_idx'] != -1:
        gs = fig.add_gridspec(3, 1)
        ax1 = fig.add_subplot(gs[:2, 0])
        ax2 = fig.add_subplot(gs[2,0])
    else:
        gs = fig.add_gridspec(4, 1)
        ax1 = fig.add_subplot(gs[:2, 0])
        ax2 = fig.add_subplot(gs[2,0])
        ax3 = fig.add_subplot(gs[3,0])

    if isinstance(data, list):
        ax1 = plot_oxygen_hmap(data[0], fig, ax1)
    elif isinstance(data, np.ndarray):
        ax1 = plot_oxygen_hmap(data, fig, ax1)
        data = [data]

    ax2 = plot_oxygen_hist(data, params, ax2)

    plt.suptitle(title, fontsize=12)

    fig.subplots_adjust(hspace=.6)
    plt.gcf().subplots_adjust(right=0.95, left=.22, top=.85)

    if params['layer_idx'] == -1:
        rKs = plot_K_on_diag(datas, ax3)
        ax3.semilogy()
        ax2.semilogy()
        #ax2.set_xlim([5, 50])# axis()

    plt.savefig(fig_name, dpi=250)
    if show: plt.show()

def plot_D99(save_folder, params, figw=4, figh=3):

    from run_all import active_vessels_arr, K0s
    from input_utils import dose_data

    fig, ax = plt.subplots(figsize=(figw, figh))
    for K0 in K0s[::2]:
        params['K_0'] = K0
        doses_arr = np.load(dose_data.format(**params))
        ax.errorbar(active_vessels_arr,
                    doses_arr.mean(axis=1),
                    yerr=doses_arr.std(axis=1),
                    label=r'$K_0={}$mmHg'.format(K0),
                    **error_bar_dict)

    ax.set_xlabel('Active vessels [%]')
    ax.set_ylabel('D99-Dose [Gy]')

    ax.set_title(r'$\alpha={alpha:.03f} ~ \beta={beta:.03f} ~ K_{rm}={K_rm:.1f}$'.format(**params))
    plt.gcf().subplots_adjust(right=0.95, left=.22, bottom=.20)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(save_folder + 'D99_D={D}_Km={K_m}_C0={C_0}_alpha={alpha}_beta={beta}_oereq={oereq}_Krm={K_rm}.png'.format(**params), dpi=250)

    plt.show()

if __name__ == '__main__':


    parser = get_cmd_args_parser()
    args = parser.parse_args()

    params = vars(args)

    print('Run arguments:')
    print(vars(args))


    datas = []
    for run_id in range(20):
        params['run_idx'] = run_id
        #for layer_idx in range(12):
        #params['layer_idx'] = layer_idx
        res_dir_f = res_dir.format(**params)

        try:
            datas.append(np.load(res_dir_f + K_names(terminated_flag)[1]))
        except FileNotFoundError:
            pass

    if len(datas) == 0:
        raise FileNotFoundError('Not a single run found.')


    folder_name = figure_folder_name.format(**params)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)



    plot_K_arr(datas, params, folder_name, figw, figh)
