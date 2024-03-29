import subprocess

layer_idxs = [3] #,4,5,7,9,10]
npixels_arr = [100]
active_vessels_arr = [15, 30, 50, 75, 100]
K0s = [10, 20, 30, 40, 50]
K_m_arr = [1, 3] #1
Ds = [2000]
rids = range(10)

if __name__ == '__main__':
    for layer_idx in layer_idxs:
        for npixels in npixels_arr:
            for K_m in K_m_arr:
                for D in Ds:
                    for active_vessels in active_vessels_arr:
                        rids_run = [0] if active_vessels==100 else rids
                        for K0 in K0s:
                            for rid in rids_run:
                                params = ['-layer_idx', str(layer_idx),
                                          '-npixels', str(npixels),
                                          '-active_vessels', str(active_vessels),
                                          '-K_0', str(K0),
                                          '-D', str(D),
                                          '-K_m', str(K_m)] #,
                                subprocess.run(['python', 'model.py'] + params + ['-run_id', str(rid)])
                                #, str(layer_idx), str(npixels), str(active_vessels), str(K0), str(rid)])

                            print('Plot res:')
                            subprocess.run(['python', 'plot_utils.py'] + params)

                    subprocess.run(['python', 'radiation_utils.py'] + params)
