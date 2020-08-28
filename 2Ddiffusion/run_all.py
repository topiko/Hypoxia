import subprocess

layer_idxs = [3, 4] #,4,5,7,9,10]
npixels_arr = [100]
active_vessels_arr = [30, 50, 75, 100]
K0s = [10, 20, 30, 40, 50]
Ds = [1500, 2000, 2500]
rids = range(5)

for layer_idx in layer_idxs:
    for npixels in npixels_arr:
        for active_vessels in active_vessels_arr:
            rids_run = [0] if active_vessels==100 else rids
            for K0 in K0s:
                for D in Ds:
                    for rid in rids_run:
                        params = ['-layer_idx', str(layer_idx),
                                  '-npixels', str(npixels),
                                  '-active_vessels', str(active_vessels),
                                  '-K_0', str(K0),
                                  '-D', str(D)] #,
                        subprocess.run(['python', 'model.py'] + params + ['-run_id', str(rid)])
                        #, str(layer_idx), str(npixels), str(active_vessels), str(K0), str(rid)])

                    print('Plot res:')
                    subprocess.run(['python', 'plot_utils.py'] + params)
