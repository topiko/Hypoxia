import numpy as np

d = np.load('vessels_dict_np=100_iz=10.npy', allow_pickle=True).item()
d_new = {}
for key, val in d.items():
    if not key.startswith('vessel'):
        print(key)
        d_new[key] = val

#print(d)
middle = np.array([d_new['W']/2, d_new['H']/2]) #+ d_new['origin']

n = 20
vessel_arr = np.zeros((n-1, 2))
r_vessel = 1e-5
i = 0
for phi in np.linspace(0, np.pi*2, n)[:-1]:
    vessel_arr[i, :] = np.array([r_vessel*np.cos(phi), r_vessel*np.sin(phi)]) + middle
    print(vessel_arr[i])
    i+=1
    d_new['vessel0000'] = {'edge': vessel_arr}

print('W = {} H = {}'.format(d_new['W'], d_new['H']))
print(vessel_arr.max(axis=0), vessel_arr.min(axis=0))
np.save('vessels_dict_np=100_iz=-1.npy', d_new)
