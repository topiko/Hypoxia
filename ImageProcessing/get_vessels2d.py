import matplotlib.pyplot as plt
import sys
import numpy as np
#from scipy.misc import imresize
from PIL import Image
from skimage import filters, measure

PATH = '../VesselData/AlignedImages/'
REMAKE = True #True #False #True

if len(sys.argv) > 1:
    NPIXELS = int(sys.argv[1])
else:
    NPIXELS = 100 #None #168 #84
# orig paper voxel edge l = 12e-6m
# 167*167*167 voxel grid where tumor is mirrored.
# --> roughly 90 pixels for tumor area required.
# Tumor is roughly half of the image --> scale so that
# roughly 200 pixels for width...
VESSEL_THRESHOLD = 64
N_MIN_CONTOUR = 12

IMG_W = 1176 #[mum]
IMG_H = 882 #[mum] 12mum layer pixel size in aligned images == layer thickness --> 98, 74 pixels

# Crop this area from the original images:
area = (600, 600)
# X crop from here
x_crop_left = 500
# Y crop from here:
y_crop_low = 300

OUTPUT_PATH = '../VesselData/'

def get_init_vessels_dict():
    """
    Make initial vessels dict containing information of the bounding box and
    pixel scales.
    """
    im_data = Image.open(PATH + '1,{:02d} merge.png'.format(1))

    xpixels_orig, ypixels_orig = im_data.size

    pixel_w = IMG_W/xpixels_orig
    pixel_h = IMG_H/ypixels_orig

    if pixel_w != pixel_h:
        raise ValueError('Pixels ought to be squares: {} != {}'.format(pixel_w, pixel_h))

    cropped_fig_w = area[0]*pixel_w
    cropped_fig_h = area[1]*pixel_h

    origin = np.array([x_crop_left*pixel_w, y_crop_low*pixel_h])

    return {'origin': origin,
            'W': cropped_fig_w * 1e-6,
            'H': cropped_fig_h * 1e-6,
            'units': ('mum', 1e-6),
            'pixel_size': pixel_w*area[1]/NPIXELS}

def get_vessels(im_dat, vessel_thres=VESSEL_THRESHOLD):
    """
    Retunr the boolean vessel array given the vessel image data.
    """
    vessel_channel = im_dat[:, :, 1]
    return vessel_channel > vessel_thres


def load_image(i, show=False):

    im_data = Image.open(PATH + '1,{:02d} merge.png'.format(i))
    xpixels_orig, ypixels_orig = im_data.size

    if area[0] != area[1]:
        raise ValueError('Only square images supported')

    x_crop_right = x_crop_left + area[0]
    y_crop_top = y_crop_low + area[1]
    im_data = im_data.crop((x_crop_left, y_crop_low, x_crop_right, y_crop_top))
    if (i == 1) and show:
        im_data.show()

    if NPIXELS is not None:
        im_data = im_data.resize((NPIXELS, NPIXELS))
        if (i == 1) and show:
            im_data.show()

    return  np.array(im_data)

def get_vessels3d(remake=False):
    """
    Build the 3D vessel array givven the vessel images.
    """
    try:
        vessels3d = np.load(OUTPUT_PATH + 'vessels3d_bool_NP={}.npy'.format(NPIXELS))
        if not remake:
            print('Found data rdy formatted 3d vessels.')
            return vessels3d
    except FileNotFoundError:
        pass

    vessels3d = [] #None
    for i in range(1, 15):
        im_data = load_image(i)
        vessels = get_vessels(im_data)
        vessels3d.append(vessels)

    vessels3d = np.stack(vessels3d, axis=2)
    np.save(OUTPUT_PATH + 'vessels3d_bool_NP={}.npy'.format(NPIXELS), vessels3d)

    return vessels3d

def get_neighbors(i, j, nx, ny, k=None, nz=None):

    imin = max(0, i-1)
    imax = min(i+1, nx-1)
    jmin = max(0, j-1)
    jmax = min(j+1, ny-1)

    prop_neighbors = np.array([[imax, j],# [imax, jmax],
                               [i, jmax],# [imin, jmax],
                               [imin, j],# [imin, jmin],
                               [i, jmin]]) #,# [imax, jmin]])

    mask = (prop_neighbors == np.array([i, j])).all(axis=1)
    return prop_neighbors[~mask]

def pixels_to_lenghts(contour, pixel_size, origin, unit):

    contour_physical = contour*pixel_size
    contour_physical *= unit[1]

    return contour_physical

def make_2d_vessel_structure():

    vessels3d = get_vessels3d(REMAKE)
    nx, ny, nz = vessels3d.shape

    if nx != ny:
        raise ValueError('Only square pictures supported!')

    for iz in range(nz):

        _, (ax1, ax2, ax3) = plt.subplots(1, 3,
                                          figsize=(6, 3),
                                          sharey=True,
                                          sharex=True)

        # PIxels and coordinates... Need to mirror
        image = vessels3d[::-1, :, iz]
        edge_image = filters.roberts(image)

        im_data = load_image(iz+1)

        # PIxels and coordinates... Need to mirror
        ax1.imshow(im_data[::-1, :, :])
        ax2.imshow(edge_image, cmap='gray')

        # Mirror the y-axis since pixels and coordinates are handled differently
        contours = measure.find_contours(vessels3d[::-1, :, iz], .95)
        vessels_dict = get_init_vessels_dict()
        i = 0
        for n, contour in enumerate(contours):
            contour = contour[:, ::-1]
            c = 'C0'
            if (contour == 0).any() or \
               (contour == nx-1).any():
                #print('Touches edge!!')
                c = 'red'
            if len(contour) < N_MIN_CONTOUR:
                #print('Discarding contour')
                continue
            if c == 'red':
                contour = np.append(contour, contour[0].reshape(-1, 2), axis = 0)

                ax3.scatter(contour[0, 0], contour[0, 1], c=c, s=5)
            ax3.plot(contour[:, 0], contour[:, 1], linewidth=1, c=c)

            contour_physical = pixels_to_lenghts(contour,
                                                 vessels_dict['pixel_size'],
                                                 vessels_dict['origin'],
                                                 vessels_dict['units'])

            # Mirror the yaxis to match coordinate and pixel representations
            contour_physical[:, 1] = ny*vessels_dict['pixel_size']*vessels_dict['units'][1] \
                    - contour_physical[:, 1]
            vessels_dict['vessel{:04d}'.format(i)] = {'edge': contour_physical[:-1]}
            i += 1

        # TODO: Figure out pixel size! Further scale with that pixel size to get correct x-y coord
        ax3.imshow(image, cmap='gray')
        plt.savefig(OUTPUT_PATH + 'Figures/layer_np={}_{}_vessels.pdf'.format(NPIXELS, iz+1))
        np.save(OUTPUT_PATH + 'vessels_dict_np={}_iz={}.npy'.format(NPIXELS, iz+1), vessels_dict)
        plt.title('Layer idx = {}'.format(iz+1))
        plt.show()

if __name__ == '__main__':
    make_2d_vessel_structure()
