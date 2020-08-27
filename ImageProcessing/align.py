import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.optimize import minimize
from scipy.misc import imresize
import numpy as np

path_to_aligned = '../VesselData/AlignedImages/'

def shift_rotate(angle, shift_x, shift_y, fig_dat):

    nx, ny, nz = fig_dat.shape
    shift_x_pixels = shift_x / (IMG_W / nx)
    shift_y_pixels = shift_y / (IMG_H / ny)
    print('Shift in pixes:', [shift_x_pixels, shift_y_pixels])

    rot_top = ndimage.rotate(dat_top, angle, reshape=False)
    top_rot_shift = ndimage.shift(rot_top,
                                  [shift_x_pixels,
                                   shift_y_pixels,
                                   0])

    return top_rot_shift


def get_distance(x,
                 dat_bot=None,
                 dat_top=None,
                 return_images=False):

    angle = x[0]
    shift_x = x[1]
    shift_y = x[2]

    n = 6
    nx, ny, nz = dat_bot.shape
    mask = np.zeros((nx, ny, nz), dtype=bool)
    mask[nx//n:nx*(n-1)//n, ny//n:ny*(n-1)//n, :] = True

    top_rot_shift = shift_rotate(angle,
                                 shift_x,
                                 shift_y,
                                 dat_top)

    score = np.sqrt(np.sum((top_rot_shift[mask] - dat_bot[mask])**2))/mask.sum()
    if return_images:
        return score, dat_bot, top_rot_shift

    return score, None, None

if __name__ == '__main__':

    IMG_W = 1176 #[mum]
    IMG_H = 888 #[mum] 12mum layer pixel size in aligned images == layer thickness --> 98, 74 pixels
    N = 40
    nfigs = 14
    angle_max = 7 #[deg]
    step_angle = .5
    shift_max = 100 #[mum]
    step_shift = 10

    resize = 20
    plot = False
    rot_shift_dict = {}

    fname = '1,{:02d} merge.tif'
    dat_bot = plt.imread(fname.format(1)).copy()
    # Save the most bottom image
    # Save the new bottom to aligned::
    plt.imsave(path_to_aligned + fname.format(1).replace('tif', 'png'), dat_bot)

    # Loop over the remaining images:
    for i in range(1, nfigs):
        dat_top = plt.imread(fname.format(i+1)).copy()

        # resize for speedup:
        # TODO: fix to use PIL and its resize instead of scipy.
        dat_bot_resized = imresize(dat_bot, resize)[:, :, 1:] #Drop the red channel
        dat_top_resized = imresize(dat_top, resize)[:, :, 1:]

        print('Resized image size:', dat_top_resized.shape)
        min_score = 1000000
        # Rotate and shift:
        for angle in np.arange(-angle_max, angle_max+step_angle, step_angle):
            for shift_x in np.arange(-shift_max, shift_max+step_shift, step_shift):
                for shift_y in np.arange(-shift_max, shift_max+step_shift, step_shift):

                    print('Angle {:.2f}, sx={:.1f}, sy={:.1f}'.format(angle, shift_x, shift_y))
                    score, dat_bot_, dat_top_ = get_distance([angle, shift_x, shift_y],
                                                             dat_bot_resized,
                                                             dat_top_resized,
                                                             plot)
                    if score < min_score:
                        print('Score improved {:.5f} --> {:.5f}'.format(min_score, score))
                        min_score = score
                        rot_shift_dict[i+1] = [angle, shift_x, shift_y]
                        if plot:
                            _, _, shaped_top = get_distance([angle, shift_x, shift_y],
                                                            dat_bot,
                                                            dat_top,
                                                            True)
                            _, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
                            ax1.imshow(dat_bot, aspect='equal')
                            ax2.imshow(shaped_top, aspect='equal')

                            ax2.set_title(fname.format(i+1))
                            ax1.set_axis_off()
                            ax2.set_axis_off()
                            plt.suptitle('angle={:.1f},'.format(angle) \
                                         + 'shift_x={:.1f}, '.format(shift_x) \
                                         + 'shift_y={:.1f} - Score {:.5f}'.format(shift_y, score))
                            plt.show()
        # The new bottom image is the current top one:
        _, _, dat_bot = get_distance(rot_shift_dict[i+1],
                                     dat_bot,
                                     dat_top,
                                     True)
        # Save the new bottom to aligned:
        plt.imsave(path_to_aligned + fname.format(i+1).replace('tif', 'png'), dat_bot)
        np.save(path_to_aligned + 'rot_shift_dict.npy', rot_shift_dict)


