# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import sys
import multiprocessing
import numpy as np

from config import cfg
from python import common
from python import masking
from python import pointing_accuracy
from python import visualisation
from python import rpc_model
from python import rpc_utils


def minmax_color_on_tile(tile_info):
    """
    Compute min and max intensities on a given tile and save them to a txt file.

    Args:
        tile_info: dictionary containing all the information needed to process
            the tile
    """
    # read info
    img = cfg['images'][0]['img']
    coords = tile_info['coordinates']
    tile_dir = tile_info['directory']
    z = cfg['subsampling_factor']

    # output files
    crop_ref = os.path.join(tile_dir, 'roi_ref.tif')
    local_minmax = os.path.join(tile_dir, 'local_minmax.txt')

    # do the job
    if not (os.path.isfile(crop_ref) and cfg['skip_existing']):
        common.cropImage(img, crop_ref, *coords, zoom=z)
    else:
        print 'roi_ref.tif for tile %s already generated, skip' % tile_dir
    if os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):
        print 'tile %s is masked, skip' % tile_dir
    elif os.path.isfile(os.path.join(tile_dir, 'local_minmax.txt')) and cfg['skip_existing']:
        print 'extrema intensities on tile %s already computed, skip' % tile_dir
    else:
        common.image_getminmax(crop_ref, local_minmax)


def pointing_correction(tile_info):
    """
    Compute the translations that correct the pointing errors on a tile.

    There is one correction per pair of images.

    Args:
        tile_info: dictionary containing all the information needed to process
            the tile
    """
    tile_dir = tile_info['directory']
    x, y, w, h = tile_info['coordinates']
    
    img1, rpc1 = cfg['images'][0]['img'], cfg['images'][0]['rpc']
    roi_msk = cfg['images'][0]['roi']
    cld_msk = cfg['images'][0]['cld']
    wat_msk = cfg['images'][0]['wat']
    cwid_msk = '%s/cloud_water_image_domain_mask.png' % (tile_dir)
    
    # check if the ROI is masked by water or out of image domain
    H = np.array([[1, 0, -x], [0, 1, -y], [0, 0, 1]])
    if masking.cloud_water_image_domain(cwid_msk, w, h, H, rpc1, roi_msk, cld_msk, wat_msk):
        print "Tile masked by water or outside definition domain, skip"
        open("%s/this_tile_is_masked.txt" % tile_dir, 'a').close()

    img1, rpc1 = cfg['images'][0]['img'], cfg['images'][0]['rpc']
    roi_msk = cfg['images'][0]['roi']
    cld_msk = cfg['images'][0]['cld']
    wat_msk = cfg['images'][0]['wat']
    cwid_msk = '%s/cloud_water_image_domain_mask.png' % (tile_dir)

    # check if the ROI is masked by water or out of image domain
    H = np.array([[1, 0, -x], [0, 1, -y], [0, 0, 1]])
    if masking.cloud_water_image_domain(cwid_msk, w, h, H, rpc1, roi_msk, cld_msk, wat_msk):
        print "Tile masked by water or outside definition domain, skip"
        open("%s/this_tile_is_masked.txt" % tile_dir, 'a').close()


    for i in range(1, tile_info['number_of_pairs'] + 1):
        paired_tile_dir = os.path.join(tile_dir, 'pair_%d' % i)

        # create a directory for the experiment
        if not os.path.exists(paired_tile_dir):
            os.makedirs(paired_tile_dir)
            
        if not os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):
            
            img2, rpc2 = cfg['images'][i]['img'], cfg['images'][i]['rpc']
            
            # output files
            pointing = '%s/pointing.txt' % paired_tile_dir
            center = '%s/center_keypts_sec.txt' % paired_tile_dir
            sift_matches = '%s/sift_matches.txt' % paired_tile_dir

        if not os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):

            img2, rpc2 = cfg['images'][i]['img'], cfg['images'][i]['rpc']

            # output files
            pointing = '%s/pointing.txt' % paired_tile_dir
            center = '%s/center_keypts_sec.txt' % paired_tile_dir
            sift_matches = '%s/sift_matches.txt' % paired_tile_dir

            # check if the tile is already done
            if os.path.isfile('%s/pointing.txt' % paired_tile_dir) and cfg['skip_existing']:
                print "pointing correction on tile %d %d (pair %d) already done, skip" % (x, y, i)
            else:
                # check if the ROI is out of the secondary image domain
                r1 = rpc_model.RPCModel(rpc1)
                r2 = rpc_model.RPCModel(rpc2)
                x2, y2, w2, h2 = rpc_utils.corresponding_roi(rpc1, rpc2, x, y, w, h)
                ww2, hh2 = common.image_size(img2)

                dont_process_this_pair = ((max(0,x2)>min(ww2,w2+x2)) or (max(0,y2)>min(hh2,y2+h2)))

                if dont_process_this_pair:
                    print "The ROI is out of the secondary image domain"
                    open("%s/dont_process_this_pair.txt" % paired_tile_dir, 'a').close()

                else:
                    # correct pointing error
                    # A is the correction matrix and m is the list of sift matches
                    A, m = pointing_accuracy.compute_correction(img1, rpc1, img2,
                                                            rpc2, x, y, w, h)
                    if A is not None:
                        np.savetxt(pointing, A, fmt='%6.3f')
                    if m is not None:
                        np.savetxt(sift_matches, m, fmt='%9.3f')
                        np.savetxt(center, np.mean(m[:, 2:4], 0), fmt='%9.3f')
                        if cfg['debug']:
                            png = '%s/sift_matches_plot.png' % paired_tile_dir
                            visualisation.plot_matches_pleiades(img1, img2, rpc1,
                                                                rpc2, m, x, y, w, h,
                                                                png)
