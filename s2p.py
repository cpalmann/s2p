#!/usr/bin/env python

# s2p - Satellite Stereo Pipeline
# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@polytechnique.org>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import shutil
import os.path
import datetime
import traceback
import numpy as np
import multiprocessing

from python.config import cfg
from python import common
from python import initialization
from python import preprocess
from python import globalvalues
from python import process
from python import globalfinalization


def show_progress(a):
    """
    Print the number of tiles that have been processed.

    Args:
        a: useless argument, but since this function is used as a callback by
            apply_async, it has to take one argument.
    """
    show_progress.counter += 1
    status = "done {:{fill}{width}} / {} tiles".format(show_progress.counter,
                                                       show_progress.total,
                                                       fill='',
                                                       width=len(str(show_progress.total)))
    if show_progress.counter < show_progress.total:
        status += chr(8) * len(status)
    else:
        status += '\n'
    sys.stdout.write(status)
    sys.stdout.flush()


def print_elapsed_time(since_first_call=False):
    """
    Print the elapsed time since the last call or since the first call.

    Args:
        since_first_call:
    """
    t2 = datetime.datetime.now()
    if since_first_call:
        print "Total elapsed time:", t2 - print_elapsed_time.t0
    else:
        try:
            print "Elapsed time:", t2 - print_elapsed_time.t1
        except AttributeError:
            print t2 - print_elapsed_time.t0
    print_elapsed_time.t1 = t2


def preprocess_tile(tile_info):
    """
    Compute pointing corrections and extrema intensities for a single tile.

    Args:
        tile_info: dictionary containing all the information needed to process a
            tile.
    """
    # create output directory for the tile
    tile_dir = tile_info['directory']
    if not os.path.exists(tile_dir):
        os.makedirs(tile_dir)

    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open(os.path.join(tile_dir, 'stdout.log'), 'w', 0)
        # the last arg '0' is for no buffering
        sys.stdout = fout
        sys.stderr = fout
        print 'preprocess_tiles on %s ...' % tile_info

    try:
        preprocess.pointing_correction(tile_info)
        preprocess.minmax_color_on_tile(tile_info)
    except Exception:
        print("Exception in preprocessing tile:")
        traceback.print_exc()
        raise

    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()


def global_values(tiles_full_info):
    """
    Compute the global pointing correction and extrema intensities for the ROI.
    """
    globalvalues.pointing_correction(tiles_full_info)
    globalvalues.minmax_intensities(tiles_full_info)


def get_disparity_maps(tile_info, pair_id):
    """
    Process a pair of images on a given tile.

    It includes rectification, disparity estimation and triangulation.

    Args:
        tile_info: dictionary containing all the information needed to process a
            tile.
        pair_id: index of the pair to process
    """
    # read all the information
    tile_dir = tile_info['directory']
    col, row, tw, th = tile_info['coordinates']
    images = cfg['images']

    img1, rpc1 = images[0]['img'], images[0]['rpc']
    img2, rpc2 = images[pair_id]['img'], images[pair_id]['rpc']

    A_global = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_%d.txt' % pair_id)

    # check whether the pair must be processed
    pair_dir = os.path.join(tile_dir, 'pair_%d' % (pair_id))
    if os.path.isfile(os.path.join(pair_dir, 'dont_process_this_pair.txt')):
        print 'Pair %s will not be processed, skip' % pair_dir
        return

    print 'processing tile %d %d...' % (col, row)

    print 'recovering neighborhood'
    neighboring_tile_dir = list()
    for neighbor in tile_info['neighborhood']:
        neighboring_pair_dir = os.path.join(neighbor, 'pair_%d' % (pair_id))
        if os.path.isfile(os.path.join(neighboring_pair_dir, 'dont_process_this_pair.txt')) == False:
            neighboring_tile_dir.append(neighboring_pair_dir)

    # rectification
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(pair_dir, 'disp_min_max.txt')) and
        os.path.isfile(os.path.join(pair_dir, 'rectified_ref.tif')) and
        os.path.isfile(os.path.join(pair_dir, 'rectified_sec.tif'))):
        print '\trectification on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\trectifying tile %d %d (pair %d)...' % (col, row, pair_id)
        process.rectify(pair_dir, np.loadtxt(A_global), img1, rpc1,
                        img2, rpc2, col, row, tw, th, None, neighboring_tile_dir )

    # disparity estimation
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(pair_dir, 'rectified_mask.png')) and
        os.path.isfile(os.path.join(pair_dir, 'rectified_disp.tif'))):
        print '\tdisparity estimation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\testimating disparity on tile %d %d (pair %d)...' % (col, row, pair_id)
        process.disparity(tile_dir, pair_id, img1, rpc1, img2, rpc2, col, row,
                          tw, th, None)

    if(cfg['clean_intermediate']):
        common.remove_if_exists(os.path.join(pair_dir,'disp_min_max.txt'))
        common.remove_if_exists(os.path.join(pair_dir,'rectified_ref.tif'))
        common.remove_if_exists(os.path.join(pair_dir,'rectified_sec.tif'))
        common.remove_if_exists(os.path.join(pair_dir,'cloud_water_image_domain_mask.png'))
        common.remove_if_exists(os.path.join(pair_dir,'rectified_cloud_water_image_domain_mask.png'))


def process_tile(tile_info):
    """
    Process a tile in order to obtain a height map.

    Args:
        tile_info: a dictionary that provides all you need to process a tile
    """

    # read all the information
    tile_dir = tile_info['directory']
    col, row, tw, th = tile_info['coordinates']
    nb_pairs = tile_info['number_of_pairs']


    #A_global = os.path.join(cfg['out_dir'],
    #                        'global_pointing_pair_%d.txt' % pair_id)

    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open('%s/stdout.log' % tile_dir, 'a', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout
        print 'process_tiles on %s ...' % tile_info

    try:
        # check that the tile is not masked
        if os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):
            print 'tile %s already masked, skip' % tile_dir
            return

        # process each pair to get the disparity maps
        for pair_id in xrange(1, nb_pairs + 1):
            get_disparity_maps(tile_info, pair_id)

        # triangulation
        if (cfg['skip_existing'] and
            os.path.isfile(os.path.join(tile_dir, 'height_map.tif'))):
            print '\ttriangulation on tile %d %d already done, skip' % (col, row)
        else:
            print '\ttriangulating tile %d %d...' % (col, row)
            process.triangulate(tile_info, None)

        # finalization
        process.finalize_tile(tile_info, cfg['utm_zone'])

        # ply extrema
        common.run("plyextrema {} {}".format(tile_dir, os.path.join(tile_dir, 'plyextrema.txt')))

    except Exception:
        print("Exception in processing tile:")
        traceback.print_exc()
        raise

    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    if cfg['clean_intermediate']:
        for i in xrange(nb_pairs):
           shutil.rmtree(os.path.join(tile_dir,'pair_%d' %(i+1)))
        common.remove_if_exists(os.path.join(tile_dir,'roi_ref.tif'))
        common.remove_if_exists(os.path.join(tile_dir,'rpc_err.tif'))
        common.remove_if_exists(os.path.join(tile_dir,'height_map.tif'))
        common.remove_if_exists(os.path.join(tile_dir,'nb_views.tif'))
        common.remove_if_exists(os.path.join(tile_dir,'local_merged_height_map.tif'))
        common.remove_if_exists(os.path.join(tile_dir,'roi_color_ref.tif'))
        common.remove_if_exists(os.path.join(tile_dir,'local_minmax.txt'))
        common.remove_if_exists(os.path.join(tile_dir,'applied_minmax.txt'))
        common.remove_if_exists(os.path.join(tile_dir,'roi_color_ref.tif'))
        common.remove_if_exists(os.path.join(tile_dir,'roi_ref_crop.tif'))

def global_extent(tiles_full_info):
    """
    Compute the global extent from the extrema of each ply file
    """
    xmin, xmax, ymin, ymax = float('inf'), -float('inf'), float('inf'), -float('inf')

    for tile in tiles_full_info:
        plyextrema_file = os.path.join(tile['directory'], 'plyextrema.txt')

        if (os.path.exists(plyextrema_file)):
            extremaxy = np.loadtxt(plyextrema_file)
            xmin = min(xmin, extremaxy[0])
            xmax = max(xmax, extremaxy[1])
            ymin = min(ymin, extremaxy[2])
            ymax = max(ymax, extremaxy[3])

    global_extent = [xmin, xmax, ymin, ymax]
    np.savetxt(os.path.join(cfg['out_dir'], 'global_extent.txt'), global_extent,
               fmt='%6.3f')


def compute_dsm(args):
    """
    Compute the DSMs

    Args:
         - args  ( <==> [config_file,number_of_tiles,current_tile])
    """
    list_of_tiles_dir = os.path.join(cfg['out_dir'],'list_of_tiles.txt')

    config_file,number_of_tiles,current_tile = args

    dsm_dir = os.path.join(cfg['out_dir'],'dsm')
    out_dsm = os.path.join(dsm_dir,'dsm_%d.tif' % (current_tile) )

    if cfg['global_extent'] is None:
        extremaxy = np.loadtxt(os.path.join(cfg['out_dir'], 'global_extent.txt'))
    else :
        extremaxy = cfg['global_extent']

    global_xmin,global_xmax,global_ymin,global_ymax = extremaxy

    global_y_diff = global_ymax-global_ymin
    tile_y_size = (global_y_diff)/(number_of_tiles)

    # horizontal cuts
    ymin = global_ymin + current_tile*tile_y_size
    ymax = ymin + tile_y_size

    # cutting info
    x, y, w, h, z, ov, tw, th, nb_pairs = initialization.cutting(config_file)
    range_y = np.arange(y, y + h - ov, th - ov)
    range_x = np.arange(x, x + w - ov, tw - ov)
    colmin, rowmin, tw, th = common.round_roi_to_nearest_multiple(z, range_x[0], range_y[0], tw, th)
    colmax, rowmax, tw, th = common.round_roi_to_nearest_multiple(z, range_x[-1], range_y[-1], tw, th)
    cutsinf = '%d %d %d %d %d %d %d %d' % (rowmin, th - ov, rowmax, colmin, tw - ov, colmax, tw, th)

    flags={}
    flags['average-orig']=0
    flags['average']=1
    flags['variance']=2
    flags['min']=3
    flags['max']=4
    flags['median']=5
    flags['interpol-asympt']=6
    flags['interpol-gauss']=7
    flags['interpol-sigmoid']=8
    flag = "-flag %d" % ( flags.get(cfg['dsm_option'],0) )
    radius = "-radius %d" % ( cfg['dsm_radius'] )
    pinterp = "-pinterp %d" % ( cfg['dsm_pinterp'] )
    minnonan = "-minnonan %d" % ( cfg['dsm_min_nonan'] )

    if (ymax <= global_ymax):
        common.run("plytodsm %s %s %s %s %f %s %f %f %f %f %s %s" % (
                                                 flag,    #%s
                                                 radius,  #%s
                                                 pinterp, #%s
                                                 minnonan, #%s
                                                 cfg['dsm_resolution'], #%f
                                                 out_dsm, #%s
                                                 global_xmin, #%f
                                                 global_xmax, #%f
                                                 ymin, #%f
                                                 ymax, #%f
                                                 cutsinf, #%s
                                                 cfg['out_dir'])) #%s

def global_finalization(tiles_full_info):
    """
    Produce a single height map, DSM and point cloud for the whole ROI.

    The height maps associated to each pair, as well as the height map obtained
    by merging all the pairs, are stored as VRT files. The final DSM is obtained
    by projecting the 3D points from the ply files obtained on each tile. The
    final point cloud is obtained as the union of all the locally merged point
    clouds, in the LidarViewer format.

    Args:
        tiles_full_info: dictionary providing all the information about the
            processed tiles
    """
    globalfinalization.write_vrt_files(tiles_full_info)
    globalfinalization.write_dsm()

    # whole point cloud (LidarViewer format)
    if common.which('LidarPreprocessor'):
        out = os.path.join(cfg['out_dir'], 'cloud.lidar_viewer')
        plys = []
        for tile_info in tiles_full_info:
            plys.append(os.path.join(os.path.abspath(tile_info['directory']),
                                     'cloud.ply'))
        globalfinalization.lidar_preprocessor(out, plys)

    # copy RPC xml files in the output directory
    for img in cfg['images']:
        shutil.copy2(img['rpc'], cfg['out_dir'])


def launch_parallel_calls(fun, list_of_args, nb_workers, extra_args=None):
    """
    Run a function several times in parallel with different given inputs.

    Args:
        fun: function to be called several times in parallel.
        list_of_args: list of (first positional) arguments passed to fun, one
            per call
        nb_workers: number of calls run simultaneously
        extra_args (optional, default is None): tuple containing extra arguments
            to be passed to fun (same value for all calls)
    """
    results = []
    show_progress.counter = 0
    pool = multiprocessing.Pool(nb_workers)
    for x in list_of_args:
        args = (x,) + extra_args if extra_args else (x,)
        results.append(pool.apply_async(fun, args=args, callback=show_progress))

    for r in results:
        try:
            r.get(3600)  # wait at most one hour per call
        except multiprocessing.TimeoutError:
            print "Timeout while running %s" % str(r)
        except common.RunFailure as e:
            print "FAILED call: ", e.args[0]["command"]
            print "\toutput: ", e.args[0]["output"]
        except ValueError as e:
            print traceback.format_exc()
            print str(r)
            pass
        except KeyboardInterrupt:
            pool.terminate()
            sys.exit(1)


    pool.close()
    pool.join()


def execute_job(config_file,params):
    """
    Execute a job.

    Args:
         - json config file
         - params  ( <==> [tile_dir,step,...])
    """
    tile_dir = params[0]
    step = int(params[1])

    tiles_full_info = initialization.init_tiles_full_info(config_file)

    if not (tile_dir == 'all_tiles' or 'dsm' in tile_dir ):
        for tile in tiles_full_info:
            if tile_dir == tile['directory']:
                tile_to_process = tile
                break

    try:

        if step == 2:#"preprocess_tiles":
            print 'preprocess_tiles on %s ...' % tile_to_process
            preprocess_tile(tile_to_process)

        if step == 3:#"global_values":
            print 'global values ...'
            global_values(tiles_full_info)

        if step == 4:#"process_tiles" :
            print 'process_tiles on %s ...' % tile_to_process
            process_tile(tile_to_process)

        if step == 5:#"global extent" :
            print 'global extent ...'
            global_extent(tiles_full_info)

        if step == 6:#"compute_dsm" :
            print 'compute_dsm ...'
            current_tile=int(tile_dir.split('_')[1]) # for instance, dsm_2 becomes 2
            compute_dsm([config_file,cfg['dsm_nb_tiles'],current_tile])

        if step == 7:#"global_finalization":
            print 'global finalization...'
            global_finalization(tiles_full_info)

    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]


def list_jobs(config_file, step):

    tiles_full_info = initialization.init_tiles_full_info(config_file)
    filename = str(step) + ".jobs"

    if not (os.path.exists(cfg['out_dir'])):
        os.mkdir(cfg['out_dir'])

    if step in [2,4]:           #preprocessing, processing
        f = open(os.path.join(cfg['out_dir'],filename),'w')
        for tile in tiles_full_info:
            tile_dir = tile['directory']
            f.write(tile_dir + ' ' + str(step) + '\n')
        f.close()
    elif step in [3,5,7]:       # global values, global extent, finalization
        f = open(os.path.join(cfg['out_dir'],filename),'w')
        f.write('all_tiles ' + str(step) + '\n')
        f.close()
    elif step ==6 :             # compute dsm
        f = open(os.path.join(cfg['out_dir'],filename),'w')
        for i in range(cfg['dsm_nb_tiles']):
            f.write('dsm_'+ str(i) + ' ' + str(step) + '\n')
        f.close()
    else:
        print "Unkown step required: %s" % str(step)


def main(config_file, step=None, clusterMode=None, misc=None):
    """
    Launch the entire s2p pipeline with the parameters given in a json file.

    It is a succession of six steps:
        initialization
        preprocessing
        global_values
        processing
        compute dsms
        global_finalization

    Args:
        config_file: path to a json configuration file
        step: integer between 1 and 5 specifying which step to run. Default
        value is None. In that case all the steps are run.
    """
    print_elapsed_time.t0 = datetime.datetime.now()

    if clusterMode == 'list_jobs':
        list_jobs(config_file, step)
    elif clusterMode == 'job':
        cfg['omp_num_threads'] = 1
        execute_job(config_file,misc)
    else:
        # determine which steps to run
        steps = [step] if step else [1, 2, 3, 4, 5, 6, 7]

        # initialization (has to be done whatever the queried steps)
        initialization.init_dirs_srtm(config_file)
        tiles_full_info = initialization.init_tiles_full_info(config_file)

        # multiprocessing setup
        nb_workers = multiprocessing.cpu_count()  # nb of available cores
        if cfg['max_nb_threads']:
            nb_workers = min(nb_workers, cfg['max_nb_threads'])

        # omp_num_threads: should not exceed nb_workers when multiplied by the
        # number of tiles
        cfg['omp_num_threads'] = max(1, int(nb_workers / len(tiles_full_info)))

        # do the job
        if 2 in steps:
            print '\npreprocessing tiles...'
            show_progress.total = len(tiles_full_info)
            launch_parallel_calls(preprocess_tile, tiles_full_info, nb_workers)
            print_elapsed_time()

        if 3 in steps:
            print '\ncomputing global values...'
            global_values(tiles_full_info)
            print_elapsed_time()

        if 4 in steps:
            print '\nprocessing tiles...'
            show_progress.total = len(tiles_full_info)
            launch_parallel_calls(process_tile, tiles_full_info, nb_workers)
            print_elapsed_time()

        if 5 in steps:
            print '\ncomputing global extent...'
            global_extent(tiles_full_info)
            print_elapsed_time()

        if 6 in steps:
            print '\ncompute dsm...'
            args = []
            for i in range(cfg['dsm_nb_tiles']):
                args.append([config_file, cfg['dsm_nb_tiles'], i])
            show_progress.total = cfg['dsm_nb_tiles']
            launch_parallel_calls(compute_dsm, args, nb_workers)
            print_elapsed_time()

        if 7 in steps:
            print '\nglobal finalization...'
            global_finalization(tiles_full_info)
            print_elapsed_time()

    # cleanup
    print_elapsed_time(since_first_call=True)
    common.garbage_cleanup()


if __name__ == '__main__':

    error = False
    steps=[1,2,3,4,5,6,7]

    if len(sys.argv) < 2:
        error = True

    elif sys.argv[1].endswith(".json"):
        if len(sys.argv) == 2:
            main(sys.argv[1])
        elif len(sys.argv) == 3 and int(sys.argv[2]) in steps:
            main(sys.argv[1], int(sys.argv[2]))
        else:
            error = True
    else:  # cluster modes
        if sys.argv[1] not in ['list_jobs', 'job']:
            error = True
        else:
            if sys.argv[1] == 'list_jobs':
                if len(sys.argv) == 4 and int(sys.argv[3]) in steps:
                    main(sys.argv[2], int(sys.argv[3]), 'list_jobs')
                else:
                    error = True

            if sys.argv[1] == 'job':
                if len(sys.argv) >= 5 and int(sys.argv[4]) in steps:
                    main(sys.argv[2], None, 'job', sys.argv[3:])
                else:
                    error = True
    if error:
        print """
        Incorrect syntax, use:
          > %s config.json [step (integer between 1 and 7)]
            1: initialization
            2: preprocessing (tilewise sift, local pointing correction)
            3: global-pointing
            4: processing (tilewise rectification, matching and triangulation)
            5: global-extent
            6: compute dsm from ply files (one per tile)
            7: finalization
            Launches the s2p pipeline.

          > %s list_jobs config.json step (integer between 2 and 7)
            Return the list of jobs for a specific step.

          > %s job config.json tile_dir step (integer between 2 and 7)
            Run a specific job defined by a json string. This mode allows to run jobs returned
            by the list_jobs running mode.


          All the parameters, paths to input and output files, are defined in
          the json configuration file.

        """ % (sys.argv[0], sys.argv[0], sys.argv[0])
        sys.exit(1)
