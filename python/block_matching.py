# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import shutil
import common
from config import cfg


# define paths of various bm binaries
b = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 'bin')
hirschmuller02 = os.path.join(b, 'subpix.sh')
hirschmuller08 = os.path.join(b, 'callSGBM.sh')
hirschmuller08_laplacian = os.path.join(b, 'callSGBM_lap.sh')
hirschmuller08_cauchy = os.path.join(b, 'callSGBM_cauchy.sh')
sgbm = os.path.join(b, 'call_sgbm.sh')
mgm = os.path.join(b, 'mgm')
msmw = os.path.join(b, 'iip_stereo_correlation_multi_win2')
msmw2 = os.path.join(b, 'iip_stereo_correlation_multi_win2_newversion')
tvl1 = os.path.join(b, 'callTVL1.sh')

def compute_disparity_map(im1, im2, out_disp, out_mask, algo, disp_min, disp_max, extra_params=''):
    """
    Runs a block-matching binary on a pair of stereo-rectified images.

    Args:
        im1, im2: rectified stereo pair
        out_disp: path to the output diparity map
        out_mask: path to the output rejection mask
        algo: string used to indicate the desired binary. Currently it can be
            one among 'hirschmuller02', 'hirschmuller08', 'hirschmuller08_laplacian', 'hirschmuller08_cauchy', 'sgbm', and 'msmw'
        disp_min : smallest disparity to consider
        disp_max : biggest disparity to consider
        extra_params: optional string with algorithm-dependent parameters
    """
    # call the block_matching binary
    if (algo == 'hirschmuller02'):
        bm_binary = hirschmuller02
        common.run("%s %s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            out_mask, disp_min, disp_max, extra_params))
        # extra_params: LoG(0) regionRadius(3)
        #    LoG: Laplacian of Gaussian preprocess 1:enabled 0:disabled
        #    regionRadius: radius of the window

    if (algo == 'hirschmuller08'):
        bm_binary = hirschmuller08
        common.run("%s %s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            out_mask, disp_min, disp_max, extra_params))
        # extra_params: regionRadius(3) P1(default) P2(default) LRdiff(1)
        #    regionRadius: radius of the window
        #    P1, P2 : regularization parameters
        #    LRdiff: maximum difference between left and right disparity maps

    if (algo == 'hirschmuller08_laplacian'):
        bm_binary = hirschmuller08_laplacian
        common.run("%s %s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            out_mask, disp_min, disp_max, extra_params))
        # extra_params: regionRadius(3) P1(default) P2(default) LRdiff(1)
        #    regionRadius: radius of the window
        #    P1, P2 : regularization parameters
        #    LRdiff: maximum difference between left and right disparity maps

    if (algo == 'hirschmuller08_cauchy'):
        bm_binary = hirschmuller08_cauchy
        common.run("%s %s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            out_mask, disp_min, disp_max, extra_params))
        # extra_params: regionRadius(3) P1(default) P2(default) LRdiff(1)
        #    regionRadius: radius of the window
        #    P1, P2 : regularization parameters
        #    LRdiff: maximum difference between left and right disparity maps

    if (algo == 'sgbm'):
        bm_binary = sgbm
        out_cost = common.tmpfile('.tif')
        common.run("%s %s %s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            out_cost, out_mask, disp_min, disp_max, extra_params))

    if (algo == 'tvl1'):
        bm_binary = tvl1
        common.run("%s %s %s %s %s" %(bm_binary, im1, im2, out_disp,
            out_mask))

    if (algo == 'msmw'):
        bm_binary = msmw
        common.run("%s -i 1 -n 4 -p 4 -W 5 -x 9 -y 9 -r 1 -d 1 -t -1 -s 0 -b 0 -o 0.25 -f 0 -P 32 -m %d -M %d %s %s %s %s" %(bm_binary,
            disp_min, disp_max, im1, im2, out_disp, out_mask))

    if (algo == 'msmw2'):
        bm_binary = msmw2
        common.run("%s -i 1 -n 4 -p 4 -W 5 -x 9 -y 9 -r 1 -d 1 -t -1 -s 0 -b 0 -o -0.25 -f 0 -P 32 -D 0 -O 25 -c 0 -m %d -M %d %s %s %s %s" % (bm_binary,
            disp_min, disp_max, im1, im2, out_disp, out_mask))

    if (algo == 'mgm'):
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(cfg['omp_num_threads'])
        env['MEDIAN'] = '1'
        env['CENSUS_NCC_WIN'] = '5'
        env['TSGM'] = '3'
        common.run("%s -r %d -R %d -s vfit -t census -O 8 %s %s %s" % (mgm,
                                                                       disp_min,
                                                                       disp_max,
                                                                       im1, im2,
                                                                       out_disp),
                  env)

        # produce the mask: rejected pixels are marked with nan of inf in disp
        # map
        common.run('plambda %s "isfinite" -o %s' % (out_disp, out_mask))

    if (algo == 'micmac'):
        # add micmac binaries to the PATH environment variable
        s2p_dir = os.path.dirname(os.path.dirname(os.path.realpath(os.path.abspath(__file__))))
        micmac_bin = os.path.join(s2p_dir, '3rdparty', 'micmac', 'bin')
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + micmac_bin

        # prepare micmac xml params file
        micmac_params = os.path.join(s2p_dir, 'data', 'micmac_params.xml')
        work_dir = os.path.dirname(os.path.abspath(im1))
        common.run('cp %s %s' % (micmac_params, work_dir))

        # run MICMAC
        common.run("MICMAC %s" % os.path.join(work_dir, 'micmac_params.xml'))

        # copy output disp map
        micmac_disp = os.path.join(work_dir, 'MEC-EPI',
                                 'Px1_Num6_DeZoom1_LeChantier.tif')
        out_disp = os.path.join(work_dir, 'rectified_disp.tif')
        common.run('cp %s %s' % (micmac_disp, out_disp))

        # compute mask by rejecting the 10% of pixels with lowest correlation score
        micmac_cost = os.path.join(work_dir, 'MEC-EPI',
                                 'Correl_LeChantier_Num_5.tif')
        out_mask = os.path.join(work_dir, 'rectified_mask.png')
        common.run('plambda %s "x x%%q10 < 0 255 if" -o %s' % (micmac_cost,
                                                               out_mask))

    if (algo == 'asp'):
        # get disp_min/disp_max
        if hasattr(disp_min,'__len__') and hasattr(disp_min,'__len__'): # epipolar == False
            hmin, vmin = disp_min
            hmax, vmax = disp_max
        else:
            hmin, hmax = int(disp_min), int(disp_max)  # epipolar == True
            vmin, vmax = 0, 1

        corr_search = "--corr-search %d %d %d %d" %(min(hmin,hmax),min(vmin,vmax),max(hmin,hmax),max(vmin,vmax))

        # temporary files
        work_dir = os.path.dirname(os.path.abspath(im1))
        tmp_asp = os.path.join(work_dir,'asp')
        output_file_prefix = os.path.join(tmp_asp,"tile")

        # create tmp directory for asp
        if not os.path.exists(tmp_asp):
            os.makedirs(tmp_asp)

        # left and right images (and subsampled versions)
        left_im = output_file_prefix+"-L.tif"
        right_im = output_file_prefix+"-R.tif"
        left_im_sub = output_file_prefix+"-L_sub.tif"
        right_im_sub = output_file_prefix+"-R_sub.tif"

        # masks of left and right images (and subsampled versions)
        left_mask = output_file_prefix+"-lMask.tif"
        right_mask = output_file_prefix+"-rMask.tif"
        left_mask_sub = output_file_prefix+"-lMask_sub.tif"
        right_mask_sub = output_file_prefix+"-rMask_sub.tif"

        # filtered version
        filtered = output_file_prefix+"-F.tif"

        # output before formatting
        out_disp_vrt= output_file_prefix+".vrt"
        out_disp_b1 = output_file_prefix+"_b1.tif"
        out_disp_b2 = output_file_prefix+"_b2.tif"

        # prepare asp block matching
        common.run('cp %s %s' % (im1, left_im))
        common.run('cp %s %s' % (im2, right_im))

        # create subsampling tiles
        common.image_safe_zoom_fft(im1, 1.6, left_im_sub)
        common.image_safe_zoom_fft(im2, 1.6, right_im_sub)

        # create mask and subsampling mask
        for mask, im in zip([left_mask, right_mask, left_mask_sub, right_mask_sub],
                            [left_im, right_im, left_im_sub, right_im_sub]):
            (nc, nr) = common.image_size_gdal(im)
            common.run('plambda zero:%dx%d "x 255 +" -o %s' % (nc, nr, mask))

        # run asp
        io = " ".join([im1,im2,"-t isis",output_file_prefix])
        exec_lambda = lambda exe,args: " ".join([exe,io,args])
        exec_list = map(exec_lambda,
                        ["stereo_corr",
                         "stereo_rfne",
                         "stereo_fltr"],
                        ["--cost-mode=2 "+corr_search,
                         "--subpixel-kernel 11 11 --subpixel-pyramid-levels 5",""])

        for cmd in exec_list:
            common.run(cmd)

        # get disparities (row, col)
        common.run('gdal_translate -b 1 %s %s' % (filtered,out_disp_b1))
        common.run('gdal_translate -b 2 %s %s' % (filtered,out_disp_b2))

        common.run('gdalbuildvrt -separate %s %s %s' %(out_disp_vrt, out_disp_b1,out_disp_b2))
        common.run('gdal_translate -of GTiff %s %s' %(out_disp_vrt, out_disp))

        # out_mask
        goodPixelMap = output_file_prefix+'-GoodPixelMap.tif'
        common.run('plambda %s "x[1] 0 = 0 255 if" -o %s' % (goodPixelMap,out_mask))

        # delete temporary directory
        if os.path.exists(tmp_asp):
            shutil.rmtree(tmp_asp, ignore_errors=True)
