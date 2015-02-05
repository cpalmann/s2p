# Copyright (C) 2013, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2013, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2013, Julien Michel <julien.michel@cnes.fr>

import numpy as np
import common
from config import cfg

def image_apply_pleiades_unsharpening_filter(im):
    """
    Returns the image convolved by the unsharpening MTF idata_0009_MTF_89x89.tif
    This filter specifically undoes the sharpening applied to the sensor perfect
    Pleiades images
    """
    unsharpening_mtf_small = common.image_pleiades_unsharpening_mtf()
    tmp_mtf_large = common.image_zeropadding_from_image_with_target_size(
            unsharpening_mtf_small, im)
    return common.image_fftconvolve(im, tmp_mtf_large)


def crop_and_apply_homography(im_out, im_in, H, w, h, subsampling_factor=1,
        convert_to_gray=False):
    """
    Warps a piece of a Pleiades (panchro or ms) image with a homography.

    Args:
        im_out: path to the output image
        im_in: path to the input (tif) full Pleiades image
        H: numpy array containing the 3x3 homography matrix
        w, h: size of the output image
        subsampling_factor (optional, default=1): when set to z>1,
            will result in the application of the homography Z*H where Z =
            diag(1/z, 1/z, 1), so the output will be zoomed out by a factor z.
            The output image will be (w/z, h/z)
        convert_to_gray (optional, default False): it set to True, and if the
            input image has 4 channels, it is converted to gray before applying
            zoom and homographies.

    Returns:
        nothing

    The homography has to be used as: coord_out = H coord_in. The produced
    output image corresponds to coord_out in [0, w] x [0, h]. The warp is made
    by Pascal Monasse's binary named 'homography'.
    """

    # crop a piece of the big input image, to which the homography will be
    # applied
    # warning: as the crop uses integer coordinates, be careful to round off
    # (x0, y0) before modifying the homograpy. You want the crop and the
    # translation representing it do exactly the same thing.
    pts = [[0, 0], [w, 0], [w, h], [0, h]]
    inv_H_pts = common.points_apply_homography(np.linalg.inv(H), pts)
    x0, y0, w0, h0 = common.bounding_box2D(inv_H_pts)
    x0, y0 = np.floor([x0, y0])
    w0, h0 = np.ceil([w0, h0])
    crop_fullres = common.image_crop_LARGE(im_in, x0, y0, w0, h0)

    # This filter is needed (for panchro images) because the original PLEAIDES
    # SENSOR PERFECT images are aliased
    if (common.image_pix_dim(crop_fullres) == 1 and subsampling_factor == 1 and
            cfg['use_pleiades_unsharpening']):
        tmp = image_apply_pleiades_unsharpening_filter(crop_fullres)
        common.run('rm %s' % crop_fullres)
        crop_fullres = tmp

    # convert to gray
    if common.image_pix_dim(crop_fullres) == 4:
        if convert_to_gray:
            crop_fullres = common.pansharpened_to_panchro(crop_fullres)

    # compensate the homography with the translation induced by the preliminary
    # crop, then apply the homography and crop.
    H = np.dot(H, common.matrix_translation(x0, y0))

    # Since the objective is to compute a zoomed out homographic transformation
    # of the input image, to save computations we zoom out the image before
    # applying the homography. If Z is the matrix representing the zoom out and
    # H the homography matrix, this trick consists in applying Z*H*Z^{-1} to
    # the zoomed image Z*Im instead of applying Z*H to the original image Im.
    if subsampling_factor == 1:
        common.image_apply_homography(im_out, crop_fullres, H, w, h)
        return

    else:
        assert(subsampling_factor >= 1)

        # H becomes Z*H*Z^{-1}
        Z = np.eye(3);
        Z[0,0] = Z[1,1] = 1 / float(subsampling_factor)
        H = np.dot(Z, H)
        H = np.dot(H, np.linalg.inv(Z))

        # w, and h are updated accordingly
        w = int(w / subsampling_factor)
        h = int(h / subsampling_factor)

        # the DCT zoom is NOT SAFE when the input image size is not a multiple
        # of the zoom factor
        tmpw, tmph = common.image_size(crop_fullres)
        tmpw, tmph = int(tmpw / subsampling_factor), int(tmph / subsampling_factor)
        crop_fullres_safe = common.image_crop_TIFF(crop_fullres, 0, 0, tmpw *
                subsampling_factor, tmph * subsampling_factor)
        common.run('rm %s' % crop_fullres)

        # zoom out the input image (crop_fullres)
        crop_zoom_out = common.image_safe_zoom_fft(crop_fullres_safe,
                subsampling_factor)
        common.run('rm %s' % crop_fullres_safe)

        # apply the homography to the zoomed out crop
        common.image_apply_homography(im_out, crop_zoom_out, H, w, h)
        return

        # TODO DCT IS STILL NOT SAFE THE position 0,0 is translated half pixel !
