# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import sys
import numpy as np

import homography_cropper
import rpc_model
import rpc_utils
import estimation
import evaluation
import common
import sift
from config import cfg


def center_2d_points(pts):
    """
    Translates 2D points.

    The input points are translated such that the output points are centered at
    origin.

    Args:
        pts: 2D array of dimension Nx2 containing the coordinates of the input
            points, one point per line

    Returns:
        new_x, new_y, T: coordinates of the transformed points, together with
            the similarity (translation) matrix. This transformation takes the
            input points on the output points.
    """
    # centroid
    cx = np.mean(pts[:, 0])
    cy = np.mean(pts[:, 1])

    # shift origin to centroid
    new_x = pts[:, 0] - cx
    new_y = pts[:, 1] - cy

    # translation matrix
    T = np.eye(3)     #              1     0    -cx
    T[0, 2] = -cx     # matrix T  =  0     1    -cy
    T[1, 2] = -cy     #              0     0     1

    return np.vstack([new_x, new_y]).T, T


def filter_matches_epipolar_constraint(F, matches, thresh):
    """
    Discards matches that are not consistent with the epipolar constraint.

    Args:
        F: fundamental matrix
        matches: list of pairs of 2D points, stored as a Nx4 numpy array
        thresh: maximum accepted distance between a point and its matched
            epipolar line

    Returns:
        the list of matches that satisfy the constraint. It is a sub-list of
        the input list.
    """
    out = []
    for match in matches:
        x = np.array([match[0], match[1], 1])
        xx = np.array([match[2], match[3], 1])
        d1 = evaluation.distance_point_to_line(x, np.dot(F.T, xx))
        d2 = evaluation.distance_point_to_line(xx, np.dot(F, x))
        if max(d1, d2) < thresh:
            out.append(match)

    return np.array(out)


def register_horizontally(matches, H1, H2, do_shear=False, flag='center'):
    """
    Adjust rectifying homographies to modify the disparity range.

    Args:
        matches: list of pairs of 2D points, stored as a Nx4 numpy array
        H1, H2: two homographies, stored as numpy 3x3 matrices
        do_shear: boolean flag indicating wheter to minimize the shear on im2
            or not.
        flag: option needed to control how to modify the disparity range:
            'center': move the barycenter of disparities of matches to zero
            'positive': make all the disparities positive
            'negative': make all the disparities negative. Required for
                Hirshmuller stereo (java)

    Returns:
        H2: corrected homography H2

    The matches are provided in the original images coordinate system. By
    transforming these coordinates with the provided homographies, we obtain
    matches whose disparity is only along the x-axis. The second homography H2
    is corrected with a horizontal translation to obtain the desired property
    on the disparity range.
    """
    # transform the matches according to the homographies
    p1 = common.points_apply_homography(H1, matches[:, 0:2])
    x1 = p1[:, 0]
    y1 = p1[:, 1]
    p2 = common.points_apply_homography(H2, matches[:, 2:4])
    x2 = p2[:, 0]
    y2 = p2[:, 1]

    # for debug, print the vertical disparities. Should be zero.
    print "Residual vertical disparities: max, min, mean. Should be zero ------"
    print np.max(y2 - y1), np.min(y2 - y1), np.mean(y2 - y1)

    # shear correction
    # we search the (s, b) vector that minimises \sum (x1 - (x2+s*y2+b))^2
    # it is a least squares minimisation problem
    if do_shear:
        A = np.vstack((y2, y2*0+1)).T
        B = x1 - x2
        z = np.linalg.lstsq(A, B)[0]
        s = z[0]
        b = z[1]
        H2 = np.dot(np.array([[1, s, b], [0, 1, 0], [0, 0, 1]]), H2)
        x2 = x2 + s*y2 + b

    # compute the disparity offset according to selected option
    t = 0
    if (flag == 'center'):
        t = np.mean(x2 - x1)
    if (flag == 'positive'):
        t = np.min(x2 - x1)
    if (flag == 'negative'):
        t = np.max(x2 - x1)

    # correct H2 with a translation
    return np.dot(common.matrix_translation(-t, 0), H2)


def disparity_range(rpc1, rpc2, x, y, w, h, H1, H2, matches, A=None):
    """
    Compute the disparity range of a ROI from a list of point matches.

    The estimation is based on the extrapolation of the affine registration
    estimated from the matches. The extrapolation is done on the whole region of
    interest.

    Args:
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers defining the rectangular ROI in the first
            image.  (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        H1, H2: two rectifying homographies, stored as numpy 3x3 matrices
        matches: Nx4 numpy array containing a list of sift matches, in the full
            image coordinates frame
        A (optional): 3x3 numpy array containing the pointing error correction
            for im2. This matrix is usually estimated with the pointing_accuracy
            module.

    Returns:
        disp_min, disp_max: horizontal disparity range
    """
    # transform the matches according to the homographies
    p1 = common.points_apply_homography(H1, matches[:, :2])
    x1 = p1[:, 0]
    p2 = common.points_apply_homography(H2, matches[:, 2:])
    x2 = p2[:, 0]
    y2 = p2[:, 1]

    # estimate an affine transformation (tilt, shear and bias)
    # that maps p1 on p2
    z = np.linalg.lstsq(np.vstack((x2, y2, y2*0+1)).T, x1)[0]
    t, s, b = z[0:3]

    # corners of ROI
    xx2 = np.array([0, w, 0, w])
    yy2 = np.array([0, 0, h, h])

    # compute the max and min disparity values according to the estimated
    # model. The min/max disp values are necessarily obtained at the ROI
    # corners
    roi_disparities_by_the_affine_model = (xx2*t + yy2*s + b) - xx2
    max_roi = np.max(roi_disparities_by_the_affine_model)
    min_roi = np.min(roi_disparities_by_the_affine_model)

    # min/max disparities according to the keypoints
    max_kpt = np.max(x2 - x1)
    min_kpt = np.min(x2 - x1)

    # compute the range with the extracted min and max disparities
    disp_m = np.floor(min(min_roi, min_kpt))
    disp_M = np.ceil(max(max_roi, max_kpt))

    # add a security margin to the disp range
    disp_m *= (1 - np.sign(disp_m) * cfg['disp_range_extra_margin'])
    disp_M *= (1 + np.sign(disp_M) * cfg['disp_range_extra_margin'])

    print "SIFT disparity range: [%f, %f]" % (disp_m, disp_M)

    # expand disparity range with srtm according to cfg params
    if (cfg['disp_range_method'] == "srtm") or (matches is None) or (len(matches) < 2):
        disp_m, disp_M = rpc_utils.srtm_disp_range_estimation(rpc1, rpc2, x, y,
                                                              w, h, H1, H2, A,
                                                              cfg['disp_range_srtm_high_margin'],
                                                              cfg['disp_range_srtm_low_margin'])
        print "SRTM disparity range: [%f, %f]" % (disp_m, disp_M)
    if cfg['disp_range_method'] == "wider_sift_srtm" and (matches is not None) and (len(matches) >= 2):
        d_m, d_M = rpc_utils.srtm_disp_range_estimation(rpc1, rpc2, x, y, w, h,
                                                        H1, H2, A,
                                                        cfg['disp_range_srtm_high_margin'],
                                                        cfg['disp_range_srtm_low_margin'])
        print "SRTM disparity range: [%f, %f]" % (d_m, d_M)
        disp_m = min(disp_m, d_m)
        disp_M = max(disp_M, d_M)

    print "Final disparity range: [%f, %f]" % (disp_m, disp_M)
    return disp_m, disp_M


def rectification_homographies(matches, x, y, w, h):
    """
    Computes rectifying homographies from point matches for a given ROI.

    The affine fundamental matrix F is estimated with the gold-standard
    algorithm, then two rectifying similarities (rotation, zoom, translation)
    are computed directly from F.

    Args:
        matches: numpy array of shape (n, 4) containing a list of 2D point
            correspondences between the two images.
        x, y, w, h: four integers definig the rectangular ROI in the first
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.

    Returns:
        S1, S2, F: three numpy arrays of shape (3, 3) representing the
        two rectifying similarities to be applied to the two images and the
        corresponding affine fundamental matrix.
    """
    # estimate the affine fundamental matrix with the Gold standard algorithm
    F = estimation.affine_fundamental_matrix(matches)

    # compute rectifying similarities
    S1, S2 = estimation.rectifying_similarities_from_affine_fundamental_matrix(F, True)

    if cfg['debug']:
        y1 = common.points_apply_homography(S1, matches[:, :2])[:, 1]
        y2 = common.points_apply_homography(S2, matches[:, 2:])[:, 1]
        err = np.abs(y1 - y2)
        print "max, min, mean rectification error on point matches: ",
        print np.max(err), np.min(err), np.mean(err)

    # pull back top-left corner of the ROI to the origin
    pts = common.points_apply_homography(S1, [[x, y], [x+w, y], [x+w, y+h], [x, y+h]])
    x0, y0 = common.bounding_box2D(pts)[:2]
    T = common.matrix_translation(-x0, -y0)
    return np.dot(T, S1), np.dot(T, S2), F


def rectify_pair(im1, im2, rpc1, rpc2, x, y, w, h, out1, out2, A=None,
                 sift_matches=None, method='rpc'):
    """
    Rectify a ROI in a pair of images.

    Args:
        im1, im2: paths to two image files
        rpc1, rpc2: paths to the two xml files containing RPC data
        x, y, w, h: four integers defining the rectangular ROI in the first
            image.  (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        out1, out2: paths to the output rectified crops
        A (optional): 3x3 numpy array containing the pointing error correction
            for im2. This matrix is usually estimated with the pointing_accuracy
            module.
        sift_matches (optional): Nx4 numpy array containing a list of sift
            matches, in the full image coordinates frame
        method (default: 'rpc'): option to decide wether to use rpc of sift
            matches for the fundamental matrix estimation.

        This function uses the parameter subsampling_factor from the
        config module. If the factor z > 1 then the output images will
        be subsampled by a factor z. The output matrices H1, H2, and the
        ranges are also updated accordingly:
        Hi = Z * Hi with Z = diag(1/z, 1/z, 1) and
        disp_min = disp_min / z  (resp _max)

    Returns:
        H1, H2: Two 3x3 matrices representing the rectifying homographies that
        have been applied to the two original (large) images.
        disp_min, disp_max: horizontal disparity range
    """
    # read RPC data
    rpc1 = rpc_model.RPCModel(rpc1)
    rpc2 = rpc_model.RPCModel(rpc2)

    # compute real or virtual matches
    if method == 'rpc':
        # find virtual matches from RPC camera models
        matches = rpc_utils.matches_from_rpc(rpc1, rpc2, x, y, w, h,
                                             cfg['n_gcp_per_axis'])

        # correct second image coordinates with the pointing correction matrix
        if A is not None:
            matches[:, 2:] = common.points_apply_homography(np.linalg.inv(A),
                                                            matches[:, 2:])
    else:
        matches = sift_matches

    # compute rectifying homographies
    H1, H2, F = rectification_homographies(matches, x, y, w, h)

    # compose H2 with a horizontal translation to center disp range around 0
    if sift_matches is not None:
        sift_matches = filter_matches_epipolar_constraint(F, sift_matches,
                                                          cfg['epipolar_thresh'])
        if len(sift_matches) < 10:
            print 'WARNING: no registration with less than 10 matches'
        else:
            H2 = register_horizontally(sift_matches, H1, H2)

    # compute disparity range
    disp_m, disp_M = disparity_range(rpc1, rpc2, x, y, w, h, H1, H2,
                                     sift_matches, A)

    # compute rectifying homographies for non-epipolar mode (rectify the secondary tile only)
    if cfg['epipolar_rectification'] == False:
        H1_inv = np.linalg.inv(H1)
        H1 = np.dot(H1_inv,H1)
        H1 = np.dot(common.matrix_translation(-x, -y),H1)
        H2 = np.dot(H1_inv,H2)
        H2 = np.dot(common.matrix_translation(-x, -y),H2)

    # compute output images size
    roi = [[x, y], [x+w, y], [x+w, y+h], [x, y+h]]
    pts1 = common.points_apply_homography(H1, roi)
    x0, y0, w0, h0 = common.bounding_box2D(pts1)
    # check that the first homography maps the ROI in the positive quadrant
    np.testing.assert_allclose(np.round([x0, y0]), 0, atol=.01)

    # apply homographies and do the crops TODO XXX FIXME cleanup here
    #homography_cropper.crop_and_apply_homography(out1, im1, H1, w0, h0, cfg['subsampling_factor'], True)
    #homography_cropper.crop_and_apply_homography(out2, im2, H2, w0, h0, cfg['subsampling_factor'], True)
    common.image_apply_homography(out1, im1, H1, w0, h0)
    common.image_apply_homography(out2, im2, H2, w0, h0)

    #  if subsampling_factor'] the homographies are altered to reflect the zoom
    if cfg['subsampling_factor'] != 1:
        Z = np.eye(3)
        Z[0, 0] = Z[1, 1] = 1.0 / cfg['subsampling_factor']

        H1 = np.dot(Z, H1)
        H2 = np.dot(Z, H2)
        disp_m = np.floor(disp_m / cfg['subsampling_factor'])
        disp_M = np.ceil(disp_M / cfg['subsampling_factor'])

    if cfg['epipolar_rectification'] == False:
        if cfg['disp_min'] is not None:
            disp_m = cfg['disp_min']
        if cfg['disp_max'] is not None:
            disp_M = cfg['disp_max']

        pts_in = [[0, 0], [disp_m, 0], [disp_M, 0]]
        pts_out = common.points_apply_homography(H1_inv,
                                                 pts_in)
        disp_m = pts_out[1,:] - pts_out[0,:]
        disp_M = pts_out[2,:] - pts_out[0,:]

    return H1, H2, disp_m, disp_M
