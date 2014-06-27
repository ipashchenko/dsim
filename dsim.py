#!/usr/bin python
# -*- coding: utf-8 -*-

import numpy as np
from math import sin, cos, acos, pi, sqrt


# R_{i} = g_{i} * ||M_{i}*exp(2i*\chi^{M}_{i})*exp(-2i*\fi_{i}^{GRT}) +
#                   D_{i}^{GRT}*exp(i*\fi_{i}^{D_GRT}) +
#                   D_{i}^{RA}*exp(-i*\fi_{i}^{D_RA})*exp(2i*(\fi_{i}^{GRT} -
#                                                           \fi_{i}^{RA}))
#                 ||
#
#   M_{i} - amplitude of source frac. polarization for current source
#   D_{i}^{RA} - RA D-term amplitude, generated from some narrow distribution.
#   D_{i}^{GRT} - D-term amplitude for current GRT, generated from some unknown
# distribution.
#   g_{i} - amplitude ratio of gains for RA & ground telescope (we got some
# systematic difference for RA R & L for some band!)
#  \fi_{i}^{D_RA} - RA D-term phase, generated from some narrow distribution.
#  \fi_{i}^{D_GRT} - D-term phase for current GRT, generated from some unknown
# distribution.
#  \fi_{i}^{GRT/RA} - parallactic angle of RA/GRT telescope (calculated).
# Parallactic angle of RA is constant (zero), but one needs to correct for
# rotation of feed (rotation ox z-axis around x-axis).


def get_z_rotation_in_image_plane(z_old, z_new, x_new):
    """
    Function that calculates rotation of Radioastron feed in image plane.

    :param z_old:
        Tuple (right ascenction, declination) of z-axis before rotation [deg].
    :param z_new:
        Tuple (right ascenction, declination) of z-axis after rotation [deg].
    :param x_new:
        Tuple (right ascenction, declination) of x-axis after rotation [deg].
    :return:
        The angle between new z-axis and projection of old z-axis to the new
        y,z-plane [rad].
    """
    # degrees -> rads
    z_old = np.asarray(z_old) * pi / 180.
    z_new = np.asarray(z_new) * pi / 180.
    x_new = np.asarray(x_new) * pi / 180.
    # Unit vectors of axes in coordinate system with unit vectors x=(\alpha=0,
    # \delta=0), y=(\alpha=+90deg, \delta=0), z along \delta=+90deg
    n_z_old = np.asarray([cos(z_old[0]) * cos(z_old[1]),
                          sin(z_old[0]) * cos(z_old[1]),
                          sin(z_old[1])])
    n_z_new = np.asarray([cos(z_new[0]) * cos(z_new[1]),
                          sin(z_new[0]) * cos(z_new[1]),
                          sin(z_new[1])])
    n_x_new = np.asarray([cos(x_new[0]) * cos(x_new[1]),
                          sin(x_new[0]) * cos(x_new[1]),
                          sin(x_new[1])])

    # Projecting ``n_z_old`` on ``yz_new`` plane
    n_z_old_2_yz_new = n_z_old - np.dot(n_z_old, n_x_new) * n_x_new
    # Normalization
    n_z_old_2_yz_new = n_z_old_2_yz_new / np.linalg.norm(n_z_old_2_yz_new)
    # Find angle between two unit vectors
    theta = acos(np.dot(n_z_old_2_yz_new, n_z_new))

    return theta


def func1(x):
    """
    Bijection (0, 1) => R
    """
    return (2. * x - 1.) / (x - x ** 2.)


class Simulation(object):
    """
    Class that represents data (cross-to-parallel hands amplitude ratios)
    generating process.
    """
    def __init__(self, mean_ampD_ra, sigma_ampD_ra, mean_phD_ra, sigma_phD_ra,
                 mean_ampD_grt, sigma_ampD_grt):
        pass

    def generate(self, m, n=10000):
        """
        Generate ``n`` values of cross-to-parallel hands ratios amplitudes
        using parameters of Radioastron given in class constructor and
        ``func1``-transform of the amplitude of source frac. LP.
        """

        M = m - 2. + math.sqrt(4. + m ** 2.)
        results = list()
        for i in range(n):
            result =
            results.append()

