#!/usr/bin python
# -*- coding: utf-8 -*-

import numpy as np
import math
from scipy.stats import norm
from PA import GRT_coordinates, PA


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
    Function that calculates rotation of RA feed in image plane. That
    is rotation of z axis along x-axis (image plane == yz-plane & x is the
    electric axis of RA)

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
    z_old = np.asarray(z_old) * math.pi / 180.
    z_new = np.asarray(z_new) * math.pi / 180.
    x_new = np.asarray(x_new) * math.pi / 180.
    # Unit vectors of axes in coordinate system with unit vectors x=(\alpha=0,
    # \delta=0), y=(\alpha=+90deg, \delta=0), z along \delta=+90deg
    n_z_old = np.asarray([math.cos(z_old[0]) * math.cos(z_old[1]),
                          math.sin(z_old[0]) * math.cos(z_old[1]),
                          math.sin(z_old[1])])
    n_z_new = np.asarray([math.cos(z_new[0]) * math.cos(z_new[1]),
                          math.sin(z_new[0]) * math.cos(z_new[1]),
                          math.sin(z_new[1])])
    n_x_new = np.asarray([math.cos(x_new[0]) * math.cos(x_new[1]),
                          math.sin(x_new[0]) * math.cos(x_new[1]),
                          math.sin(x_new[1])])

    # Projecting ``n_z_old`` on ``yz_new`` plane
    n_z_old_2_yz_new = n_z_old - np.dot(n_z_old, n_x_new) * n_x_new
    # Normalization
    #n_z_old_2_yz_new = n_z_old_2_yz_new / np.linalg.norm(n_z_old_2_yz_new)
    n_z_old_2_yz_new /= np.linalg.norm(n_z_old_2_yz_new)
    # Find angle between two unit vectors
    theta = math.acos(np.dot(n_z_old_2_yz_new, n_z_new))

    return theta


def func1(x):
    """
    Bijection (0, 1) => R
    """
    return (2. * x - 1.) / (x - x ** 2.)


class Simulation(object):
    """
    Class that represents data (cross-to-parallel hands amplitude ratios for
    RA baseline with some ground radio telescope) generating process.

    :param mu_ampD_ra:
        Mean of RA D-term amplitude prior distribution.

    :param sigma_ampD_ra:
        Std of RA D-term amplitude prior distribution.

    :param mu_phD_ra:
        Mean of RA D-term phase prior distribution.

    :param sigma_phD_ra:
        std of RA D-term phase prior distribution.

    :param mu_ampD_grt:
        Mean of GRT D-term amplitude prior distribution.

    :param sigma_ampD_grt:
        Std of GRT D-term amplitude prior distribution.

    :param mu_phD_grt:
        Mean of GRT D-term phase prior distribution.

    :param sigma_phD_grt:
        Std of GRT D-term phase prior distribution.

    :param mean_r_G:
        Group mean of amplitude ratio for R&L gain.

    :param sigma_r_G:
        Group std of amplitude ratio for R&L gain
    """
    def __init__(self, mu_ampD_ra=0.075, sigma_ampD_ra=0.01, mu_phD_ra=0.,
                 sigma_phD_ra=0.2, mu_ampD_grt=0.02, sigma_ampD_grt=0.005,
                 mu_phD_grt=-math.pi, sigma_phD_grt=0.2, mean_r_G=0.,
                 sigma_r_G=0.0125, grt_coordinates=GRT_coordinates):
        self.mu_ampD_ra = mu_ampD_ra
        self.sigma_ampD_ra = sigma_ampD_ra
        self.mu_phD_ra = mu_phD_ra
        self.sigma_phD_ra = sigma_phD_ra
        self.mu_ampD_grt = mu_ampD_grt
        self.sigma_ampD_grt = sigma_ampD_grt
        self.mu_phD_grt = mu_phD_grt
        self.sigma_phD_grt = sigma_phD_grt
        self.mean_r_G = mean_r_G
        self.sigma_r_G = sigma_r_G
        self.grt_coordinates = grt_coordinates

    def generate(self, ampM, phiM, n, grt='EFF', source_coords=(110., 71.),
                 jds=[2456003.083333], z_old=(0., 100.), z_new=(-0.9, 7.0),
                 x_new=(10.5, 1.4)):
        """
        Generate value of cross-to-parallel hands ratios amplitude using
        parameters given in class constructor and amplitude&phase of source
        complex fractional LP.

        :param ampM:
            Amplitude of source complex fractional LP.

        :param phM:
            Phase of source complex fractional LP.

        :param n:
            Number of simulated values to return.

        :param grt:
            GRT on RA-Earth baseline.

        :param source_coords:
            Iterable of source coordinates (ra, dec).

        :param jds:
            Iterable of Julian Dates of observations.

        :param z_old:
            Orientation of RA z-axis while measuring D_{RA}.

        :param z_new:
            Orientation of RA z-axis at the moments of observations. Iterable of
            n tuples - (ra_z, dec_z).

        :param x_new:
            Orientation of RA x-axis at the moments of observations. Iterable of
            n tuples - (ra_x, dec_x).

        :return:
            List of n generated values of cross-hand ratios (amplitude).

        """
        if not grt in self.grt_coordinates:
            raise Exception("Uknown GRT!")
        else:
            grt_lat, grt_long = self.grt_coordinates[grt]

        if z_new is None:
            # generate randomly `z_new`
            raise NotImplementedError
        elif x_new is None:
            # generate randomly `x_new`
            raise NotImplementedError
        elif jds is None:
            # generate randomly `jds` in some range
            raise NotImplementedError
        else:
            assert(len(z_new) == n)
            assert(len(x_new) == n)
            assert (len(jds) == n)

        pa_grt = np.asarray(PA(jds, source_coords[0], source_coords[1], grt_lat,
                               grt_long))
        delta_pa_ra = np.asarray(get_z_rotation_in_image_plane(z_old, z_new,
                                                               x_new))
        n_pa_ra_diff_grt = delta_pa_ra - pa_grt

        n_ampD_ra = get_samples_from_truncated_normal(self.mu_ampD_ra,
                                                      self.sigma_ampD_ra, 0.,
                                                      1., n)
        n_phD_ra = get_samples_from_folded_normal(self.mu_phD_ra,
                                                  self.sigma_phD_ra,
                                                  -math.pi, math.pi, 1., n)
        n_ampD_grt = get_samples_from_truncated_normal(self.mu_ampD_grt,
                                                        self.sigma_ampD_grt,
                                                        0., 1.)
        n_phD_grt = get_samples_from_folded_normal(self.mu_phD_grt,
                                                   self.sigma_phD_grt,
                                                   -math.pi, math.pi, n)
        n_r = np.random.lognormal(mean=1., sigma=self.sigma_r_G, size=n)

        return n_r * abs(ampM * np.exp(2j * phiM) * np.exp(-2j * pa_grt) +
                         n_ampD_grt * np.exp(1j * n_phD_grt) +
                         n_ampD_ra * np.exp(-1j * n_phD_ra) *
                         np.exp(-2j * n_pa_ra_diff_grt))


def get_samples_from_truncated_normal(mean, sigma, a, b, size):
    """
    Function that returns ``size`` number of samples from truncated normal
    distribution with ``mu``, ``sigma`` and truncated at [``a``, ``b``].
    """

    kwargs = {'loc': mean, 'scale': sigma}
    us = np.random.uniform(norm.cdf(a, **kwargs), norm.cdf(b, **kwargs),
                           size=size)
    return norm.ppf(us, **kwargs)

def get_samples_from_folded_normal(mean, sigma, a, b, size):
    """
    Function that returns samples from normal distribution folded to [a,b]
    interval (e.g. [-pi, +pi])
    :param mean:
    :param sigma:
    :param a:
    :param b:
    :param size:
    :return:
    """
    if mean < a:
        dx = (a - mean) % (b - a)
        mean = b - dx
    elif mean > b:
        dx = (mean - b) % (b - a)
        mean = a + dx

    us = np.random.normal(mean, sigma, size)
    _du = (a - us) % (b - a)
    du_ = (us - b) % (b - a)
    _du_ixs = np.where(_du > 0)
    du_ixs = np.where(du_ > 0)
    us[_du_ixs] = b - _du
    us[du_ixs] = a + du_

    return us


if __name__ == '__main__':
    # 0716+714
    source_coords=(110., 71.)
