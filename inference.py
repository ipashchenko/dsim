#!/usr/bin python
# -*- coding: utf-8 -*-

# Modified Bessel function of order 0
from scipy.special import i0
from scipy.misc import logsumexp
import math


class LnPost(object):
    """
    Class that represents posterior pdf of parameters given data.
    """
    def __init__(self, obs, lnpr, prargs=[], prkwargs={}, liargs=[],
                 likwargs={}):
        self._lnpr = lnpr
        self._lnlike = LnLike(obs)
        self.prargs = prargs
        self.prkwargs = prkwargs
        self.liargs = liargs
        self.likwargs = likwargs

    def lnpr(self, p):
        return self._lnpr(p, *self.prargs, **self.prkwargs)

    def lnlik(self, p):
        return self._lnlik(p, *self.liargs, **self.likwargs)

    def __call__(self, p):
        return self.lnlike(p) + self.lnpr(p)


class LnLike(object):
    """
    Class that represents likelihood function.
    """
    def __init__(self, obs, model, *args, **kwargs):
        """
        :param obs:
            Observed data. Nested iterable where outer iter corresponds to
            groups and inner iter corresponds to current group's observations.
            So len(obs) = #groups, len(obs[i]) = #obs in i-th group.
        :param model:
            Callable that given parameter vector `p` and some other arguments
            returns amplitude of cross-to-parallel hands ratio.
        :param args:
            Positional arguments to ``model`` callable.
        :param kwargs:
            Keyword arguments to ``model`` callable.
        """
        self.obs = obs
        self.args = args
        self.kwargs = kwargs

    def __call__(self, p, *args, **kwargs):
        """
        Returns ln of likelihood function of data ``obs`` given parameters
        ``p``.
        """
        pass


class LnPr(object):
    """
    Class that represents prior pdf of parameters.
    """
    def __init__(self):
        pass

    def __call__(self, p, *args, **kwargs):
        pass


def model(p):
    """
    Model of cross-to-parallel hands ratio.

    :param p:
        Parameter vector.
    :return:
        Amplitude of cross-to-parallel hands ratio for given parameters ``p``.
    """
    pass


def lnpdf_von_mises(x, mu, k):
    """
    Returns pdf of Von Mises distribution with parameters ``mu`` & ``k``
    (~precision of noraml distribution).
    """
    return k * math.cos(x - mu) - math.log(2. * math.pi * i0(k))
    

def lnpdf_wrapped_normal(x, mu, tau, T=None):
    """
    Returns pdf of Wrapped Normal distribution with period ``T`` & interval
    [0, T] or [-T/2, T/2].
    """
    k = np.arange(-5,5)
    return -0.5 * math.log(2. * math.pi) + 0.5 * math.log(tau) + logsumexp(-tau * (x - mu + T * k)
