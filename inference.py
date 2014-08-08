#!/usr/bin python
# -*- coding: utf-8 -*-


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
        ``p``, where vector of parameters:
        p = [mu_amp_M_G, tau_amp_M_G, mu_ph_M_G, tau_ph_M_G, mu_amp_M_i,
        ta_amp_M_i, mu_amp_D_grt_G, tau_amp_D_grt_G, mu_ph_D_grt_G,
        tau_ph_D_grt_G, mu_amp_D_grt_i, tau_amp_D_grt_i, mu_ph_D_grt_i,
        tau_ph_D_grt_i, mu_r_G, tau_r_G, mu_r_i, tau_r_i, mu_amp_D_ra,
        tau_amp_D_ra, mu_ph_D_ra, tau_ph_D_ra]

        ``mu_amp_M_G`` - group mean of amplitude of frac. pol.
        ``tau_amp_M_G`` - group precision of amplitude of frac. pol.
        ``mu_ph_M_G`` - group mean of phase of frac. pol.
        ``tau_ph_M_G`` - group precision of phase of frac. pol.
        ``mu_amp_M_i`` - mean of amplitude of frac. pol. for i-th group.
        ``ta_amp_M_i`` - precision of amplitude of frac. pol. for i-th group.
        ``mu_amp_D_grt_G`` - group mean of amplitude of GRT D-term.
        ``tau_amp_D_grt_G`` - group precision of amplitude of GRT D-term.
        ``mu_ph_D_grt_G`` - group mean of phase of GRT D-term.
        ``tau_ph_D_grt_G`` - group precision of phase of GRT D-term.
        ``mu_amp_D_grt_i`` - mean of amplitude of GRT D-term for i-th
            experiment.
        ``tau_amp_D_grt_i`` - precision of amplitude of GRT D-term for i-th
            group.
        ``mu_ph_D_grt_i`` - mean of phase of GRT D-term for i-th group.
        ``tau_ph_D_grt_i`` - precision of phase of GRT D-term for i-th group.
        ``mu_r_G`` - group mean of amplitude of L,R gains ratio.
        ``tau_r_G`` - group precision of amplitude L,R gains ratio.
        ``mu_r_i`` - mean of amplitude of L,R gains ration for i-th experiment.
        ``tau_r_i`` - precision of amplitude of L,R gains ration for i-th
            experiment.
        ``mu_amp_D_ra`` - mean of amplitude of RA D-term.
        ``tau_amp_D_ra`` - precision of amplitude of RA D-term.
        ``mu_ph_D_ra`` - mean of phase of RA D-term.
        ``tau_ph_D_ra`` - precision of phase of RA D-term.
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
