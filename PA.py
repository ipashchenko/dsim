#!/usr/bin/env python
#-*- coding: utf-8 -*-
import math

# Dictionary of GRT latitudes and longitudes
GRT_coordinates = {'AR': (18.344167, -66.752778),
                   'GBT': (38.433056, -79.839722),
                   'EFF': (50.524722, 6.882778),
                   'WSRT': (52.914722, 6.603333)}

def JD_to_LST(JDs, llong):
    """Returns LST [frac. of day] using:
        JD - [iterable] - values of Julian data
        llong - local longitude [+/-degrees].
    """

    LSTs = list()

    llong_h = llong / 15.

    for JD in list(JDs):

        D = JD - 2451545.0
        GMST = 18.697374558 + 24.06570982441908 * D
        GMST = GMST % 24

        if GMST < 0:
            GMST += 24.
        elif GMST >= 24.0:
            GMST -= 24.0

        LST = GMST + llong_h
        #convert to fraction of day
        LST = LST / 24.

        LSTs.append(LST)

    return LSTs


def LST_to_HA(LSTs, RA):
    """
    RA [degrees] - right ascenction
    LST [fr. of days] - local sidireal time
    Returns Hour Angle [rads]
    """

    HAs = list()

    RA_rad = RA * math.pi / 180.

    for LST in list(LSTs):
        HA = 2. * math.pi * LST - RA_rad
        HAs.append(HA)

    return HAs


def PA(JDs, ra, dec, latitude, longitude):
    """Function returns parallactic angles of source with given coordinates
    observed at given Julian Dates at given geographic location.
    'east' or 'west' of Greenwitch => +/- sign for longitude.

        :param jds:
            Iterable of Julian Dates.
        :param ra:
            Right assention of source. Dergees.
        :param dec:
            Declination of source. Dergees.
        :param longitude:
            Geographical longitude of observatory. If west of Greenwitch => <
            0. Degrees.
        :param latitude:
            Geographical latitude of observatory. Degrees.
    """

    PAs = list()

    LSTs = JD_to_LST(JDs, longitude)

    HAs = LST_to_HA(LSTs, ra)

    for HA in list(HAs):

        PA = math.atan(math.sin(HA) / (math.tan(latitude) * math.cos(dec) -\
            math.sin(dec) * math.cos(HA)))
        PAs.append(PA)

    return PAs
