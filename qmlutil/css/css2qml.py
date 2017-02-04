# -*- coding: utf-8 -*-
#
# Copyright 2016 University of Nevada, Reno
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""
qmlutil.css.css2qml

    Mark C. Williams (2015)
    Nevada Seismological Laboratory

    Converter classes to map CSS3.0 to QuakeML schema


Classes
=======
CSSToQMLConverter : methods to convert CSS to QuakeML schema
---------------------
USE:
>>> c = CSSToQuakeMLConverter(
...    agency='NN',
...    rid_factory=ResourceURIGenerator('quakeml', 'org.nvseismolab'),
...    utc_factory=timestamp2isostr
...    )


NOTES
-----
This is a shot at implementing a schema converter in pure python. The goal is
to enable dict->dict conversion between formats. There are some caveats, based
on the schemas:

1) CSS3.0 input dicts can contain joined, or "view" records for completeness.
2) CSS3.0 view keys are namespaced by table for uniqueness (only if necessary)
3) This QML schema uses 'xmltodict' style (for keys, dicts, lists). This means
   it's maybe incompatible with JSON-LD. Hopefully falls under YAGNI.

"""
import math
import datetime
import uuid

from qmlutil import Dict, Root, anss_params
# from qmlutil import ResourceURIGenerator

# Default weight to use based on timedef
TIMEDEF_WEIGHT = dict(d=1.0, n=0.0)

# Default CSS3.0 etypes to QML event types
# TODO: update to support IASPEI 2-letter coding of type & certainty
ETYPE_MAP = {
    'qb': "quarry blast",
    'eq': "earthquake",
    'me': "rock burst",
    'ex': "explosion",
    'o': "other",
    'l': "earthquake",
    'r': "earthquake",
    't': "earthquake",
    'f': "earthquake",
}


def _dt(timestamp):
    """Returns the UTC dateTime"""
    try:
        return datetime.datetime.utcfromtimestamp(timestamp)
    except ValueError as ex:  # analysis:ignore
        return None


def _ts(datetime_input):
    """
    Return timestamp from datetime object
    """
    datetime_reference = datetime.datetime(1970, 1, 1, 0, 0, 0)
    return (datetime_input - datetime_reference).total_seconds()


def _str(item):
    """Return a string no matter what"""
    if item is not None:
        return str(item)
    else:
        return ''


def _dict(*args, **kwargs):
    """
    Return a dict only if at least one value is not None
    """
    dict_ = Dict(*args, **kwargs)
    if dict_.values() == [None] * len(dict_):
        return None
    return dict_


def _quan(*args, **kwargs):
    """
    Return a dict only if the value for key "value" is not None
    """
    dict_ = Dict(*args, **kwargs)
    if dict_.get('value') is None:
        return None
    return dict_


def _km2m(dist):
    """Convert from km to m only if dist is not None"""
    if dist is not None:
        return dist * 1000.
    else:
        return None


def _m2deg_lat(dist):
    return dist / 110600.


def _m2deg_lon(dist, lat=0.):
    average_earth_radius = 6367449.
    return (dist / (math.pi / 180.) / average_earth_radius /
            math.cos(math.radians(lat)))


def _eval_ellipse(semi_major, semi_minor, angle):
    return semi_major*semi_minor/(math.sqrt(
        (semi_minor*math.cos(math.radians(angle)))**2 +
        (semi_major*math.sin(math.radians(angle)))**2))


def _get_ne_on_ellipse(semi_major, semi_minor, strike):
    """
    Return the solution for points N and E on an ellipse

    A : float of semi major axis
    B : float of semi minor axis
    strike : angle of major axis from North

    Returns
    -------
    n, e : floats of ellipse solution at north and east
    """
    north = _eval_ellipse(semi_major, semi_minor, strike)
    east = _eval_ellipse(semi_major, semi_minor, strike - 90)
    return north, east


def extract_etype(origin):
    """Return a CSS3.0 etype flag stored in an origin"""
    if 'css:etype' in origin:
        return origin['css:etype']
    else:
        if origin.get('comment'):
            for comm in origin.get('comment'):
                # This will throw error if comment is a dict not list:
                # TODO: check?
                if comm.get('@id', '').endswith('etype'):
                    return comm.get('text')


class CSSToQMLConverter(Root):
    """
    Converter to QuakeML schema from CSS3.0 schema

    Attributes
    ----------
    agency  : str to use when agency can't be inferred from author
    default_network  : str to use when waveform identifier is missing network
    automatic_authors : list of authors to mark as "auto"
    rid_factory : ResourceUIDGenerator function which returns ID's
    utc_factory : function that converts a float timestamp

    Methods
    -------
    get_event_type : static class method to convert CSS origin type flag

    """
    etype_map = dict(ETYPE_MAP)
    automatic_authors = []  # list of authors to mark as automatic

    def __init__(self, *args, **kwargs):
        """
        Set event
        """
        # self.event = Dict()

        # Allow setting of map at class level by noclobber update
        if 'etype_map' in kwargs:
            etype_map = kwargs.pop('etype_map')
            _etypemap = dict(self.etype_map)
            _etypemap.update(etype_map)
            self.etype_map = _etypemap

        # TODO: inherit from Root class and call super for this stuff
        super(CSSToQMLConverter, self).__init__(*args, **kwargs)

    def get_event_type(self, etype):
        """
        Map a CSS3.0 etype origin flag to a QuakeML event type

        Inputs
        ------
        etype : str of a valid etype
        """
        return self.etype_map.get(etype, "not reported")

    def origin_event_type(self, origin):
        """
        Return a proper event_type from a CSS3.0 etype flag stored in an origin
        """
        etype = extract_etype(origin)
        return self.get_event_type(etype)

    def get_event_status(self, posted_author):
        """
        Return mode and status based on author
        """
        for auto_author in self.automatic_authors:
            if not posted_author or auto_author in posted_author:
                mode = "automatic"
                status = "preliminary"
                return mode, status
        mode = "manual"
        status = "reviewed"
        return mode, status

    def get_method_model(self, algorithm):
        """
        Return method and model based on algorithm
        """
        if algorithm is None:
            return None, None, None

        if ':' in algorithm:
            method, model = algorithm.split(':', 1)
        else:
            method, model = 'grassoc', algorithm

        if '(' in model and ')' in model:
            model, quality = model.split('(')
            quality = quality.split(')')[0]
        else:
            quality = ''

        method_rid = "{0}/{1}".format('method', method)

        return method_rid, model, quality

    def split_auth(self, auth,
                   recognized_magnitudes=('mb', 'ml', 'mn', 'mw', 'mwb',
                                          'mwp', 'mwr', 'mww', ' m')):
        '''
        Splits an Antelope 'auth' string into its component parts.
        Anything prior to a colon is considered the agency.
        An attempt is made to strip trailing magnitude types.

        :returns: agency, author, magnitude_types, method, info
        '''
        if auth is None:
            return None, None, None, None, None

        # usage for picks in arrival table
        if 'dbp' in auth:
            method, author, info = auth.split(':', 2)
            return None, author, None, method, info

        # usage in amplitudes in arrival table, netmag and stamag table
        if 'dbevproc' in auth:
            if ':' in auth:
                method, info = auth.split(':', 1)
            else:
                method, info = auth, None
            return None, None, None, method, info

        # usage in origin table and magnitudes from foreign agencies
        if ':' in auth:
            agency, author = auth.split(':', 1)
        else:
            agency, author = '', auth
        agency = agency.strip()
        author = author.strip()

        # deal with silly concatenation of magnitudes onto auth field
        magnitude_types = []
        while any([author[-len(recognized_magnitude):].lower() ==
                   recognized_magnitude.lower()
                   for recognized_magnitude in recognized_magnitudes]):
            for recognized_magnitude in recognized_magnitudes:
                magnitude_type = author[-len(recognized_magnitude):]
                if magnitude_type.lower() == recognized_magnitude.lower():
                    magnitude_types += [magnitude_type]
                    author = author[:len(author) - len(recognized_magnitude)]
                    author = author.strip()
        magnitude_types = []
        return agency, author, magnitude_types, None, None

    def map_origin2origin(self, db):
        """
        Return a dict of QuakeML origin from a dict of CSS key/values

        Inputs
        ======
        db : dict of key/values of CSS fields related to the origin (see Join)

        Returns
        =======
        dict of key/values of QuakeML fields of "origin"

        Notes
        =====
        Any object that supports the dict 'get' method can be passed as
        input, e.g. OrderedDict, custom classes, etc.

        Join
        ----
        origin <- origerr [orid] (outer)
        """
        css_etype = _str(db.get('etype'))
        agency, author, _, _, _ = self.split_auth(_str(db.get('auth')))
        method_rid, model, _ = self.get_method_model(_str(db.get('algorithm')))
        model_rid = "{0}/{1}".format('vmodel',
                                     db.get('vmodel') or model or uuid.uuid4())

        mode, status = self.get_event_status(author)
        agency if agency is not '' else self.agency
        origin_rid = "{0}/{1}".format('origin',
                                      db.get('orid') or uuid.uuid4())

        # -- Solution Uncertainties ----------------------------------
        # in CSS the ellipse is projected onto the horizontal plane
        # using the covariance matrix
        semi_major = _km2m(db.get('smajax'))
        semi_minor = _km2m(db.get('sminax'))
        strike = db.get('strike')

        if all([semi_major, semi_minor, strike]):
            north, east = _get_ne_on_ellipse(semi_major, semi_minor, strike)
            lat_u = _m2deg_lat(north)
            lon_u = _m2deg_lon(east, lat=db.get('lat') or 0.0)

            uncertainty = Dict([
                ('preferredDescription', "uncertainty ellipse"),
                ('maxHorizontalUncertainty', semi_major),
                ('minHorizontalUncertainty', semi_minor),
                ('azimuthMaxHorizontalUncertainty', strike),
            ])
            if db.get('conf') is not None:
                uncertainty['confidenceLevel'] = db.get('conf') * 100.
        else:
            lat_u = None
            lon_u = None
            uncertainty = None

        # -- Basic Hypocenter ----------------------------------------
        origin = Dict([
            ('@publicID', self._uri(origin_rid)),
            ('latitude', _quan([
                ('value', db.get('lat')),
                ('uncertainty', lat_u),
                ])),
            ('longitude', _quan([
                ('value', db.get('lon')),
                ('uncertainty', lon_u),
                ])),
            ('depth', _quan([
                ('value', _km2m(db.get('depth'))),
                ('uncertainty', _km2m(db.get('sdepth'))),
                ])),
            ('time', _quan([
                ('value', self._utc(db.get('time'))),
                ('uncertainty', db.get('stime')),
                ])),
            ('quality', Dict([
                ('standardError', db.get('sdobs')),
                ('usedPhaseCount', db.get('ndef')),
                ('associatedPhaseCount', db.get('nass')),
                ])),
            ('originUncertainty', uncertainty),
            ('methodID', self._uri(method_rid, schema="smi")),
            ('earthModelID', self._uri(model_rid, schema="smi")),
            ('evaluationMode', mode),
            ('evaluationStatus', status),
            ('creationInfo', Dict([
                ('creationTime', self._utc(db.get('lddate'))),
                ('agencyID', agency),
                ('author', author),
                # ('version', db.get('orid')),
                ])),
            ('comment', [Dict([  # TODO: move this to extra info
                ('@id', self._uri(origin_rid, local_id="etype")),
                ('text', css_etype),
                ])]),
            ('arrival', []),
            # ('css:etype', css_etype),
        ])
        return origin

    def map_stamag2stationmagnitude(self, db):
        """
        Map stamag record to StationMagnitude
        """
        agency, author, _, method, info = self.split_auth(_str(db.get('auth')))

        css_sta = db.get('sta')
        css_chan = db.get('chan')
        wf_rid = "{0}/{1}-{2}-{3}".format(
            'wfdisc',
            css_sta,
            css_chan,
            int(db.get('time') * 10**6),
        )
        origin_rid = "{0}/{1}".format('origin',
                                      db.get('orid') or uuid.uuid4())
        stamag_rid = "{0}/{1}-{2}-{3}-{4}".format(
            'stamag',
            db.get('sta'),
            db.get('magtype'),
            db.get('orid') or uuid.uuid4(),
            db.get('magid') or uuid.uuid4(),
        )

        stationmagnitude = Dict([
            ('@publicID', self._uri(stamag_rid)),
            ('mag', Dict([
                ('value', db.get('magnitude')),
                ('uncertainty', db.get('uncertainty'))])),
            ('waveformID', Dict([
                ('@stationCode', db.get('fsta') or css_sta),
                ('@channelCode', db.get('fchan') or css_chan),
                ('@networkCode', db.get('snet') or self.default_network),
                ('@locationCode', db.get('loc') or ""),
                ('#text', self._uri(wf_rid, schema="smi")),
                # 'resourceURI' in schema
                ])),
            ('type', db.get('magtype')),
            ('creationInfo', Dict([
                ('creationTime', self._utc(db.get('lddate'))),
                ('agencyID', agency if agency != '' else self.agency),
                ('author', author),
                # ('version', db.get('magid')),
                ])),
            ('originID', self._uri(origin_rid)),
        ])
        return stationmagnitude

    def map_netmag2magnitude(self, db):
        """
        Return a dict of QuakeML magnitude from a dict of CSS key/values
        corresponding to one record.

        Inputs
        ======
        db : dict of key/values of CSS fields from the 'netmag' table

        Returns
        =======
        dict of key/values of QuakeML fields of "magnitude"

        Notes
        =====
        Any object that supports the dict 'get' method can be passed as
        input, e.g. OrderedDict, custom classes, etc.
        """
        agency, author, _, method, info = self.split_auth(_str(db.get('auth')))
        mode, status = self.get_event_status(author)

        origin_rid = "{0}/{1}".format('origin',
                                      db.get('orid') or uuid.uuid4())
        netmag_rid = "{0}/{1}".format('netmag',
                                      db.get('magid') or uuid.uuid4())

        magnitude = Dict([
            ('@publicID', self._uri(netmag_rid)),
            ('mag', _quan([
                ('value', db.get('magnitude')),
                ('uncertainty', db.get('uncertainty')),
                ])),
            ('type', db.get('magtype')),
            ('stationCount', db.get('nsta')),
            ('originID', self._uri(origin_rid)),
            ('evaluationMode', mode),
            ('evaluationStatus', status),
            ('creationInfo', Dict([
                ('creationTime', self._utc(db.get('lddate'))),
                ('agencyID', agency or self.agency),
                ('author', author),
                # ('version', db.get('magid')),
                ])),
        ])
        return magnitude

    def map_origin2magnitude(self, db, mtype='ml'):
        """
        Return a dict of magnitude from an dict of CSS key/values
        corresponding to one record.

        Inputs
        ======
        db : dict of key/values of CSS fields from the 'origin' table
        mtype : str of valid field from origin to use as mag ('ml', 'mb') etc

        Returns
        =======
        dict of key/values of QuakeML fields of "magnitude"

        Notes
        =====
        Any object that supports the dict 'get' method can be passed as
        input, e.g. OrderedDict, custom classes, etc.
        """
        agency, author, _, _, _ = self.split_auth(_str(db.get('auth')))
        mode, status = self.get_event_status(author)

        origin_rid = "{0}/{1}".format('origin',
                                      db.get('orid') or uuid.uuid4())

        # If foreign key to netmag table exists, use it as a unique id,
        # otherwise unique id is unique origin + local field
        netmagid = "{0}id".format(mtype)
        if db.get(netmagid):
            origmag_rid = "{0}/{1}".format('netmag', db.get(netmagid))
            public_uri = self._uri(origmag_rid)
        else:
            public_uri = self._uri(origin_rid, local_id=mtype)

        magnitude = Dict([
            ('@publicID', public_uri),
            ('mag', _quan(value=db.get(mtype))),
            ('type', mtype),
            ('originID', self._uri(origin_rid)),
            ('evaluationMode', mode),
            ('evaluationStatus', status),
            ('creationInfo', Dict([
                ('creationTime', self._utc(db.get('lddate'))),
                ('agencyID', agency or self.agency),
                ('author', author),
                # ('version', db.get('orid')),
                ])),
        ])
        return magnitude

    def map_arrival2pick(self, db):
        """
        Experimental map of just CSS arrival to QML pick.

        IF snetsta and schanloc are joined, will use those for SEED SNCL.
        Otherwise, will just use your converter agencyID for net and the
        sta/chan recorded with the pick.

        Inputs
        ======
        db : dict of key/values of CSS fields related to the phases (see Join)

        Returns
        =======
        dict of QuakeML schema for Pick type

        Notes
        =====
        Any object that supports the dict 'get' method can be passed as
        input, e.g. OrderedDict, custom classes, etc.

        Join
        ----
        arrival <- snetsta [sta] (outer) <- schanloc [sta chan] (outer)
        """
        agency, author, _, method, info = self.split_auth(_str(db.get('auth')))

        css_sta = db.get('sta')
        css_chan = db.get('chan')
        wf_rid = "{0}/{1}-{2}-{3}".format(
            'wfdisc',
            css_sta,
            css_chan,
            int(db.get('time') * 10**6),
        )
        pick_rid = "{0}/{1}".format('arrival',
                                    db.get('arid') or uuid.uuid4())

        on_qual = _str(db.get('qual')).lower()
        if 'i' in on_qual:
            onset = "impulsive"
        elif 'e' in on_qual:
            onset = "emergent"
        elif 'w' in on_qual:
            onset = "questionable"
        else:
            onset = None

        pol = _str(db.get('fm')).lower()
        if 'c' in pol or 'u' in pol:
            polarity = "positive"
        elif 'd' in pol or 'r' in pol:
            polarity = "negative"
        elif '.' in pol:
            polarity = "undecidable"
        else:
            polarity = None

        mode = "automatic"
        if 'orbassoc' not in _str(db.get('auth')):
            mode = "manual"

        status = "preliminary"
        if mode is "manual":
            status = "reviewed"

        pick = Dict([
            ('@publicID', self._uri(pick_rid)),
            ('time', _quan([
                ('value', self._utc(db.get('time'))),
                ('uncertainty', db.get('deltim')),
                ])),
            ('waveformID', Dict([
                ('@stationCode', db.get('fsta') or css_sta),
                ('@channelCode', db.get('fchan') or css_chan),
                ('@networkCode', db.get('snet') or self.default_network),
                ('@locationCode', db.get('loc') or ""),
                ('#text', self._uri(wf_rid, schema="smi")),
                # 'resourceURI' in schema
                ])),
            ('phaseHint', db.get('iphase')),  # 'code' in schema
            ('polarity', polarity),
            ('onset', onset),
            ('backazimuth', _quan([
                ('value', db.get('azimuth')),
                ('uncertainty', db.get('delaz'))
                ])),
            ('horizontalSlowness', _quan([
                ('value', db.get('slow')),
                ('uncertainty', db.get('delslo'))
                ])),
            ('creationInfo', Dict([
                ('creationTime', self._utc(db.get('arrival.lddate') or
                                           db.get('lddate'))),
                ('agencyID', agency or self.agency),
                ('author', author),
                # ('version', db.get('arid')),
                ])),
            ('evaluationMode', mode),
            ('evaluationStatus', status),
        ])
        return pick

    def map_assoc2arrival(self, db):
        """
        Experimental to map CSS assoc just to QML arrival

        Inputs
        ======
        db : dict of key/values of CSS fields related to the phases (see Join)

        Returns
        =======
        dict of QuakeML Arrival type

        Notes
        =====
        Any object that supports the dict 'get' method can be passed as
        input, e.g. OrderedDict, custom classes, etc.

        Join
        ----
        assoc
        """
        css_timedef = _str(db.get('timedef'))
        pick_rid = "{0}/{1}".format('arrival',
                                    db.get('arid') or uuid.uuid4())
        model_rid = "{0}/{1}".format('vmodel',
                                     db.get('vmodel') or uuid.uuid4())
        assoc_rid = "{0}/{1}-{2}".format('assoc',
                                         db.get('orid') or uuid.uuid4(),
                                         db.get('arid') or uuid.uuid4())

        arrival = Dict([
            ('@publicID', self._uri(assoc_rid)),
            ('pickID', self._uri(pick_rid)),
            ('phase', db.get('phase')),
            ('azimuth', db.get('esaz')),
            ('distance', db.get('delta')),
            ('timeResidual', db.get('timeres')),
            ('timeWeight', db.get('wgt')),
            ('earthModelID', self._uri(model_rid, schema="smi")),
            ('creationInfo', Dict([
                ('creationTime', self._utc(db.get('lddate'))),
                ('agencyID', self.agency),
                # ('version', db.get('arid')),
                ])),
            # ('css:timedef', css_timedef),
        ])

        # Assign a default weight based on timedef if none in db
        if arrival.get('timeWeight') is None:
            arrival['timeWeight'] = TIMEDEF_WEIGHT.get(css_timedef)

        return arrival

    def map_assocarrival2pickarrival(self, db):
        """
        Return tuple of quakeML (pick, arrival) from a dict of CSS key/values
        corresponding to one record. See the 'Join' section for the implied
        database table join expected.

        Inputs
        ======
        db : dict of key/values of CSS fields related to the phases (see Join)

        Returns
        =======
        tuple of dicts of key/values of QuakeML: ("pick", "arrival")

        Notes
        =====
        Any object that supports the dict 'get' method can be passed as
        input, e.g. OrderedDict, custom classes, etc.

        Join
        ----
        assoc <- arrival [arid] <- snetsta [sta] (outer)
            <- schanloc [sta chan] (outer)
        """
        pick = self.map_arrival2pick(db)
        arrival = self.map_assoc2arrival(db)
        return (pick, arrival)

    # TODO: amplitudes & station magnitude contributions

    def map_fplane2focalmech(self, db):
        """
        Return a dict of focalMechanism from an dict of CSS key/values
        corresponding to one record. See the 'Join' section for the implied
        database join expected.

        Inputs
        ======
        db : dict of key/values of CSS fields from the 'fplane' table

        Returns
        =======
        dict of key/values of QuakeML "focalMechansim"

        Notes
        =====
        Any object that supports the dict 'get' method can be passed as
        input, e.g. OrderedDict, custom classes, etc.
        """
        #
        # NOTE: Antelope schema for this is wrong, no nulls defined
        #
        origin_rid = "{0}/{1}".format('origin',
                                      db.get('orid') or uuid.uuid4())
        fplane_rid = "{0}/{1}".format('fplane',
                                      db.get('mechid') or uuid.uuid4())

        agency, author, _, _, _ = self.split_auth(_str(db.get('auth')))
        mode, status = self.get_event_status(author)
        method_rid = self.get_method_model(_str(db.get('algorithm')))[0]

        nodal_planes = Dict([
            ('nodalPlane1', Dict([
                ('strike', Dict(value=db.get('str1'))),
                ('dip', Dict(value=db.get('dip1'))),
                ('rake', Dict(value=db.get('rake1'))),
                ])),
            ('nodalPlane2', Dict([
                ('strike', Dict(value=db.get('str2'))),
                ('dip', Dict(value=db.get('dip2'))),
                ('rake', Dict(value=db.get('rake2'))),
                ])),
            ('@preferredPlane', 1),
        ])

        principal_axes = Dict([
            ('tAxis', Dict([
                ('azimuth', Dict(value=db.get('taxazm'))),
                ('plunge', Dict(value=db.get('taxplg'))),
                ])),
            ('pAxis', Dict([
                ('azimuth', Dict(value=db.get('paxazm'))),
                ('plunge', Dict(value=db.get('paxplg'))),
                ])),
        ])

        focal_mechanism = Dict([
            ('@publicID', self._uri(fplane_rid)),
            ('triggeringOriginID', self._uri(origin_rid)),
            ('nodalPlanes', nodal_planes),
            ('principalAxes', principal_axes),
            ('creationInfo', Dict([
                ('creationTime', self._utc(db.get('lddate'))),
                ('agencyID', agency if agency != '' else self.agency),
                ('author', author),
                # ('version', db.get('mtid')),
                ])),
            ('evaluationMode', mode),
            ('evaluationStatus', status),
            ('methodID', self._uri(method_rid)),
        ])
        return focal_mechanism

    # TODO: also do 'moment'???
    def map_moment2focalmech(self, db):
        """
        Map moment record to a FocalMechanism
        """
        raise NotImplementedError("No moment tensor support yet")

    def map_mt2focalmech(self, db):
        """
        Map BRTT CSS table 'mt' record to a FocalMechanism

        Notes
        =====
        1) This should not be first choice, mt table lacks many attributes of a
        moment tensor solution, only use if nothing else is available.

        2) This treats derived parameters weirdly, there can be an orid, but
        also derived lat/lon etc, which should be in the origin table? So this
        method uses orid as the triggering origin and ignores the derived
        origin parameters. A more comprehensive one would build
        origin/magnitude and make the necessary ID's but leaves that for other
        methods, i.e. map_mt2origin, map_mt2magnitude. There should be an "Mw"
        in the netmag table, not here.
        """
        origin_rid = "{0}/{1}".format('origin',
                                      db.get('orid') or uuid.uuid4())
        mt_rid = "{0}/{1}".format('mt',
                                  db.get('mtid') or uuid.uuid4())

        agency, author, _, _, _ = self.split_auth(_str(db.get('auth')))

        moment_tensor = Dict([
            ('@publicID', self._uri(mt_rid, local_id="tensor")),
            ('derivedOriginID', self._uri(origin_rid)),
            ('scalarMoment', _quan(value=db.get('scm'))),
            ('doubleCouple', db.get('pdc')),
            ('tensor', Dict([
                ('Mrr', Dict(value=db.get('tmrr'))),
                ('Mtt', Dict(value=db.get('tmtt'))),
                ('Mpp', Dict(value=db.get('tmpp'))),
                ('Mrt', Dict(value=db.get('tmrt'))),
                ('Mrp', Dict(value=db.get('tmrp'))),
                ('Mtp', Dict(value=db.get('tmtp'))),
                ])),
            ('creationInfo', Dict([
                ('creationTime', self._utc(db['lddate'])),
                ('agencyID', self.agency),
                ('author', author),
                # ('version', db.get('mtid')),
                ])),
        ])

        nodal_planes = Dict([
            ('nodalPlane1', Dict([
                ('strike', Dict(value=db.get('str1'))),
                ('dip', Dict(value=db.get('dip1'))),
                ('rake', Dict(value=db.get('rake1'))),
                ])),
            ('nodalPlane2', Dict([
                ('strike', Dict(value=db.get('str2'))),
                ('dip', Dict(value=db.get('dip2'))),
                ('rake', Dict(value=db.get('rake2'))),
                ])),
            ('@preferredPlane', 1),
        ])

        principal_axes = Dict([
            ('tAxis', Dict([
                ('azimuth', Dict(value=db.get('taxazm'))),
                ('plunge', Dict(value=db.get('taxplg'))),
                ('length', Dict(value=db.get('taxlength'))),
                ])),
            ('pAxis', Dict([
                ('azimuth', Dict(value=db.get('paxazm'))),
                ('plunge', Dict(value=db.get('paxplg'))),
                ('length', Dict(value=db.get('paxlength'))),
                ])),
            ('nAxis', Dict([
                ('azimuth', Dict(value=db.get('naxazm'))),
                ('plunge', Dict(value=db.get('naxplg'))),
                ('length', Dict(value=db.get('naxlength'))),
                ])),
        ])

        focal_mechanism = Dict([
            ('@publicID', self._uri(mt_rid, local_id="focalmech")),
            ('triggeringOriginID', self._uri(origin_rid)),
            ('nodalPlanes', nodal_planes),
            ('principalAxes', principal_axes),
            ('momentTensor', moment_tensor),
            ('creationInfo', Dict([
                ('creationTime', self._utc(db.get('lddate'))),
                ('agencyID', agency if agency != '' else self.agency),
                ('author', author),
                # ('version', db.get('mtid')),
                ])),
            # These determined from auth?? or weird incorrect mt fields...
            ('evaluationMode', db.get('rstatus')),
            ('evaluationStatus', db.get('estatus')),
        ])
        return focal_mechanism

    def convert_origins(self, records):
        """
        Origin converter

        Inputs
        ======
        records : iterable sequence of dict-like database records

        Returns
        =======
        list of dict
        """
        return [self.map_origin2origin(row) for row in records]

    def convert_phases(self, records):
        """
        Phase (Pick + Arrival) converter

        Inputs
        ======
        records : iterable sequence of dict-like database records

        Returns : tuple of (picks, arrivals) where
        =======
        picks    : list of dict
        arrivals :  list of dict
        """
        pick_arr_pairs = [self.map_assocarrival2pickarrival(row)
                          for row in records]
        return map(list, zip(*pick_arr_pairs))

    def convert_focalmechs(self, records, schema='fplane'):
        """
        FocalMechanism converter

        Inputs
        ======
        records : iterable sequence of dict-like database records
        schema : str name of CSS table indicating input type

        Returns
        =======
        list of dict

        Notes
        =====
        FocalMechanisms can be from first-motion or moment tensor inversions,
        the CSS3.0 extensions are denormalized, so you need the schema of the
        table, which can be specified in the function call.
        """
        if schema == "fplane":
            return [self.map_fplane2focalmech(row) for row in records]
        elif schema == "mt":
            return [self.map_mt2focalmech(row) for row in records]
        else:
            pass  # schema == "moment" case too?

    @staticmethod
    def description(nearest_string, type="nearest cities"):
        """
        Return a dict eventDescription of type 'nearest cities'

        Inputs
        ======
        nearest_string : str of decription for the text field
        """
        return Dict(text=nearest_string, type=type)

    def map_event(self, db, anss=False):
        """
        Create a QML Event from a CSS event
        (will also accept a CSS origin row dict)
        """
        evid = db.get('evid')
        lddate = db.get('lddate') or _ts(datetime.datetime.utcnow())
        prefor = db.get('prefor') or db.get('orid')
        event_rid = "{0}/{1}".format('event', evid)

        event = Dict([
            ('@publicID', self._uri(event_rid)),
            ('type', "not reported"),
            ('creationInfo', Dict([
                ('creationTime', self._utc(lddate)),
                ('agencyID', self.agency),
                # ('version', str(evid)),
                ])),
        ])
        # Set the prefor if you gave on origin or the event table has one
        if prefor:
            origin_rid = "{0}/{1}".format('origin', prefor)
            event['preferredOriginID'] = self._uri(origin_rid)
        #
        # Add the ANSS NS parameters automatically
        #
        if anss:
            cat_event_attr = anss_params(self.agency, evid)
            event.update(cat_event_attr)
        return event
