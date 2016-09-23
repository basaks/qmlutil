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
qmlutil.xml

Mark Williams (2015) Nevada Seismological Laboratory

Module to convert QML (QuakeML Modeling Language) structures
to XML (QuakeML) and vice versa

Dumping is easy, assuming you keep the xmltodict syntax in your object.
Loading from file/string is easy for the non-validating parser, but if you
want types, we need to have access to the QuakeML XML schema. This module
attempts to do this.

Note: Typing is thankfully NOT needed if the application is just creating SQL
statements or anything in HTTP request/responses (except JSON, of course).

"""
from qmlutil.lib import xmltodict

# Module-level dict for XSD -> python type mapping
# - Change xs:dataTime to datetime.datetime.strptime e.g.
# - Maybe make class attribute for flexibility
# TODO: add other legit ones, float, real, etc
XSTYPES = {
    'xs:string': str,
    'xs:integer': int,
    'xs:double': float,
    'xs:boolean': bool,
    'xs:dateTime': str,
    'xs:anyURI': str,
}


def dt(datetime_string):  # pylint: disable=unused-argument
    """
    Parse datetime string to python datetime.datetime
    """
    raise NotImplementedError("Not done yet")


class TypeExtractor(object):
    """Object to validate/entype XML"""
    XSDtypes = None  # holds flat map of nested elements/XML types
    PYtypes = None   # holds flat map of nested keys/python types

    delim = '|'  # Nested key delimiter in type maps
    ns = "bed"   # namespace of XSDtypes keys

    def flatten(self, node, name=""):
        """
        Craete a flat map of XML types from xsd node
        """
        if isinstance(node, dict):
            if '@name' in node:
                # print "Name: {0}, Type: {1}".format(node['@name'],
                #                                     node.get('@type'))
                name += self.delim + node['@name']
                name = name.strip(self.delim)
                if '@type' in node:
                    self.XSDtypes[name] = node['@type']
            if '@base' in node:
                # print "Name: {0}, Base: {1}".format(node.get('@name'),
                #                                     node.get('@base'))
                self.XSDtypes[name] = node['@base']

            for subnode in node:
                # print "Key: {0}".format(subnode)
                if not subnode.startswith('@'):
                    self.flatten(node[subnode], name)
                # else:
                #    print "Attribute: {}, STOPPING".format(subnode)
        elif isinstance(node, list):
            for subnode in node:
                self.flatten(subnode, name)

    def __init__(self, qml):
        self.qml = qml
        self.XSDtypes = dict()
        self.PYtypes = dict()

    def entype(self, node, name=""):
        """
        Entype the whole dict/list struct under "node" given previously built
        types in self.PYtypes

        todo, make deepcopy??
        """
        if isinstance(node, dict):
            for subnode in node:
                rname = self.delim.join([name, subnode]).strip(self.delim)
                type_ = self.PYtypes.get(rname)
                if isinstance(node[subnode], list):
                    for subsubnode in node[subnode]:
                        self.entype(subsubnode, rname)
                elif isinstance(node[subnode], dict):
                    self.entype(node[subnode], rname)
                elif type_:
                    node[subnode] = type_(node[subnode])
        return self.qml

    def gentypes(self, node, name="", realname=""):
        """
        Generate flat map of python types for every node in the tree

        Map contains python types for given nodes
        Recursively try nested nodes in dict or list
        """
        if isinstance(node, dict):
            # Get new node names based on key
            for subnode in node:
                keyname = self.delim.join([name, subnode.lstrip('@')])
                keyname = keyname.strip(self.delim)
                rname = self.delim.join([realname, subnode]).strip(self.delim)
                self.gentypes(node[subnode], keyname, rname)
        elif isinstance(node, list):
            # Just pass on node names
            for subnode in node:
                self.gentypes(subnode, name, realname)
        else:
            # Try to get a type
            type_ = self._gettype(name)
            if isinstance(type_, str) or isinstance(type_, unicode):
                if type_ in XSTYPES:
                    type_ = XSTYPES[type_]
                elif type_.startswith("bed:"):
                    type_ = self._gettype(type_.lstrip("bed:"))
                    if type_ in XSTYPES:
                        type_ = XSTYPES[type_]
                # Got a code, add to map of python types
                if isinstance(type_, type):
                    self.PYtypes[realname] = type_

    # TODO: Ugly -- clean this up
    # TODO: use generic settable Ns instead of bed diectly
    def _gettype(self, key):
        """
        Follow the types through linked keys to get a basic type
        """
        type_ = None
        key_parts = key.split(self.delim)
        if len(key_parts) <= 2:
            if key in self.XSDtypes:
                value = self.XSDtypes[key]
                if value.startswith("bed:"):
                    value = self.XSDtypes[value.lstrip("bed:")]
                    type_ = self._gettype(value)
                type_ = value
            else:
                if key_parts[0] in self.XSDtypes:
                    value = self.XSDtypes[key_parts[0]].lstrip("bed:")
                    key_parts[0] = value
                    newkey = self.delim.join(key_parts)
                    type_ = self._gettype(newkey)
        else:
            if key_parts[0] in self.XSDtypes:
                value = self.XSDtypes[key_parts[0]].lstrip("bed:")
                key_parts[0] = value
                newkey = self.delim.join(key_parts)
                type_ = self._gettype(newkey)
            elif self.delim.join(key_parts[:2]) in self.XSDtypes:
                value = self.XSDtypes[self.delim.join(key_parts[:2])].lstrip("bed:")
                key_parts = [value] + key_parts[2:]
                newkey = self.delim.join(key_parts)
                type_ = self._gettype(newkey)
        return type_

    def extract_typed(self):
        """
        Build type map from data and convert types
        """
        self.gentypes(self.qml)
        return self.entype(self.qml)


class Rounder(object):  # pylint: disable=too-few-public-methods
    """
    Rounder is an xmltodict preprocessor function for generating NSL QuakeML

    Notes
    -----
    Rounds specified values for objects in an Event because the Client doesn't
    understand precision in computing vs. precision in measurement, or the
    general need for reproducibility in science.
    """
    @staticmethod
    def _round(dict_, key, num_digits):
        """
        Round a number given a dict, a key, and # of places.
        """
        value = dict_.get(key)
        if value is not None:
            value = round(value, num_digits)
            dict_[key] = value

    def __init__(self):
        pass

    def __call__(self, key, dict_):
        # Case of integer attribute
        if key == "nodalPlanes" and dict_.get("@preferredPlane"):
            dict_['@preferredPlane'] = str(dict_['@preferredPlane'])

        # USGS can't handle ID in content yet
        if key == "waveformID":
            dict_.pop('#text')

        # Don't serialize empty stuff
        if dict_ is None:
            return None
        # Caveat to that is, have to enforce QuakeML rules:
        #
        # arrival: must have phase
        if key == "arrival" and isinstance(dict_, list):
            dict_ = [p for p in dict_ if p.get('phase') is not None]

        # Round stuff TODO: move to decorator/method
        if key == "depth":
            self._round(dict_, 'value', -2)
            self._round(dict_, 'uncertainty', -2)
            # TODO: lowerUncertainty, upperUncertainty, confidenceLevel??
        elif key == "latitude":
            self._round(dict_, 'uncertainty', 4)
        elif key == "longitude":
            self._round(dict_, 'uncertainty', 4)
        elif key == "time":
            self._round(dict_, 'uncertainty', 6)
        elif key == "time":
            self._round(dict_, 'uncertainty', 6)
        elif key == "originUncertainty":
            self._round(dict_, 'horizontalUncertainty', -1)
            self._round(dict_, 'minHorizontalUncertainty', -1)
            self._round(dict_, 'maxHorizontalUncertainty', -1)
        elif key == "mag":
            self._round(dict_, 'value', 2)
            self._round(dict_, 'uncertainty', 3)
        return key, dict_


def ignore_null(key, value):
    """
    Preprocessor for xmltodict.unpasre that ignores keys with None value
    """
    if value is None:
        return None
    return key, value


def dumps(input_dict, *args, **kwargs):
    """
    Dump QML dict object to XML string
    """
    return xmltodict.unparse(input_dict, *args, **kwargs)


# TODO: add kwargs for typing/conversions???
def loads(xml_input, *args, **kwargs):
    """Load QML dict object from XML"""
    return xmltodict.parse(xml_input, *args, **kwargs)


# Testing
def main():
    with open('data/QuakeML-BED-1.2.xsd') as f:
        schema = loads(f)

    with open('tests/quakeml.xml') as f:
        qroot = loads(f)
    qml = qroot['q:quakeml']

    # extractor set with qml pulled from xmltodict
    type_extractor = TypeExtractor(qml)
    # build type schema map from xsd pulled from xmltodict
    type_extractor.flatten(schema)

    qml = type_extractor.extract_typed()
    return qml, type_extractor, schema

# Run test
# q, te, schema = main()
