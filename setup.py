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
setup.py for qmlutil

- Mark Williams, Nevada Seismological Laboratory
- University of Nevada, Reno, (2015)
- Sudipta Basak, Geoscience Australia (2017)

"""
from setuptools import setup

s_args = {
    'name': 'qmlutil',
    'version': '0.5.0',
    'description': 'QuakeML/Seiscomp3ML utils for python',
    'author': 'Mark Williams, Sudipta Basak',
    'maintainer': 'started by Nevada Seismological Laboratory, '
                  'Now manintained by Geoscience Australia',
    'url': 'https//github.com/basaks/qmlutil',
    'setup_requires': ['numpy >= 1.9.2', 'Cython'],  # required due to scipy
    'install_requires': [
        'scipy >= 0.15.1',
        'matplotlib >= 1.5.1',
        'obspy'
    ],
    'dependency_links': [
        'https://github.com/obspy/obspy/archive/'
        'master.zip#egg={package}'.format(package='obspy')
      ],
    'extras_require': {
        'demos': [
            'matplotlib'
        ],
        'kmz': [
            'simplekml',
            'pillow'
        ],
        'dev': [
            'sphinx',
            'ghp-import',
            'sphinxcontrib-programoutput'
        ]
    },
    'packages': [
        'qmlutil',
        'qmlutil.plugins',
        'qmlutil.css',
        'qmlutil.data',  # TODO: add static files IMPORTANT!!
        'qmlutil.lib',
        'qmlutil.ichinose',
    ],
    'package_dir': {'qmlutil': 'qmlutil'},
    'package_data': {
        'qmlutil' : [
            'data/*.rng',
            'lib/README',
        ],
    },
}

setup(**s_args)
