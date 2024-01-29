#!/usr/bin/env python

# ========================================================================
#  Copyright Het Nederlands Kanker Instituut - Antoni van Leeuwenhoek
#
#  Licensed under the 3-clause BSD License
# ========================================================================

from setuptools import setup

with open('requirements.txt', 'r') as fp:
  requirements = list(filter(bool, (line.strip() for line in fp)))

setup(
  name='nrrdify',

  author='Joost van Griethuysen',
  author_email='j.v.griethuysen@nki.nl',

  version='0.3',

  packages=['nrrdify'],
  zip_safe=False,

  entry_points={
        'console_scripts': [
            'nrrdify=nrrdify.__main__:main'
        ]},

  description='DICOM conversion script to check and sort DICOM slices prior to conversion using SimpleITK',
  license='BSD License',

  classifiers=[
    'Development Status :: 1 - Planning',
    'Environment :: Console',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Operating System :: Microsoft :: Windows'
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Utilities'
  ],

  install_requires=requirements,

  keywords='cancerimaging medicalresearch computationalimaging'
)
