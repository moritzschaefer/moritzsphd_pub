#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Setup file for moritzsphd.
    Use setup.cfg to configure your project.

    This file was generated with PyScaffold 3.1.
    PyScaffold helps you to put up the scaffold of your new Python project.
    Learn more under: https://pyscaffold.org/
"""
import sys

from pkg_resources import VersionConflict, require
from setuptools import find_packages, setup

try:
    require('setuptools>=38.3')
except VersionConflict:
    print("Error: version of setuptools is too old (<38.3)!")
    sys.exit(1)


if __name__ == "__main__":
    setup(use_pyscaffold=False,
        version="0.3",
        package_dir={"": "src"},
        packages=find_packages('src/'),
        install_requires=[
            'pandas==1.2.3',
            'numpy==1.21.0',
            'scipy==1.6.3',
            'requests==2.25.1',
            'biopython==1.76',
            'importmagic==0.1.7',
            'matplotlib==3.3.2',
            'tensorboard==2.5.0',
            # 'pycrypto',
            'pyproj==3.1.0',
            'scikit-learn==0.24.2',
            'seaborn==0.11.1',
            'lxml==4.6.3',
            'pyensembl==1.9.1',
            'gffutils==0.10.1',
            'pytz==2021.1',
            'venn==0.1.3',
            'xlrd==2.0.1',
            'gseapy==0.10.4',
            'pyperclip==1.8.2',
            'scikit-plot @ git+https://github.com/moritzschaefer/scikit-plot@feature/label-dots#egg=scikitplot',
            'tqdm==4.61.1'
            # 'pybedtools'
        ]
    )
