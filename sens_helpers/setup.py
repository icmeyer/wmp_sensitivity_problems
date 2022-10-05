#!/usr/bin/env python
import glob
from setuptools import setup, find_packages

with open('sens_helpers/__init__.py', 'r') as f:
    version = f.readlines()[-1].split()[-1].strip("'")

kwargs = {
    'name': 'sens_helpers',
    'version': version,
    'packages': find_packages(),

    # Metadata
    'author': 'Isaac Meyer',
    'author_email': 'icmeyer@mit.edu',
    'description': 'package for running OpenMC CLUTCH sensitivity',
}

setup(**kwargs)
