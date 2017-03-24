#!/usr/bin/env python

from setuptools import setup

setup(name='cnstools',
    version='0.2.1',
    description='',
    author='David Lyon',
    author_email='dlyon@fandm.edu',
    url='https://github.com/solgenomics/cnstools',
#    install_requires=['pillow'],
    packages=['cnstools'],
    entry_points = {
        'console_scripts': ['cnstools = cnstools.__main__:main'],
    }
)