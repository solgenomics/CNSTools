from setuptools import setup
from __consts import dependencies

setup(name='cnstools',
      version='0.1',
      description='Conserved noncoding sequence identification tools. Developed for the Mueller lab at Boyce Thompson Institute.',
      author='David Lyon',
      author_email='dal333@cornell.edu',
      license='?',
      packages=['cnstools'],
      install_requires= dependencies,
      zip_safe=True)