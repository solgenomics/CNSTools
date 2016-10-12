from setuptools import setup
import sys
if sys.version_info < (2,7) or not sys.version_info < (2,8):
    sys.exit('cnstools requires python 2.7')
setup(name='cnstools',
      version='0.1',
      description='Conserved noncoding sequence identification tools. Developed for the Mueller lab at Boyce Thompson Institute.',
      author='David Lyon',
      author_email='dal333@cornell.edu',
      license='?',
      packages=['cnstools'],
      install_requires= [],
      zip_safe=True)