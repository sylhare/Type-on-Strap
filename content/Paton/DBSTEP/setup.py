from setuptools import setup
import io

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name='dbstep',
  packages=['dbstep'],
  version='1.1.0',
  description='DFT Based Steric Parameters',
  long_description=long_description,
  long_description_content_type='text/markdown',
  author='Guilian Luchini, Toby Patterson, Robert Paton',
  author_email='patonlab@colostate.edu',
  url = 'https://github.com/bobbypaton/DBSTEP',
  download_url = 'https://github.com/patonlab/DBSTEP/archive/refs/tags/1.1.0.zip',
  keywords=['compchem', 'steric', 'sterimol', 'informatics'],
  license="MIT",
  classifiers=['License :: OSI Approved :: MIT License',],
  install_requires=["numpy",'numba','scipy','cclib'],
  python_requires='>=3.6',
  include_package_data=True,
  package_data={'': ['*.csv']},
)
