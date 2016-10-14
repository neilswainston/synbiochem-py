'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from setuptools import find_packages, setup


setup(name='synbiochem-py',
      version='0.1.0',
      description='synbiochem-py: Core python modules for SYNBIOCHEM',
      long_description='synbiochem-py: Core python modules for SYNBIOCHEM',
      url='https://github.com/synbiochem/synbiochem-py',
      author='Neil Swainston',
      author_email='neil.swainston@manchester.ac.uk',
      license='MIT',
      classifiers=[
                   'Development Status :: 4 - Beta',
                   'Intended Audience :: Developers',
                   'Topic :: Software Development :: Build Tools',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 2.7'
      ],
      keywords='synbio synthetic biology',
      packages=find_packages(),
      test_suite='synbiochem.utils.test',
      install_requires=['regex', 'numpy', 'scipy', 'matplotlib', 'biopython'])
