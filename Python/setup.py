"""Setup file to Not Another Neuroimaging Slicer

"""

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.abspath(path.join(here, '../README.md')), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='QUIT',
    version='0.1.0',
    description='nipype interfaces to QUantitative Imaging Tools',
    url='https://github.com/spinicist/quit',
    author='Tobias Wood',
    author_email='tobias@spinicist.org.uk',
    py_modules=['QUIT'],
    install_requires=['nipype',
                      'nibabel'],
    license='MPL',
    classifiers=['Development Status :: 4 - Beta',
                 'Topic :: Scientific/Engineering :: Visualization',
                 'Programming Language :: Python :: 3',
                 ],
    keywords='neuroimaging nifti',
    packages=find_packages()
)
