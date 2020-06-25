"""Setup file to Not Another Neuroimaging Slicer

"""

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.abspath(path.join(here, '../README.md')), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='qipype',
    version='1.0',
    description='nipype interfaces to QUantitative Imaging Tools',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/spinicist/quit',
    author='Tobias Wood',
    author_email='tobias@spinicist.org.uk',
    py_modules=['qipype'],
    install_requires=['nipype>=1.2.3',
                      'nibabel>=2.5.1'],
    python_requires='>=3',
    license='MPL',
    classifiers=['Topic :: Scientific/Engineering :: Physics',
                 'Programming Language :: Python :: 3',
                 ],
    keywords='neuroimaging mri',
    packages=find_packages()
)
