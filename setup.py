#! /usr/bin/env python
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext
import setuptools

with open("src/README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="matplotlib_user_interface", 
    version="2021.03.19",
    author="Kevin Volk",
    author_email="kvolk@stsci.edu",
    description="A Matplotlib GUI program reminiscent of xmgrace.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KevinVolkSTScI/matplotlib-gui-code",
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'astropy>=4.0',
        'matplotlib>=3.0.0',
        'numpy',
        'scipy>=1.1.0',
        'matplotlib'
    ],                                                                 
)
