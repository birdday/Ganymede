#!/usr/bin/env python

from distutils.core import setup

setup(
    name="lcparser",
    version="0.0.1",
    description="Liquid Chromatography data parser for Ganymede interview.",
    packages=["lcparser"],
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "mpl_interactions",
        "BaselineRemoval",
        "ipympl",
        "ipykernel",
        "pytest",
    ],
)
