#!/usr/bin/env python3
"""
ACMG Variant Classification System
Version: 1.1.0

Setup script for package installation.
"""

from setuptools import setup, find_packages
import os

# Read version from _version.py
version = {}
with open(os.path.join(".", "_version.py"), "r") as f:
    exec(f.read(), version)

# Read requirements
def get_requirements():
    with open("requirements.txt", "r") as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="acmg-variant-classification",
    version=version["__version__"],
    description="ACMG-based variant classification system with Bayesian scoring",
    long_description=open("README.md", "r", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Variant Classification Team",
    author_email="",
    url="",
    packages=find_packages(exclude=["test*", "archives*", "ClinGen*"]),
    install_requires=get_requirements(),
    python_requires=">=3.9",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    entry_points={},
    include_package_data=True,
    zip_safe=False,
)
