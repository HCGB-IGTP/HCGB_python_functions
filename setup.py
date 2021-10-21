#!/usr/bin/env python3
import setuptools
from devel import setup_module

long_description_text = ""
with open("README.md", "r") as fh:
    long_description_text = fh.read()

setuptools.setup(
    name="HCGB",
    version=setup_module.get_version("./VERSION"),

    author="Jose F. Sanchez-Herrero",
    description="Useful python functions",

    author_email="jfbioinformatics@gmail.com",

    long_description_content_type="text/markdown",
    long_description=long_description_text,

    url="https://github.com/HCGB-IGTP/HCGB_python_functions/",

    packages=setuptools.find_packages(),
    license='MIT License',

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=setup_module.get_require_modules("./HCGB/config/python_requirements.txt"),
    #['pandas', 'termcolor', 'biopython', 'wget', 'xlsxwriter', 'patool'],
)
