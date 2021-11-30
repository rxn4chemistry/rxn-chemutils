# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

import io
import os
import re

from setuptools import setup, find_packages

# Get the version
match = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
    io.open('rxn_chemutils/__init__.py', encoding='utf_8_sig').read()
)
if match is None:
    raise SystemExit('Version number not found.')
__version__ = match.group(1)

# Get the long description from the README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='rxn_chemutils',
    version=__version__,
    author='IBM RXN team',
    description='RXN chemistry-related utilities',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    package_data={'rxn_chemutils': ['py.typed']},
    scripts=[
        'scripts/chemutils-canonicalize',
        'scripts/chemutils-detokenize',
        'scripts/chemutils-tokenize',
        'scripts/chemutils-depictready',
        'scripts/chemutils-remove-mapping',
    ],
    install_requires=['attrs>=21.2.0', 'click>=7.0', 'loguru>=0.5.3', 'rxn-utilities>=0.8.1'],
    extras_require={
        'dev': [
            'flake8>=3.8.4',
            'mypy>=0.910',
            'pytest>=5.3.4',
            'yapf>=0.31.0',
        ],
    },
)
