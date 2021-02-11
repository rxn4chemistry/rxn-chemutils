# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

import io
import os
import re

from setuptools import setup, find_packages

match = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
    io.open('rxn_chemutils/__init__.py', encoding='utf_8_sig').read()
)
if match is None:
    raise SystemExit('Version number not found.')
__version__ = match.group(1)

setup(
    name='rxn_chemutils',
    version=__version__,
    author='IBM RXN team',
    packages=find_packages(),
    package_data={'rxn_chemutils': ['py.typed']},
    scripts=[
        'scripts/chemutils-canonicalize',
        'scripts/chemutils-detokenize',
        'scripts/chemutils-tokenize',
        'scripts/chemutils-depictready',
    ],
    install_requires=[
        'attrs>=19.1.0',
        'click>=7.0',
        'loguru>=0.5.3',
        'rxn_utilities '
        '@ git+https://{}@github.ibm.com/rxn/rxn_utilities@latest'.format(os.environ['GHE_TOKEN']),
    ]
)
