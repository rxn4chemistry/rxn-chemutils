import io
import os
import re

from setuptools import setup

__version__ = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
    io.open('rxn_chemutils/__init__.py', encoding='utf_8_sig').read()
).group(1)

setup(
    name='rxn_chemutils',
    version=__version__,
    author='IBM RXN team',
    packages=['rxn_chemutils'],
    package_data={'rxn_chemutils': ['py.typed']},
    install_requires=[
        'attrs>=19.1.0',
        'click>=7.0',
        'rxn_utilities '
        '@ git+https://{}@github.ibm.com/rxn/rxn_utilities@latest'.format(os.environ['GHE_TOKEN']),
    ]
)
