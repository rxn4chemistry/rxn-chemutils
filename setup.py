import io
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
    install_requires=[]
)
