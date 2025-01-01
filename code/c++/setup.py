from setuptools import setup
from setuptools.extension import Extension

setup(
    name='SGEN_Py',
    version='0.1',
    packages=[],
    ext_modules=[
        Extension('SGEN_Py', sources=[])
    ],
    data_files=[('', ['./build/SGEN_Py.cpython-38-x86_64-linux-gnu.so'])]
)
