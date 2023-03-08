import sys

from glob import glob
from pybind11 import get_cmake_dir
# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import intree_extensions, build_ext
from setuptools import setup, find_packages

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

__version__ = "0.1.1"

#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

ext_modules = intree_extensions(sorted(glob("n_link_sim/*.cpp")))

setup(
    name="n-link-simulator",
    version=__version__,
    author="Sebastian Markgraf",
    author_email="sebastian-markgraf@t-online.de",
    url="https://github.com/sebimarkgraf/n_link_simulator",
    description="Simulators for double and quad link dynamics, taken from pypost toolbox",
    long_description=long_description,
    ext_modules=ext_modules,

    extras_require={"test": "pytest"},
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=["numpy"],
    packages=find_packages(include=["n_link_sim"])
)

