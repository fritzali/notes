"""

PLATO Stellar Light-curve Simulator (PSLS)

Copyright (c) 2014, October 2017, R. Samadi (LESIA - Observatoire de Paris)

This is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this code.  If not, see <http://www.gnu.org/licenses/>.
"""
from setuptools import setup # , find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.txt").read_text(encoding="utf-8")

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(

    name="psls",
    version="1.9",
    description='PLATO Stellar Light-curve Simulator (SLS): Simulate stochastically-excited oscillations together with other stellar and instrumental components',
    long_description=long_description,  # Optional
    long_description_content_type="text/plain",  # Optional (see note above)
    url="https://sites.lesia.obspm.fr/psls/",  # Optional
    author="R. Samadi",  # Optional
    author_email="reza.samadi@obspm.fr",  # Optional
    # package_dir={"": "src"},  # Optional
    scripts=['psls.py'],
    py_modules=['sls', 'universal_pattern', 'FortranIO', 'transit', 'spotintime','flares'],
    python_requires=">=3.5",
    install_requires=['numpy', 'scipy', 'pyyaml', 'h5py', 'packaging','astropy','matplotlib']
)