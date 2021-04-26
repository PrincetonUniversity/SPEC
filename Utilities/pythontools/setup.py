import setuptools
from py_spec import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="py_spec",
    version=__version__,
    description="SPEC(Stepped-Pressure Equilibrium Code) python utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
    ],
    url="https://princetonuniversity.github.io/SPEC/",
    author="SPEC developers",
    license="GNU 3.0",
    packages=setuptools.find_packages(),
)
