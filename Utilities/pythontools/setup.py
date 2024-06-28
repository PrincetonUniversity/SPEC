import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="py_spec",
    version="3.3.5",
    install_requires=["numpy>=1.21.1", 
                      "f90nml",
                      "h5py", 
                      "matplotlib",
                      "coilpy; python_version<'3.12'",
                      "scipy>=1.7.0"],
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
    packages=['py_spec', 'py_spec.input', 'py_spec.output', 'py_spec.ci']
)
