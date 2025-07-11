import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyCHNOSZ",
    version="0.9.1",
    author="Grayson Boyer",
    author_email="gmboyer@asu.edu",
    description="Python wrapper for the R package CHNOSZ.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={},
    packages=['pyCHNOSZ'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    install_requires=['rpy2', 'pandas', 'ipython', 'plotly', 'simplegeneric', 'chemparse', 'WORMutils', 'matplotlib', 'numpy', 'wormutils_r', 'statistics'],
    package_data={'': ['*.r']},
    include_package_data=True,
    zip_safe=False
)

