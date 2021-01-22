# pyCHNOSZ

Author: Dr. Grayson Boyer, GEOPIG Lab, Arizona State University

This is a python wrapper for certain functions in the R package [CHNOSZ](https://www.chnosz.net/) by Dr. Jeffrey Dick.

## Features

Currently, the following CHNOSZ functions are supported in pyCHNOSZ:

- [subcrt](https://chnosz.net/manual/subcrt.html) - Calculate standard state partial molal thermodynamic properties of reactions at elevated temperatures and pressures.

## Requirements

This package must be installed into an environment that has an [R](https://www.r-project.org/) installation. See [these instructions](https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/) for installing R with Anaconda. Additionally, the CHNOSZ package for R must be installed (see instructions below).

## Installation

### Installing CHNOSZ

Open an R session. Install the CHNOSZ package with:

```r
install.packages("CHNOSZ")
```

Once CHNOSZ is installed you may quit the R session.

### Installing pyCHNOSZ

Install pyCHNOSZ using pip:

```
$ pip install pyCHNOSZ
```

## Usage

### Get thermodynamic properties of reactions

See the original documentation for [subcrt](https://chnosz.net/manual/subcrt.html). Useage in pyCHNOSZ is the same, except python lists are used in place of R's vectors. The function produces a dictionary of results stored in pandas dataframes.

```r
subcrt(species=["H2", "O2", "H2O"], coeff=[-1.0, -0.5, 1.0],
       state=["aq", "gas", "liq"], T=[30, 50], P=100)
```
```
subcrt: 3 species at 2 values of T (ºC) and P (bar) (wet) [energy units: cal]

{'reaction':       coeff    name formula state  ispecies
 62     -1.0      H2      H2    aq      62.0
 2612   -0.5  oxygen      O2   gas    2612.0
 1       1.0   water     H2O   liq       1.0,
 'out':       T    P       rho       logK             G             H         S  \
 1  30.0  100  1.000017  43.855086 -60832.380282 -67420.887872 -21.89070   
 2  50.0  100  0.992305  40.834419 -60379.262657 -67882.530994 -23.36663   
 
           V         Cp  
 1 -7.494052 -24.126268  
 2 -8.259704 -20.941879  }
```