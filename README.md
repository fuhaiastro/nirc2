## An IDL package to prepare and reduce Keck NIRC2 observations
*written by H.Fu around 2012* 

## References
- [2012ApJ...745...67F](https://ui.adsabs.harvard.edu/abs/2012ApJ...745...67F): *The Nature of Double-peaked [O III] Active Galactic Nuclei*
- [2012ApJ...753..134F](https://ui.adsabs.harvard.edu/abs/2012ApJ...753..134F): *A Comprehensive View of a Strongly Lensed Planck-Associated Submillimeter Galaxy*

## Demo Example

An example with raw data, the reduced image, and the reduction script
are provided for in the [example
folder](https://github.com/fuhaiastro/nirc2/tree/main/example). Once you
have followed the steps below to install the code, you can do the
following commands to test your installation:
```shell
cd ~/idl/nirc2/example/
gunzip cal/*.gz
gunzip cosmos01_ks/*.gz
nirc2
NIRC2> .r redux.idlsrc
```
The script takes about 5 min to complete the full reduction of the data
set and generates a new output FITS file. 

## Installation Guide

1. Download source code
```shell
cd ~/idl
git clone https://github.com/fuhaiastro/nirc2.git
```
2. Set environmental variables and define start-up command to launch `spfit`.
```shell
vi ~/.bash_profile
alias nirc2='source ~/idl/nirc2/setup.sh; /Applications/harris/idl/bin/idl -IDL_PROMPT "SPFIT>"'
```
3. Make sure the appropriate paths are set correctly for your host.
```shell
vi ~/idl/nirc2/setup.sh
```
Below is the content of the setup file:
```shell
# Bash shell commands to define IDL environment variables and aliases.
. /Applications/harris/idl/bin/idl_setup.bash

# start up IDL path
export IDL_PATH=+$IDL_DIR/examples:+$IDL_DIR/lib

# local user's IDL folder
export IDL=$HOME/idl

# NIRC2_IDL - NIRC2 data reduction & analysis
export NIRC2_DIR=$IDL/pipelines/nirc2
export IDL_PATH=${IDL_PATH}:+$NIRC2_DIR/pro

# If preferred, use latest Astro, Coyote, MPFIT packages
#export IDL_PATH=${IDL_PATH}:+$IDL/astron/pro:+$IDL/coyote:+$IDL/mpfit

# IDLUTILS (which includes Astrolib, Coyote, & MPFIT)
export IDLUTILS_DIR=$IDL/idlutils
export IDL_PATH=${IDL_PATH}:+$IDLUTILS_DIR/goddard/pro:+$IDLUTILS_DIR/pro
export PATH=$IDLUTILS_DIR/bin:$PATH

# Other useful IDL programs
export IDL_PATH=${IDL_PATH}:+$NIRC2_DIR/utilities
```
4. Install IDLUTILS
```shell
cd ~/idl
mkdir idlutils
cd idlutils
svn co https://svn.sdss.org/public/repo/sdss/idlutils/trunk .
export IDL_DIR=/Applications/harris/idl
export IDL=$HOME/idl
export IDLUTILS_DIR=$IDL/idlutils
$IDLUTILS_DIR/bin/evilmake clean
$IDLUTILS_DIR/bin/evilmake
```

## Routines for Observation Preparation

- `nirc2findtt`
	- finds NIRC2 tip-tilt stars for a list of targets given RA & Dec 
	- uses both SDSS and USNO-B from QueryVizier
- `nirc2fchart`
	- Makes finder charts to prepare for NIRC2 observations.
	- Uses SDSS color image as the background
	- Overlay tip-tilt reference star patrol field

- `nirc2_obs_manual.txt`
	- some notes for observers

## Routines for Imaging Data Reduction

- `nirc2dark`
	- median combine dark fields

- `nirc2flat`
	- generates flat-field using lamp-on and lamp-off images

- `nirc2redux`
	- A simple data reduction pipeline to reduce NIRC2 and OSIRIS
	  near-IR imaging data

- `nirc2_wide_X_distortion.fits`, `nirc2_wide_X_distortion.fits`
	- NIRC2 wide camera distortion solution derived from M92. The solution maps the image to a grid defined by the nominal pixel size and the NIRC2 PA = ROTPOSN - INSTANGL - 0.252 degrees. 
	
- `nirc2distort`
	- uses the distortion solution to correct the images. To apply within IRAF drizzle, follow the same procedure on [AstroBetter](https://www.astrobetter.com/wiki/tiki-index.php?page=NIRC2+Distortion+Correction).

- `nirc2ast`
	- Calibrate the astrometry of reduced NIRC2 images against SDSS
	  or DSS catalogs.


