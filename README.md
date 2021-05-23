## An IDL package to prepare and reduce Keck NIRC2 observations
*written by H.Fu around 2012* 

## References
- [2012ApJ...745...67F](https://ui.adsabs.harvard.edu/abs/2012ApJ...745...67F): *The Nature of Double-peaked [O III] Active Galactic Nuclei*
- [2012ApJ...753..134F](https://ui.adsabs.harvard.edu/abs/2012ApJ...753..134F): *A Comprehensive View of a Strongly Lensed Planck-Associated Submillimeter Galaxy*

## Demo Example

An example with raw data, the reduced image, and the reduction script are provided for in the [example folder](https://github.com/fuhaiastro/nirc2/tree/main/example).

## Observation Preparation

- `nirc2findtt`
	- finds NIRC2 tip-tilt stars for a list of targets given RA & Dec 
	- uses both SDSS and USNO-B from QueryVizier
- `nirc2fchart`
	- Makes finder charts to prepare for NIRC2 observations.
	- Uses SDSS color image as the background
	- Overlay tip-tilt reference star patrol field

- `nirc2_obs_manual.txt`
	- some notes for observers

## Reduction of Imaging Data 

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


