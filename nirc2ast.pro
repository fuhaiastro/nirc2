pro nirc2ast,infile,catalog=catalog,reset=reset
;+
; NAME:
;       NIRC2AST
; PURPOSE:
;       updated with FIXAST.pro
;       Calibrate astrometry for NIRC2 images against SDSS or DSS. 
;
; INPUTS:
;       infile - input file to be calibrated
;       catalog - default SDSS-DR8, could also use 'II/319/las9' for UKIDSS
;
; REVISION HISTORY:
;       Written - Hai Fu - April 2012
;       
;-

; load image, extension 0 only
img = mrdfits(infile+'.fits',0,h0)
if keyword_set(reset) then begin
	camera = strtrim(sxpar(h0,'CAMNAME'),2)
	if (camera EQ 'narrow') then $ 
   	pixscale = 0.009942d $ ; arcsec/pix narrow camera
	else $
   	pixscale = 0.039686d ; arcsec/pix wide camera
	par = (sxpar(h0,'ROTDEST')-sxpar(h0,'INSTANGL')-0.252)*!dtor ; PA in radian
	putast,h0,[[-1.*cos(par),sin(par)],[sin(par),cos(par)]]*pixscale/3600d,$
		[512+15+0.5,512+30+0.5],[sxpar(h0,'ra'),sxpar(h0,'dec')],$
		['RA---TAN','DEC--TAN']
endif
; fix astrometry
fixast,img,h0,fwhm=5,extendbox=10,catalog=catalog
; update FITS file header only
modfits,infile+'.fits',0,h0,exten_no=0

end
