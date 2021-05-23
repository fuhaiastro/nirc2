function nirc2distort,img
;+
; NAME
;	NIRC2DISTORT
;
; PURPOSE
;	Distortion correction for NIRC2 Wide Camera based on 15 dithered
;	images of M92. Images were taken by A. Stockton on 2011 Aug 28 with
;	300 coadds of 0.2 s. The solution was obtained by fitting 5th order 
;	2D Polynomials with SFIT.
;
; HISTORY
;	H. Fu  - Written - 11/02/2012
;
;-
; load distortion solutions from M92 (data courtesy of A. Stockton)
; following Yelda+12, pixel values are the required shift (in pixels) to undistort an image.
which,'nirc2distort',files=files,/silent
dir = file_dirname(files[0])
xdistort = mrdfits(dir+'/nirc2_wide_X_distortion.fits',/silent)
ydistort = mrdfits(dir+'/nirc2_wide_Y_distortion.fits',/silent)
; generate pixel grids in observed frame
nx = 1024l
ny = 1024l
xo = findgen(nx) # replicate(1., ny)
yo = replicate(1.,nx) # findgen(ny)
; distortion corrected positions
xt = xo + xdistort
yt = yo + ydistort
; undistort by interpolation
img2 = warp_tri(reform(xt,nx*ny),reform(yt,nx*ny),reform(xo,nx*ny),reform(yo,nx*ny),img)
return,img2

end
