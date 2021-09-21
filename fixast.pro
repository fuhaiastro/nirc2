pro fixast,img,h,imgout=imgout,fwhm=fwhm,extendbox=extendbox,snrcut=snrcut,$
	catalog=catalog,fullsol=fullsol
;+
; Purpose:
;	Fix astrometry offset by comparing image with SDSS/GSC2.3/USNO-B1 catalogs
;	Simple correction will only apply offset to the WCS header, while
;	/fullsol will distort and rotate the image
; 
; Example: 
;	fixast,img,h,imgout=imgout,fwhm=5,snrcut=snrcut,extendbox=10,catalog=catalog
;
; History:
;	Written - Hai Fu - July 2013
;	Added fixast keywords in FITS header - Hai Fu - Sep 2013  
; 
;-

if ~keyword_set(catalog) then catalog = 'SDSS-DR9' ; GSC2.3, USNO-B1
if ~keyword_set(fwhm) then fwhm = 10 ; for centroiding
if ~keyword_set(extendbox) then extendbox = 0 ; for centroiding
if ~keyword_set(snrcut) then snrcut = 5 ; above which considered good sources

; add Equinox if not present
if sxpar(h,'Equinox') eq 0 then sxaddpar,h,'Equinox',2000.0

; load reduced image
nx = (size(img))[1]
ny = (size(img))[2]
b = 10 ; pixels from edge
; central RA & Dec
xyad,h,nx/2,ny/2,ra,dec
getrot,h,pa,cdelt ; get PA and plate scale
radius = max([nx,ny])*mean(abs(cdelt))*60. ; arcmin
; download catalog and image and overplot
tmp = QueryVizier(catalog,[ra,dec],radius)
; select only good photometry
;if catalog eq 'SDSS-DR9' then tmp = tmp[where(strtrim(tmp.q_mode,2) eq '+')]
; select primary sample to avoid duplicated sources
if catalog eq 'SDSS-DR9' then tmp = tmp[where(tmp.mode eq 1)]

print,' running Astronomy calibration '
;window,0,xsize=1000,ysize=1000

; 1. crude astrometry solution based on a single source
; show input image
loadct,0
astrim,asinh_scale(img,img,max=10),tit=sxpar(h,'object')

; overlay catalog
adxy,h,tmp.ra_icrs,tmp.de_icrs,x,y
plots,x,y,psym=8,syms=2,noclip=0,color=cgcolor('green')
; mark a source to align catalog to
print,' Click on a source ... '
print,' If nothing is visible then click to the left of the image '
cursor,xcur,ycur,/wait,/up
; if there is no visible source
if xcur lt 0 then begin
	print,' Error: none of the catalog sources is visible, do nothing'
	imgout = img
	goto,theend
endif
CNTRD,img,xcur,ycur,xcen0,ycen0,fwhm,extendbox=extendbox,/silent
if sqrt((xcur-xcen0)^2+(ycur-ycen0)^2) gt 2*fwhm then begin ; drifted too far
	xcen0 = xcur
	ycen0 = ycur
endif
plots,xcen0,ycen0,psym=6,symsize=2,color=cgcolor('red')
print,' Click on the corresponding object in the catalog ... '
cursor,xcur,ycur,/wait,/up
mindis = min(sqrt((x-xcur)^2+(y-ycur)^2),idx)
plots,x[idx],y[idx],psym=8,syms=2,color=cgcolor('red')
; update Astrometry
sxaddpar,h,'crpix1',xcen0+.5
sxaddpar,h,'crpix2',ycen0+.5
sxaddpar,h,'crval1',tmp[idx].ra_icrs
sxaddpar,h,'crval2',tmp[idx].de_icrs

; remove catalog objects outside of the image
adxy,h,tmp.ra_icrs,tmp.de_icrs,x,y ; predicted positions
s = where(x gt b and x lt nx-b and y gt b and y lt ny-b)
tmp = tmp[s]

; 2. refine astrometry with the S/N weighted mean offset from all objects
loadct,0
astrim,asinh_scale(img,img,max=10),tit=sxpar(h,'object')
; show all sources in the catalog within the image
adxy,h,tmp.ra_icrs,tmp.de_icrs,x,y ; predicted positions
plots,x,y,psym=8,syms=2,color=cgcolor('blue')
; recenter
CNTRD,img,x,y,xcen,ycen,fwhm,extendbox=extendbox,/silent
; test if drifted too far
for i=0,n_elements(x)-1 do begin
	if sqrt((x[i]-xcen[i])^2+(y[i]-ycen[i])^2) gt 2*fwhm then begin
		xcen[i] = x[i]
		ycen[i] = y[i]
	endif
endfor
; aperture photometry to get S/N for weighting the average offsets
snr = x*0
for i=0,n_elements(x)-1 do begin
	if xcen[i] eq 0 or ycen[i] eq 0 then continue
	phot,img,xcen[i],ycen[i],fwhm,5*fwhm,flux=flux
	snr[i] = flux[0]/flux[1]
endfor
forprint,x,y,xcen,ycen,snr
; include only successful centroids and above certain S/N
s = where(xcen gt 0 and ycen gt 0 and snr gt snrcut)
plots,x[s],y[s],psym=8,syms=2,color=cgcolor('red')
plots,xcen[s],ycen[s],psym=1,syms=2,color=cgcolor('green')
; weighted mean based on S/N
dx = total( (xcen[s]-x[s]) * (snr[s]/total(snr[s])) )
dy = total( (ycen[s]-y[s]) * (snr[s]/total(snr[s])) )
; update WCS header with median difference
sxaddpar,h,'crpix1',xcen0+.5+dx
sxaddpar,h,'crpix2',ycen0+.5+dy
; update fixast related keywords
sxaddpar,h,'fixast_catalog',catalog,$
	' Catalog used for astrometry calibration'
sxaddpar,h,'fixast_nstars',n_elements(s),$
	' Number of sources used for astrometry calibration'
sxaddpar,h,'fixast_method','Simple Shift',$
	' Method used for astrometry calibration'

; pause a bit to examine the result
print,' Click to proceed ...'
cursor,aa,bb,/up

fullsolution:
; 3. full non-linear astrometry solution
;    this will distort the origional image
if ~keyword_set(fullsol) then begin
	imgout = img ; if not full solution then output=input, except that
					 ; the header has been updated
endif else begin
	; plot image again
	loadct,0
	astrim,asinh_scale(img,img,max=10),tit=sxpar(h,'object')
	; show predicted positions
	adxy,h,tmp.ra_icrs,tmp.de_icrs,x,y
	; show observed positions
	CNTRD,img,x,y,xcen,ycen,fwhm,extendbox=extendbox,/silent
	; aperture photometry to get S/N for weighting the average offsets
	snr = x*0
	for i=0,n_elements(x)-1 do begin
		if xcen[i] eq 0 or ycen[i] eq 0 then continue
		phot,img,xcen[i],ycen[i],fwhm,5*fwhm,flux=flux
		snr[i] = flux[0]/flux[1]
	endfor
	; include only successful centroids
	s = where(xcen gt b and xcen lt nx-b and ycen gt b and ycen lt ny-b $
				and snr gt snrcut)
	cat = tmp[s] & xcen = xcen[s] & ycen = ycen[s] & x = x[s] & y = y[s]
	plots,x,y,psym=8,syms=2,color=cgcolor('red')
	plots,xcen,ycen,psym=1,syms=2,color=cgcolor('green')
	; print differences
	print,'Initial: Astrometry Predicted vs. Centroid Positions (pixels):'
	print,mean(xcen-x),median(xcen-x),stddev(xcen-x)
	print,mean(ycen-y),median(ycen-y),stddev(ycen-y)
	
	; input RA & Dec
	ra_opt  = cat.ra_icrs
	dec_opt = cat.de_icrs
	; reference coordinate
	xcref = nx/2.0
	ycref = ny/2.0
	xyad,h,xcref,ycref,raref,decref
	; remove bad sources
	bad = bytarr(n_elements(xcen))
	print,'Click on bad sources:'
	cursor,xc,yc,/wait,/up
	while xc gt 0 do begin
		mindis = min(sqrt((xc-x)^2+(yc-y)^2),idx)
		bad[idx] = 1b
		plots,x[idx],y[idx],psym=7,syms=2,color=cgcolor('blue'),thick=4
		cursor,xc,yc,/wait,/up
	endwhile
	; if more than 4 good sources left, solve for 5 terms
	if n_elements(where(bad eq 0b)) gt 4 then terms = ['CONST','X','Y','XX','YY'] $
		                   else terms = ['CONST','X','Y']
	print,terms
	; update fixast related keywords
	sxaddpar,h,'fixast_nstars',n_elements(where(bad eq 0b)),$
		' Number of sources used for astrometry calibration'
	sxaddpar,h,'fixast_method',strjoin(terms,','),$
		' Method used for astrometry calibration'
	; solve for astrometry
	renorm = max([xcref,ycref])
	astrd2sn,ra_opt*!dtor,dec_opt*!dtor,raref*!dtor,decref*!dtor,xi,eta,/arcsec
	astsolve,(xcen-xcref)/renorm,(ycen-ycref)/renorm,xi,eta,terms,renorm,bad,cxi,ceta
	ininfo = {renormfac:renorm,cxi:cxi,ceta:ceta,terms:terms,$
			  xcref:xcref,ycref:ycref,raref:raref*!dtor,decref:decref*!dtor}
	;; check how good the fit is
	;astrd2xy,ra_opt*!dtor,dec_opt*!dtor,ininfo,x,y,/full
	;print,'Fitted Position vs. Centroid Positions (pixels):'
	;print,mean(xcen-x),median(xcen-x),stddev(xcen-x)
	;print,mean(ycen-y),median(ycen-y),stddev(ycen-y)
	; dewarp image
	pscale = min([sqrt(cxi[1]^2+cxi[2]^2),sqrt(ceta[1]^2+ceta[2]^2)])
	outinfo = {renormfac:renorm,cxi:[cxi[0],-1*pscale,0],ceta:[ceta[0],0,pscale],$
		terms:['CONST','X', 'Y'], prot:0.0,$
		xcref:xcref,ycref:ycref,raref:raref*!dtor,decref:decref*!dtor}
	nx = sxpar(h,'naxis1')
	ny = sxpar(h,'naxis2')
	; Create two arrays that carry the native output pixel coordinates
	outy = indgen(long(nx)*long(ny),/LONG)
	outx = outy mod nx
	outy = temporary(outy)/nx
	; Convert output pixel coordinates to position on sky
	astxy2rd,outx,outy,outinfo,outra,outdec,/FULL
	; Convert sky positions to input pixel coordinates
	astrd2xy,outra,outdec,ininfo,inx,iny,/FULL
	; Create output image
	imgout = fltarr(nx,ny)
	; interpolate output pixels that map onto input array
	zg = where(inx ge 0 and inx lt nx and $
	           iny ge 0 and iny lt ny, countg)
	imgout[zg] = interpolate(img,inx[zg],iny[zg]) ;,cubic=-0.5)
	; insert Astrometry with original reference point
	xcref = nx/2.0
	ycref = ny/2.0
	astxy2rd,xcref,ycref,outinfo,raref,decref,/FULL
	raref /= !dtor
	decref /= !dtor
	putast,h,[[-1,0],[0,1]]*pscale/renorm/3600d, [xcref+1.0,ycref+1.0], $
		[raref,decref], ['RA---TAN','DEC--TAN']
endelse

; 4. check result
loadct,0
astrim,asinh_scale(imgout,imgout,max=10),tit=sxpar(h,'object')
; show predicted positions
adxy,h,tmp.ra_icrs,tmp.de_icrs,x,y
plots,x,y,psym=8,syms=2,color=cgcolor('blue')

theend:
end
