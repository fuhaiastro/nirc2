pro phot,img,xc,yc,apr,dis,flux=flux,sky=sky,debug=debug
;+
; 
; aperture photometry
; first do photometry on central source
; then do photometry of the same aperture on an annulus at dis
;
; Output
; 	flux - sky-subtracted flux of the source and its error based on the
;		standard deviation of the apertures placed around the source
;	sky - fluxes of the sky apertures, no background subtraction
;
; History
;	Written by Hai Fu, Aug 2013
;-

;window,0,xsize=500,ysize=500
;multiplot,[1,1],/square
;plot,[0,200],[0,200],/nodata,/xs,/ys

apr *= 1.0 ; make them real numbers
dis *= 1.0 
if keyword_set(debug) then begin
	astrim,asinh_scale(img,img)
	tvcircle,apr,xc,yc,/data
	tvcircle,dis,xc,yc,/data,lines=1
endif
aper,img,xc,yc,flux,ferr,sky,skyerr,1.0,apr,dis+[-0.5,0.5]*apr,[-1,-1],$
	/exact,/flux,/nan,/silent

ang = atan(apr/dis)*2.1 ; in radian
for deg = 0.,2*!PI-ang,ang do begin
	x = xc + dis*cos(deg)
	y = yc + dis*sin(deg)
	if keyword_set(debug) then tvcircle,apr,x,y,/data
	aper,img,x,y,f,df,s,serr,1.0,apr,-1,[-1,-1],$
		/exact,/flux,/nan,/silent,setsky=0.0
	sky = [sky,f]
endfor
sky = sky[1:*]
flux = [flux,robust_sigma(sky)]
if flux[1] lt 0 then flux[1] = stddev(sky)

end
