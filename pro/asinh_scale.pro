FUNCTION ASINH_SCALE,IM,SCALE_IM,MINV=minv,MAXV=maxv,LINEAR=linear,$
	sclpars=sclpars,ncolors=ncolors, skyval=skyval, skynoi=skynoi
;+
; NAME
;   ASINH_SCALE
;
; PURPOSE
;   Scale an 2D image to a byte array with asinh or linear scale. 
;   Asinh scaling from ATV.pro
; 
; SYNTAX 
;   scaled = ASINH_SCALE(img, SCALE_img, MINV=, MAXV=, /LINEAR, $
;	   sclpars=, ncolors=, skyval=, skynoi=)
;
;-

if ~keyword_set(ncolors) then $
	ncolors = !d.table_size-1 

; scaling parameters
; note: skynoi could be zero, so use n_elements instead of keyword_set
if n_elements(skynoi) eq 0 or n_elements(skyval) eq 0 then $
	sky,SCALE_IM,skyval,skynoi,/nan,/silent
if ~keyword_set(minv) then min_value = skyval-(2.0*skynoi) $
                      else min_value = skyval+minv*skynoi
if ~keyword_set(maxv) then max_value = max(SCALE_IM,/nan) $
                      else max_value = skyval+maxv*skynoi

if ~keyword_set(linear) then begin
	asinh_beta = skynoi
	scaled_image = bytscl(asinh((im-min_value)/asinh_beta),min=0, $
   		max=asinh((max_value-min_value)/asinh_beta),/nan,top=ncolors-1)
	sclpars = [skyval,skynoi,min_value,max_value,asinh_beta]
endif else begin
   	scaled_image = bytscl(im,min=min_value,max=max_value,/nan,top=ncolors-1) 
	sclpars = [skyval,skynoi,min_value,max_value]
endelse

return,scaled_image

END
