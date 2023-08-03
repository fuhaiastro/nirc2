function shift_sub, image, x0, y0, cubic=cubic, missing=missing
;+
; NAME: 
;	SHIFT_SUB
;
; PURPOSE:
;     	Shift an image with subpixel accuracies with INTERPOLATE
; 	By default, it bilinear interpolates the images;
; 	Set CUBIC to a value between -1 and 0 to use the cubic
;	convolution interpolation method.
;
; CALLING SEQUENCE:
;      	Result = shift_sub(image, x0, y0, [cubic=, missing=])
; 
;	The same as shift(), positive shifts are to the right and up while 
;	left and down shifts are expressed as a negative number.
; 
; EXAMPLE
;	; generate an image with a Gaussian star near center
;	gauss2d,100,100,50,50,5,img
;	; shift the star to the upper right
;	img1 = shift_sub(img,5,5,cubic=-0.5,missing=0)
;	img2 = shift(img,5,5)
; 
; HISTORY
;       9-Sep-2010 - H. Fu - written
;	     1-Apr-2022 - added documentation
;-

;if fix(x0)-x0 eq 0. and fix(y0)-y0 eq 0. then $
;	return, shift(image, x0, y0)

x1 = findgen((size(image))[1])-x0[0]
y1 = findgen((size(image))[2])-y0[0]

if n_elements(cubic) ne 0 then begin
	if n_elements(missing) ne 0 then begin
	   result = interpolate(image,x1,y1,/grid,cubic=cubic,missing=missing)
	endif else begin 
	   result = interpolate(image,x1,y1,/grid,cubic=cubic) 
	endelse
endif else begin 
	if n_elements(missing) ne 0 then begin
	   result = interpolate(image,x1,y1,/grid,missing=missing) 
	endif else begin
	   result = interpolate(image,x1,y1,/grid)
	endelse
endelse 

return,result

end 

