pro nirc2_link,dir,n1,n2,outdir,suffix=suffix,clear=clear
;+
; NAME
; 	NIRC2_LINK
;
; PURPOSE
; 	create symbolic links from frames n1 to n2 in dir to outdir
;
; SYNTAX
; 	nirc2_link,'../raw_lin',200,201,'1445+3032',suffix='luci2.20160122.',/clear
;
; HISTORY
;	written by HF to organize NIRC2 images
;	expanded to be able to organize FITS file from other instruments
;
;-
if ~file_test(outdir) then spawn,'mkdir '+outdir
if keyword_set(clear) then spawn,'del '+outdir+'/*.fits'
if ~keyword_set(suffix) then suffix = 'n'
cd,outdir
for i=n1,n2,1 do begin
	file = '../'+dir+'/'+suffix+string(i,f='(i04)')+'.fits'
	if ~file_test(file) then begin
		print,file+' does not exist!'
		continue
	endif
	spawn,'ln -s '+file+' .'
endfor
cd,'..'

end
