;+
; NAME
;  nirc2findtt
; 
; PURPOSE:
;  Find NIRC2 tip-tilt stars for a list of targets given RA & Dec 
;  Uses both SDSS and USNO-B from QueryVizier
;
; CALLING SEQUENCE:
;  nirc2findtt,id, RA, Dec, prefix=prefix,rmin=rmin,rmax=rmax,dmin=dmin,dmax=dmax
;
; INPUTS:
;  ID  - Designation (String Array)
;  RA  - RA in degrees (Float Array) 
;  DEC - Dec in degrees (Float Array)
; 
; OPTIONAL INPUTS:
;  prefix - prefix for output files
;  rmin - brightest R magnitudes to consider, default: 8.0
;  rmax - faintest R magnitudes to consider, default: 19.0
;  dmin - minimum distance to target, default:  0.0"
;  dmax - maximum distance to target, default: 80.0". Max distance from
;         NIRC2 chip center is ~62". So anything greater than ~59" requires 
;         an offset (leaving 3" margin).
; 
; OUTPUT:
;  1. Keck starlist for all targets with tip-tilt stars
;  2. IDL structure FITS file with necessary tags
;
; REVISION HISTORY:
;  Sep 2011 - Written, H.Fu
;  Apr 2021 - updated to SDSS DR12, H.Fu
;-

pro nirc2findtt,id,ra,dec,prefix=prefix,rmin=rmin,rmax=rmax,dmin=dmin,dmax=dmax

if n_params() lt 3 then begin
	print,' syntax: nirc2findtt,id(string),ra(deg),dec(deg)'
	print,'   [,prefix=prefix,rmin=rmin,rmax=rmax,dmin=dmin,dmax=dmax]'
	return
endif

if ~keyword_set(prefix) then prefix = 'targets'
; set limits for TTref
;      R < 17, Sep < 60" => Low Risk LGSAO
; 17 < R < 19, Sep < 60" => Medium Risk LGSAO
; NIRC2/Wide is able to do Sep < 80" (leaving 5" to edge on each side)
; OSIRIS is able to do Sep < 72"
if ~keyword_set(rmax) then rmax = 19.0 ; R-band magnitudes
if ~keyword_set(rmin) then rmin =  8.0
if ~keyword_set(dmax) then dmax = 80.0 ; = 60+15*sqrt(2) ; arcsec
if ~keyword_set(dmin) then dmin =  0.0

coord = strtrim(repstr(adstring(ra,dec,2),'  ',' '),2)
forprint,id,coord,f='(a-16,1x,a25,1x,"2000.0 lgs=1")',$
	text=prefix+'.lst',/nocomment

; read in the starlist
in0 = readstarlist(prefix+'.lst',/lgs)
spawn,'rm -f '+prefix+'.lst'
remove_tags,in0,['rah','ram','ras','design','ded','dem','des','equ','vmag'],in
; build TTref structure 
tags =['ra','dec','Rmag','NIRC2_PA','ttrmag','ttstrehl','ttsep','ttPA','ttsource']
values = ['0d','0d','0.','0.','0.','0.','999.9','999.9','""']
ttr = struct_combine(in,mrd_struct(tags,values,n_elements(in)))
; search radius in arcmin
radius = dmax/60*1.1  
; Query SDSS/USNO photometry database 
iso = 40.0 ; LGSAO isokinetic angle, in arcseconds
effwave = 2.124 ; Effective central wavelength
kstrehl = [[10, 0.35], [11, 0.35], [12, 0.35], $ ; K Strehl vs. R mag.
           [13, 0.34], [14, 0.32], [15, 0.30], $ ;  (07/04-09/04 Eng.)
           [16, 0.27], [17, 0.23], [18, 0.18], $
           [19, 0.12]]
rms = kstrehl
rms[1,*] = 2124/(2*!pi)*SQRT(-1*ALOG(kstrehl[1,*])) ; RMS WF error (nm)
; NIRC2 TT ref patrol field
which,'nirc2distort',files=files,/silent
readcol,file_dirname(files[0])+'/NIRC2_ttref_field.txt',xb,yb

for i=0,n_elements(in)-1 do begin
	tmp = QueryVizier('V/147/sdss12',[in[i].raj2000,in[i].dej2000],radius)
	if size(tmp,/type) eq 8 then begin
		; SDSS r --> Rc in Vega, Fukugita+96
		usno_rmag = tmp.rmag-0.154254-0.0889*(tmp.gmag-tmp.rmag)
		b_v = (tmp.gmag-tmp.rmag+0.23)/1.05
		ttr[i].ttsource = 'SDSS'
	endif else begin ; if not in SDSS, query USNO
		usno = QUERYVIZIER('USNO-B1',[in[i].raj2000,in[i].dej2000],radius)
		if size(usno,/type) eq 8 then begin
			; require at least one measurement in B and R
			s = where((finite(usno.r1mag) or finite(usno.r2mag)) $
			      and (finite(usno.b1mag) or finite(usno.b2mag)),ct)
			if ct ge 1 then usno = usno[s] else continue
	   	; record magnitudes
			usno_rmag = usno.r1mag
			usno_bmag = usno.b1mag
			for j=0,n_elements(usno)-1 do begin
				usno_rmag[j] = finite(usno[j].r1mag)?usno[j].r1mag:usno[j].r2mag
				usno_bmag[j] = finite(usno[j].b1mag)?usno[j].b1mag:usno[j].b2mag
			endfor
			s = where(finite(usno.r1mag) and finite(usno.r2mag),ct)
			if ct ge 1 then usno_rmag[s] = 0.5*(usno[s].r1mag + usno[s].r2mag)
			s = where(finite(usno.b1mag) and finite(usno.b2mag),ct)
			if ct ge 1 then usno_bmag[s] = 0.5*(usno[s].b1mag + usno[s].b2mag)
			;forprint,usno_rmag,usno.r1mag,usno.r2mag
      	b_v = 0.556*(usno_bmag - usno_rmag) ; B-V color magnitude 
			tmp = usno
			ttr[i].ttsource = 'USNO'
		endif else continue ; if not in USNO then jump to next object
	endelse
	; Target Rc magnitude in Vega 
	if min(sphdist(tmp.ra_ICRS,tmp.de_ICRS,$
		in[i].raj2000,in[i].dej2000,/deg)*3600,idx) lt 1 $
	then ttr[i].rmag = usno_rmag[idx]
	; compute Strehl ratios
	usno_sep = sphdist(tmp.ra_ICRS,tmp.de_ICRS,$
		ttr[i].raj2000,ttr[i].dej2000,/deg)*3600
	rms_mag = INTERPOL(rms[1,*], rms[0,*], usno_rmag)
	rms_iso = ((effwave*1e3*usno_sep)/(2*!pi*iso))^(5/6.) ; Hardy Eq7.61
	rms_tot = SQRT(rms_mag^2 + rms_iso^2)
	strehl = EXP(-1*(2*!pi*rms_tot/(effwave*1e3))^2) ; Est. Strehl
	; for non-stellar sources and others, set strehl to zero
	if ttr[i].ttsource eq 'SDSS' then begin
		idx = where(usno_sep gt dmax or tmp.class ne 6 or usno_rmag gt rmax)
		if idx[0] ge 0 then strehl[idx] = 0 
	endif else begin
		idx = where(usno_sep gt dmax or usno_rmag gt rmax)
		if idx[0] ge 0 then strehl[idx] = 0
	endelse
	; if there are TT refs within 59", invalid TT refs beyond 59"
	idx1 = where(usno_sep gt 59. and strehl gt 0,ct1)
	idx2 = where(usno_sep le 59. and strehl gt 0,ct2)
	if ct2 ge 1 and ct1 ge 1 then strehl[idx1] = 0
	; save best TT ref
	junk = max(strehl,idx)
	if junk gt 0 then begin
		ttr[i].ttraj2000 = tmp[idx].ra_ICRS
		ttr[i].ttdej2000 = tmp[idx].de_ICRS
		ttr[i].ttrmag = usno_rmag[idx]
		ttr[i].ttstrehl = strehl[idx]
		ttr[i].ttcomment = $
			'R='+strtrim(string(usno_rmag[idx],f='(f5.1)'),2)+$
			' Sep='+strtrim(string(usno_sep[idx],f='(f5.1)'),2)+$
			' B-V='+strtrim(string(b_v[idx],f='(f5.2)'),2)+$
			' S='+strtrim(string(strehl[idx],f='(f5.2)'),2)
	endif
	if i mod 5 eq 0 then print,i+1,' / ',n_elements(in) ; show status
endfor
; fill in additional tags so that we can use starstring.pro
ttr.ra = ttr.raj2000
ttr.dec = ttr.dej2000
; compute separation (arcsec) and PA (deg)
ttsep=sphdist(ttr.raj2000,ttr.dej2000,ttr.ttraj2000,ttr.ttdej2000,/deg)*3600
posang,1,ttr.raj2000/15d,ttr.dej2000,ttr.ttraj2000/15d,ttr.ttdej2000,ttpa
s = where(ttsep lt 100)
ttr[s].ttsep = ttsep[s]
ttr[s].ttpa = ttpa[s]
ttr[s].ttref = 'tt'+in[s].targ
; determine NIRC2 PA to fit TTref w/i the patrol field 
; (62.4" max distance from chip center)
for i=0,n_elements(ttr)-1 do begin
	x = ttr[i].TTsep*sin(ttr[i].TTpa*!dtor)*(-1)
	y = ttr[i].TTsep*cos(ttr[i].TTpa*!dtor)
   if ttr[i].ttsep gt 100 then continue
	if inside(x,y,xb,yb) then pa=0 else $
		if ttr[i].TTsep gt 55 $
			then pa=ttr[i].TTpa-144 $ ; lowerleft corner
			else pa=ttr[i].TTpa-50    ; upperleft corner
	if pa lt 0 then pa+=360
	if pa lt 0 then pa+=360
	ttr[i].NIRC2_PA = pa
endfor
; add NIRC2 PA to comments
ttr.comment += ' R='+strtrim(string(ttr.rmag,f='(f6.1)'),2)+$
	       ' PA='+strtrim(string(ttr.nirc2_pa,f='(f6.1)'),2)

; Save IDL structure
mwrfits,ttr,prefix+'.fits',/create
; Save starlist
s = where(ttr.ttrmag ne 0)
str1 = starstring(ttr[s],str2)
str3 = ''
for i=0,n_elements(str1)-1 do str3 = [str3,str1[i],str2[i]]
forprint,str3[1:*],textout=prefix+'.lst',/nocomment

return

end
