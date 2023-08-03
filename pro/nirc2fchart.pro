;+
; Finder Charts for NIRC2 Observations
; Requires Keck starlist format target list
; Also checks magnitudes from 2MASS
;-
pro nirc2fchart,starlist,scale=scale,npix=npix,overwrite=overwrite

; load starlist 
targs = readstarlist(starlist,/lgs)
targs.targ = strtrim(targs.targ,2)

dir_sdss = 'sdss_fchart/'
if ~file_test(dir_sdss) then spawn,'mkdir '+dir_sdss
dir = 'nirc2_fchart/' ; save Finder Charts to
if ~file_test(dir) then spawn,'mkdir '+dir

if ~keyword_set(scale) then scale = 0.15 ; arcsec/pixel
if ~keyword_set(npix) then   npix = 1001 ; # of pixels
width = scale*npix/60.0 ; arcmin per side
xr = [-1,1]*npix*scale/2 ; size of the downloaded SDSS image

; TT ref patrol field
which,'nirc2distort',files=files,/silent
readcol,file_dirname(files[0])+'/NIRC2_ttref_field.txt',xb,yb

for i=0,n_elements(targs)-1 do begin
	; skip if set not to overwite
	if file_test(dir+targs[i].targ+'.png') and ~keyword_set(overwrite) then continue
	; take PA, offset (from TT to chip center) from starlist
	PA = float(strmid(targs[i].comment,strpos(targs[i].comment,'PA=')+3,5))
	offset = -1 ; default is auto-determine offset
	if strpos(targs[i].comment,'off=') ge 0 then $
		offset = float(strmid(targs[i].comment,strpos(targs[i].comment,'off=')+4,4))
	; setup window
	window,1,xsize=800,ysize=800 
	multiplot,[1,1],/square
	coord = repstr(strtrim(adstring(targs[i].raj2000,targs[i].dej2000),2),' ',':')
	tit = targs[i].targ+' '+repstr(coord,'::',' ');+$
	;	' PA='+strtrim(string(PA,f='(f6.1)'),2)+$
	;	' '+targs[i].ttcomment
	; download SDSS image
	if ~file_test(dir_sdss+targs[i].targ+'.jpg') then begin
		QueryURL='http://skyservice.pha.jhu.edu/DR10/ImgCutout/'+$
			'getjpeg.aspx?ra='+strtrim(string(targs[i].raj2000,f='(f12.6)'),2)+$
			'&dec='+strtrim(string(targs[i].dej2000,f='(f12.6)'),2)+$
			'&scale='+strc(scale)+'&width=1001&height=1001&opt=G&query='
		Result = webget(QueryURL,copyfile=dir_sdss+targs[i].targ+'.jpeg')
		; convert BMP format to JPG
		spawn,'convert '+dir_sdss+targs[i].targ+'.jpeg '+dir_sdss+targs[i].targ+'.jpg'
		spawn,'rm -f '+dir_sdss+targs[i].targ+'.jpeg'
	endif
	read_jpeg,dir_sdss+targs[i].targ+'.jpg',im
	; download DSS images
	if (size(im))[0] eq 2 then begin
		if ~file_test(dir_sdss+targs[i].targ+'.gif') then begin
			catalog = 'phase2_gsc2'
			QueryURL='http://archive.stsci.edu/cgi-bin/dss_search?v='+catalog+$
				'&r='+strtrim(string(targs[i].raj2000,f='(f12.6)'),2)+$
				'&d='+strtrim(string(targs[i].dej2000,f='(f12.6)'),2)+$
				'&e=J2000&h='+strtrim(string(width,f='(f12.3)'),2)+$
				'&w='+strtrim(string(width,f='(f12.3)'),2)+'&f=gif&c=none&fov=NONE&v3='
			Result = webget(QueryURL,copyfile=dir_sdss+targs[i].targ+'.gif')
		endif
		spawn,'convert '+dir_sdss+targs[i].targ+'.gif '+dir_sdss+targs[i].targ+'.jpg'
		read_jpeg,dir_sdss+targs[i].targ+'.jpg',im
	endif
	; display image
	tvimage,im,true=((size(im))[0] eq 3)?1:0,position=!p.position
	plot,[1,1],[1,1],/nodata,xr=xr,yr=xr,/xs,/ys,/noerase,$
		xtickinterval=20,ytickinterval=20,xminor=4,yminor=4,$
		xticklen=0.01,yticklen=0.01,xthick=2,ythick=2,$
		position=!p.position,tit=tit
	; mark TTref star
	TTsep = sphdist(targs[i].raj2000,targs[i].dej2000,targs[i].ttraj2000,targs[i].ttdej2000,/deg)*3600d
	posang,1,targs[i].raj2000/15d,targs[i].dej2000,targs[i].ttraj2000/15d,targs[i].ttdej2000,ttpa
	x = TTsep*sin(TTpa*!dtor)*(-1)
	y = TTsep*cos(TTpa*!dtor)
	tvbox,6,x,y,/data,lines=0,color=cgcolor('green')
	; label TTref
	ttstr = targs[i].ttcomment ;+$ ;'Sep='+strtrim(string(ttsep,f='(f4.1)'),2)+$
		;' PA='+strtrim(string(ttpa,f='(f6.1)'),2) ;+$
		;' R='+string(targs[i].ttrmag,f='(f4.1)')
	;ttstr = 'TTPA='+strtrim(string(ttpa,f='(f6.1)'),2)
	xyouts,x,y-8,ttstr,/data,align=0.5,chars=1.0
	; dX & dY are offset from TTref to chip center in instrument coordinates
	if offset lt 0 then begin
		if TTsep gt 59 then begin
			; if TTsep > 59", then by default move TT to lower left corner
			offset = 59.0
		endif else begin
			; if TTsep < 59", then put target at chip center
			offset = TTsep
		endelse
	endif
	dX = offset*sin((ttPA-pa)*!dtor)
	dY = -1*offset*cos((ttPA-pa)*!dtor)
	; xoff & yoff are offset from target to chip center in sky coordinates
	xoff = -1.0*(TTsep-sqrt(dX^2+dY^2))*sin(TTpa*!dtor)
	yoff = (TTsep-sqrt(dX^2+dY^2))*cos(TTpa*!dtor)
	;; mark target
	;tvcircle,5,0,0,/data,lines=1
	;; xx, yy are in pixel coordinates
	;xx = 512.+(TTsep-sqrt(dX^2+dY^2))/0.04*sin((TTpa-pa)*!dtor)
	;yy = 512.-(TTsep-sqrt(dX^2+dY^2))/0.04*cos((TTpa-pa)*!dtor)
	;xyouts,0,5,string(round(xx),f='(i4)')+string(round(yy),f='(i4)'),align=.5,chars=1.5
	; label the offset amounts from the TTref to the chip center
	xyouts,0,60,'dX='+string(dX,f='(f5.1)')+' dY='+string(dY,f='(f5.1)'),align=.5,chars=1.5
	; show vignetted patrol field
	; counter-clockwise rotate TTref patrol field in degrees
	tmp = [[xb],[yb]]
	for j=0,n_elements(xb)-1 do $
		tmp[j,*] = [[cos(pa*!dtor),sin(pa*!dtor)],[-1*sin(pa*!dtor),cos(pa*!dtor)]]#[xb[j],yb[j]]
	oplot,[tmp[*,0],tmp[0,0]]+xoff,[tmp[*,1],tmp[0,1]]+yoff,lines=3
	; patrol field is slightly off center
	x0y0 = [[cos(pa*!dtor),sin(pa*!dtor)],[-1*sin(pa*!dtor),cos(pa*!dtor)]]#[-1.01,-2.58]
	tvbox,[96.22,95.30],x0y0[0]+xoff,x0y0[1]+yoff,lines=1,/data,angle=-1*pa
	; detector
	tvbox,40,xoff,yoff,/data,lines=2,angle=-1*pa,color=cgcolor('royal blue')
	;plots,xoff,yoff,/data,psym=1
	; save screen to PNG file
	write_png,dir+targs[i].targ+'.png',tvrd(true=1)
endfor

end

