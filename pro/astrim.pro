;+
; Name
;	astrim
;
; Purpose
;	plot images with the option of showing axes of equatorial coordinate
;
; 	compared to tvimage, astrim has the advantage of giving the pixel grids
; 	so that overplotting becomes easier
;-

PRO TAG_SPEC,STRUCT,INDEX,TAG_STRING,TAG_DESCRIPTOR
;
code = ['','B','I','J','F','D','','A','']
tags = strlowcase(tag_names(struct))
nind = n_elements(index)
tag_string = ''
tag_descriptor = ''
for i=0,nind-1 do begin 
    tag_string = [tag_string,tags(index(i))]
    descrn = strtrim(string(n_elements(struct.(index(i)))),2)
    descrt = code(datatype(struct.(index(i)),2))
    tag_descriptor = tag_descriptor+descrn+descrt
    if i lt nind-1 then tag_descriptor = tag_descriptor+','
endfor
tag_string = tag_string(1:*)

END

PRO ASTRIM,IM,HDR,GRID=GRID,NOSCALE=NOSCALE,_EXTRA=EXTRA
;
astro = 1
rgb = 0
if n_elements(grid) eq 0 then grid = 1
case n_params() of

0: begin
    print,' usage: ASTRIM,IM_DISP,HDR_DISP'
    print,'        optional keywords: GRID,NOSCALE,<most graphics keywords>'
    return
   end

1: begin
    astro = 0
    imim = im
   end

2: begin
    imim = im
    hdrim = hdr
    check_fits,imim,hdrim,/update   ;force header to be consistent with image
   end 

else: begin
    print,' usage: ASTRIM,IM_DISP,HDR_DISP'
    print,'        optional keywords: GRID=GRID,<most graphics keywords>'
    return
   end
endcase

hdrdim = size(hdrim)
if hdrdim(0) gt 1 then begin
   print,' ERROR: 2nd parameter must be a FITS header'
   return
endif
imdim = size(imim)
if imdim(0) ne 2 and imdim(0) ne 3 then begin
   print,' ERROR: 1st parameter must be a 2d image'
   return
endif

xdim = float(imdim(1))
ydim = float(imdim(2))

if astro then begin
get_astro:
   gsss = sxpar(hdrim,'PPO1',COUNT=N_ppo1)
   if N_ppo1 EQ 1 then gsss_stdast,hdrim

   ctyp = sxpar(hdrim,'ctype*',count=n_ctype)
   if n_ctype gt 0 and strmid(ctyp(0),0,7) eq 'sky_pix' then begin
      sxaddpar,hdrim,'ctype1','RA---TAN'
      sxaddpar,hdrim,'ctype2','DEC--TAN'
      print,' WARNING: header 1 contains bizarre CTYPE specifications - astrometry may be screwed up!'
   endif
   equi = sxpar(hdrim,'equinox',count=n_equi)
   if n_equi le 0 then begin
      equi = sxpar(hdrim,'epoch',count=n_equi)
      if n_equi le 0 then begin
         print,' WARNING: header 1 contains neither EQUINOX nor EPOCH '+$
               'keywords - astrometry may be screwed up!'
         equinox = '(equinox unknown)'
      endif else begin
         print,' WARNING: header 1 contains no EQUINOX keyword '+$
               '- using EPOCH instead'
         equinox = '(J'+strtrim(string(equi,form='(i4)'),2)+')'
      endelse
   endif else equinox = '(J'+strtrim(string(equi,form='(i4)'),2)+')'

   extast,hdrim,astrim    ; extract astrometry parameters for displayed image
 
   xyad,hdrim,(indgen(3)*(xdim-1)/2)#replicate(1,3),$
              replicate(1,3)#((indgen(3)*(ydim-1)/2)),ra,dec

   ramin = min(ra,iramin) & ramax = max(ra,iramax)
   decmin = min(dec,idecmin) & decmax = max(dec,idecmax)

;   find P.A. of North: 

   adxy,hdrim,(ramax+ramin)/2*[1,1],(decmax+decmin)/2+[-0.1,0.1],x,y
   if x(0) gt x(1) then begin
      pa = atan((y(1)-y(0))/(x(1)-x(0))) + !pi/2
   endif else if x(0) lt x(1) then begin
      pa = atan((y(1)-y(0))/(x(1)-x(0))) - !pi/2
   endif else begin
      pa = 0 + !pi*(y(0) gt y(1))
   endelse

   pa = 180/!pi*pa

   if abs(pa) gt 20 then begin
      print,' image misoriented (north is not up) - please rotate, i.e., type:'
      print,'              HROT,IM,HDR,-1,-1,'+strtrim(string(pa,form='(f7.1)'),2)+',-1,-1,2'
      print,'                   then run ASTRIM again'
      return
   endif
   if (iramin mod 3) lt (iramax mod 3) then ramax = ramax+360.
   ntics = 6 ; Hai Fu 20080205
   if grid gt 0 then ntics = fix(ntics*grid)
   case xdim gt ydim of
   1: begin
        nxtics = ntics
        nytics = 3>fix(nxtics*ydim/xdim)
      end
   0: begin
        nytics = ntics
        nxtics = 3>fix(nytics*xdim/ydim)
      end
   endcase
   raticmin = min(ra(*,0)) & raticmax = max(ra(*,0))
   decticmin = min(dec(0,*)) & decticmax = max(dec(0,*))
   xticsize = xdim/nxtics
   tics,raticmin,raticmax,xdim-1,xticsize,ra_incr,/ra
   nxtics = round(xdim/xticsize)+2
   tic_one,raticmin,xticsize,ra_incr,ra_1,ratic_1,/ra
   ticlabels,ra_1,nxtics,ra_incr,xticlab,/ra,delta=1
   yticsize = ydim/nytics
   tics,decticmin,decticmax,ydim-1,yticsize,dec_incr
   nytics = round(ydim/yticsize)+2
   tic_one,decticmin,yticsize,dec_incr,dec_1,dectic_1
   ticlabels,dec_1,nytics,dec_incr,yticlab,delta=1
   if nxtics gt 30 or nytics gt 30 then begin
      print,' ERROR: number of tick marks to large; reset GRID keyword to smaller value'
      return
   endif

endif

x = findgen(101)*xdim/100 
y = findgen(101)*ydim/100 

imimmax = float(max(imim,min=imimmin))

psave = !p
xsave = !x
ysave = !y

!x.tickname = replicate(' ',30)
!y.tickname = replicate(' ',30)
!x.ticks = 1
!y.ticks = 1
!x.minor = 1
!y.minor = 1
!x.tickv = minmax(x)
!y.tickv = minmax(y)
!p.ticklen = 0.01 ;0.0001
!x.title = 'X'
!y.title = 'Y'

if astro then begin
   !x.title = 'Right Ascension '+equinox
   !y.title = 'Declination '+equinox
   xtickv = xdim-1 - (ratic_1+indgen(nxtics)*xticsize)
   ix = where(xtickv ge 0 and xtickv le xdim)
   !x.tickv = xtickv(ix)   
   !x.tickname = xticlab(ix)
   !x.ticks = n_elements(ix) - 1
   ytickv = dectic_1+indgen(nytics)*yticsize
   iy = where(ytickv ge 0 and ytickv le ydim)
   !y.tickv = ytickv(iy)
   !y.tickname = yticlab(iy)
   !y.ticks = n_elements(iy) - 1
endif   

; if !d.name eq 'PS' then device,xsize=18,xoff=2,ysize=27.5,yoff=1,bit=8
; if !d.name eq 'PS' then device,xsize=18,xoff=2.5,ysize=29.5,yoff=1,bit=8
; the following assumes A4 paper!
;if !d.name eq 'PS' then $     
;	if xdim/ydim gt 20/29. then device,bit=8,xsize=18,xoff=1,$
;                     		       ysize=18.*ydim/xdim,$
;				       yoff=14.5-9*ydim/xdim $
;	else $
;;				    device,bit=8,ysize=27,yoff=1,$
;				    device,bit=8,ysize=25,yoff=1,$
;				       xsize=25*ydim/xdim,$
;				       xoff=10-12.5*ydim/xdim
if !d.name eq 'X' and !d.window eq -1 then $
                  window,xsize=600*sqrt(xdim/ydim),ysi=600*sqrt(ydim/xdim)

nxplot = !p.multi(1)>1
nyplot = !p.multi(2)>1
chscale = 1
if !p.charsize gt 0 then chscale = !p.charsize
if (nxplot>nyplot) gt 2 then chscale = 0.5*chscale

xoff = (yoff = 0.)
xrsize = !d.x_vsize/nxplot
yrsizemax = !d.y_vsize/nyplot
if total(psave.region) ne 0 then begin
   xy = convert_coord(psave.region[[0,2]],psave.region[[1,3]],/to_device,/normal)
   xrsize = xy[0,1]-xy[0,0]
   yrsizemax = xy[1,1]-xy[1,0]
endif
xpsize = xrsize-total(!x.margin)*!d.x_ch_size*chscale 
ypsize = xpsize*ydim/xdim
yrsize = ypsize+total(!y.margin)*!d.y_ch_size*chscale
yoff = 0.5*(!d.y_vsize-nyplot*yrsize)
if yrsize gt yrsizemax then begin 
   yrsize = yrsizemax                       
   ypsize = yrsize-total(!y.margin)*!d.y_ch_size*chscale
   xpsize = ypsize*xdim/ydim
   xrsize = xpsize+total(!x.margin)*!d.x_ch_size*chscale
   xoff = 0.5*(!d.x_vsize-nxplot*xrsize)
   yoff = 0.
endif
psize = convert_coord(xpsize,ypsize,/device,/to_normal)
rsize = convert_coord(xrsize,yrsize,/device,/to_normal)
offset = convert_coord(xoff,yoff,/device,/to_normal)
mplot = shift(reverse(indgen(nxplot*nyplot)),1)
iplot = mplot(!p.multi(0))
xr = ((iplot mod nxplot)+[0,1])*rsize(0) + offset(0) 
yr = 1-(iplot/nxplot+[1,0])*rsize(1)     - offset(1)
if total(psave.region) ne 0 then begin
   xr = psave.region[0]+[0,xr[1]-xr[0]]
   yr = psave.region[1]+[0,yr[1]-yr[0]]
endif
!p.region = [xr(0),yr(0),xr(1),yr(1)]
;if total(psave.region) eq 0 then !p.region = [xr(0),yr(0),xr(1),yr(1)]

plot,x,y,/xst,/yst,/nodata,_extra=extra

tvlct,r,g,b,/get

if !d.name eq 'PS' THEN BEGIN
   if imdim[0] eq 2 then begin    ; simple 2d image
      r[0] = 0                  ;reset lowest color index to BLACK
      g[0] = 0
      b[0] = 0
      tvlct,r,g,b
      ;top = !d.n_colors-2
		top = !d.n_colors-13 ; Hai Fu, reserve colors for lines
      if keyword_set(noscale) then top = max(imim)
      ;tv,!d.n_colors-1-bytscl(imim,top=top),!x.window(0),!y.window(0),$
      tv,!d.n_colors-12-bytscl(imim,top=top),!x.window(0),!y.window(0),$
                             xsize=!x.window(1)-!x.window(0),$
                             ysize=!y.window(1)-!y.window(0),/norm
   endif else begin               ; image is taken to be r,g,b cube
      loadct,0
      device,/color
      if keyword_set(noscale) then top = max(imim)
      tv,bytscl(imim,top=top),!x.window(0),!y.window(0),$
                             xsize=!x.window(1)-!x.window(0),$
                             ysize=!y.window(1)-!y.window(0),/norm,true=3
   endelse
endif else begin
   px = !x.window * !d.x_vsize
   py = !y.window * !d.y_vsize
   sx = px(1)-px(0)+1
   sy = py(1)-py(0)+1
   top = !d.n_colors-2
   if keyword_set(noscale) then top = max(imim)
   if imdim[0] eq 2 then begin    ; simple 2d image
      imax = (!d.n_colors-1)<255
      if r(!p.background)+g(!p.background)+b(!p.background) eq 0 then begin
         r(imax) = 255          ;background is black
         g(imax) = 255          ;reset highest color index to White
         b(imax) = 255
      endif else begin
         r(imax) = 0       ;reset highest color index to Black
         g(imax) = 0
         b(imax) = 0
      endelse
      tvlct,r,g,b
      tv,congrid(bytscl(imim,top=top),sx,sy),px(0),py(0)
   endif else begin               ; image is taken to be r,g,b cube
      loadct,0
      tv,congrid(bytscl(imim,top=top),sx,sy,3),px(0),py(0),true=3
   endelse
endelse

if astro and grid gt 0 then begin
   ra = ra_1
   while ra gt ramin do ra = ra - 15*ra_incr/60.
   while ra lt ramax do begin
       ra = ra+15*ra_incr/60.
       yg = findgen(ydim+1)
       xg = cons_ra(ra,yg,astrim)                      
       if !d.name eq 'PS' then col=!d.n_colors/2 else col=cgcolor('light gray')
       oplot,xg,yg,color=col
   endwhile

   dec = dec_1
   while dec gt decmin do dec = dec - dec_incr/60.
   while dec lt decmax do begin
       dec = dec+dec_incr/60.
       xg = findgen(xdim+1)
       yg = cons_dec(dec,xg,astrim)
       if !d.name eq 'PS' then col=!d.n_colors/2 else col=cgcolor('light gray')
       oplot,xg,yg,color=col
   endwhile
endif

plot,x,y,/xst,/yst,/nodata,/noerase,_extra=extra

!x.tickname = xsave.tickname
!y.tickname = ysave.tickname
!x.ticks = xsave.ticks
!y.ticks = ysave.ticks
!x.minor = xsave.minor
!y.minor = ysave.minor
!x.tickv = xsave.tickv
!y.tickv = ysave.tickv
!p.ticklen = psave.ticklen
!p.region = psave.region
!p.position = psave.position
!x.title = xsave.title
!y.title = ysave.title

END
