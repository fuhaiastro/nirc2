PRO NIRC2REDUX,rawfiles,output=output,caldir=caldir,$
	fltfile=fltfile,bpmfile=bpmfile,drkfile=drkfile,skyfile=skyfile,bottom=bottom,$
	fwhm=fwhm,rmin=rmin,rmax=rmax,snr_cut=snr_cut,rejsig=rejsig,niter=niter,$
	initial=initial,recenter=recenter,ngs=ngs,exam=exam,undistort=undistort,reset=reset
;+
; NAME:
;   NIRC2REDUX
; PURPOSE:
;   Reduce NIRC2/OSIRIS Imaging data. 
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;   nirc2redux,'lock01_kp/n*.fits',output=,caldir=,fltfile='flat_Kp_W_10jun',$
;            bpmfile='badpix_12oct',drkfile='drk_16_4_2_10jun',$
;	          /initial,/recenter,/ngs,/exam,/undistort,/reset
;
; INPUTS:
;   rawfiles - e.g., 'lock01_kp/n*.fits'
;
; OUTPUT:
;   Reduced image - e.g., 'lock01_kp.fits'
;      Extensions - 1. intensity map w/ bspline subtraction (e-/s)
;                   2. noise map (e-/s)
;                   3. exposure map (s)
;                   4. intensity map w/o bspline subtraction (e-/s) 
;
; OPTIONAL INPUTS:
;   output - output file name, if unspecified, use directory name
;   fltfile - flat field, doesn't have to be normalized
;   bpmfile - bad pixel map
;   drkfile - dark current map; default is zero
;   skyfile - separate background only map
;   fwhm - default 5 pixels
;   rmin - default 2xFWHM
;   rmax - default 10xFWHM
;   snr_cut - default 1.2, S/N ratio above which a pixel can be excluded
;   rejsig - default 2, a frame would be discarded if its seeing is
;            rejsig*sigma above average (only works in /recenter mode)
;   niter - number of iterations for sky subtraction in Step 4
;   /bottom - if set and in wide camera, subtract a constant sky level for the bottom 1/3
;             of the image instead of using the bspline-fitted 2D sky model 
;   /initial - if set, do not go through the full reduction, just save (x,y) poisitions
;   /recenter - remeasure centroids on the primary source, instead of trusting the header
;               turn on this option when the alignment object is bright enough
;   /ngs - observed in NGS mode or the AOTSX/Y keywords do not make sense. 
;          In this mode, one has to click through the images and mark
;          the primary alignment source. 
;          The script determines if it's in NGS mode from header
;   /exam - examine the centroids of the primary source, click through images
;   /undistort - undistort images based on M92 solutions
;   /reset - if set, delete all intermediate files in the data folder before proceeding
;
; OPERATIONAL NOTES:
;
; HISTORY:
;   Apr 2012 - Written - H. Fu
;   Jun 2012 - Added blind registration based on LGS header - H. Fu
;   Nov 2012 - Added NGS-mode, on-sky distortion solution, and rough WCS - H. Fu
;   Aug 2013 - Enabled discarding frames in /exam mode - H. Fu
;            - Removed /background option, always save it as the 4th extension
;            - Added /initial option
;   Feb 2015 - Automatically look for the directory of distortion
;              solution FITS files - H. Fu
;            - Added /bottom option
;-

if ( N_params() EQ 0 ) then begin
	print,'Syntax - NIRC2REDUX,rawfiles,[output=,caldir= fltfile=,bpmfile=,drkfile=,skyfile=,'
   print,'                  fwhm=,rmin=,rmax=,snr_cut=,rejsig=,catalog=,'
   print,'                  /recenter,/ngs,/exam,/undistort,/reset'
   return
endif

; default parameters
dir = strmid(rawfiles,0,strpos(rawfiles,'/'))
if ~keyword_set(output) then output = dir
if ~keyword_set(caldir) then caldir = '../cal/'
if ~keyword_set(bpmfile) then bpmfile = 'badpix_12oct' 
if ~keyword_set(fltfile) then fltfile = 'flat_Ks_W_11apr'
if ~keyword_set(drkfile) then drkfile = 'zero' ; default dark is zero image, i.e., no dark subtraction
if ~keyword_set(fwhm)    then fwhm = 5
if ~keyword_set(rmin)    then rmin = 2*fwhm
if ~keyword_set(rmax)    then rmax = 10*fwhm
if ~keyword_set(snr_cut) then snr_cut = 1.2
if ~keyword_set(bkspace) then bkspace = 50
if ~keyword_set(rejsig)  then rejsig = 2.0 ; 2-sigma rejection
if ~keyword_set(niter)  then niter = 3 ; 3 iterations to update mask
sepsky = keyword_set(skyfile)
if keyword_set(reset) then spawn,'rm -f '+dir+'/*.txt'

spawn,'date',begin_time
loadct,0
;==============================
; Step 0:
; LOAD ALL IMAGES INTO ARRAYS
;==============================
; badpix mask, 0 - bad pixel, non-0 - good pixel
bpm = 1b-mrdfits(caldir+bpmfile+'.fits',0,/silent)
; flat field
flat0 = mrdfits(caldir+fltfile+'.fits',0,/silent)
fixpix,flat0,bpm,flat1,/silent ; fix bad pixels
flat = flat1/median(flat1) ; re-normalize
; dark image
dark = mrdfits(caldir+drkfile+'.fits',0,/silent)
; load images in DN (IMGS) and subtract dark and readin offsets
ins = file_search(rawfiles,count=ct)
print,'Input: ',ct,' frames'
imgs = mrdfits(ins[0],0,h0,/silent)*1.0-dark
aotsx = sxpar(h0,'AOTSX')
aotsy = sxpar(h0,'AOTSY')
camera = strtrim(sxpar(h0,'CAMNAME'),2)
if (camera EQ 'narrow') then $ 
   pixscale = 0.009942d $ ; arcsec/pix narrow camera
else $
   pixscale = 0.039686d ; arcsec/pix wide camera
for i=1,ct-1 do begin
	imgs = [[[imgs]],[[mrdfits(ins[i],0,h,/silent)-dark]]]
   aotsx = [aotsx,sxpar(h,'AOTSX')]
   aotsy = [aotsy,sxpar(h,'AOTSY')]
endfor
; test if it is observed in NGS mode
if strtrim(sxpar(h0,'AOFCLGCT')) eq 'OFF' then ngs = 1
if keyword_set(ngs) then recenter=1 ; if NGS then always recenter on bright sources
; convert offset in header to pixels
if ~keyword_set(ngs) then begin
	; coordinate transformation because the detector is 
	; rotated by 0.7+0.252 degrees (PA = rotposn-instangl-0.252)
	angle = (0.7+0.252)/180.0*!pi
	aotsx = aotsx*cos(angle) - aotsy*sin(angle)
	aotsy = aotsx*sin(angle) + aotsy*cos(angle)
	platescale = 0.727 ; mm/arcsec
	xref = aotsx[0] ; relative the first frame
	yref = aotsy[0]
	for i=0, ct-1 do begin    
	  aotsx[i] = (xref - aotsx[i]) *(-1.0)
	  aotsy[i] = (yref - aotsy[i]) *(-1.0)
	  ;convert to pixels
	  aotsx[i] =  aotsx[i] / platescale / pixscale 
	  aotsy[i] =  aotsy[i] / platescale / pixscale 
	  ;print, aotsx[i], aotsy[i] , ' pixels'
	endfor
endif else begin
	aotsx *= 0
	aotsy *= 0
endelse
; load distortion solution
if keyword_set(undistort) then begin
	findpro,'nirc2redux',/noprint,dirlist=dirlist ; find this procedure's directory
	xdistort = mrdfits(dirlist[0]+'nirc2_'+camera+'_X_distortion.fits',/silent)
	ydistort = mrdfits(dirlist[0]+'nirc2_'+camera+'_Y_distortion.fits',/silent)
	; generate pixel grids in observed frame
	nx = 1024l
	ny = 1024l
	xorig = findgen(nx) # replicate(1., ny)
	yorig = replicate(1.,nx) # findgen(ny)
	; distortion corrected positions
	xtrue = reform(xorig + xdistort,nx*ny)
	ytrue = reform(yorig + ydistort,nx*ny)
	xorig = reform(xorig,nx*ny)
	yorig = reform(yorig,nx*ny)
endif
; CCD parameters
gain = sxpar(h,'gain') ; e/DN
no_reads =  (sxpar(h,'sampmode') eq 3)?sxpar(h,'multisam'):2
coadds = sxpar(h,'coadds')
rd_noise = 60./sqrt(no_reads) ; e-
if camera eq '0' then rd_noise = 0 ; if OSIRIS
exptime = sxpar(h,'itime')*coadds ; seconds
; fix bad pixels (IMG1)
fixpix,imgs,bpm,img1,/silent
; correct for gain
img1 = img1*gain ; DN -> e-
; generate Poisson noise cube 
noise = img1*0   ; save noise cube in e-
for i=0,ct-1 do noise[*,*,i] = sqrt(img1[*,*,i]+coadds*rd_noise^2)/flat

; Note: img1 and noise are the input cubes for subsequent steps

;==============================
; Step 1: 
; SUBTRACT MEDIAN SKY W/O MASKING
; minimum 3-dithering positions
;==============================
img2 = img1*0 ; datacube w/ direct median sky subtracted
; median combined sky
sky1set = img1 ; save masked and re-normalized sky frames
for i=0,ct-1 do begin
	; renormalize
	sky,sky1set[*,*,i],skymode,skysig,/nan,/silent
	sky1set[*,*,i] /= skymode
endfor
if ~sepsky then begin 
	sky1 = median(sky1set,dim=3,/even)
endif else begin ; if too few frames use supersky
	print,'Use supersky frame'
	sky0 = mrdfits(skyfile,0,/silent)
	fixpix,sky0,bpm,sky1,/silent
endelse 
sky,sky1,sky1mode,sky1sig,/nan,/silent
sky1 /= sky1mode ; renormalize
; subtract sky & flat-fielding (IMG2)
for i=0,ct-1 do begin
	sky,img1[*,*,i],skymode,skysig,/nan,/silent
	img2[*,*,i] = (img1[*,*,i]-sky1*skymode)/flat
endfor

;==============================
; Step 2:
; Mark objects in the frames
; either manually or automatically based on AO headers
; can decide whether or not to keep a certain frame
;==============================
setx
window,0,xsize=800,ysize=800
; select: holds whether or not to keep a frame
if ~keyword_set(exam) and file_test(output+'/discard.txt') then $
; if not in /exam mode and discard.txt exist then read it in
	readcol,output+'/discard.txt',jk1,select,f='a,i' $
else $ ; otherwise keep every frame at this point
	select = intarr(n_elements(ins))+1
for i=0,ct-1 do begin
	loadct,0
	astrim,asinh_scale(img2[*,*,i],img2[*,*,i],max=10),title=ins[i]
	if file_test(ins[i]+'.txt') then begin 
		; if txt file exist then skip the clicks on the primary source
		; still need to go through the discard process if in /exam mode
		readcol,ins[i]+'.txt',xc,yc
		plots,xc,yc,psym=6,syms=2
		if keyword_set(recenter) then begin ; if recenter, update xc/yc
			CNTRD,img2[*,*,i],xc,yc,xnew,ynew,1.5*fwhm,/silent
			; update if drift is reasonable
			ind = where(sqrt((xc-xnew)^2+(yc-ynew)^2) lt 2*fwhm,tmpct) 
			if tmpct gt 0 then begin
				xc[ind] = xnew[ind]
				yc[ind] = ynew[ind]
			endif
			plots,xc,yc,psym=8,syms=2				
		endif
		; use offsets therein (aotsx/y)
		if i eq 0 then begin
			xref = xc[0]
			yref = yc[0]
		endif else begin
			aotsx[i] = xc-xref
			aotsy[i] = yc-yref
		endelse
	endif else begin
		if i eq 0 then begin ; first frame
			xc = 0
			yc = 0
			a = 0
			print,' Identify objects in the image ... '
			print,'    - first, click on the primary alignment obj '
			print,'    - click on the left of the image to finish'
			while a ge 0 do begin
				cursor,a,b,/wait,/up
				CNTRD,img2[*,*,i],a,b,x,y,2*fwhm,/silent
				if sqrt((x-a)^2+(y-b)^2) gt 2*fwhm then begin ; drifted too far
					x = a
					y = b
				endif
				plots,x,y,psym=6,syms=2,color=cgcolor('red')
				xc = [xc,x]
				yc = [yc,y]
			endwhile
			print,' First frame finished'
			xc = xc[1:n_elements(xc)-2]
			yc = yc[1:n_elements(yc)-2]
			forprint,xc,yc,text=ins[i]+'.txt',/nocomment	
			; define xref and yref
			xref = xc[0]
			yref = yc[0]
		endif else begin ; subsequent frames
			if ~keyword_set(ngs) then begin ; if LGS then use header offset
				if keyword_set(recenter) then begin ; if recenter, update aotsx/y
					xc = xref+aotsx[i]
					yc = yref+aotsy[i]
					CNTRD,img2[*,*,i],xc,yc,xnew,ynew,1.5*fwhm,/silent
					; update if drift is reasonable
					if sqrt((xc-xnew)^2+(yc-ynew)^2) lt 2*fwhm then begin 
						aotsx[i] += xnew-xc
						aotsy[i] += ynew-yc
					endif
				endif
				xc = xref+aotsx[i]
				yc = yref+aotsy[i]
			endif else begin ; if NGS then click on primary for each frame
				print,' Click on the primary object '+ins[i]
				cursor,a,b,/wait,/up
				if a lt 0 then begin
					print,' Discarded frame '+ins[i]
					select[i] = 0
				endif
				CNTRD,img2[*,*,i],a,b,xc,yc,1.5*fwhm,/silent
				if sqrt((xc-a)^2+(yc-b)^2) gt 2*fwhm then begin ; drifted too far
					xc = a
					yc = b
				endif
				aotsx[i] = xc-xref
				aotsy[i] = yc-yref
			endelse
			; show primary object - large square
		   plots,xc,yc,psym=6,syms=3,color=cgcolor('red')
			; recentering - small circle
			if keyword_set(recenter) then begin
				CNTRD,img2[*,*,i],xc,yc,tmpx,tmpy,2*fwhm,/silent
			   if sqrt((xc-tmpx)^2+(yc-tmpy)^2) lt 2*fwhm then begin	
					xc = tmpx
					yc = tmpy
				endif
				plots,xc,yc,psym=8,syms=2,color=cgcolor('red')
			endif
			; save centroids in a text file
			forprint,xc,yc,text=ins[i]+'.txt',/nocomment	
		endelse
	endelse
	if keyword_set(exam) and ~keyword_set(ngs) then begin
		if i eq 0 then begin
			print,' Examine each frame'
			print,'    - Click anywhere inside the frame to keep it'
			print,'    - Click on the left of the image to discard it'
		endif else begin ; always include first frame
			cursor,a,b,/wait,/up
			if a lt 0 then begin
				print,' Discarded frame '+ins[i]
				select[i] = 0
			endif
		endelse
	endif
endfor
if keyword_set(exam) then begin
	print,' /exam finished ...'
	forprint,ins,select,text=output+'/discard.txt',/nocomment
endif
; if a separate sky frame is provided
if sepsky then begin ; skip subsequent steps
	sky3 = sky1
	goto,finalsteps
endif

;==============================
; Step 3:
; APPLY INITIAL OBJECT MASK TO ALL FRAMES
; THEN SUBTRACT MEDIAN SKY
; these masks are based on individual sky-subtracted images
; in step 4 the mask will be based on the stacked image
;==============================
for kk=1,2 do begin 
; run twice:
; 1st time sky-subtracted cube (img2) is based on direct median sky
; 2nd time sky-subtracted cube (img2=img5) is based on object-masked median sky
	img3 = img1*0 ; img1 w/ sources masked out as NaN
	img4 = img1*0 ; img3 renormalized by its median
	img5 = img1*0 ; img1 w/ median sky subtracted 
	for i=0,ct-1 do begin
		readcol,ins[i]+'.txt',xc,yc,f='f,f'
		if i eq 0 then begin ; 1st frame has coords for all marked sources
			x0 = xc
			y0 = yc
		endif
		; generate mask based on source location and S/N map
		snr = filter_image(reform(img2[*,*,i]/noise[*,*,i]),smooth=5)
		mask = bytarr((size(img1))[1:2])
		for j=0,n_elements(x0)-1 do begin
			; mask any pixels above snr_cut & within rmax
			dist_ellipse,im,(size(snr))[1:2],x0[j]+(xc[0]-x0[0]),y0[j]+(yc[0]-y0[0]),1,0
			ind = where(im lt rmin OR (snr gt snr_cut and im lt rmax),tmpct)
			if tmpct gt 0 then mask[ind] = 1b
		endfor
		; remove isolated islands
		nmask = remove_islands(mask,50,/all_neighbors)
		; apply masks to the cubes
		tmp = reform(img1[*,*,i]) 
		ind = where(nmask,tmpct)
		if tmpct ge 1 then tmp[ind] = !values.d_nan
		img3[*,*,i] = tmp
		img4[*,*,i] = tmp/median(tmp)
	endfor
	; median combine renomalized sky
	sky2 = median(img4,dim=3)
	sky,sky2,sky2mode,sky2sig,/nan,/silent
	sky2 /= sky2mode ; normalize
	; normal median-sky subtraction then flat fielding (IMG5)
	for i=0,ct-1 do begin
		sky,img3[*,*,i],skymode,skysig,/nan,/silent
		img5[*,*,i] = (img1[*,*,i]-sky2*skymode)/flat
	endfor
	; update img2 for S/N calculation
	img2 = img5
endfor

;==============================
; Step 4:
; GET OBJECT MASK FROM STACKED IMAGE
; REDO SKY SUBTRACTION BASED ON NEW MASK
; REMOVE RESIDUAL BACKGROUND W/ BSPLINE
;==============================
; generate a mask area for bspline subtraction
if camera eq 'wide' and keyword_set(bottom) then begin
	; define a special mask for the bottom 1/3 problematic area
	bmsk = bytarr((size(img1))[1],(size(img1))[2])
	x1 = 0. & y1 = 120. & x2 = 1023. & y2 = 450.
	for x=0,1023 do begin
		for y=0,(y2-y1)/(x2-x1)*(x-x1)+y1 do begin
			bmsk[x,y] = 1b
		endfor
	endfor
	bmidx = where(bmsk eq 0b)
endif
for kk=1,niter do begin
; run niter times:
; Initially, S/N map based on img5: sky subtracted w/ individual obj masks
; subsequently, S/N map based on img5 = bkg5: bspline+sky subtracted w/ stacked mask
	print,' Iteration '+strc(kk)+'/'+strc(niter)
	; enlarge images to get the full mosaic
	ind = where(select eq 1)
	xoff = abs(round(minmax(-1*aotsx[ind])))
	yoff = abs(round(minmax(-1*aotsy[ind])))
	zero = fltarr(1024+xoff[1]+xoff[0],1024+yoff[1]+yoff[0],ct)
	bimg = zero+!values.d_nan ; intensity map (e-)
	bvar = zero+!values.d_nan ; variance map  (e-^2)
	bexp = zero               ; exposure map  (s)
	for i=0,ct-1 do begin
		bimg[xoff[0]:xoff[0]+1023,yoff[0]:yoff[0]+1023,i] = img5[*,*,i] ; sky subtracted
		bvar[xoff[0]:xoff[0]+1023,yoff[0]:yoff[0]+1023,i] = noise[*,*,i]^2 ; variance = noise^2
		bexp[xoff[0]:xoff[0]+1023,yoff[0]:yoff[0]+1023,i] = exptime
	endfor
	; align intensity & noise maps 
	img6 = bimg
	var6 = bvar
	exp6 = bexp
	for i=0,ct-1 do begin
		img6[*,*,i] = shift_sub(img6[*,*,i],-1*aotsx[i],-1*aotsy[i],missing=!values.d_nan)
		var6[*,*,i] = shift_sub(var6[*,*,i],-1*aotsx[i],-1*aotsy[i],missing=!values.d_nan)
		exp6[*,*,i] = shift_sub(exp6[*,*,i],-1*aotsx[i],-1*aotsy[i],missing=0.0)
	endfor
	; median combine the selected, generate S/N map
	ind = where(select eq 1)
	img = median(img6[*,*,ind],dim=3)/exptime ; e-/s
	exm = total(exp6[*,*,ind],3,/nan) ; seconds
	noi = sqrt(median(var6[*,*,ind],dim=3))/exptime/sqrt(exm/exptime) ; e-/s
	snr = filter_image(reform(img/noi),smooth=5)	
	; load object list
	readcol,ins[0]+'.txt',xc,yc,f='f,f' ; read original list
	xc += xoff[0]
	yc += yoff[0]
	; modify object list
	if keyword_set(exam) then begin
		; display image w/ selected sources
		loadct,0
		astrim,asinh_scale(snr,snr,max=10),tit=output+' S/N map: Click on additional sources'
		tvcircle,rmax,xc,yc,/data
		print,' Click on additional objects in the combined image ... '
		a = 0
		while a ge 0 do begin
			cursor,a,b,/wait,/up
			CNTRD,img,a,b,x,y,2*fwhm,/silent ; measure centroid
			if sqrt((x-a)^2+(y-b)^2) gt 2*fwhm then begin ; drifted too far
				x = a
				y = b
			endif
			tvcircle,rmax,x,y,/data
			xc = [xc,x]
			yc = [yc,y]
		endwhile
		print,' Click finished'
		xc = xc[0:n_elements(xc)-2]
		yc = yc[0:n_elements(yc)-2]
		; re-display image and objects
		astrim,asinh_scale(snr,snr,max=10),tit=output+' S/N map: Click to delete'
		tvcircle,rmax,xc,yc,/data
		print,' Click on objects to delete ... '
		a = 0
		tokeep = intarr(n_elements(xc))+1 ; 
		while a ge 0 do begin
			cursor,a,b,/wait,/up
			mindis = min(sqrt((xc-a)^2+(yc-b)^2),idx)
			if a gt 0 and mindis lt rmax then begin
				tokeep[idx] = 0
				plots,xc[idx],yc[idx],psym=7,syms=5
			endif
		endwhile
		xc = xc[where(tokeep)]
		yc = yc[where(tokeep)]
		print,' Click finished'
		; re-display image and objects
		astrim,asinh_scale(snr,snr,max=10),tit=output+' S/N map'
		tvcircle,rmax,xc,yc,/data
		; update object list
		forprint,xc-xoff[0],yc-yoff[0],text=ins[0]+'.txt',/nocomment
	endif

	; re-measure centroids
	if keyword_set(recenter) then begin
		CNTRD,img,xc,yc,tmpx,tmpy,1.5*fwhm,/silent 
		dis = sqrt((tmpx-xc)^2+(tmpy-yc)^2)
		ind = where(dis lt 2*fwhm,tmpct)
		if tmpct gt 0 then begin
			xc[ind] = tmpx[ind]
			yc[ind] = tmpy[ind]
		endif
		forprint,xc-xoff[0],yc-yoff[0],text=ins[0]+'.txt',/nocomment ; update
	endif
	
	; skip the following steps if running in /initial mode
	if keyword_set(initial) then goto,theend

	; generate a mask image
	mask = bytarr((size(img))[1:2]) ; 0 - keep, 1 - masked
	for j=0,n_elements(xc)-1 do begin
		; mask any pixels above certain sigma & within rmax
		dist_ellipse,im,(size(snr))[1:2],xc[j],yc[j],1,0
		; use small masks for 1st iteration then enlarge mask w/ S/N map
		if kk eq 1 then ind = where(im lt rmin,tmpct) else $
		                ind = where(snr gt snr_cut and im lt rmax,tmpct)
		                ;ind = where(im lt rmin or (snr gt snr_cut and im lt rmax),tmpct)
		if tmpct gt 0 then mask[ind] = 1b
	endfor
	; remove isolated islands
	nmask = remove_islands(mask,50,/all_neighbors)
	; apply mask to datacube
	img3 = img1*0 ; img1 w/ sources masked out as NaN
	img4 = img1*0 ; img3 renormalized by its median
	bkg5 = img1*0 ; bspline subtracted img5
	for i=0,ct-1 do begin
		readcol,ins[i]+'.txt',xc,yc,f='f,f'
		if i eq 0 then begin
			x0 = xc
			y0 = yc
		endif
		; shift mask and apply to datacube
		mask = shift_sub(nmask,xc[0]-x0[0],yc[0]-y0[0],missing=0b)
		mask = mask[xoff[0]:xoff[0]+1023,yoff[0]:yoff[0]+1023]
		tmp = reform(img1[*,*,i]) 
		ind = where(mask,tmpct)
		if tmpct ge 1 then tmp[ind] = !values.d_nan
		img3[*,*,i] = tmp
		img4[*,*,i] = tmp/median(tmp) 
	endfor
	; median combine sky after sources being masked out
	sky3 = median(img4,dim=3)
	fixpix,sky3,sky3*0+1,sky3,/nan,/silent ; remove NaN from super sky
	sky,sky3,sky3mode,sky3sig,/nan,/silent
	sky3 /= sky3mode ; normalize
	; sky subtraction - update img5
	for i=0,ct-1 do begin
		sky,img3[*,*,i],skymode,skysig,/nan,/silent
		img5[*,*,i] = (img1[*,*,i]-sky3*skymode)/flat ; recalculate img5
		; normal median sky subtracted image
		tmp = reform(img5[*,*,i])
		; invalidate pixels under sources -> residual sky
		tmpsky0 = tmp
		ind = where(img3[*,*,i] ne img3[*,*,i],tmpct)
		if tmpct gt 0 then tmpsky0[ind] = !values.d_nan
		; fit b-splines
		tmpskyx = fltarr((size(tmp))[1:2]) ; columns
		tmpskyy = fltarr((size(tmp))[1:2]) ; rows
		xx = indgen((size(tmp))[2])
		for j=0,(size(tmp))[1]-1 do begin
			ind = where(tmpsky0[j,*] eq tmpsky0[j,*])
			sset = bspline_iterfit(xx[ind],tmpsky0[j,ind],bkspace=bkspace)
			tmpskyx[j,*] = bspline_valu(xx,sset)
		endfor
		xx = indgen((size(tmp))[1])
		for j=0,(size(tmp))[2]-1 do begin
			ind = where(tmpsky0[*,j] eq tmpsky0[*,j])
			sset = bspline_iterfit(xx[ind],tmpsky0[ind,j],bkspace=bkspace)
			tmpskyy[*,j] = bspline_valu(xx,sset)
		endfor
		tmpsky1 = (tmpskyx+tmpskyy)/2.0 ; bspline fitted residual background
		; bspline sky subtraction
		tmpbkg = tmp-tmpsky1
		; subtract a constant sky for the bottom 1/3 problematic area
		if camera eq 'wide' and keyword_set(bottom) then tmpbkg[bmidx]=tmp[bmidx]-median(tmpsky1[bmidx])
		bkg5[*,*,i] = tmpbkg
		; display background subtracted images
		; top left (tmpsky0): median-sky subtracted frame with mask applied
		; top right (tmpsky1): bspline best-fit model of tmpsky0
		; bottom left (img5[*,*,i]): median-sky-subtracted frame 
		; bottom right (bkg5[*,*,i]): median-sky-subtracted frame - bspline-fit residual sky
		astrim,asinh_scale(([[img5[*,*,i],bkg5[*,*,i]],[tmpsky0,tmpsky1]]),img5[*,*,i],max=10),$
			tit=ins[i]
		if camera eq 'wide' and keyword_set(bottom) then contour,[[bmsk,bmsk],[bmsk,bmsk]],/overplot,levels=1
	endfor
	; for every but the last iteration, replace img5 with bkg5 to compute the combined S/N map
	; after the S/N map is computed and mask updated, img5 is recomputed in each iteration,
	; img5 therefore remains as the median-sky subtracted datacube before
	; the bspline subtraction of any residual sky
	if kk lt niter then img5 = bkg5  
endfor
; save the final object mask in a PNG file, overlaid with S/N map and 
loadct,0
astrim,asinh_scale(img,img,maxv=10),tit=output+': Intensity + S/N Map + Object Mask'
contour,snr,color=cgcolor('red'),levels=snr_cut*[1,4,16],/overplot
contour,nmask,color=cgcolor('green'),levels=1,/overplot
save_screen,dir+'/mask.png'
loadct,0

finalsteps:
;==============================
; Step 5:
; DO CAMERA DISTORTION
; REJECT FRAMES THAT ARE BAD OR HAVE POOR IMAGE QUALITY
; REALIGN AND COMBINE IMAGES W/ 2D PEAK FIT
;==============================
; sky subtraction, flat fielding, distortion correction (IMG5)
noi5 = noise ; to save noise maps
; distortion correction
if keyword_set(undistort) then begin
	print,'Working on distortion correction'
	for i=0,ct-1 do begin
		;; camera distortion correction
		;; 02/10/2012: image reduced with dewarp has worse quality than w/o
		;;             dewarp, so decided to skip this step
		;if camera eq '0' then img5[*,*,i] = tmp-tmpsky1 else $
		;img5[*,*,i] = nirc2dewarp(tmp-tmpsky1,camera=camera)
		; new solution based on M92 observations
   	img5[*,*,i] = warp_tri(xtrue,ytrue,xorig,yorig,reform(img5[*,*,i]))
   	bkg5[*,*,i] = warp_tri(xtrue,ytrue,xorig,yorig,reform(bkg5[*,*,i]))
		noi5[*,*,i] = warp_tri(xtrue,ytrue,xorig,yorig,reform(noi5[*,*,i]))
	endfor
endif

; after applying camera distortion, now we need to update offsets 
if keyword_set(recenter) then begin
	for i=0,ct-1 do begin
		readcol,ins[i]+'.txt',a,b,f='f,f'	
		CNTRD,img5[*,*,i],a[0],b[0],x,y,1.5*fwhm,/silent
		if i eq 0 then begin
			if sqrt((x-a[0])^2+(y-b[0])^2) gt 2*fwhm then begin ; drifted too far
				x0 = a[0]
				y0 = b[0]
			endif else begin
				x0 = x
				y0 = y
			endelse
		endif else begin
			if sqrt((x-a)^2+(y-b)^2) gt 2*fwhm then begin ; drifted too far
				x = a
				y = b
			endif
			aotsx[i] = x-x0
			aotsy[i] = y-y0
		endelse
	endfor
endif

; enlarge images to get the full mosaic
ind = where(select eq 1)
xoff = abs(round(minmax(-1*aotsx[ind])))
yoff = abs(round(minmax(-1*aotsy[ind])))
zero = fltarr(1024+xoff[1]+xoff[0],1024+yoff[1]+yoff[0],ct)
bimg = zero+!values.d_nan ; intensity map (e-)
bbkg = zero+!values.d_nan ; bspline subtracted intensity map (e-)
bvar = zero+!values.d_nan ; variance map  (e-^2)
bexp = zero               ; exposure map  (s)
for i=0,ct-1 do begin
	bimg[xoff[0]:xoff[0]+1023,yoff[0]:yoff[0]+1023,i] = img5[*,*,i]
	bbkg[xoff[0]:xoff[0]+1023,yoff[0]:yoff[0]+1023,i] = bkg5[*,*,i]
	bvar[xoff[0]:xoff[0]+1023,yoff[0]:yoff[0]+1023,i] = noi5[*,*,i]^2 ; variance = noise^2
	bexp[xoff[0]:xoff[0]+1023,yoff[0]:yoff[0]+1023,i] = exptime
endfor
; update primary aligment source position
x0 += xoff[0]
y0 += yoff[0]

; align maps - first pass using previous determined centroids
img6 = bimg
for i=0,ct-1 do $ 
	img6[*,*,i] = shift_sub(bimg[*,*,i],-1*aotsx[i],-1*aotsy[i],cubic=-0.5,missing=!values.d_nan)

; align maps - second pass using centroid from 2D peak fit
gwidth = fltarr(ct) 
if keyword_set(recenter) then begin
	; measure FWHM on the alignment source for each frame
	x1 = (x0[0]-7*fwhm)>0
	x2 = (x0[0]+7*fwhm)<((size(img6))[1]-1)
	y1 = (y0[0]-7*fwhm)>0
	y2 = (y0[0]+7*fwhm)<((size(img6))[2]-1)
	for i=0,ct-1 do begin
		if i eq 0 then begin
			; get parameter estimates on stacked image
			ind = where(select eq 1)
			tmp = median(img6[x1:x2,y1:y2,ind],dim=3)
			result = mpfit2dpeak(tmp,a0,/LORENTZIAN,/tilt)
		endif
		result = mpfit2dpeak(img6[x1:x2,y1:y2,i],a,/LORENTZIAN,/tilt,estimates=a0)
		gwidth[i] = min(a[2:3]) ; Gaussian sigma or HWHM
		; Gaussian peak centroids, adjust offsets 
		aotsx[i] += (a[4]-a0[4])
		aotsy[i] += (a[5]-a0[5])
	endfor
	; align frames again (IMG6)
	img6 = bimg
	for i=0,ct-1 do $
		img6[*,*,i] = shift_sub(bimg[*,*,i],-1*aotsx[i],-1*aotsy[i],cubic=-0.5,missing=!values.d_nan)
endif

; align variance and exposure maps
bkg6 = bbkg
var6 = bvar 
exp6 = bexp
for i=0,ct-1 do begin 
	bkg6[*,*,i] = shift_sub(bbkg[*,*,i],-1*aotsx[i],-1*aotsy[i],cubic=-0.5,missing=!values.d_nan)
	var6[*,*,i] = shift_sub(bvar[*,*,i],-1*aotsx[i],-1*aotsy[i],cubic=-0.5,missing=!values.d_nan)
	exp6[*,*,i] = shift_sub(bexp[*,*,i],-1*aotsx[i],-1*aotsy[i],missing=0.0)
endfor

; reject images with large FWHMs: rejsigma x sigma above median
if rejsig gt 0 and keyword_set(recenter) then begin
   idx = where(gwidth gt median(gwidth)+rejsig*robust_sigma(gwidth),ct)
	if ct gt 0 then select[idx] = 0
endif
; median combine the selected
ind = where(select eq 1)
exm = total(exp6[*,*,ind],3,/nan) ; seconds
img = median(img6[*,*,ind],dim=3)/exptime ; e-/s
bkg = median(bkg6[*,*,ind],dim=3)/exptime ; e-/s
noi = sqrt(median(var6[*,*,ind],dim=3))/exptime/sqrt(exm/exptime) ; e-/s
print,' Combined: ',n_elements(ind),' frames'

; save combined frames info
forprint,ins,gwidth,select,f='a,f,i',text=dir+'/fwhm.txt',$
	comment='# File, Lorentian HWHM [pix]: median = '+strc(median(gwidth))+$
	' stddev = '+strc(robust_sigma(gwidth))+', combined?'

; insert rough WCS
h0 = headfits(ins[0],ext=0) 
par = (sxpar(h0,'ROTDEST')-sxpar(h0,'INSTANGL')-0.252)*!dtor ; PA in radian
putast,h0,[[-1.*cos(par),sin(par)],[sin(par),cos(par)]]*pixscale/3600d,$
	[512+15+xoff[0]+0.5,512+30+yoff[0]+0.5],[sxpar(h0,'ra'),sxpar(h0,'dec')],$
	['RA---TAN','DEC--TAN']

; write FITS file
sxaddpar,h0,'gain',1,' ADU to e-' 
sxaddpar,h0,'ncombine',n_elements(where(select eq 1)),' Number of combined frames'
sxaddpar,h0,'bunit','e-/s',' Intensity map unit: e- per second' ; unit in e-/s
sxaddpar,h0,'history','Reduced by Hai Fu on '+begin_time
mwrfits,bkg,output+'.fits',h0,/create ; save intensity map w/ bspline subtraction
;
mkhdr,h,noi,/image
sxaddpar,h,'bunit','e-/s',' Noise map unit: e- per second'
mwrfits,noi,output+'.fits',h ; save noise 
;
mkhdr,h,exm,/image
sxaddpar,h,'bunit','s',' Exposure map unit: seconds'
mwrfits,exm,output+'.fits',h ; save exposure map 
;
mkhdr,h,img,/image
sxaddpar,h,'bunit','e-/s',' Intensity map unit: e- per second'
mwrfits,img,output+'.fits',h ; save intensity map w/o bspline subtraction 

theend:
; show time ellapsed
print,begin_time
spawn,'date'
end
