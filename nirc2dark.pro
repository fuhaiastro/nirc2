pro nirc2dark,dir,list,output

list = 'n'+string(list,f='(i04)')+'.fits'
; load all files
print,output
for i=0,n_elements(list)-1 do begin
	spawn,'gethead '+dir+list[i]+' itime coadds multisam sampmode filter',result
	if i eq 0 then dark = mrdfits(dir+list[i],0,h,/silent) $
             else dark = [[[dark]],[[mrdfits(dir+list[i],0,/silent)]]]
	print,list[i]+' '+result+' '+strc(median(dark[*,*,i]))	
endfor
; save
mwrfits,median(dark,dim=3),output,h,/create

end

