pro nirc2flat,dir,onlist,offlist,output

onlist = 'n'+string(onlist,f='(i04)')+'.fits'
offlist = 'n'+string(offlist,f='(i04)')+'.fits'

; load all files
print,output
print,'Lamp on: '
for i=0,n_elements(onlist)-1 do begin
	spawn,'gethead '+dir+onlist[i]+' itime coadds multisam sampmode filter', result
	if i eq 0 then flaton = mrdfits(dir+onlist[i],0,h,/silent) $
	          else flaton = [[[flaton]],[[mrdfits(dir+onlist[i],0,/silent)]]]
	print,onlist[i]+' '+result+' '+strc(median(flaton[*,*,i]))
endfor

print,'Lamp off:'
for i=0,n_elements(offlist)-1 do begin
	spawn,'gethead '+dir+offlist[i]+' itime coadds multisam sampmode filter',result
	if i eq 0 then flatoff = mrdfits(dir+offlist[i],0,h,/silent) $
	          else flatoff = [[[flatoff]],[[mrdfits(dir+offlist[i],0,/silent)]]]
	print,offlist[i]+' '+result+' '+strc(median(flatoff[*,*,i]))
endfor

flat = median(flaton,dim=3)-median(flatoff,dim=3)

mwrfits,flat,output,h,/create

end
