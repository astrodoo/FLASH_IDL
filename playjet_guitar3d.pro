pro playjet_guitar3d,fname,var=var,sample=sample,startn=startn,endn=endn,step=step,png=png,xrange=xrange,yrange=yrange
device, decomposed=0

if not keyword_set(fname) then fname='PWN3d_hdf5_plt_cnt_'
if not keyword_set(var) then var='dens'
if not keyword_set(sample) then sample=0
if not keyword_set(step) then step=1

fnames = file_search(fname+'*')
nfiles = n_elements(fnames)

if keyword_set(startn) then begin
   str_start = string(startn,format='(I4.4)')
   start_j_tmp = where(strpos(fnames,str_start) ne -1)
   start_j = start_j_tmp[0]
   if start_j eq -1 then begin
      print,'out of range of startnumber'
      stop
   endif
endif else start_j = 0
if keyword_set(endn) then begin
   str_end = string(endn,format='(I4.4)')
   end_j_tmp = where(strpos(fnames,str_end) ne -1)
   end_j = end_j_tmp[0]
   if end_j eq -1 then begin
      print,'out of range of startnumber'
      stop
   endif
endif else end_j = nfiles-1
nfiles = end_j - start_j + 1

; Making PNG files
if keyword_set(png) then begin
   spawn,'mkdir png_'+var
   opt=',/pixmap' 
endif else opt=''

yc= 0
zc= 0
;maxd = -23 & mind = -28   ;10l case
;maxd = -22 & mind = -28    ;10h case

maxd = 7.e-24 & mind = 1.e-28

;if not keyword_set(xrange) then xrange=[1.e18,2.4e20]

;x1=1.e18 & x2=2.4e20
loadct,0,/sil
!p.background=255
for j=start_j,end_j,step do begin
    print,j+1,'   of',end_j+1,'  ', fnames[j]
;    d = loaddata(fnames[j],var,xra=xrange,time=time,sample=smp,xcoords=x,ycoords=y,zcoords=z)
    d = loaddata(fnames[j],var,time=time,sample=sample,xcoords=x,ycoords=y,zcoords=z)
    sd = size(d,/dimension)
   
;    y0ind_tmp = where(y ge yc) & y0ind = y0ind_tmp[0]
    z0ind_tmp = where(z ge zc) & z0ind = z0ind_tmp[0]

;    imgxz = reform(d[*,y0ind,*])
    imgxy = reform(d[*,*,z0ind])

    exestr = execute('window,0,xs='+strtrim(sd[0],1)+',ys='+strtrim(sd[1],1)+opt)
    tv,bytscl(alog(imgxy),min=alog(mind),max=alog(maxd))
;    tv,bytscl(alog(imgxz),min=alog(mind),max=alog(maxd)),0
;    xyouts,10,sd[1]*2-20,/dev,'X-Z plane',color=0
;    tv,bytscl(alog10(imgxy),min=mind,max=maxd),0,0
;    xyouts,10,sd[1]-30,/dev,'X-Y plane',color=0

;    legend,strtrim(time/60./60./24./365.,1)+ ' yr',/right,/top,textcolor=0,box=0

    if keyword_set(png) then $
       draw,'png_'+var+'/Jet_plt_'+string(j,format='(I4.4)')
endfor
stop
end
