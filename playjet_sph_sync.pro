pro playjet_sph_sync,fname,var=var,startn=startn,endn=endn,step=step,png=png
device, decomposed=0

if not keyword_set(fname) then fname='JetSet_hdf5_plt_cnt_'
if not keyword_set(var) then var='pres'
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
endif else  begin
   start_j = 0 & startn = 0
endelse

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
   spawn,'mkdir png_sync'
   opt=',/pixmap' 
endif else opt=''

xc0 = 0. & yc0 = 0. & zc0 = 0.
;maxd = 2.2d-13 & mind = 1.1d-15
;maxd = 2.76d-13 & mind = 2.9d-16
;maxd = 2.76d-13 & mind = 1.d-16
;maxd = 8.e-14 & mind = 3.d-18
maxp = 5.e5 & minp = 100

m1 = 4.d34
m2 = 2.d34
mtot = m1+m2
peri = 3.d12

xc1 = -peri*m1/mtot

jeton = 9.76e4

!p.background=255
for j=start_j,end_j,step do begin
    print,j+1,'   of',end_j+1,'  ', fnames[j]
    if not keyword_set(sample) then begin
       read_amr,fnames[j],var='pres',tree=tree
       smp = max(tree.lrefine) - 7
    endif else smp = sample
;    p = dload(j,var=var,xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=smp,time) 
    p = dload(j+startn,var=var,xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=smp,time) 

    pSum = total(p^1.75,2)

    img0 = congrid(pSum,768,768)
    sd = size(img0)

    loadct,0,/sil
    exestr = execute('window,0,xs=768,ys=768'+opt)
;    plot,x,y,/xst,/yst,/nodata,xr=[x[0],x[sd[1]-1]],yr=[y[0],y[sd[2]-1]],/iso,position=[0,0,sd[1],sd[2]],/dev
    tv,bytscl(alog10(img0),max=alog10(maxp),min=alog10(minp))
;    overorbit

;    legend,strtrim(time,1)+ ' s',/right,/top,text=fsc_color('yellow'),box=0
     legend,string(time/60./60./24.,format='(f4.2)')+' days',box=0,textcolor=fsc_color('yellow') $
           ,position=[600,750],/dev

    if keyword_set(png) then $
;       draw,'png_sync/Jet_plt_'+string(j,format='(I4.4)')
       draw,'png_sync/Jet_plt_'+string(j+startn,format='(I4.4)')
endfor
stop
end
;

