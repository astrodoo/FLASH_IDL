pro playjet_recoll,fname,var=var,sample=sample,startn=startn,endn=endn,step=step,png=png,zoom=zoom
device, decomposed=0

if not keyword_set(fname) then fname='JetSet_hdf5_plt_cnt_'
if not keyword_set(var) then var='dens'
;if not keyword_set(sample) then sample=0
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

xc0 = -2.e12 & yc0 = 0. & zc0 = 0.
;maxd = 2.2d-13 & mind = 1.1d-15
;maxd = 2.76d-13 & mind = 2.9d-16
;maxd = 2.76d-13 & mind = 1.d-16
maxd = 3.d-11 & mind = 1.85d-16
sample=4

m1 = 4.d34
m2 = 2.d34
mtot = m1+m2
peri = 3.d12

zoom=1
for j=start_j,end_j,step do begin
    print,j+1,'   of',end_j+1,'  ', fnames[j]
    jind = fix(strmid(fnames[j],3,4,/rev))
    if not keyword_set(zoom) then begin
       if not keyword_set(sample) then begin
          read_amr,fnames[j],var='dens',tree=tree
          smp = max(tree.lrefine) - 7
       endif else smp = sample
       d = dload(jind,var=var,xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=smp,time) 
    endif else begin
       if not keyword_set(sample) then smp=0 else smp=sample
       smp=3
       xrange = [-2.5e12,-1.5e12]
       yrange = [-8.e11,8.e11]
       zrange = [-1.e11,2.e12]
       d = loaddata(fnames[j],var,xra=xrange,yra=yrange $
           ,zra=zrange,time=time,sample=smp,xcoords=x,ycoords=y,zcoords=z)
    endelse
    sd = size(d,/dimension)
    x0 = 100 & y0 = 60
;maxd = max(d)   & mind = min(d)

    dxz = reform(d[*,sd[1]/2,*])
    xind = where(x ge -2.e12) & xind = reform(xind[0])
    dyz = reform(d[xind,*,*])

    loadct,0,/sil
    exestr = execute('window,0,xs='+strtrim(sd[0]+sd[1]+x0+20,1)+',ys='+strtrim(sd[2]+y0+20,1)+opt)
;    plot,x,z,/xst,/yst,/nodata,xr=[x[0],x[sd[0]-1]],yr=[z[0],z[sd[2]-1]],/iso,position=[x0,y0,x0+sd[0],y0+sd[2]],/dev
;    tv,bytscl(alog10(d[*,*,256]),max=alog10(maxd),min=alog10(mind))
    tv,bytscl(alog10(dxz),max=alog10(maxd),min=alog10(mind)),x0,y0
    plot,x,z,/xst,/yst,/nodata,xr=[x[0],x[sd[0]-1]],yr=[z[0],z[sd[2]-1]],/iso,position=[x0,y0,x0+sd[0],y0+sd[2]],/dev,/noerase,xtitle='x [cm]',ytitle='z [cm]',xtickinterval=4.e11
    

;    plot,y,z,/xst,/yst,/nodata,xr=[y[1],y[sd[1]-1]],yr=[z[0],z[sd[2]-1]],/iso,position=[x0+sd[0],y0,x0+sd[0]+sd[1],y0+sd[2]],/dev,/noerase
;    tv,bytscl(alog10(d[*,*,256]),max=alog10(maxd),min=alog10(mind))
    tv,bytscl(alog10(dyz),max=alog10(maxd),min=alog10(mind)),x0+sd[0],y0
    plot,y,z,/xst,/yst,/nodata,xr=[y[0],y[sd[1]-1]],yr=[z[0],z[sd[2]-1]],/iso,position=[x0+sd[0],y0,x0+sd[0]+sd[1],y0+sd[2]],/dev,/noerase, ytickformat='(a1)',xtitle='y [cm]',xtickinterval=4.e11

    tvlct,r,g,b,/get
    legend,'time = '+strtrim(time/60./60.,2)+' hrs',/right,/top,color=fsc_color('yellow'),box=0
    tvlct,r,g,b

;    legend,strtrim(time,1)+ ' s',/right,/top,text=fsc_color('yellow'),box=0
;     legend,string(sqrt(!unit.g * mtot / peri^3.d0)/2./!pi*time,format='(f4.2)')+' Period',/right,/top,box=0,textcolor=fsc_color('yellow')

    if keyword_set(png) then $
       snapshot,'png_'+var+'/Jet_plt_'+string(jind,format='(I4.4)')
endfor
stop
end


