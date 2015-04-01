pro clip2d,fname,xc0=xc0,yc0=yc0,xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,startn=startn,endn=endn $
          ,xcut=xcut,ycut=ycut,zcut=zcut,maxlref_sim=maxlref_sim,out=out

; for '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/sphere_1e37'
;xrange=[-1.5e13,5.e12]
;yrange=[-1.e13,1.e13]
;zrange=[-3.e11,1.9e13]
;startn=399
;endn=400
;sample=3
;maxlref_sim=11
;xcut=-2.e12
;ycut=0.

; for /d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L38_M30_lg
zrange=[-3.e11,2.e13]
;startn=2634
;endn=2634
sample=5
maxlref_sim=12
xcut=-2.e12
ycut=0.

if not keyword_set(fname) then fname='plt_cnt'
if not keyword_set(sample) then sample=3
if not keyword_set(maxlref_sim) then maxlref_sim = 11
if not keyword_set(out) then out='clip2d.sav'

fnames = file_search('*'+fname+'*')
nfiles = n_elements(fnames)

if keyword_set(startn) then begin
   str_start = string(startn,format='(I4.4)')
   start_j = (where(strpos(fnames,str_start) ne -1))[0]
   if start_j eq -1 then begin
      print,'out of range of startnumber'
      stop
   endif
endif else start_j = 0
if keyword_set(endn) then begin
   str_end = string(endn,format='(I4.4)')
   end_j = (where(strpos(fnames,str_end) ne -1))[0]
   if end_j eq -1 then begin
      print,'out of range of endnumber'
      stop
   endif
endif else end_j = nfiles-1
nfiles = end_j - start_j + 1

read_amr,fnames[start_j],param=p,tree=t,/nodata
maxlref = max(t.lrefine)
if not keyword_set(xrange) then xrange = [min(t.bndbox[0,0]), max(t.bndbox[1,0])]
if not keyword_set(yrange) then yrange = [min(t.bndbox[0,1]), max(t.bndbox[1,1])]
if not keyword_set(zrange) then zrange = [min(t.bndbox[0,2]), max(t.bndbox[1,2])]

print,'xra: ',xrange
print,'yra: ',yrange
print,'zra: ',zrange

d = loaddata(fnames[start_j],'dens',xra=xrange,yra=yrange,zra=zrange,sample=sample-(maxlref_sim-maxlref) $
            ,xcoords=x,ycoords=y,zcoords=z,time=time)
sz = size(d,/dimension)

if (n_elements(zcut) ne 0) then begin
   optxy = 'dxy:fltarr(sz[0],sz[1]),pxy:fltarr(sz[0],sz[1]),jxy:fltarr(sz[0],sz[1])'+ $
           ',v1xy:fltarr(sz[0],sz[1]),v2xy:fltarr(sz[0],sz[1]),v3xy:fltarr(sz[0],sz[1])' 
   zc_ind = (where(z ge zcut))[0]
   strzcut = 'zcut'
endif else begin
   optxy=''
   strzcut=''
endelse
if (n_elements(xcut) ne 0) then begin
   optyz = 'dyz:fltarr(sz[1],sz[2]),pyz:fltarr(sz[1],sz[2]),jyz:fltarr(sz[1],sz[2])'+ $
           ',v1yz:fltarr(sz[1],sz[2]),v2yz:fltarr(sz[1],sz[2]),v3yz:fltarr(sz[1],sz[2])' 
   xc_ind = (where(x ge xcut))[0]
   strxcut = 'xcut'
endif else begin
   optyz=''
   strxcut=''
endelse
if (n_elements(ycut) ne 0) then begin
   optxz = 'dxz:fltarr(sz[0],sz[2]),pxz:fltarr(sz[0],sz[2]),jxz:fltarr(sz[0],sz[2])'+ $
           ',v1xz:fltarr(sz[0],sz[2]),v2xz:fltarr(sz[0],sz[2]),v3xz:fltarr(sz[0],sz[2])' 
   yc_ind = (where(y ge ycut))[0]
   strycut='ycut'
endif else begin
   optxz=''
   strycut=''
endelse

if ((optxy ne '') and (optyz ne '')) then begin
   optyz=','+optyz
   strxcut=','+strxcut
endif
if (((optxy ne '') or (optyz ne '')) and (optxz ne '')) then begin
   optxz=','+optxz
   strycut=','+strycut
endif

strexe = execute('data = replicate({'+optxy+optyz+optxz+ $
       ',time:0.},nfiles)')
k=0
for jj=start_j,end_j do begin
    print,jj+1,'   of',end_j+1,'  ', fnames[jj]   

    read_amr,fnames[jj],tree=t,/nodata
    maxlref = max(t.lrefine)
    
    if (n_elements(zcut) ne 0) then begin
       data[k].dxy = reform(loaddata(fnames[jj],'dens',xra=xrange,yra=yrange,zra=[zcut,zcut] $
                 ,sample=sample-(maxlref_sim-maxlref), time=time))
       data[k].pxy = reform(loaddata(fnames[jj],'pres',xra=xrange,yra=yrange,zra=[zcut,zcut] $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].jxy = reform(loaddata(fnames[jj],'jet ',xra=xrange,yra=yrange,zra=[zcut,zcut] $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].v1xy = reform(loaddata(fnames[jj],'velx',xra=xrange,yra=yrange,zra=[zcut,zcut] $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].v2xy = reform(loaddata(fnames[jj],'vely',xra=xrange,yra=yrange,zra=[zcut,zcut] $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].v3xy = reform(loaddata(fnames[jj],'velz',xra=xrange,yra=yrange,zra=[zcut,zcut] $
                 ,sample=sample-(maxlref_sim-maxlref)))
    endif
    if (n_elements(xcut) ne 0) then begin
       data[k].dyz = reform(loaddata(fnames[jj],'dens',xra=[xcut,xcut],yra=yrange,zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref), time=time))
       data[k].pyz = reform(loaddata(fnames[jj],'pres',xra=[xcut,xcut],yra=yrange,zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].jyz = reform(loaddata(fnames[jj],'jet ',xra=[xcut,xcut],yra=yrange,zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].v1yz = reform(loaddata(fnames[jj],'velx',xra=[xcut,xcut],yra=yrange,zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].v2yz = reform(loaddata(fnames[jj],'vely',xra=[xcut,xcut],yra=yrange,zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].v3yz = reform(loaddata(fnames[jj],'velz',xra=[xcut,xcut],yra=yrange,zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
    endif
    if (n_elements(ycut) ne 0) then begin
       data[k].dxz = reform(loaddata(fnames[jj],'dens',xra=xrange,yra=[ycut,ycut],zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref), time=time))
       data[k].pxz = reform(loaddata(fnames[jj],'pres',xra=xrange,yra=[ycut,ycut],zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].jxz = reform(loaddata(fnames[jj],'jet ',xra=xrange,yra=[ycut,ycut],zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].v1xz = reform(loaddata(fnames[jj],'velx',xra=xrange,yra=[ycut,ycut],zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].v2xz = reform(loaddata(fnames[jj],'vely',xra=xrange,yra=[ycut,ycut],zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
       data[k].v3xz = reform(loaddata(fnames[jj],'velz',xra=xrange,yra=[ycut,ycut],zra=zrange $
                 ,sample=sample-(maxlref_sim-maxlref)))
    endif

    data[k].time = time
    k++
endfor

strexe2 = execute("save,file='"+out+"',x,y,z,data,"+strzcut+strxcut+strycut)

stop
end


pro playclip2d,file=file
device,decomposed=0

;file='clip2d_tmp.sav'
if not keyword_set(file) then file='clip2d.sav'

restore,file
print,'complete to read data'

ndata = (size(data,/dimension))[0]
plotsz = size(data[0].dyz,/dimension)

spawn,'mkdir tmp_playclip2d'

xp0=100. & yp0=80.

loadct,0,/sil
nx = xp0+plotsz[0]*4+10 & ny = yp0+plotsz[1]*2+200

mind=1.e-16 & maxd=1.e-13
minjv= 0. & maxjv=1.
minp=1. & maxp = 650.

vj = 3.d9

; loop
for i=0,ndata-1 do begin
print,i, ' of ', ndata-1
;swindow,xs=nx,ys=ny
window,xs=nx,ys=ny,/pixmap
jvyz = data[i].jyz*data[i].v3yz/vj
jvxz = data[i].jxz*data[i].v3xz/vj

rampyz = data[i].dyz*(data[i].v1yz*data[i].v1yz + data[i].v2yz*data[i].v2yz)
rampxz = data[i].dxz*(data[i].v1xz*data[i].v1xz + data[i].v2xz*data[i].v2xz)

loadct,0,/sil
tvcoord,bytscl(alog10(data[i].dyz),min=alog10(mind),max=alog10(maxd)),y,z,/axes,pos=[xp0,yp0+plotsz[1]+80,xp0+plotsz[0],yp0+2*plotsz[1]+80], ytitle='z [cm]'
color_bar,lim=[mind,maxd],/log,/top,bartitle='density',titlegap=0.03
loadct,3,/sil
tvcoord,bytscl(jvyz,min=minjv,max=maxjv),y,z,/axes,pos=[xp0+plotsz[0],yp0+plotsz[1]+80,xp0+2*plotsz[0],yp0+2*plotsz[1]+80],ytickformat='(a1)'
color_bar,lim=[minjv,maxjv],/top,bartitle='j*vz/vj',titlegap=0.03
loadct,13,/sil
tvcoord,bytscl(alog10(data[i].pyz),min=alog10(minp),max=alog10(maxp)),y,z,/axes,pos=[xp0+2*plotsz[0],yp0+plotsz[1]+80,xp0+3*plotsz[0],yp0+2*plotsz[1]+80],ytickformat='(a1)'
color_bar,lim=[minp,maxp],/top,bartitle='pressure',titlegap=0.03
tvcoord,bytscl(alog10(rampyz),min=alog10(minp),max=alog10(maxp)),y,z,/axes,pos=[xp0+3*plotsz[0],yp0+plotsz[1]+80,xp0+4*plotsz[0],yp0+2*plotsz[1]+80],ytickformat='(a1)'
color_bar,lim=[minp,maxp],/top,bartitle='ram pressure',titlegap=0.03

loadct,0,/sil
xyouts,nx/2.,yp0+plotsz[1]+30,/dev,'y [cm]'
xyouts,nx/2.,30,/dev,'x [cm]'
tvcoord,bytscl(alog10(data[i].dxz),min=alog10(mind),max=alog10(maxd)),x,z,/axes,pos=[xp0,yp0,xp0+plotsz[0],yp0+plotsz[1]], ytitle='z [cm]'
loadct,3,/sil
tvcoord,bytscl(jvxz,min=minjv,max=maxjv),x,z,/axes,pos=[xp0+plotsz[0],yp0,xp0+plotsz[0],yp0+2*plotsz[1]],ytickformat='(a1)'
loadct,13,/sil
tvcoord,bytscl(alog10(data[i].pxz),min=alog10(minp),max=alog10(maxp)),x,z,/axes,pos=[xp0+2*plotsz[0],yp0,xp0+3*plotsz[0],yp0+plotsz[1]],ytickformat='(a1)'
tvcoord,bytscl(alog10(rampxz),min=alog10(minp),max=alog10(maxp)),x,z,/axes,pos=[xp0+3*plotsz[0],yp0,xp0+4*plotsz[0],yp0+plotsz[1]],ytickformat='(a1)'

xyouts,30,ny-30,/dev,string(data[i].time/60./60.,format='(f6.2)')+' hrs',color=fsc_color('yellow'),charsize=3.

snapshot,'tmp_playclip2d/playclip2d_'+string(i,format='(I4.4)')
endfor
stop
end
