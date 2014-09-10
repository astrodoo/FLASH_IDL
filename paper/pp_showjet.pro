pro pp_showjet,n,sample=sample,block=block,zoom=zoom,dx=dx,lineread=lineread,var=var $
   ,noorbit=noorbit,nolog=nolog,veloxy=veloxy,veloxz=veloxz,bh=bh,c0=xc0,yc0=yc0,zc0=zc0 $
   ,maskoff=maskoff,ps=ps,out=out,bar=bar,xzoom=xzoom,yzoom=yzoom,zzoom=zzoom,copy=copy
device,decomposed=0

bh=1
sample=3
zoom=1
xzoom=[-1.9158e13,3.e12]
yzoom=[-5.e12,5.e12]
zzoom=[-1.9158e13,1.925e13]

x0 = 280 & y0 = 150
x1 = 10  & y1 = 20
barx=300
; for sph_1e36 970
;cd,'/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_high'
!p.charsize=1.5
;n=971
;copy=1
;bar=1

; for sph_1e36 50
;cd,'/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36'
;!p.charsize=2.
;n=50   ; sph_1e36 earlier case
;for sph_1e36 300
;n=300   ; sph_1e36 earlier case
;copy=1

; for sph_1e36 30/60 deg
!p.charsize=2.
n=530  ; 30deg
;n=462  ; 60deg

; for sph_1e36 75 deg
;!p.charsize=1.5
;n=700
;bar=1

; for sph_1e36 M10 636
;!p.charsize=2.
;n=636
;bar=1
;sample=2
;xzoom=[-1.4e13,3.e12]
;zzoom=[-7.46e12,7.46e12]
;x0 = 170
;barx = 250

; for sph_1e35
;!p.charsize=2.
;cd,'/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e35/'
;n=700

; for sph_1e37
;cd,'/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37/'
;!p.charsize=2.
;n=491
;bar=1

if not keyword_set(out) then out='pp_showjet_1e36_'+string(n,format='(I4.4)')+'.eps'

if not keyword_set(var) then var='dens'
fname = 'JetSet_hdf5_plt_cnt_'+string(n,format='(I4.4)')
read_amr,fname,var='dens',parameters=params,tree=tree
time = params.time

jetOn = 9.76e4
time = time - jetOn

print,'time = ', time

; find position of BB
; temporary calculation to circular orbit. (in case of ellipse, need more tasks)
m1 = 4.d34
m2 = 2.d34
mtot = m1+m2
peri = 3.d12
th0 = !pi
th = th0 + sqrt(!unit.g * mtot / peri^3.d0) * time

th = th0
; center of mass, position of BB
cm2x = peri*m1/mtot*cos(th)
cm2y = peri*m1/mtot*sin(th)

; center of mass, position of Star
cm2xSt = peri*m2/mtot*cos(th+!pi)
cm2ySt = peri*m2/mtot*sin(th+!pi)

if keyword_set(bh) then begin
   xc0 = cm2x
   yc0 = cm2y
endif else begin
   if not keyword_set(xc0) then xc0 = 0.
   if not keyword_set(yc0) then yc0 = 0.
endelse
zc0 = 0.

if keyword_set(zoom) then begin
   if not keyword_set(dx) then dx=3.e12
   if not keyword_set(sample) then sample=0
;   d = loaddata(fname,var,xra=[xc0-dx,xc0+dx],yra=[yc0-dx,yc0+dx] $
;      ,zra=[zc0-dx,zc0+dx],time=time,sample=sample,xcoords=x,ycoords=y,zcoords=z)
   d = loaddata(fname,var,xra=xzoom,yra=yzoom $
      ,zra=zzoom,time=time,sample=sample,xcoords=x,ycoords=y,zcoords=z)
endif else begin
   if not keyword_set(sample) then sample=0
   ;d = dload(n,var=var,xc=xc0-0.2d12,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
   d = dload(n,var=var,xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
endelse

time = time - jetOn

sd = size(d)

if (var eq 'dens') then begin
;maxd = 1.e-13 & mind = 1.d-16
;maxd = 2.2e-13 & mind = 1.1d-15
maxd = 8.e-14 & mind = 3.d-18
endif else begin
maxd = max(d) & mind = min(d)
endelse

;;indx_tmp = where(x ge cm2x) & indx = indx_tmp[0]
;indy_tmp = where(y ge cm2y) & indy = indy_tmp[0]
indy_tmp = where(y ge yc0) & indy = indy_tmp[0]
indz_tmp = where(z ge zc0) & indz = indz_tmp[0]

ddxy = reform(d[*,*,indz])
ddxz = reform(d[*,indy,*])

if keyword_set(maskoff) then begin
print,'start Masking off the star'
;cord3 = replicate({x:0.,y:0.,z:0.},sd[1],sd[2],sd[3])
;for i=0,sd[1]-1 do cord3[i,*,*].x = x[i]
;for j=0,sd[2]-1 do cord3[*,j,*].y = y[j]
;for k=0,sd[3]-1 do cord3[*,*,k].z = z[k]

cordxy = replicate({x:fltarr(sd[1])},sd[2])
cordyx = replicate({y:fltarr(sd[2])},sd[1])

cordxy.x = x
xy_x = cordxy.x
cordyx.y = y
xy_y = cordyx.y & xy_y = transpose(xy_y)

Str_xy = sqrt((xy_x-cm2xSt)^2.+(xy_y-cm2ySt)^2.)

ddxy[where(Str_xy le 1.4e12)] = !values.f_nan

cordxz = replicate({x:fltarr(sd[1])},sd[3])
cordzx = replicate({z:fltarr(sd[3])},sd[1])

cordxz.x = x
xz_x = cordxz.x
cordzx.z = z
xz_z = cordzx.z & xz_z = transpose(xz_z)

Str_xz = sqrt((xz_x-cm2xSt)^2.+xz_z^2.)

ddxz[where(Str_xz le 1.4e12)] = !values.f_nan
print,'finish Masking off'
endif

if keyword_set(zoom) then begin
   if (sd[1] gt 512) then begin
      xs = sd[1] & ys=sd[2] & zs=sd[3]
      ddxy2 = ddxy
      ddxz2 = ddxz
      xcoord = x & ycoord=y & zcoord=z
   endif else begin
      xs = 512 & ys = 512 & zs = 512
      ddxy2 = congrid(ddxy,512,512)
      ddxz2 = congrid(ddxz,512,512)
      xcoord = findgen(512)/512.*(x[sd[1]-1]-x[0]) + x[0]
      ycoord = findgen(512)/512.*(y[sd[2]-1]-y[0]) + y[0]
      zcoord = findgen(512)/512.*(z[sd[3]-1]-z[0]) + z[0]
   endelse
endif else begin
   xs = sd[1] & ys = sd[2] & zs = sd[3]
   ddxy2 = ddxy
   ddxz2 = ddxz
   xcoord = x & ycoord = y & zcoord=z
endelse

loadct,0,/sil
if keyword_set(bar) then xs2 = xs + x0+x1+barx $
                    else xs2 = xs + x0+x1
zs2 = zs + y0+y1

if keyword_set(ps) then begin
;   psxs = 10. & psys=10.*float(zs2)/xs2
   psys = 15.3803 & psxs= psys*float(xs2)/float(zs2) 
   mkeps,out,xs=psxs,ys=psys
endif else window,1,xs=xs2,ys=zs2

plot,x,z,/xst,/yst,/nodata,position=posnorm([x0,y0,x0+xs,y0+zs],nx=xs2,ny=zs2),/norm $
    ,xtitle='x [cm]', ytitle='z [cm]',/iso $
    ,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]],xtickinterval=1.e13
if keyword_set(bar) then $
   color_bar, pos=posnorm([x0+xs+10,y0,x0+xs+40,y0+zs],nx=xs2,ny=zs2),/norm,/right,lim=[mind,maxd] $
   ,bartitle=textoidl('Density [g cm^{-3}]'),/log,titlegap=0.15,/minor

if keyword_set(maskoff) then begin
   maxd = max(ddxz) & mind = min(ddxz)
endif
;if keyword_set(nolog) then $
;   tv,bytscl(ddxz2,max=maxd,min=mind) $
; else $

if keyword_set(copy) then begin
   ddxz2[*,0:510] = reverse(ddxz2[*,512:1022],2)
endif 

if keyword_set(ps) then $
   tv,bytscl(alog10(ddxz2),max=alog10(maxd),min=alog10(mind)),float(x0)/xs2 *psxs,float(y0)/zs2*psys,/centimeter $
   ,xsize=float(xs)/xs2*psxs, ysize=float(zs)/zs2*psys $
  else $ 
   tv,bytscl(alog10(ddxz2),max=alog10(maxd),min=alog10(mind)),float(x0)/xs2,float(y0)/zs2,/norm 

plot,x,z,/xst,/yst,/nodata,position=posnorm([x0,y0,x0+xs,y0+zs],nx=xs2,ny=zs2),/norm $
    ,xtitle='x [cm]', ytitle='z [cm]',/iso ,/noerase $
    ,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]],xtickinterval=1.e13
plots,ring(1.e12,0.,1.4e12),/data,color=0
;legend,'X-Z',/top,/left,textcolor=255,box=0

;legend,string(time/60./60.,format='(f6.2)')+' hrs',/right,/top,box=0,textcolor=255

;if keyword_set(maskoff) then $
;   plots,ring(cm2xSt/!unit.au,0.,1.4e12/!unit.au),/data,line=1,color=0

; draw asymptotic lines
;inc   = 1.7
;inc   = 0.33
;inc   = 10.
;asymy = [inc*(x[0]-cm2xSt) + cm2ySt,-inc*(x[0]-cm2xSt) + cm2ySt]
;oplot, [x[0],cm2xSt],[asymy[0],cm2ySt],line=2,color=fsc_color('cyan')
;oplot, [x[0],cm2xSt],[asymy[1],cm2ySt],line=2,color=fsc_color('cyan')

;---------------------------------------------------------------------------------------
; show blocks
;---------------------------------------------------------------------------------------
if keyword_set(block) then begin
   for i=0L,params.totblocks-1 do begin
; find the box in which yc0 is located.   
       if ((tree[i].bndbox[0,1] le yc0) and (tree[i].bndbox[1,1] ge yc0)) then begin
          x1 = tree[i].bndbox[0,0] & x2 = tree[i].bndbox[1,0]
          z1 = tree[i].bndbox[0,2] & z2 = tree[i].bndbox[1,2]
          plots,[x1,x2],[z1,z1],/data,color=0,thick=1
          plots,[x2,x2],[z1,z2],/data,color=0,thick=1
          plots,[x1,x2],[z2,z2],/data,color=0,thick=1
          plots,[x1,x1],[z1,z2],/data,color=0,thick=1
       endif
   endfor
endif

if keyword_set(lineread) then begin
n=100
ang=45
ddx = 2.e11
x0 = xc0 & z0=0.
xx = x0 + findgen(n)/float(n)*ddx - ddx/2.
zz = tan(float(ang) *!pi/180.)*(xx-x0) + z0

;oplot,xx,zz,color=fsc_color('cyan')

xind = intarr(n) & zind = intarr(n) & linedata = fltarr(n) & ll = fltarr(n)
for i=0,n-1 do begin
    xind_tmp = where(xcoord ge xx[i]) & xind[i] = xind_tmp[0]
    zind_tmp = where(zcoord ge zz[i]) & zind[i] = zind_tmp[0]
    if ((xcoord[xind[i]]-x0) le 0) then ll[i] = -sqrt((xcoord[xind[i]]-x0)^2.+zcoord[zind[i]]^2.) $
                                   else ll[i] = sqrt((xcoord[xind[i]]-x0)^2.+zcoord[zind[i]]^2.)
    linedata[i] = ddxz2[xind[i],zind[i]]
endfor

plots,xcoord[xind],zcoord[zind],color=fsc_color('cyan'),line=0

window,2
plot,ll,linedata,/xst,xtitle='L [cm]',ytitle='density',xmargin=13.
oplot,[0.,0.],!y.crange,line=2

endif

if keyword_set(veloxz) then begin
   if keyword_set(zoom) then begin
      velx = loaddata(fname,'velx',xra=[xc0-dx,xc0+dx],yra=[yc0-dx,yc0+dx] $
          ,zra=[zc0-dx,zc0+dx],sample=sample)
      velz = loaddata(fname,'velz',xra=[xc0-dx,xc0+dx],yra=[yc0-dx,yc0+dx] $
          ,zra=[zc0-dx,zc0+dx],sample=sample)
   endif else begin
      ;d = dload(n,var=var,xc=xc0-0.2d12,yc=yc0,zc=zc0,sample=sample,time)
      velx = dload(n,var='velx',xc=xc0,yc=yc0,zc=zc0,sample=sample)
      velz = dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,sample=sample)
   endelse

   velx_xz = reform(velx[*,indy,*])
   velz_xz = reform(velz[*,indy,*])
   if keyword_set(zoom) then begin
      if (sd[1] gt 512) then begin
         velx_xz2 = velx_xz
         velz_xz2 = velz_xz
      endif else begin
         velx_xz2 = congrid(velx_xz,512,512)
         velz_xz2 = congrid(velz_xz,512,512)
      endelse
   endif else begin
      velx_xz2 = velx_xz
      velz_xz2 = velz_xz
   endelse

   vect_2d, velx_xz2, velz_xz2, xcoord, zcoord,vmax=6.e8

endif

if keyword_set(ps) then epsfree
stop
end
