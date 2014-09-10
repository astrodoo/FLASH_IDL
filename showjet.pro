pro showjet,n,sample=sample,block=block,bcolor=bcolor,zoom=zoom,dx=dx,lineread=lineread,var=var $
   ,noorbit=noorbit,nolog=nolog,veloxy=veloxy,veloxz=veloxz,bh=bh,star=star,xc0=xc0,yc0=yc0,zc0=zc0 $
   ,maskoff=maskoff, out=out, totv=totv, ct=ct, chk=chk, rotback=rotback, ionparam=ionparam
device,decomposed=0

if (n_elements(chk) eq 0) then chk=0

if keyword_set(bcolor) then tvlct,r,g,b,/get
if not keyword_set(ct) then ct=0

if not keyword_set(var) then var='dens'
if (where(strmatch(['velx','vely','velz','totv'],var) eq 1) ne -1) then nolog=1 

if keyword_set(chk) then fname = 'JetSet_hdf5_chk_'+string(n,format='(I4.4)') $
   else fname = 'JetSet_hdf5_plt_cnt_'+string(n,format='(I4.4)')

if keyword_Set(chk) then print,'check point file'
read_amr,fname,var='dens',parameters=params,tree=tree ,/nodata
time = params.time

;time=time-103197.44

print,'time = ', time

; find position of BB
; temporary calculation to circular orbit. (in case of ellipse, need more tasks)
m1 = 4.d34
m2 = 2.d34
mtot = m1+m2
peri = 3.d12
th0 = !pi
th = th0 + sqrt(!unit.g * mtot / peri^3.d0) * time
;th = th0

if keyword_Set(rotback) then begin
   dth = (th-th0) / !dtor
   bh=0 
endif

; center of mass, position of BB
cm2x = peri*m1/mtot*cos(th)
cm2y = peri*m1/mtot*sin(th)

; center of mass, position of Star
cm2xSt = peri*m2/mtot*cos(th+!pi)
cm2ySt = peri*m2/mtot*sin(th+!pi)

;xc0 = -1992463065262.43
;yc0 = -173467384729.911
if keyword_set(bh) then begin
   xc0 = cm2x
   yc0 = cm2y
   zc0 = 0.
endif else if keyword_set(star) then begin
   xc0 = cm2xSt
   yc0 = cm2ySt
   zc0 = 0.
endif else begin
   if not keyword_set(xc0) then xc0 = 0.
   if not keyword_set(yc0) then yc0 = 0.
   if not keyword_set(zc0) then zc0 = 0.
endelse

if (var eq 'totv') then begin
   vx=dload(n,var='velx',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time,lref=lref,chk=chk)
   vy=dload(n,var='vely',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time,chk=chk)
   vz=dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time,chk=chk)
   d = sqrt(vx*vx + vy*vy + vz*vz)
   nolog=1
endif else if (var eq 'velxy') then begin
   vx=dload(n,var='velx',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time,lref=lref,chk=chk)
   vy=dload(n,var='vely',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time,chk=chk)
   d = sqrt(vx*vx + vy*vy)
   nolog=1
endif else begin

if keyword_set(zoom) then begin
   if not keyword_set(dx) then dx=3.e12
   if not keyword_set(sample) then sample=0
   d = loaddata(fname,var,xra=[xc0-dx,xc0+dx],yra=[yc0-dx,yc0+dx] $
      ,zra=[zc0-dx,zc0+dx],time=time,sample=sample,xcoords=x,ycoords=y,zcoords=z,lref=lref)
endif else begin
   if not keyword_set(sample) then sample=0
   ;d = dload(n,var=var,xc=xc0-0.2d12,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
   d = dload(n,var=var,xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time,lref=lref,chk=chk)
endelse

endelse

if keyword_set(rotback) then begin
   rotd = rot_3d(d,rotaxis=3,degree=dth,/interp)
   d= rotd
endif

sd = size(d)

;if (var eq 'dens') then begin
;maxd = 2.2e-13 & mind = 1.1d-15
;endif else begin
maxd = max(d) & mind = min(d)
;endelse

;;indx_tmp = where(x ge cm2x) & indx = indx_tmp[0]
;indy_tmp = where(y ge cm2y) & indy = indy_tmp[0]
indy_tmp = where(y ge yc0) & indy = indy_tmp[0]
indz_tmp = where(z ge zc0) & indz = indz_tmp[0]

ddxy = reform(d[*,*,indz])
ddxz = reform(d[*,indy,*])


; coordinate into 2d
cordxy = replicate({x:fltarr(sd[1])},sd[2])
cordyx = replicate({y:fltarr(sd[2])},sd[1])

cordxy.x = x
xy_x = cordxy.x
cordyx.y = y
xy_y = cordyx.y & xy_y = transpose(xy_y)

if keyword_set(maskoff) then begin
print,'start Masking off the star'
;cord3 = replicate({x:0.,y:0.,z:0.},sd[1],sd[2],sd[3])
;for i=0,sd[1]-1 do cord3[i,*,*].x = x[i]
;for j=0,sd[2]-1 do cord3[*,j,*].y = y[j]
;for k=0,sd[3]-1 do cord3[*,*,k].z = z[k]

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

loadct,ct,/sil
window,0,xs=xs,ys=ys
plot,x,y,/xst,/yst,/nodata,xr=[x[0],x[sd[1]-1]],yr=[y[0],y[sd[2]-1]],/iso,position=[0,0,xs,ys],/dev
if keyword_set(maskoff) then begin 
   maxd = max(ddxy) & mind = min(ddxy)
endif
if keyword_set(nolog) then $
   tv,bytscl(ddxy2,max=maxd,min=mind) $
 else $
   tv,bytscl(alog10(ddxy2),max=alog10(maxd),min=alog10(mind))
legend,'X-Y',/top,/left,textcolor=255,box=0

timestr = string(sqrt(!unit.g * mtot / peri^3.d0)/2./!pi*time,format='(f4.2)')+' Period'
legend, timestr,/right,/top,box=0
;legend,string(sqrt(!unit.g * mtot / peri^3.d0)/2./!pi*(time-9.76e4),format='(f4.2)')+' Period',/right,/top,box=0
;legend,string(sqrt(!unit.g * mtot / peri^3.d0)/2./!pi*(time-515987.19),format='(f4.2)')+' Period',/right,/top,box=0
;plots,ring(0.,0.,1.4e12),color=0,/data,thick=1
;plots,ring(0.,0.,1.4e12*1.25),color=0,/data,line=2,thick=1
if not keyword_set(noorbit) then overorbit
if keyword_set(maskoff) then $
   plots,ring(cm2xSt,cm2ySt,1.4e12),color=fsc_color('cyan'),/data

if keyword_set(block) then showblock,tree,params,zc=zc0

if keyword_set(veloxy) then begin
   if keyword_set(zoom) then begin
      velx = loaddata(fname,'velx',xra=[xc0-dx,xc0+dx],yra=[yc0-dx,yc0+dx] $
          ,zra=[zc0-dx,zc0+dx],sample=sample)
      vely = loaddata(fname,'vely',xra=[xc0-dx,xc0+dx],yra=[yc0-dx,yc0+dx] $
          ,zra=[zc0-dx,zc0+dx],sample=sample)
   endif else begin
      ;d = dload(n,var=var,xc=xc0-0.2d12,yc=yc0,zc=zc0,sample=sample,time)
      velx = dload(n,var='velx',xc=xc0,yc=yc0,zc=zc0,sample=sample,chk=chk)
      vely = dload(n,var='vely',xc=xc0,yc=yc0,zc=zc0,sample=sample,chk=chk)
   endelse

   velx_xy = reform(velx[*,*,indz])
   vely_xy = reform(vely[*,*,indz])
   if keyword_set(zoom) then begin
      if (sd[1] gt 512) then begin
         velx_xy2 = velx_xy
         vely_xy2 = vely_xy
      endif else begin
         vely_xy2 = congrid(velx_xy,512,512)
         vely_xy2 = congrid(vely_xy,512,512)
      endelse
   endif else begin
      velx_xy2 = velx_xy
      vely_xy2 = vely_xy
   endelse

   vect_2d, velx_xy2, vely_xy2, xcoord, ycoord,vmax=6.e8

endif

if keyword_set(ionparam) then begin

  Jet_xy2 = (xy_x-cm2x)^2.+(xy_y-cm2y)^2.
  Lx = 1.d36

  mu=1.
  ionparams = Lx * mu*!unit.mh / ddxy/ Jet_xy2
 
  dsize=size(ddxy,/dimension)
  
  ioncrit100 = fltarr(dsize[0],dsize[1])
  ioncrit100[where(ionparams ge 100)]=1

  tvlct,r,g,b,/get
  yellow = transpose(fsc_color('yellow',/triple))
  tvlct,yellow[0],yellow[1],yellow[2],255
  
  contour,ionparams,x,y,levels=100,color=fsc_color('yellow'),/overplot, thick=1.5
  window,10,xs=dsize[0],ys=dsize[1],/pixmap
  tvscl,ioncrit100
  dumbImg = tvrd(true=1)
  dumbImg2 = cgTransparentImage(dumbImg, missing_value='black', transparent=80)
  tvlct,r,g,b
  wset,0
  cgImage,dumbImg2

endif ; ionparam


;oplot,!x.crange, [cm2y,cm2y],line=2,color=0
;plots,cm2x,cm2y,psym=1,color=0,/data
;jetrad = 1.25e10
;plots,ring(cm2x,cm2y,3.*jetrad),color=255

loadct,ct,/sil
window,1,xs=xs,ys=zs
;mkeps,'tmp_vel',xs=20.,ys=20.
plot,x,z,/xst,/yst,/nodata,xr=[x[0],x[sd[1]-1]],yr=[z[0],z[sd[3]-1]],/iso,position=[0,0,xs,zs],/dev
;plot,x,z,/xst,/yst,/nodata,xr=[x[0],x[sd[1]-1]],yr=[z[0],z[sd[3]-1]],/iso,position=[0,0,1,1],/norm
;tvscl,alog10(d[*,indy,*])

if keyword_set(maskoff) then begin
   maxd = max(ddxz) & mind = min(ddxz)
endif

maxd=max(ddxz2) & mind = min(ddxz2)

if keyword_set(nolog) then $
   tv,bytscl(ddxz2,max=maxd,min=mind) $
 else $
   tv,bytscl(alog10(ddxz2),max=alog10(maxd),min=alog10(mind))
legend,'X-Z',/top,/left,textcolor=255,box=0
legend, timestr,/right,/top,box=0
;legend,string(time/60./60./24.,format='(f4.2)')+' days',/right,/top,box=0,color=0
if keyword_set(maskoff) then $
   plots,ring(cm2xSt,0.,1.4e12),color=fsc_color('cyan'),/data

if keyword_set(out) then begin
   id='smp0_2'
   save,filename=var+'_'+string(n,format='(I4.4)')+'_'+id+'.sav',x,z,ddxz
   stop
endif   

if keyword_set(block) then showblock,tree,params,yc=yc0

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
      velx = dload(n,var='velx',xc=xc0,yc=yc0,zc=zc0,sample=sample,chk=chk)
      velz = dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,sample=sample,chk=chk)
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
;epsfree
stop
end

pro showblock, tree, params, yc=yc, zc=zc, color=color
if not keyword_set(color) then color=0
for i=0L,params.totblocks-1 do begin
    if (n_elements(yc) ne 0) then begin
       if ((tree[i].bndbox[0,1] le yc) and (tree[i].bndbox[1,1] ge yc)) then begin
          x1 = tree[i].bndbox[0,0] & x2 = tree[i].bndbox[1,0]
          z1 = tree[i].bndbox[0,2] & z2 = tree[i].bndbox[1,2]
          plots,[x1,x2,x2,x1,x1],[z1,z1,z2,z2,z1],/data,color=color,thick=1
       endif
    endif else if (n_elements(zc) ne 0) then begin
       if ((tree[i].bndbox[0,2] le zc) and (tree[i].bndbox[1,2] ge zc)) then begin
          x1 = tree[i].bndbox[0,0] & x2 = tree[i].bndbox[1,0]
          y1 = tree[i].bndbox[0,1] & y2 = tree[i].bndbox[1,1]
          plots,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],/data,color=color,thick=1
       endif
    endif
endfor
end
