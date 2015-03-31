pro pp_jet2

dir1 = '/d/d2/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/sphere_1e35/'
dir2 = '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/sphere_1e37/'

;restore,file=dir1+'ddxz_0700.sav'
restore,file=dir1+'mkdataSnap_0700.sav'
d35 = dxz & x35 = x & z35 = z & t35 = time & j35 = jxz & vz35 = vz_xz
;restore,file=dir2+'ddxz_0463.sav'
restore,file=dir2+'mkdataSnap_0463.sav'
d37 = dxz & x37 = x & z37 = z & t37 = time & j37 = jxz & vz37 = vz_xz

sz = size(d35,/dimension)

; find jet materials
crit=1.e8
vj35 = abs(vz35)*j35
jind35 = where(vj35 ge crit)

xjin35 = fltarr(sz[1]/2)
xjot35 = fltarr(sz[1]/2)
for i=0,sz[1]/2-1 do begin
   jjind = where(vj35[*,i] ge crit,count)
   if (count ne 0) then begin 
      jdummy = x[jjind]  
      xjin35[i] = jdummy[0]
      xjot35[i] = jdummy[n_elements(jdummy)-1]
   endif else begin
      xjin35[i] = !values.f_nan
      xjot35[i] = !values.f_nan
   endelse
endfor 
dummy = xjin35
zj35 = z35[0:sz[1]/2-1]
zj35 = zj35[where(finite(dummy))]
xjin35 = xjin35[where(finite(dummy))]
xjot35 = xjot35[where(finite(dummy))]


crit=2.9e9
vj37 = abs(vz37)*j37
jind37 = where(vj37 ge crit)

xjin37 = fltarr(sz[1]/2)
xjot37 = fltarr(sz[1]/2)
for i=0,sz[1]/2-1 do begin
   jdummy = x[where(vj37[*,i] ge crit)]  
   xjin37[i] = jdummy[0]
   xjot37[i] = jdummy[n_elements(jdummy)-1]
endfor 
zj37 = z37[0:sz[1]/2-1]



x0 = 280 & y0 = 150 
barx = 280
xbump = 200

xs = x0 + sz[0]*2 + xbump + barx
ys = y0 + sz[1] + 20 


maxd = 1.d-13 & mind = min(d37)
;maxd = 1.d-12 & mind = min(d971)
;====================================================================================================
; for transparent mark of jet materials

loadct,0,/sil
window,xs=sz[0],ys=sz[1],/pixmap
tvcoord,bytscl(alog10(d35),max=alog10(maxd),min=alog10(mind)),x35,z35,/scale
win1 = tvrd(/true)

tvlct,r,g,b,/get
oplot,xjin35,zj35,color=fsc_color('magenta')
oplot,xjot35,zj35,color=fsc_color('magenta')
polyfill,[xjin35,reverse(xjot35)],[zj35,reverse(zj35)],color=fsc_color('magenta')
tvlct,r,g,b
win2 = tvrd(/true)

alpha = 0.5
img_35  = alpha*win2 + (1.-alpha)*win1


window,xs=sz[0],ys=sz[1],/pixmap
tvcoord,bytscl(alog10(d37),max=alog10(maxd),min=alog10(mind)),x37,z37,/scale
win1 = tvrd(/true)

tvlct,r,g,b,/get
oplot,xjin37,zj37,color=fsc_color('magenta')
oplot,xjot37,zj37,color=fsc_color('magenta')
polyfill,[xjin37,reverse(xjot37)],[zj37,reverse(zj37)],color=fsc_color('magenta')
tvlct,r,g,b
win2 = tvrd(/true)

alpha = 0.2
img_37  = alpha*win2 + (1.-alpha)*win1

loadct,0,/sil
psxs=30. & psys=psxs*float(ys)/float(xs)
mkeps,'pp_jet2_1e35_37',xs=psxs, ys=psys


;tv,bytscl(alog10(d35),max=alog10(maxd),min=alog10(mind)),float(x0)/xs*psxs, float(y0)/ys*psys,/centimeter $
;  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
;plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
;    ,/nodata,position=posnorm([x0,y0,x0+sz[0],y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
;    ,xtitle='x [cm]', ytitle='z [cm]',/noerase

tvcoord,img_35,/true,x,z,position=[float(x0)/xs,float(y0)/ys],psx=float(sz[0])/xs $
    ,xtitle='x [cm]',/axes, xtickinterval=1.e13 , ytitle='z [cm]'

tvlct,r,g,b,/get

oplot,xjin35,zj35,color=fsc_color('magenta')
oplot,xjot35,zj35,color=fsc_color('magenta')
plots,ring(1.e12,0.,1.4e12),/data,color=0

; draw asymptotic lines
cm2xSt = 1.e12 & cm2ySt=0.
inc   = 0.33
asymy = [inc*(x[0]-cm2xSt) + cm2ySt,-inc*(x[0]-cm2xSt) + cm2ySt]
oplot, [x[0],cm2xSt],[asymy[0],cm2ySt],line=2,color=fsc_color('cyan')
oplot, [x[0],cm2xSt],[asymy[1],cm2ySt],line=2,color=fsc_color('cyan')


tvlct,r,g,b



;tv,bytscl(alog10(d37),max=alog10(maxd),min=alog10(mind)),float(x0+sz[0]+xbump)/xs*psxs, float(y0)/ys*psys,/centimeter $
;  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 

;tv,win3,/true,float(x0+sz[0]+xbump)/xs*psxs, float(y0)/ys*psys,/centimeter $
;  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
;plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
;    ,/nodata,position=posnorm([x0+sz[0]+xbump,y0,x0+2*sz[0]+xbump,y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
;    ,xtitle='x [cm]',/noerase

tvlct,r,g,b,/get
tvcoord,img_37,/true,x,z,position=[float(x0+sz[0]+xbump)/xs,float(y0)/ys],psx=float(sz[0])/xs $
    ,xtitle='x [cm]',/axes, xtickinterval=1.e13 

oplot,xjin37,zj37,color=fsc_color('magenta')
oplot,xjot37,zj37,color=fsc_color('magenta')

plots,ring(1.e12,0.,1.4e12),/data,color=0

; draw asymptotic lines
cm2xSt = 1.e12 & cm2ySt=0.
inc   = 10.
asymy = [inc*(x[0]-cm2xSt) + cm2ySt,-inc*(x[0]-cm2xSt) + cm2ySt]
oplot, [x[0],cm2xSt],[asymy[0],cm2ySt],line=2,color=fsc_color('cyan')
oplot, [x[0],cm2xSt],[asymy[1],cm2ySt],line=2,color=fsc_color('cyan')
tvlct,r,g,b


color_bar, pos=posnorm([x0+2*sz[0]+10+xbump,y0,x0+2*sz[0]+40+xbump,y0+sz[1]],nx=xs,ny=ys),/norm,/right,lim=[mind,maxd] $
         , bartitle=textoidl('Density [g cm^{-3}]'),/log,titlegap=0.08,/minor



;tvcoord,alog(d50),x50,z50,/scale

epsfree
stop
end
