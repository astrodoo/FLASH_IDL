pro pp_jet3

dir = '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/sphere_1e36/'

;restore,dir+'ddxz_0330.sav'
restore,dir+'mkdataSnap_0330.sav'
ddxz1 = dxz & t1 = time
;restore,dir+'ddxz_0370.sav'
restore,dir+'mkdataSnap_0370.sav'
ddxz2 = dxz & t2 = time
;restore,dir+'ddxz_0480.sav'
restore,dir+'mkdataSnap_0480.sav'
ddxz3 = dxz & t3 = time & vz3 = vz_xz & j3 = jxz &x3 = x & z3 = z

sz = size(dxz,/dimension)

maxd = max(ddxz3) & mind = min(ddxz3)

;====================================================================================================
; for transparent mark of jet materials
; find jet materials
crit=1.e9
vj3 = abs(vz3)*j3
jind = where(vj3 ge crit)

xjin = fltarr(sz[1]/2)
xjot = fltarr(sz[1]/2)
for i=0,sz[1]/2-1 do begin
   jjind = where(vj3[*,i] ge crit,count)
   if (count ne 0) then begin 
      jdummy = x3[jjind]  
      xjin[i] = jdummy[0]
      xjot[i] = jdummy[n_elements(jdummy)-1]
   endif else begin
      xjin[i] = !values.f_nan
      xjot[i] = !values.f_nan
   endelse
endfor 
dummy = xjin
zj = z3[0:sz[1]/2-1]
zj = zj[where(finite(dummy))]
xjin = xjin[where(finite(dummy))]
xjot = xjot[where(finite(dummy))]

loadct,0,/sil
window,xs=sz[0],ys=sz[1],/pixmap
tvcoord,bytscl(alog10(ddxz3),max=alog10(maxd),min=alog10(mind)),x3,z3,/scale
win1 = tvrd(/true)

tvlct,r,g,b,/get
oplot,xjin,zj,color=fsc_color('magenta')
oplot,xjot,zj,color=fsc_color('magenta')
polyfill,[xjin,reverse(xjot)],[zj,reverse(zj)],color=fsc_color('magenta')
tvlct,r,g,b
win2 = tvrd(/true)

alpha = 0.4
win3  = alpha*win2 + (1.-alpha)*win1


x0 = 280 & y0 = 150
barx = 280

xs = x0 + sz[0]*3 + barx
ys = y0 + sz[1] + 20 

psxs=30. & psys=30.*float(ys)/float(xs)
mkeps,'pp_jet3_1e36',xs=psxs, ys=psys


;tv,bytscl(alog10(ddxz1),max=alog10(maxd),min=alog10(mind)),float(x0)/xs*psxs, float(y0)/ys*psys,/centimeter $
;  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
;plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
;    ,/nodata,position=posnorm([x0,y0,x0+sz[0],y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
;    , ytitle='z [cm]',/noerase

tvcoord,bytscl(alog10(ddxz1),max=alog10(maxd),min=alog10(mind)),x,z $
    ,position=[float(x0)/xs,float(y0)/ys],psx=float(sz[0])/xs $
    ,/axes, xtickinterval=1.e13, ytitle='z [cm]'

plots,ring(1.e12,0.,1.4e12),/data,color=0
tjon=9.76e4
legend,string((t1-tjon)/60./60.,format='(f6.2)')+' hrs',box=0,textcolor=255,pos=[-6.e12,1.1e13],/data

;tv,bytscl(alog10(ddxz2),max=alog10(maxd),min=alog10(mind)),float(x0+sz[0])/xs*psxs, float(y0)/ys*psys,/centimeter $
;  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
;plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
;    ,/nodata,position=posnorm([x0+sz[0],y0,x0+2*sz[0],y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
;    ,xtitle='x [cm]',/noerase, ytickformat='(a1)'
tvcoord,bytscl(alog10(ddxz2),max=alog10(maxd),min=alog10(mind)),x,z $
    ,position=[(float(x0)+sz[0])/xs,float(y0)/ys],psx=float(sz[0])/xs $
    ,xtitle='x [cm]',/axes, xtickinterval=1.e13, ytickformat='(a1)'

plots,ring(1.e12,0.,1.4e12),/data,color=0
legend,string((t2-tjon)/60./60.,format='(f6.2)')+' hrs',box=0,textcolor=255,pos=[-6.e12,1.1e13],/data

;tv,bytscl(alog10(ddxz3),max=alog10(maxd),min=alog10(mind)),float(x0+2*sz[0])/xs*psxs, float(y0)/ys*psys,/centimeter $
;  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
;plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
;    ,/nodata,position=posnorm([x0+2*sz[0],y0,x0+3*sz[0],y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
;    ,/noerase, ytickformat='(a1)'

tvcoord,win3,/true,x,z,position=[(float(x0)+2*sz[0])/xs,float(y0)/ys],psx=float(sz[0])/xs $
    ,/axes, xtickinterval=1.e13, ytickformat='(a1)' 

plots,ring(1.e12,0.,1.4e12),/data,color=0
legend,string((t3-tjon)/60./60.,format='(f6.2)')+' hrs',box=0,textcolor=255,pos=[-6.e12,1.1e13],/data

;
tvlct,r,g,b,/get
; draw asymptotic lines
cm2xSt = 1.e12 & cm2ySt=0.
inc   = tan(50.*!dtor)
asymy = [inc*(x[0]-cm2xSt) + cm2ySt,-inc*(x[0]-cm2xSt) + cm2ySt]
oplot, [x[0],cm2xSt],[asymy[0],cm2ySt],line=2,color=fsc_color('cyan')
oplot, [x[0],cm2xSt],[asymy[1],cm2ySt],line=2,color=fsc_color('cyan')
tvlct,r,g,b


color_bar, pos=posnorm([x0+3*sz[0]+10,y0,x0+3*sz[0]+40,y0+sz[1]],nx=xs,ny=ys),/norm,/right,lim=[mind,maxd] $
         , bartitle=textoidl('Density [g cm^{-3}]'),/log,titlegap=0.08,/minor



;tvcoord,alog(d50),x50,z50,/scale

epsfree
stop
end
