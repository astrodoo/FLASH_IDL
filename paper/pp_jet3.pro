pro pp_jet3

dir1 = '/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36/'
dir2 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_high/'

restore,dir1+'sph1e36_0050_ddxz.sav'
d50 = ddxz & x50 = x & z50 = z & t50 = time
restore,dir1+'sph1e36_0300_ddxz.sav'
d300 = ddxz & x300 = x & z300 = z & t300 = time
restore,dir2+'sph1e36_0971_ddxz.sav'
d971 = ddxz & x971 = x & z971 = z & t971 = time

sz = size(d50,/dimension)

x0 = 280 & y0 = 150
barx = 280

xs = x0 + sz[0]*3 + barx
ys = y0 + sz[1] + 20 


psxs=30. & psys=30.*float(ys)/float(xs)
mkeps,'pp_jet3_1e36',xs=psxs, ys=psys

maxd = max(d971) & mind = min(d971)
;maxd = 1.d-12 & mind = min(d971)

tv,bytscl(alog10(d50),max=alog10(maxd),min=alog10(mind)),float(x0)/xs*psxs, float(y0)/ys*psys,/centimeter $
  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
    ,/nodata,position=posnorm([x0,y0,x0+sz[0],y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
    ,xtitle='x [cm]', ytitle='z [cm]',/noerase
plots,ring(1.e12,0.,1.4e12),/data,color=0
legend,string(t50/60./60.,format='(f6.2)')+' hrs',box=0,textcolor=255,pos=[-8.e12,1.8e13],/data

tv,bytscl(alog10(d300),max=alog10(maxd),min=alog10(mind)),float(x0+sz[0])/xs*psxs, float(y0)/ys*psys,/centimeter $
  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
    ,/nodata,position=posnorm([x0+sz[0],y0,x0+2*sz[0],y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
    ,xtitle='x [cm]',/noerase, ytickformat='(a1)'
plots,ring(1.e12,0.,1.4e12),/data,color=0
legend,string(t300/60./60.,format='(f6.2)')+' hrs',box=0,textcolor=255,pos=[-8.e12,1.8e13],/data

tv,bytscl(alog10(d971),max=alog10(maxd),min=alog10(mind)),float(x0+2*sz[0])/xs*psxs, float(y0)/ys*psys,/centimeter $
  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
    ,/nodata,position=posnorm([x0+2*sz[0],y0,x0+3*sz[0],y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
    ,xtitle='x [cm]',/noerase, ytickformat='(a1)'
plots,ring(1.e12,0.,1.4e12),/data,color=0
legend,string(t971/60./60.,format='(f6.2)')+' hrs',box=0,textcolor=255,pos=[-8.e12,1.8e13],/data

tvlct,r,g,b,/get
; draw asymptotic lines
cm2xSt = 1.e12 & cm2ySt=0.
inc   = 1.7
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
