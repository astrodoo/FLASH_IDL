pro pp_jet2

dir1 = '/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e35/'
dir2 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37/'

restore,file=dir1+'sph1e35_0700_ddxz.sav'
d35 = ddxz2 & x35 = x & z35 = z & t35 = time
restore,file=dir2+'sph1e37_0491_ddxz.sav'
d37 = ddxz2 & x37 = x & z37 = z & t37 = time

sz = size(d35,/dimension)

x0 = 280 & y0 = 150 
barx = 280
xbump = 200

xs = x0 + sz[0]*2 + xbump + barx
ys = y0 + sz[1] + 20 


psxs=30. & psys=psxs*float(ys)/float(xs)
mkeps,'pp_jet3_1e35_37',xs=psxs, ys=psys

maxd = 1.d-13 & mind = min(d37)
;maxd = 1.d-12 & mind = min(d971)

tv,bytscl(alog10(d35),max=alog10(maxd),min=alog10(mind)),float(x0)/xs*psxs, float(y0)/ys*psys,/centimeter $
  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
    ,/nodata,position=posnorm([x0,y0,x0+sz[0],y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
    ,xtitle='x [cm]', ytitle='z [cm]',/noerase
plots,ring(1.e12,0.,1.4e12),/data,color=0

tvlct,r,g,b,/get
; draw asymptotic lines
cm2xSt = 1.e12 & cm2ySt=0.
inc   = 0.33
asymy = [inc*(x[0]-cm2xSt) + cm2ySt,-inc*(x[0]-cm2xSt) + cm2ySt]
oplot, [x[0],cm2xSt],[asymy[0],cm2ySt],line=2,color=fsc_color('cyan')
oplot, [x[0],cm2xSt],[asymy[1],cm2ySt],line=2,color=fsc_color('cyan')
tvlct,r,g,b

tv,bytscl(alog10(d37),max=alog10(maxd),min=alog10(mind)),float(x0+sz[0]+xbump)/xs*psxs, float(y0)/ys*psys,/centimeter $
  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
    ,/nodata,position=posnorm([x0+sz[0]+xbump,y0,x0+2*sz[0]+xbump,y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
    ,xtitle='x [cm]',/noerase
plots,ring(1.e12,0.,1.4e12),/data,color=0

tvlct,r,g,b,/get
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
