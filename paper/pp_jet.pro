pro pp_jet

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_M10/'

restore,file=dir+'sph1e36_M10_636_ddxz.sav'
d35 = ddxz & x35 = x & z35 = z & t35 = time

sz = size(d35,/dimension)

x0 = 300 & y0 = 150 
barx = 280

xs = x0 + sz[0] + barx
ys = y0 + sz[1] + 20 


psxs=25. & psys=psxs*float(ys)/float(xs)
mkeps,'pp_jet_M10',xs=psxs, ys=psys

maxd = 1.d-13 & mind = min(d35)
;maxd = 1.d-12 & mind = min(d971)

tv,bytscl(alog10(d35),max=alog10(maxd),min=alog10(mind)),float(x0)/xs*psxs, float(y0)/ys*psys,/centimeter $
  ,xsize=float(sz[0])/xs*psxs, ysize=float(sz[1])/ys*psys 
plot,x,z,/xst,/yst,xr=[x[0],x[n_elements(x)-1]],yr=[z[0],z[n_elements(z)-1]] $
    ,/nodata,position=posnorm([x0,y0,x0+sz[0],y0+sz[1]],nx=xs,ny=ys),/norm, xtickinterval=1.e13 $
    ,xtitle='x [cm]', ytitle='z [cm]',/noerase, ytickinterval=4.e12
plots,ring(1.e12,0.,1.4e12),/data,color=0

color_bar, pos=posnorm([x0+sz[0]+10,y0,x0+sz[0]+40,y0+sz[1]],nx=xs,ny=ys),/norm,/right,lim=[mind,maxd] $
         , bartitle=textoidl('Density [g cm^{-3}]'),/log,titlegap=0.12,/minor

epsfree
stop
end
