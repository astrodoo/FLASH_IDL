pro guitar_head

xra=[0.,1.e16]
yra=[-2.e15,2.e16]
sample=0
fname='PWN2d_hdf5_plt_cnt_0180'

d = loaddata(fname,'dens',sample=sample,xCoord=yy,yCoord=xx,xra=xra,yra=yra)
spE = 1e33
d0 = 1.67d-24
vs = 1.5d8
vp = 1.d10
Rst1 = sqrt(spE / (4.d0*!dpi*d0*vs^2.*vp))
stR0 = sqrt(spE / (4.*!pi*d0*vs*vs*vp))
xoff1 = 4.e14

d = transpose(d)

d2 = loaddata(fname+'_d0.7','dens',sample=sample,xCoord=yy,yCoord=xx,xra=xra,yra=yra)
d2 = transpose(d2)
d02 = d0*0.7
stR02 = sqrt(spE / (4.*!pi*d02*vs*vs*vp))
Rst2 = sqrt(spE / (4.d0*!dpi*d02*vs^2.*vp))
xoff2 = 5.e14

dsz = size(d,/dimension)


restore,file='analBow_0180.dat'
x1 = x & w11 = w1 & w21 = w2
restore,file='analBow_0180_d0.7.dat'
x2 = x & w12 = w1 & w22 = w2

; drawing
x0 = 200 & y0 = 100

xs = x0 + dsz[0]+20
ys = y0 + dsz[1]*2 + 20
psxs=20. & psys=psxs*float(ys)/float(xs)

loadct,0,/sil
mkeps,'guitar_head',xs=psxs,ys=psys

maxd = max(d) & mind = min(d)
tv,bytscl(alog10(d),max=alog10(maxd),min=alog10(mind)),float(x0)/xs*psxs, float(y0+dsz[1])/ys*psys,/centimeter $
  ,xsize=float(dsz[0])/xs*psxs, ysize=float(dsz[1])/ys*psys
plot,xx,yy,/xst,/yst,xr=[xx[0],xx[n_elements(xx)-1]],yr=[yy[0],yy[n_elements(yy)-1]] $
    ,/nodata,position=posnorm([x0,y0+dsz[1],x0+dsz[0],y0+dsz[1]*2],nx=xs,ny=ys),/norm $
    , ytitle='y [cm]',/noerase, xtickformat='(a1)'
oplot,[-stR0,-stR0],!y.crange
oplot,x1-Rst1-xoff1,w11
oplot,x1-Rst1-xoff1,w21
xyouts,1.e15,8.5e15,/data,textoidl('\rho_{1}')


tv,bytscl(alog10(reverse(d2,2)),max=alog10(maxd),min=alog10(mind)),float(x0)/xs*psxs, float(y0)/ys*psys,/centimeter $
  ,xsize=float(dsz[0])/xs*psxs, ysize=float(dsz[1])/ys*psys
plot,xx,-yy,/xst,/yst,xr=[xx[0],xx[n_elements(xx)-1]],yr=[-yy[n_elements(yy)-1],-yy[0]] $
    ,/nodata,position=posnorm([x0,y0,x0+dsz[0],y0+dsz[1]],nx=xs,ny=ys),/norm $
    , xtitle='x [cm]', ytitle='y [cm]',/noerase
oplot,[-stR02,-stR02],!y.crange
oplot,x2-Rst2-xoff2,-w12
oplot,x2-Rst2-xoff2,-w22
xyouts,1.e15,-8.5e15,/data,textoidl('\rho_{2}=0.7\rho_{1}')

epsfree

loadct,39,/sil
mkeps,'guitar_head2',xs=psxs,ys=psxs*4./8.

plot,x1-Rst1-xoff1,w11,/iso,xtitle='x [cm]',ytitle='y [cm]',xtickinterval=1.e16,xra=[-5.e15,2.e16],/xst,/yst
oplot,x1-Rst1-xoff1,w21
oplot,x2-Rst2-xoff2,w12,color=250
oplot,x2-Rst2-xoff2,w22,color=250

legend,[textoidl('\rho_{1}'),textoidl('\rho_{2}')],/left,/top,color=[0,250],textcolor=[0,250],box=0,line=0

epsfree

stop
end
