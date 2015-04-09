forward_function Rst

pro guitar_head2

xra=[0.,3.e16]
yra=[-5.e15,7e16]
sample=1
fname='PWN2d_hdf5_plt_cnt_0500'

d = loaddata(fname+'_dl0.1','dens',sample=sample,xCoord=yy,yCoord=xx,xra=xra,yra=yra,time=tt)
spE = 1e33
d0 = 1.67d-24
d1 = d0*0.7
vs = 1.5d8
vp = 1.d10
Rst1 = Rst(d0=d0)
xoff1 = 8.e14

d = transpose(d)

d2 = loaddata(fname+'_dh0.1','dens',sample=sample,xCoord=yy,yCoord=xx,xra=xra,yra=yra)
d2 = transpose(d2)
d02 = d0/0.7
Rst2 = Rst(d0=d02)
xoff2 = 5.e14

dsz = size(d,/dimension)

xgr1 = -3.75d16 + vs*tt
xgr2 = -7.06d16 + vs*tt


;restore,file='analBow_0180.dat'
;x1 = x & w11 = w1 & w21 = w2
;restore,file='analBow_0180_d0.7.dat'
;x2 = x & w12 = w1 & w22 = w2

; drawing
x0 = 450 & y0 = 180

xs = x0 + dsz[0]+420
ys = y0 + dsz[1]*2 + 20
psxs=20. & psys=psxs*float(ys)/float(xs)

;analytic lines
restore,file='analBow2_0500.sav'
xan1 = x & w1an1 = w1
restore,file='analBow2_0500_0.7.sav'
xan07 = x & w1an07 = w1
restore,file='analBow2_0500_0.5.sav'
xan05 = x & w1an05 = w1
restore,file='analBow2_0500_0.3.sav'
xan03 = x & w1an03 = w1
restore,file='analBow2_0500_0.1.sav'
xan01 = x & w1an01 = w1

restore,file='analBow2_0500_3.sav'
xan3 = x & w1an3 = w1
restore,file='analBow2_0500_5.sav'
xan5 = x & w1an5 = w1
restore,file='analBow2_0500_7.sav'
xan7 = x & w1an7 = w1
restore,file='analBow2_0500_10.sav'
xan10 = x & w1an10 = w1

loadct,0,/sil
mkeps,'guitar_head2',xs=psxs,ys=psys

maxd = max(d) & mind = min(d)
tv,bytscl(alog10(d),max=alog10(maxd),min=alog10(mind)),float(x0)/xs*psxs, float(y0+dsz[1])/ys*psys,/centimeter $
  ,xsize=float(dsz[0])/xs*psxs, ysize=float(dsz[1])/ys*psys
plot,xx,yy,/xst,/yst,xr=[xx[0],xx[n_elements(xx)-1]],yr=[yy[0],yy[n_elements(yy)-1]] $
    ,/nodata,position=posnorm([x0,y0+dsz[1],x0+dsz[0],y0+dsz[1]*2],nx=xs,ny=ys),/norm $
    ,/noerase, xtickformat='(a1)', ytickinterval=1.e16
oplot,[xgr1,xgr1],!y.crange,line=1
oplot,[xgr2,xgr2],!y.crange,line=1
anacl = 'magenta'
oplot,xan1[where(xan1 le 6.e16)]-Rst(d0=d0)-xoff1,w1an1[where(xan1 le 6.e16)],color=fsc_color(anacl)
oplot,xan07[where(xan07 le 5.e16)]-Rst(d0=d0*0.7)-xoff1,w1an07[where(xan07 le 5.e16)],color=fsc_color(anacl)
oplot,xan05[where(xan05 le 4.e16)]-Rst(d0=d0*0.5)-xoff1,w1an05[where(xan05 le 4.e16)],color=fsc_color(anacl)
oplot,xan03[where(xan03 le 3.e16)]-Rst(d0=d0*0.3)-xoff1,w1an03[where(xan03 le 3.e16)],color=fsc_color(anacl)
oplot,xan01[where(xan01 le 2.e16)]-Rst(d0=d0*0.1)-xoff1,w1an01[where(xan01 le 2.e16)],color=fsc_color(anacl)
arrow,3.e16,2.e16,1.5e16,2.e16,/data
xyouts,1.8e16,2.2e16,/data,textoidl('\nabla\rho_{0} < 0')
;oplot,x1-Rst1-xoff1,w11
;oplot,x1-Rst1-xoff1,w21
;xyouts,1.e15,8.5e15,/data,textoidl('\rho_{1}')

loadct,0,/sil
tv,bytscl(alog10(reverse(d2,2)),max=alog10(max(d2)),min=alog10(min(d2))),float(x0)/xs*psxs, float(y0)/ys*psys,/centimeter $
  ,xsize=float(dsz[0])/xs*psxs, ysize=float(dsz[1])/ys*psys
plot,xx,-yy,/xst,/yst,xr=[xx[0],xx[n_elements(xx)-1]],yr=[-yy[n_elements(yy)-1],-yy[0]] $
    ,/nodata,position=posnorm([x0,y0,x0+dsz[0],y0+dsz[1]],nx=xs,ny=ys),/norm $
    , xtitle='x [cm]',/noerase, ytickinterval=1.e16
oplot,[xgr1,xgr1],!y.crange,line=1
oplot,[xgr2,xgr2],!y.crange,line=1
oplot,xan1[where(xan1 le 6.e16)]-Rst(d0=d0)-xoff1,-w1an1[where(xan1 le 6.e16)],color=fsc_color(anacl)
oplot,xan3[where(xan3 le 5.e16)]-Rst(d0=d0*3)-xoff1,-w1an3[where(xan3 le 5.e16)],color=fsc_color(anacl)
oplot,xan5[where(xan5 le 4.e16)]-Rst(d0=d0*5)-xoff1,-w1an5[where(xan5 le 4.e16)],color=fsc_color(anacl)
oplot,xan7[where(xan7 le 3.e16)]-Rst(d0=d0*7)-xoff1,-w1an7[where(xan7 le 3.e16)],color=fsc_color(anacl)
oplot,xan10[where(xan10 le 2.e16)]-Rst(d0=d0*10)-xoff1,-w1an10[where(xan10 le 2.e16)],color=fsc_color(anacl)
arrow,3.e16,-2.e16,1.5e16,-2.e16,/data
xyouts,1.8e16,-2.4e16,/data,textoidl('\nabla\rho_{0} > 0')
;
;oplot,x2-Rst2-xoff2,-w12
;oplot,x2-Rst2-xoff2,-w22
;xyouts,1.e15,-8.5e15,/data,textoidl('\rho_{2}=0.7\rho_{1}')

xyouts,0.04,0.45,/norm,'y [cm]',orientation=90

loadct,0,/sil

color_bar,/right,position=posnorm([x0+dsz[0]+20,y0,x0+dsz[0]+70,y0+dsz[1]*2],nx=xs,ny=ys),/norm,lim=[min(d2),max(d2)],/log, bartitle=textoidl('\rho [g cm^{-3}]'),titlegap=0.12


epsfree

stop
end

function Rst,spE=spE,d0=d0,vs=vs,vp=vp
if not keyword_set(spE) then spE = 1e33
if not keyword_set(d0) then d0 = 1.67d-24
if not keyword_set(vs) then vs = 1.5d8
if not keyword_set(vp) then vp = 1.d10

return, sqrt(spE / (4.d0*!dpi*d0*vs^2.*vp))
end
