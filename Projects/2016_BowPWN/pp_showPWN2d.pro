pro pp_showPWN2d,ps=ps,velo=velo

nth = 300
th0 = 0.1 * !dtor
th1 = 180. * !dtor
th  = findgen(nth) / nth * (th1-th0) + th0

Esp = 6.66d35
dam = 1.67d-24
vw  = 6.d7
vp  = 1.d10

r0 = sqrt(Esp / (4.*!pi*dam*vw*vw*vp))
rth = r0 / sin(th) * sqrt(3.*(1.-th/tan(th)))

xx = - rth * cos(th)
yy = rth * sin(th)

xrange = [-2.e17,5.e17] & yrange = [0.,4.e17]
dir='/d/d7/yoon/out_FLASH3.3_mhd/out_PWN2d/out_PWN2d_600kms/'
d = loaddata(dir+'PWN2d_hdf5_plt_cnt_0150','dens',xCoord=y,yCoord=x,lref=lref,time=t,yra=xrange,xra=yrange)
d2 = transpose(d)

if keyword_set(velo) then begin
vx = loaddata(dir+'PWN2d_hdf5_plt_cnt_0150','vely',yra=xrange,xra=yrange)
vy = loaddata(dir+'PWN2d_hdf5_plt_cnt_0150','velx',yra=xrange,xra=yrange)
endif 


loadct,0,/sil
if not keyword_set(ps) then window,0,xs=800,ys=650 $
 else mkeps,'pp_showPWN2d_600kms_0150',xs=20.,ys=20.*40./80.

if keyword_set(ps) then $
   tvcoord,alog(d2),x,y,/scale,/on,xtitle='x [cm]',ytitle='y [cm]',charsize=1.5,position=[0.17,0.2],imgsize=0.63 $
else $
   tvcoord,alog(d2),x,y,/scale,/on,xtitle='x [cm]',ytitle='y [cm]' 

plots,ring(0,0,1.24e16),/data,color=255
oplot,xx,yy,color=fsc_color('magenta'),line=2

r0 = 7.e16
rth = r0 / sin(th) * sqrt(3.*(1.-th/tan(th)))
xx = - rth * cos(th)
yy = rth * sin(th)
oplot,xx,yy,color=fsc_color('magenta'),line=2

if keyword_set(velo) then vect_2d,vx,vy,x,y,sample=10.,vmax=1.e9

loadct,0,/sil
color_bar,lim=[min(d2),max(d2)],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),titlegap=0.12,charsize=1.5







if keyword_set(ps) then epsfree
stop
end
