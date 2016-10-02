pro pp_showPWN2d_guitar,ps=ps

dir='/d/d7/yoon/out_FLASH3.3_mhd/out_PWN2d/out_PWN2d_guitar2/'
d = loaddata(dir+'PWN2d_hdf5_plt_cnt_0117','dens',xCoord=x,yCoord=y,lref=lref,time=t,sample=3)

dd = transpose(d[0:199,200:899])
xx = y[200:899]
yy = x[0:199]

dsz = size(dd,/dimension)
d = fltarr(dsz[0],dsz[1]*2)
d[*,0:dsz[1]-1] = reverse(dd,2)
d[*,dsz[1]:dsz[1]*2-1] = dd

x = xx
y = [-reverse(yy),yy]

;window,xs=700,ys=400
;tvcoord,alog(d),x,y,/scale

loadct,0,/sil
if not keyword_set(ps) then window,0,xs=700+250,ys=400+150 $
 else mkeps,'pp_showPWN2d_guitar',xs=20.,ys=20.*55./95.

;contour,alog(transpose(d)),y,x,nlevel=256,/fill,/iso,xtitle='x [cm]',ytitle='y [cm]'

if keyword_set(ps) then $
   tvcoord,alog(d),x,y,/scale,/on,xtitle='x [cm]',ytitle='y [cm]',charsize=1.5,position=[0.18,0.22],imgsize=0.65,xtickinterval=1.e18  $
else $
   tvcoord,alog(d),x,y,/scale,/on,xtitle='x [cm]',ytitle='y [cm]' 

plots, 0.,0.,psym=7,color=0,/data,symsize=1.5

;t2 = t - 1.10167e+11 ; time for second bubble exploding (9.35e10 + 1.e18/6.e7)
t2 = t - (6.45e17/1.5d8) ; time for second bubble exploding (9.35e10 + 1.e18/6.e7)

Esp = 1d33
damb1 = 2.76d-26
vpwn = 1.5d8
rr1 = (125./154./!pi)^0.2*(Esp/damb1)^0.2*t^0.6 
plots,ring(t*vpwn,0.,rr1),/data
damb2 = 1.23e-25
rr2 = (125./154./!pi)^0.2*(Esp/damb2)^0.2*t2^0.6 
plots,ring(t2*vpwn,0.,rr2),/data

loadct,0,/sil
color_bar,lim=[min(d),max(d)],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),titlegap=0.11,charsize=1.5

if keyword_set(ps) then epsfree
stop
end
