pro pp_showPWN2d2,ps=ps

dir='/d/d7/yoon/out_FLASH3.3_mhd/out_PWN2d/out_PWN2d_600kms/l10_2/'
d = loaddata(dir+'PWN2d_hdf5_plt_cnt_0250','dens',xCoord=x,yCoord=y,lref=lref,time=t,sample=3)

dsz = size(d,/dimension)
d2 = fltarr(dsz[0]*2,dsz[1])

d2[0:dsz[0]-1,*] = reverse(d,1)
d2[dsz[0]:dsz[0]*2-1,*] = d

x2 = [-reverse(x),x]

loadct,0,/sil
if not keyword_set(ps) then window,0,xs=1000,ys=650 $
 else mkeps,'pp_shoPWN2d_600kms_0250',xs=20.,ys=20.*43./80.

;contour,alog(transpose(d)),y,x,nlevel=256,/fill,/iso,xtitle='x [cm]',ytitle='y [cm]'

if keyword_set(ps) then $
   tvcoord,alog(transpose(d2)),y,x2,/scale,/on,xtitle='x [cm]',ytitle='y [cm]',charsize=1.5,position=[0.17,0.14],imgsize=0.65 $
else $
   tvcoord,alog(transpose(d2)),y,x2,/scale,/on,xtitle='x [cm]',ytitle='y [cm]' 

;plots,ring(0,0,1.24e16),/data,color=255
plots, 0.,0.,psym=7,color=255,/data,symsize=1.5

t2 = t - 6.86333e+10 ; time for second bubble exploding (6.03e10 + 5.e17/6.e7)

Esp = 6.66d35
dam = 1.67d-24
fac = 0.9
rr1 = (25./14./!pi)^0.2*(Esp/dam)^0.2*t^0.6 *fac
plots,ring(t*6.e7,0.,rr1),/data
dam2 = 1.67d-25
rr2 = (25./14./!pi)^0.2*(Esp/dam2)^0.2*t2^0.6 *fac
plots,ring(t2*6.e7,0.,rr2),/data

loadct,0,/sil
color_bar,lim=[min(d),max(d)],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),titlegap=0.12,charsize=1.5

if keyword_set(ps) then epsfree
stop
end
