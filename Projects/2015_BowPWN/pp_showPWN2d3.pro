pro pp_showPWN2d3,ps=ps

;dir='/d/d7/yoon/out_FLASH3.3_mhd/out_PWN2d/out_PWN2d_600kms_wall_2/'
dir='/d/d7/yoon/out_FLASH3.3_mhd/out_PWN2d/out_PWN2d_600kms_wall/'
;d = loaddata(dir+'PWN2d_hdf5_plt_cnt_0692','dens',xCoord=x,yCoord=y,lref=lref,time=t,sample=4)
d = loaddata(dir+'PWN2d_hdf5_plt_cnt_0550','dens',xCoord=x,yCoord=y,lref=lref,time=t,sample=4)

d = loaddata('/d/d7/yoon/out_FLASH3.3_mhd/out_PWN3d/tst_PWN3d/PWN3d_hdf5_plt_cnt_0071','dens' $
             ,xra=[-1.e18,1.e18],yra=[-5.e17,5.e17],zra=[-1.e18,1.e18],xCoord=x,yCoord=y,lref=lref,time=t)

dsz = size(d,/dimension)
d2 = fltarr(dsz[0]*2,dsz[1])

d2[0:dsz[0]-1,*] = reverse(d,1)
d2[dsz[0]:dsz[0]*2-1,*] = d

x2 = [-reverse(x),x]

loadct,0,/sil
if not keyword_set(ps) then window,0,xs=1300,ys=650 $
; else mkeps,'pp_shoPWN2d_wall2_0691',xs=20.,ys=20.*60./130
; else mkeps,'pp_shoPWN2d_wall2_0550',xs=20.,ys=20.*60./100
 else mkeps,'pp_showPWN2d_45deg',xs=20.,ys=20.*60./100

;contour,alog(transpose(d)),y,x,nlevel=256,/fill,/iso,xtitle='x [cm]',ytitle='y [cm]'

if keyword_set(ps) then $
   tvcoord,alog(transpose(d2)),y,x2,/scale,/on,xtitle='x [cm]',ytitle='y [cm]',charsize=1.5,position=[0.18,0.22],imgsize=0.65,xtickinterval=1.e19  $
else $
   tvcoord,alog(transpose(d2)),y,x2,/scale,/on,xtitle='x [cm]',ytitle='y [cm]' 

;plots,ring(0,0,1.24e16),/data,color=255
plots, 0.,0.,psym=7,color=255,/data,symsize=1.5

if keyword_set(anal) then begin
;t2 = t - 1.10167e+11 ; time for second bubble exploding (9.35e10 + 1.e18/6.e7)
t2 = t - 7.31667e+10 ; time for second bubble exploding (9.35e10 + 1.e18/6.e7)

Esp = 6.66d35
dam = 1.67d-24
rr1 = (125./154./!pi)^0.2*(Esp/dam)^0.2*t^0.6 
plots,ring(t*6.e7,0.,rr1),/data
dam2 = 1.67d-24
rr2 = (125./154./!pi)^0.2*(Esp/dam2)^0.2*t2^0.6 
plots,ring(t2*6.e7,0.,rr2),/data

endif

loadct,0,/sil
color_bar,lim=[min(d),max(d)],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),titlegap=0.12,charsize=1.5

if keyword_set(ps) then epsfree
stop
end
