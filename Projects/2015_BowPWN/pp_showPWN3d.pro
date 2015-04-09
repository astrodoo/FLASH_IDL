pro pp_showPWN3d,ps=ps

d = loaddata('/d/d7/yoon/out_FLASH3.3_mhd/out_PWN3d/tst_PWN3d/PWN3d_hdf5_plt_cnt_0071','dens' $
             ,xra=[-1.e18,1.e18],yra=[-5.e17,5.e17],zra=[-1.e18,1.e18],xCoord=x,yCoord=y,lref=lref,time=t)

sd = size(d,/dimension)

dxy = reform(d[*,*,sd[2]/2])

smd = smooth(dxy,8)

dxy=smd
loadct,0,/sil
if not keyword_set(ps) then window,0,xs=1000,ys=500 $
 else mkeps,'pp_showPWN3d_45deg',xs=20.,ys=20.*5./10.


;contour,alog(transpose(d)),y,x,nlevel=256,/fill,/iso,xtitle='x [cm]',ytitle='y [cm]'

levs = [1.e-25,2.e-25,3.e-25,4.e-25,5.e-25,1.e-24,2.e-24]
if keyword_set(ps) then begin
   tvcoord,alog(dxy),x,y,/scale,/on,xtitle='x [cm]',ytitle='y [cm]',charsize=1.5,position=[0.18,0.22],imgsize=0.65
   contour,dxy,x,y,levels=levs,/overplot,color=0
endif else begin
   tvcoord,alog(dxy),x,y,/scale,/on,xtitle='x [cm]',ytitle='y [cm]' 
   contour,dxy,x,y,levels=levs,/overplot,color=0
endelse

plots, 0.,0.,psym=7,color=255,/data,symsize=1.5

arrow,5.e16,-4.e17,5.e16,-2.5e17,/data,color=0,thick=2.

loadct,0,/sil
color_bar,lim=[min(d),max(d)],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),titlegap=0.12,charsize=1.5

if keyword_set(ps) then epsfree
stop
end
