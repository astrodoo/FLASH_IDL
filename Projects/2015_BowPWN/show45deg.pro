pro show45deg

dir='/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN3d/out_PWN3d_600kms_45deg_highRes/'
file='PWN3d_hdf5_plt_cnt_0071'

d=reform(loaddata(dir+file,'dens',sample=0,xcoord=x,ycoord=y,xra=[-1.e18,2.e18],yra=[-1.e18,1.e18],zra=[0.,0.]))

d1=smooth(d,10)
d2=smooth(d,20)
dsz=size(d,/dimension)

;drawing
pltx0=250. & plty0=150.
winxs=pltx0+dsz[0]+300 & winys=plty0+dsz[1]+20
winxs_cm=25. & winys_cm=winxs_cm*winys/winxs

mkeps,'show45deg',xs=winxs_cm,ys=winys_cm
;!p.background=255 & !p.color=0
;window,xs=dsz[0]+200,ys=dsz[1]+200

tvcoord,alog10(d1),x,y,/scale,/axes,/black,xtickinterval=1.e18,pos=[pltx0/winxs,plty0/winys],imgsize=float(dsz[0])/winxs $
       ,xtitle='x [cm]', ytitle='y [cm]'
;contour,alog10(d2),x,y,level=[-24.5,-24.3,-24.2,-24.1,-24,-23.9,-23.8,-23.7],/overplot
contour,alog10(d2),x,y,level=[-24.5,-24.45,-24.4,-24.1,-23.9,-23.76,-23.75,-23.7],/overplot

plots,[0.,0.],psym=7,symsize=2,color=255,/data
arrow,5.e16,-8.e17,5.e16,-3.5e17,thick=3.,/data


color_bar,lim=[min(d),max(d)],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),bargap=0.01


epsfree

stop
end
