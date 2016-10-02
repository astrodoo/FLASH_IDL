pro guitar_chatterjee

img = read_png('guitar_chatterjee_04.png')

imgsz = size(img,/dimension)
window,xs=imgsz[1],ys=imgsz[2],/pixmap
tv,img,/true 
img = tvrd()

imgsz = size(img,/dimension)

;-------------------------------------
; rotate
deg = -39.
img_rot = rot(img,deg,missing=255) 

restore,file= '/home/jianiye/Work/Data/2015_BowPWN/out_PWN2d_guitar_60inc/phot3d_rot_60.sav'
simul = prj_rot
mins = min(simul) & maxs = max(simul)
szs = size(simul,/dimension)

xxx=2.28e18
coordarray,x,y,xout=x2d,yout=y2d
simul[where(x2d ge xxx)] = mins

pltx0=150. & plty0=80.
pltxs=670. & pltys=pltxs*(y[szs(1)-1]-y[0])/(x[szs(0)-1]-x[0])
winxs=pltx0+pltxs+30. & winys=plty0+pltys*2.+30.
winxscm = 20. & winyscm = winxscm * winys/winxs

;window,xs=winxs,ys=winys
mkeps,'guitar_chatterjee',xs=winxscm,ys=winyscm

xmv = 30 & ymv = -110    ; right & down movement of image (emperically obtained by eye inspection in window (no ps))
img2 = bytarr(pltxs,pltys)
img2 = img_rot[pltx0-xmv:pltx0-xmv+pltxs-1,plty0-ymv:plty0-ymv+pltys-1]

tv,img2,(pltx0/winxs)*winxscm,((plty0+pltys)/winys)*winyscm,/centimeter,xs=pltxs/winxs*winxscm
levs=[1.e-3,1.e-2,2.e-2,3.e-2,4.e-2,5.e-2,6.e-2,7.e-2,8.e-2]
contour,simul,x,y,levels=levs,/iso,/noerase,/xst,/yst $
    ,position=posnorm([pltx0,plty0+pltys,pltx0+pltxs,plty0+2*pltys],nx=winxs,ny=winys),/norm,/nodata $
    ,xtickformat='(a1)',ytickinterval=4.e17,xtickinterval=1.e18

;tv,img_rot,xmv,ymv
tv,img2,(pltx0/winxs)*winxscm,(plty0/winys)*winyscm,/centimeter,xs=pltxs/winxs*winxscm
;tv,img_rot,4.2,-0.15,/centimeter,xsize=18.
;plots,[0.,imgsz[0]],[imgsz[1],imgsz[1]]/2.+37,/dev

contour,simul,x,y,levels=levs,/iso,/noerase,/xst,/yst,xtitle='x [cm]' $
    ,position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm $
    ,ytickinterval=4.e17,xtickinterval=1.e18, thick=3.
;oplot,[xxx,xxx],!y.crange,line=2

xyouts,0.05,0.47,/norm,'y [cm]',orientation=90

epsfree

stop
end
