pro gaensler2_2

img = read_png('gaensler.png')

imgsz = size(img,/dimension)
window,xs=imgsz[1],ys=imgsz[2],/pixmap
tv,img,/true
img = tvrd()

img = transpose(img)
imgsz = size(img,/dimension)

x0=-1.19e17 & dx=6.26e17
y0=-2.58e17 & dy=5.85e17
x=[x0,x0+dx]
y=[y0,y0+dy]
imgxs=515. & imgys=imgxs*(y[1]-y[0])/(x[1]-x[0])

pltx0=150. & plty0=70.
winxs=pltx0+imgxs*2+30 & winys=plty0+imgys+30
loadct,0,/sil

xscm=30 & yscm=xscm*winys/winxs
mkeps,'gaensler2_2',xs=xscm,ys=yscm

deg = -29.
img_rot = rot(img,deg) 

;mvright = 120
mvright = 150
mvdown = 33

img2 = bytarr(pltx0+imgxs,imgsz[1]) & img2[*]=255
;img2 = bytarr(imgsz[0]+mvright,imgsz[1]) & img2[*]=255
;img2 = bytarr(imgsz[0],imgsz[1]) & img2[*]=255
imgsz2 = size(img2,/dimension)
img2[mvright:imgsz2[0]-1,*]=img_rot[0:imgsz2[0]-mvright-1,*]

img2[*,0:imgsz2[1]-mvdown-1]=img2[*,mvdown:imgsz2[1]-1]
img2[*,imgsz2[1]-mvdown:imgsz2[1]-1]=255

img2[*,0:plty0-1]=255
img2[*,plty0+imgys:imgsz2[1]-1]=255
img2[0:pltx0-1,*]=255
;img2[pltx0+imgxs:imgsz2[0]-1,*]=255

;img3 = bytarr(imgsz

;tv,img2
tv,img2,float(imgxs)/winxs*xscm,0.,/centimeter,xsize=imgsz2[0]/winxs*xscm,ysize=imgsz[1]/winys*yscm,true=true
tv,img2,0,0,/centimeter,xsize=imgsz2[0]/winxs*xscm,ysize=imgsz[1]/winys*yscm,true=true

lthk=5
;plot,x,y,/iso,/xst,/yst,position=[pltx0,plty0,pltx0+imgxs,plty0+imgys],/dev,/noerase,/nodata,color=0,xtitle='x [cm]',ytitle='y [cm]'
plot,x,y,/iso,/xst,/yst,position=posnorm([pltx0,plty0,pltx0+imgxs,plty0+imgys],nx=winxs,ny=winys),/norm,/noerase,/nodata,color=0 $
    ,xtitle='x [cm]',ytitle='y [cm]',xtickinterval=2.e17
arrow, -3.e16,0.,-1.e17,0.,/data,thick=10,hsize=400,color=0

asymang=38.*!dtor
asymy = 8.e16
tvlct,r,g,b,/get
oplot,!x.crange,[asymang*!x.crange[0]+asymy,asymang*!x.crange[1]+asymy],thick=lthk,color=fsc_color('cyan'),line=2
tvlct,r,g,b


plot,x,y,/iso,/xst,/yst,position=posnorm([pltx0+imgxs,plty0,pltx0+2*imgxs,plty0+imgys],nx=winxs,ny=winys),/norm,/noerase,/nodata,color=0 $
    ,xtitle='x [cm]',ytickformat='(a1)',xtickinterval=2.e17

;tvlct,r,g,b,/get
;oplot,!x.crange,[0.,0.],line=2,color=fsc_color('cyan')
;oplot,[0.,0.],!y.crange,line=2,color=fsc_color('cyan')
;tvlct,r,g,b


restore,file='multibow_obs.sav'
loadct,39,/sil

xoff = 8.e15
oplot,ds[600].x-ds[600].x[0]-xoff,-ds[600].w1,color=254,thick=lthk

oplot,ds[100].x-ds[100].x[0]-xoff,ds[100].w1,color=30,thick=lthk
oplot,ds[50].x-ds[50].x[0]-xoff,ds[50].w1,color=50,thick=lthk
oplot,ds[20].x-ds[20].x[0]-xoff,ds[20].w1,color=150,thick=lthk
oplot,ds[10].x-ds[10].x[0]-xoff,ds[10].w1,color=210,thick=lthk
;oplot,ds[5].x-ds[5].x[0]-xoff,ds[5].w1,color=220,thick=lthk

legend,[textoidl('\rho_{0,1}'),textoidl('\rho_{0,2}'),textoidl('\rho_{0,3}'),textoidl('\rho_{0,4}'),textoidl('\rho_{0,5}')],/left,/bottom $
      ,textcolor=[210,150,50,30,254],line=0,color=[210,150,50,30,254],/clear,thick=5 ,charsize=1.2

arrow, -3.e16,0.,-1.e17,0.,/data,thick=10,hsize=400,color=0

epsfree

print,'d1= ',ds[10].d0
print,'d2= ',ds[20].d0
print,'d3= ',ds[50].d0
print,'d4= ',ds[100].d0
print,'d5= ',ds[600].d0

stop
end
