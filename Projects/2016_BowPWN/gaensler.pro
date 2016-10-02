pro gaensler

img = read_png('gaensler.png')

imgsz = size(img,/dimension)
window,xs=imgsz[1],ys=imgsz[2],/pixmap
tv,img,/true
img = tvrd()

img = transpose(img)
imgsz = size(img,/dimension)

;-------------------------------------
; rotate
deg = -29.
img_rot = rot(img,deg) 

pltx0=100. & plty0=100
winxs=pltx0+imgsz[0]+50 & winys=plty0+imgsz[1]+50
loadct,0,/sil
!p.background=255
window,xs=winxs,ys=winys
tv,img_rot,pltx0,plty0

ytrack = 415

tvlct,r,g,b,/get
plots,[0.,winxs-1],[ytrack,ytrack],/dev,color=fsc_color("blue")

; y=ax+b
a1=38.*!dtor
b1=355.
plots,[0.,winxs-1],[b1,a1*(winxs-1)+b1],/dev,color=fsc_color("red")

a2=-12.*!dtor
b2=430.
plots,[0.,winxs-1],[b2,a2*(winxs-1)+b2],/dev,color=fsc_color("green")

tvlct,r,g,b

snapshot,'gaensler_ang'
stop
end
