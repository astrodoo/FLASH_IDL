pro tau1,sigma=sigma,noerase=noerase,color=color

; find the location where the optical depth (x-ray) would be 1.
; /d/d3/yoon/outputs/out_mHD_Binary_beta 
; saved data: denxy,x,y,time
;restore,'showjet_450_den.sav'

restore,'showjet_320_den.sav'
denxy=denxz
y=z

sd = size(denxy)

bh  = -3.e12 * 2./3.
st  = 3.e12 / 3.
str = 1.4e12

if not keyword_set(noerase) then begin
loadct,0,/sil
window,0,xs=sd[1],ys=sd[2]
plot,x,y,/xst,/yst,/nodata,xr=[x[0],x[sd[1]-1]],yr=[y[0],y[sd[2]-1]],/iso,position=[0,0,sd[1],sd[2]],/dev
tvscl,alog10(denxy)
plots,bh,0,/data,psym=1,thick=1,color=0
plots,st,0,/data,psym=1,thick=1,color=0
plots,ring(st,0,str),color=0,/data,thick=1
endif

;sigma = !unit.tcs  ; Thompson Cross-section
if not keyword_set(sigma) then sigma = !unit.tcs

; picking up the line-coordinate
ds = x[2]-x[1]
smax = (x[sd[1]-1]-x[0])/2.
ns = fix(smax/ds)

s = findgen(ns)/float(ns)*smax
ds = s[2]-s[1]

nth=100
thmax=130. ;!degress

th = !dtor*(findgen(nth)/float(nth)*thmax*2. - thmax)
;th = !pi-th

if not keyword_set(color) then color=0

for j=0,nth-1 do begin 
    xx = bh + s*cos(th[j])
    yy = s*sin(th[j])
    indxy = convert_coord(xx,yy,/data,/to_device)
;    plots,xx,yy,psym=3,thick=1,color=0

    den_int = reform(interpolate(denxy,indxy[0,*],indxy[1,*]))

    tau=0.
    const = sigma / !unit.mh
    for i=0,ns-1 do begin
        tau += const*ds*den_int[i] 
        if (tau ge 1.) then begin
           xt = xx[i] & yt = yy[i]
           goto, jump
        endif else if (sqrt((xx[i]-st)^2.+yy[i]^2.) le str) then begin
           xt = xx[i] & yt = yy[i]
           goto, jump
        endif   
    endfor
    xt = !values.f_nan & yt = !values.f_nan

Jump:
    plots,xt,yt,/data,psym=dsym(8,/fill),thick=1,color=color
endfor  ; j

;stop
end

pro extau1
;excute tau1

; w/ tompson cross section
tau1
loadct,39,/sil
tau1,sigma=1.e-23,/noerase,color=50
tau1,sigma=5.e-23,/noerase,color=100
tau1,sigma=1.e-22,/noerase,color=150
tau1,sigma=5.e-22,/noerase,color=190
tau1,sigma=1.e-21,/noerase,color=210

legend,textoidl('\sigma = ')+['tompsons','1.e-23','5.e-23','1.e-22','5.e-22','1.e-21'],/right,/bottom $
      ,psym=dsym(8,/fill),color=[0,60,100,150,190,210],textcolor=[0,60,100,150,190,210],box=0
stop
end


pro cross_phio;,nu=nu

; Osterbrook (Tab. 8.7)
; Hydrogen

nuT  = 1.097e5
sigT = 6.3e-18
bet  = 1.34
s    = 2.99

nuMax = 120.e3 * !unit.ev / !unit.h / !unit.c
nuMin = 120. * !unit.ev / !unit.h   / !unit.c
nuMin = nuT
nnu   = 10000
nu    = findgen(nnu)/float(nnu)*(nuMax-nuMin) + nuMin

ee    = nu*!unit.c *!unit.h/!unit.ev

crs = sigT*(bet*(nu/nuT)^(-s) + (1.-bet)*(nu/nuT)^(-s-1.))

;window,0,xs=800.,ys=700.
mkeps,'tau1_crs.eps',xs=20.,ys=20.*700./800.
plot,ee/1.e3,crs,/ylog,/xst,xtitle=textoidl('E_{ph} [keV]'), ytitle=textoidl('cross section [cm^{2}]')
plot,ee/1.e3,crs,/ylog,/xst,xr=[0.,12.],pos=[0.5,0.5,0.92,0.92],/norm,/noerase

epsfree
stop
end

