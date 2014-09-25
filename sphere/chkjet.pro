pro chkjet,n,xc0=xc0,yc0=yc0,zc0=zc0, jetgamma=jetgamma, mkdata=mkdata

; check jet kinetic power

if not keyword_set(xc0) then xc0=-2.e12
if not keyword_set(yc0) then yc0=0.
if not keyword_set(zc0) then zc0=0.

; parameters
if not keyword_set(jetgamma) then jetgamma=4./3.
gam=jetgamma

if keyword_set(mkdata) then begin
d=dload(n,var='dens',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=0,time)
p=dload(n,var='pres',xc=xc0,yc=yc0,zc=zc0,sample=0)
v=dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,sample=0)

zind_tmp = where(z ge zc0) & zind = zind_tmp[0]
print,'z index: ', zind+2

dxy = reform(d[*,*,zind+2])
pxy = reform(p[*,*,zind+2])
vxy = reform(v[*,*,zind+2])
save,file='chkjet.sav',dxy,vxy,pxy,x,y,z
endif else restore, file='chkjet.sav'

dj = max(dxy)
pj = max(pxy)
vj = max(vxy) 

yind_tmp = where(y ge yc0) & yind = yind_tmp[0]
dx = reform(dxy[*,yind])
px = reform(pxy[*,yind])
vx = reform(vxy[*,yind])
window,0,xs=800,ys=800
multiplot,[1,3]
plot,x,dx,/xst,yst=2,/ylog,ytitl='density'
oplot,!x.crange,[dj,dj],line=1
multiplot
plot,x,px,/xst,yst=2,/ylog,ytitle='pressure'
oplot,!x.crange,[pj,pj],line=1
multiplot
plot,x,vx,/xst,yst=2,ytitle='velocity',xtitle='x [cm]'
oplot,!x.crange,[vj,vj],line=1
multiplot,/reset,/init

print, "jet gamma (input): ",jetgamma
print, "jet velocity: ", vj

mach = vj*sqrt(1./jetgamma* dj/pj)
print, "jet Mach: ",mach

sz = size(dxy,/dimension)
window,1,xs=sz[0]*2,ys=sz[1]*2
tvscl,alog(dxy),0
tvscl,alog(pxy),1
tvscl,vxy,2

dx = x[2] - x[1]
dy = y[2] - y[1]

L = total(dxy*vxy^3.*dx*dy) * 2.*(0.5 + 1./(gam*(gam-1.)*mach))

print,"jet Power: ", L

stop
end

