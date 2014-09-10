pro chkpower,n,xc0=xc0,yc0=yc0,zc0=zc0

if not keyword_set(xc0) then xc0=-2.e12
if not keyword_set(yc0) then yc0=0.
if not keyword_set(zc0) then zc0=0.

; parameters
gam=4./3.
mach=30.

mkdata=1

if (mkdata eq 1) then begin
d=dload(n,var='dens',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=0,time)
v=dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=0,time)

dxy = reform(d[*,*,258])
vxy = reform(v[*,*,258])

save,file='chkpower.sav',dxy,vxy,x,y

endif else restore, file='chkpower.sav'

sz = size(dxy,/dimension)
window,0,xs=sz[0],ys=sz[1]*2
tvscl,alog(dxy),0
tvscl,vxy,1

dx = x[2] - x[1]
dy = y[2] - y[1]

L = total(dxy*vxy^3.*dx*dy) * 2.*(0.5 + 1./(gam*(gam-1.)*mach))

print,L

stop
end

