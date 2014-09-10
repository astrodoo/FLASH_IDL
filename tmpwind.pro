pro tmpwind

R0 = 1.4e12
r = findgen(1000)/1000.*5.*R0

bt = 0.8

vtm = 2.e8
mdot = 6.5e20

distf = 1.-R0/r*0.7
;distf = 1.-R0/r*0.999
distf[where(distf lt 0)]=0.
v = vtm * distf^bt
rho = mdot / (4.*!pi*r*r*v)

;v5 = vtm*(1.-1./5.*0.999)^0.8
;rho5 = mdot / (4.*!pi*25.*R0*R0*v5)
v1 = vtm*(1.-0.7)^0.8
rho1 = mdot / (4.*!pi*R0*R0*v1)
;rho[where(rho ge 5.*rho5)]=5.*rho5
;rho[where(rho ge 500.*rho5)]=500.*rho5
rho[where(rho ge rho1)]=rho1
v[where(v le v1)]=0.
window,0,xs=600,ys=700
;mkeps,'tmpwind.eps',xs=20.,ys=25.
multiplot,[1,2]
plot,r/R0,v,/xst,ytitle='wind velocity',yst=2
multiplot
plot,r/R0,rho,/xst,ytitle='density',xtitle='distance [r/R_*]',yst=2
multiplot,/reset
;epsfree
stop
end
