pro bubmodel,time=time, xrange=xrange, yrange=yrange

if not keyword_set(xrange) then xrange = [-1.e18,3.e18]
if not keyword_set(yrange) then yrange = [-1.e18,1.e18]
if not keyword_set(time) then time = 3.45e11

; analytic bubble
Esp = 1.d36 ;4.3d33 ;6.66d35
vw  = 1.d8 ;6.d6 ;6.d7
vp  = 1.d10

;Esp=1.d33
;vw=1.5d8
;vp=1.d10
;xrange=[-1.e18,3.e18]
;yrange=[0.,1.e18]

d1 = 1.67d-23
;d1 = 2.76d-25
t1 = time
r1 = (125./154./!pi)^0.2*(Esp/d1)^0.2*t1^0.6 
l1 = time*vw
bub1r = ring(l1,0.,r1)
print,'location 1st bubble: ',l1

;ldisc1 = 5.61e18
ldisc1 = 1e19
;d2 = 1.67d-24
d2 = 1.67d-25
t2 = time - (ldisc1/vw)
l2 = t2*vw
t2 = t2>0.
r2 = (125./154./!pi)^0.2*(Esp/d2)^0.2*t2^0.6 
bub2r = ring(l2,0.,r2)
print,'location 2nd bubble: ',l2

hdisc1 = 1.17e18
th1 = time - (hdisc1/vw)
lh1 = th1*vw



;nth = 300
;th0 = 0.1 * !dtor
;th1 = 180. * !dtor
;th  = findgen(nth) / nth * (th1-th0) + th0
;r0 = sqrt(Esp / (4.*!pi*d2*vw*vw*vp))
;rth = r0 / sin(th) * sqrt(3.*(1.-th/tan(th)))
;xx = - rth * cos(th)
;yy = rth * sin(th)

;window,0,xs=1000,ys=1000.*(yrange[1]-yrange[0])/(xrange[1]-xrange[0])
window,0
plot,xrange,yrange,/nodata,/xst,/yst,/iso,xtitle='x [cm]',ytitle='y [cm]'
plots,[0.,0.],psym=7,symsize=2
oplot,[l1,l1],!y.crange,line=2
oplot,bub1r[0,*],bub1r[1,*]
;oplot,[l2,l2],!y.crange,line=2
;oplot,bub2r[0,*],bub2r[1,*]

;oplot,[lh1,lh1],!y.crange,line=2

;oplot,xx,yy

legend,strtrim(time,2)+' s',/top,/right,box=0

stop
end
