pro bubmodel,time=time, xrange=xrange, yrange=yrange

if not keyword_set(xrange) then xrange = [-1.e18,3.e18]
if not keyword_set(yrange) then yrange = [-1.e18,1.e18]
if not keyword_set(time) then time = 3.45e11

; analytic bubble
Esp = 1.d36 ;4.3d33 ;6.66d35
vw  = 1.d8 ;6.d6 ;6.d7
vp  = 1.d10

eta = 125./154./!pi
d1 = 1.67d-23
;d1 = 2.76d-25
t1 = time
r1 = eta^0.2*(Esp/d1)^0.2*t1^0.6 
l1 = time*vw
bub1r = ring(l1,0.,r1)
print,'location bubble x-center: ',l1

tbr = sqrt(eta*Esp/d1)*vw^(-2.5)
print,'break time: ',tbr, ' s'

window,0
plot,xrange,yrange,/nodata,/xst,/yst,/iso,xtitle='x [cm]',ytitle='y [cm]'
plots,[0.,0.],psym=7,symsize=2
oplot,[l1,l1],!y.crange,line=1
oplot,!x.crange,[0.,0.],line=1
oplot,bub1r[0,*],bub1r[1,*]

legend,strtrim(time,2)+' s',/top,/right,box=0

stop
end
