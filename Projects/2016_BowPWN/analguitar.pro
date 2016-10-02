pro analguitar,mkdata=mkdata,d13=d13,d23=d23, distance=distance


if not keyword_set(d13) then d13=5.       ; density ratio between region 1 & 3 (neck & bubble)
if not keyword_set(d23) then d23=0.5      ; density ratio betseen region 2 & 3 (middle bubble-like & bubble)
if not keyword_set(distance) then distance=1.8   ; distance [kpc]

eta = 125.d0/154.d0/!dpi

vsag = 182.d-3 /60./60.*!dtor/60./60./24./365.  ; projected angular velocity (radian/s)
l_bag = 65. /60./60.*!dtor  ; distance to bubble center (radian)
r_bag = 16. /60./60.*!dtor  ; radius of bubble center (radian) 
Edot = 1.e33    ; pulsar spin-down loss energy (ergs s^-1)

dst=distance*1.e3*!unit.pc   ; distance (cm)

d0_b = eta*Edot*l_bag^3./(r_bag^5.*vsag^3.*dst^5.)
d1 = d13*d0_b
d2 = d23*d0_b

time = l_bag / vsag

vs = vsag * dst
l_b = l_bag * dst
r_b = r_bag * dst

if keyword_set(mkdata) then begin
   analBow, time=time, Ls=Edot, d0=d1, vs=vs, P0=1.7d-12,x,w,dwdx,out='analBow1_'+string(distance,format='(f3.1)')+'.sav'
   x1 = x & w1 = w
   analBow, time=time, Ls=Edot, d0=d2, vs=vs, P0=1.7d-12,x,w,dwdx,out='analBow2_'+string(distance,format='(f3.1)')+'.sav',/bubp
   x2 = x & w2 = w
endif else begin
   restore,file='analBow1_'+string(distance,format='(f3.1)')+'.sav'
   x1=x & w1=w
   restore,file='analBow2_'+string(distance,format='(f3.1)')+'.sav'
   x2=x & w2=w
endelse

disc1 = 30. /60./60.*!dtor *dst
disc2 = 53. /60./60.*!dtor *dst

ind1 = where(x1 le disc1) & ind1_r = ind1[n_elements(ind1)-1]
ind2 = where((x2 ge disc1) and (x2 le disc2)) & ind2_r = ind2[n_elements(ind2)-1]
ind2_l = ind2[0]
;window,0
mkeps,'analguitar_'+string(distance,format='(f3.1)'),xs=20.,ys=20.*5./8.

plot,findgen(10),/nodata,/xst,/yst,xra=[0.,1.3*l_b]/dst/!dtor*60.*60., $
    yra=1.5*[-r_b,r_b]/dst/!dtor*60.*60.,/iso,xtitle='x/distance [arcsec]',ytitle='y/distance [arcsec]'
plots,ring(l_b/dst/!dtor*60.*60.,0.,r_b/dst/!dtor*60.*60.),/data, color=0
oplot,x1[ind1]/dst/!dtor*60.*60.,w1[ind1]/dst/!dtor*60.*60.
oplot,x1[ind1]/dst/!dtor*60.*60.,-w1[ind1]/dst/!dtor*60.*60.
oplot,x2[ind2]/dst/!dtor*60.*60.,w2[ind2]/dst/!dtor*60.*60.
oplot,x2[ind2]/dst/!dtor*60.*60.,-w2[ind2]/dst/!dtor*60.*60.

oplot,[disc1,disc1]/dst/!dtor*60.*60.,!y.crange,line=1
oplot,[disc2,disc2]/dst/!dtor*60.*60.,!y.crange,line=1
oplot,!x.crange,[0.,0.],line=1

epsfree

print,'density (1,2,3):',d1,d2,d0_b
print,'density ratio (d1/d0_b) (d2/d0_b) :',d1/d0_b,d2/d0_b
print,'width ratio (w1/r_b) (w2/r_b_) :',w1[ind1_r]/r_b,w2[ind2_l]/r_b


stop
end

