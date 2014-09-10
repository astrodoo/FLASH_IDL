pro overorbit,m1=m1,m2=m2,peri=peri,ecc=ecc,nt=nt,color=color,thick=thick

if not keyword_set(m1) then m1 = 4.d34
if not keyword_set(m2) then m2 = 2.d34
mtot = m1+m2

if not keyword_set(peri) then peri = 3.d12
if not keyword_set(ecc) then ecc = 0.
if not keyword_set(nt) then nt = 100
if not keyword_set(color) then color=0
if not keyword_set(thick) then thick=1

l = sqrt(peri*(1.d0-ecc*ecc)*!unit.g*mtot)
T = 2.d0*!pi*sqrt(peri*peri*peri/!unit.g/mtot)   ; period

dt = T / nt ; *0.9
th0_e = !pi
th_e = th0_e
cm1x_e = dblarr(nt) & cm1y_e = dblarr(nt)
cm2x_e = dblarr(nt) & cm2y_e = dblarr(nt)

for i=0,nt-1 do begin
   rr = peri*(1.d0-ecc*ecc) / (1.d0 + ecc*cos(th_e))
   th_e = th_e + l/rr/rr * dt

; center of mass
   cm1x_e[i] = -rr * m2/mtot * cos(th_e)
   cm1y_e[i] = -rr * m2/mtot * sin(th_e)
   cm2x_e[i] = rr * m1/mtot * cos(th_e)
   cm2y_e[i] = rr * m1/mtot * sin(th_e)
endfor

plots,cm1x_e,cm1y_e,line=2,/data,color=color,thick=thick
plots,cm2x_e,cm2y_e,line=2,/data,color=color,thick=thick

end
