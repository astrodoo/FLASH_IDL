pro galprofile
; analytic solution for galatic density profile with rotation.

; blackhole mass
Mbh = 1.d9*!unit.msun
; schwarzschild radius
rbh = 2.*!unit.g*Mbh/!unit.c/!unit.c
rgal = 30.d3 * !unit.pc         ; size of the galaxy ~ 30 kpc
; orbital angular velocity
omg = 1.d5 / (1.d3*!unit.pc)    ; rotational velocity  1 km/s at 1 kpc
; initial density at (r=rbh, theta=0)
d0 = 1.67d-23
p0 = 1.d-2

; adiabatic index
gam = 5./3.


nr = 100
nth = 100
r1d = exp(dindgen(nr))/exp(double(nr)-1.)*(rgal-rbh) + rbh
th1d = dindgen(nth)/double(nth)*180.

coordarray,r1d, th1d, xout=r, yout=th

dnr = d0 * ( 1. + (gam-1.)/gam * d0/p0 *( !unit.g*Mbh*(1./r1d - 1./rbh)  ) )^(1./(gam-1.)) 
;d = d0 * ( 1. + (gam-1.)/gam * d0/p0 *( abs(!unit.g*Mbh*(1./r - 1./rbh) + 0.5*r*r*omg*omg*sin(th)*sin(th) ) ) )^(1./(gam-1.)) 


tm1 = !unit.g*Mbh/r1d
tm2 = !unit.g*Mbh/rbh
tm3 = 0.5*r1d*r1d*omg*omg










stop
end
