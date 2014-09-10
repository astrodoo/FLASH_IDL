pro radio,fname,ps=ps

if not keyword_set(fname) then fname='jet_xrbs_hdf5_plt_cnt_0730'  ; lowres

fname_r = 'radio_lowres_plt_0730.eps'

cutx = 6.e19
cuty = 9.874e20
cutz = 9.874e20

;xce = 2.e19           ; 300 km/s std case
yce = cuty
zce = cutz

;xsc = 485./2.*!unit.pc   ; x scale ; 100 km/s
xsc = 55*!unit.pc      ; lowres_.1t lowres lowres_cooling

;hyz = 70.*!unit.pc       ; half of the distance y or z scale
hyz = 20.*!unit.pc       ; lowres_cooling_.1t lowres lowres_cooling

;xrg = [-0.001e19,xce+xsc]
xrg = [-0.001e19,xsc]
yrg = [yce-hyz,yce+hyz]
zrg = [zce-hyz,zce+hyz]

data = loaddata(fname,'dens',sample=1,xcoords=x, ycoords=y,zcoords=z,time=time $
      ,xrange=xrg,yrange=yrg,zrange=zrg)
denxyz = data
data = loaddata(fname,'pres',sample=1, xrange=xrg,yrange=yrg,zrange=zrg)
prexyz = data
data = loaddata(fname,'Bubb',sample=1, xrange=xrg,yrange=yrg,zrange=zrg)
bubxyz = data
print,'loaddata done!!'

cutzind_tmp = where(z ge cutz) & cutzind = cutzind_tmp[0]
;--------------------------------------------------------------------------------------------
; parameters
pw = 2.5d0
gammin = 1.d0
nu = 100.d9
;nu = 1.d9
wnu = 2.d0*!pi*nu
d   = 1.d4 *!unit.pc
;--------------------------------------------------------------------------------------------
; assume equipartition (B^2 / 8 pi = 3 pres)
B = sqrt(3.d0*prexyz * 8.d0*!pi)

bubxyz[where(bubxyz ge 1.e-6)] = 1.
E_tot = prexyz *bubxyz    ;/ denxyz * !unit.mh

CC = E_tot/!unit.me/!unit.c/!unit.c *(-2.d0+pw)/gammin^(2.d0-pw) 

; Rybicki p180
sina = 1.d0 ; sin(a) = 1 ; maximum value
P_tot = sqrt(3.d0)*!unit.e^3.d0*CC*B*sina/(2.d0*!pi*!unit.me*!unit.c*!unit.c*(pw+1.d0))$
       *gamma(pw/4. + 19./12.)*gamma(pw/4. - 1./12.) $
       *(!unit.me*!unit.c*wnu/(3.d0*!unit.e*B*sina))^((1.d0-pw)/2.d0)  

; self-absorption for synchrotron
selfabsc = sqrt(3.d0)*!unit.e^3.d0/(8.d0*!pi*!unit.me)  $
          *(3.d0*!unit.e/(2.d0*!pi*!unit.me^3.d0*!unit.c^5.d0))^(pw/2.d0) $
          *gamma((3.d0*pw+2.d0)/12.d0)*gamma((3.d0*pw+22.d0)/12.d0) $
          * CC * (B*sina)^((pw+2.d0)/2.d0) * nu^(-(pw+4.d0)/2.d0)

dx = x[2]-x[1]

tauself = total(selfabsc,3)*dx  ; tau for self-absorption of synchrotron

prj_P_tot = total(P_tot,3)*dx^3.d0

sufb = prj_P_tot / 4.d0 /!pi / dx / dx

minv_sufb = alog10(min(sufb))
maxv_sufb = alog10(max(sufb))
;minv_sufb = -20.
;maxv_sufb = -17.880018
minv_sufb = -21.
maxv_sufb = -17.


minv_sufb2 = minv_sufb + 26. - 10. - alog10(4.25452)  ; convert to [mJy/arcsec^2] for log scale
maxv_sufb2 = maxv_sufb + 26. - 10. - alog10(4.25452)  ; convert to [mJy/arcsec^2] for log scale

nlev_sufb = 256
levs_sufb = fltarr(nlev_sufb+4)
levs_sufb[2:257] = findgen(nlev_sufb)/float(nlev_sufb)*(maxv_sufb - minv_sufb) + minv_sufb
levs_sufb[1] = minv_sufb - 1000. & levs_sufb[0] = levs_sufb[1] - 1000.
levs_sufb[258] = maxv_sufb + 1000. & levs_sufb[259] = levs_sufb[258] + 1000.

minv_T_b = alog10(!unit.c^2.d0/(2.d0*nu*nu*!unit.k)) + minv_sufb
maxv_T_b = alog10(!unit.c^2.d0/(2.d0*nu*nu*!unit.k)) + maxv_sufb

loadct,3,/sil
nx=800. & ny=700.
;nx=800. & ny=550.
if not keyword_set(ps) then ps=0
if (ps eq 1) then begin
psx = 20. & psy = 20.*ny/nx
mkeps,outdir+fname_r,xsize=psx,ysize=psy
;mkeps,'radio_sufb_0050_lw2.eps',xsize=psx,ysize=psy
endif else begin
!p.background=255 & !p.color=0
window,0,xs=nx,ys=ny
endelse
contour,alog10(sufb),x/!unit.pc,(y-cuty)/!unit.pc,/fill,/xst,/yst,/iso,levels=levs_sufb $
       ,xtitle=textoidl('x [pc] / (L/10^{37} ergs s^{-1})^{1/2}') $
       ,ytitle=textoidl('y [pc] / (L/10^{37} ergs s^{-1})^{1/2}')

legend,lab_r,/right,/top,textcolor=255,box=0

lab_bar2 = textoidl('Surface Brightness (Syncrotron) [mJy arcsec^{-2}] / (L/10^{37} ergs s^{-1})^{1/2} (Log scale)')
;lab_bar2 = textoidl('Surface Brightness (Syncrotron) [erg cm^{-2} s^{-1} Str^{-1} Hz^{-1}] (Log scale)')
lab_bar3 = textoidl('Brightness Temperature [K] / (L/10^{37} ergs s^{-1})^{1/2} (Log scale)')

xx0 = !p.clip[0] & xx1 = !p.clip[2]
yy0 = !p.clip[1] & yy1 = !p.clip[3]
;print,!p.clip
;xx0 = 90 & xx1 = 773
;yy0 = 60 & yy1 = 443
if (ps eq 1) then begin
  xx0 = xx0 / psx / 1000. * nx
  xx1 = xx1 / psx / 1000. * nx
  yy0 = yy0 / psy / 1000. * ny
  yy1 = yy1 / psy / 1000. * ny
endif

color_bar,pos=posnorm([xx0,yy1+50,xx1,yy1+70],nx=nx,ny=ny),/norm,/up,lim=[minv_sufb2,maxv_sufb2] $
         ,bartitle=lab_bar2,titlegap=30./ny,/white,dual_title=lab_bar3,dual_lim=[minv_T_b,maxv_T_b],dual_gap=40./ny
;color_bar,pos=posnorm([xx0,yy1+50,xx1,yy1+70],nx=nx,ny=ny),/norm,/up,lim=[minv_sufb,maxv_sufb] $
;         ,bartitle=lab_bar2,titlegap=30./ny,/white,dual_title=lab_bar3,dual_lim=[minv_T_b,maxv_T_b],dual_gap=40./ny
;color_bar,pos=posnorm([90,470,775,490],nx=nx,ny=ny),/norm,/up,lim=[minv_sufb,maxv_sufb],bartitle=lab_bar2,titlegap=30./ny,/white

if (ps eq 1) then epsfree

minv_tau = alog10(min(tauself))
maxv_tau = alog10(max(tauself))

nlev_tau = 256
levs_tau = fltarr(nlev_tau+4)
levs_tau[2:257] = findgen(nlev_sufb)/float(nlev_tau)*(maxv_tau - minv_tau) + minv_tau
levs_tau[1] = minv_tau - 1000. & levs_tau[0] = levs_tau[1] - 1000.
levs_tau[258] = maxv_tau + 1000. & levs_tau[259] = levs_tau[258] + 1000.

!p.background=255 & !p.color=0
window,1,xs=nx,ys=ny
contour,alog10(tauself),x/!unit.pc,(y-cuty)/!unit.pc,/fill,/xst,/yst,/iso,levels=levs_tau $&
       ,xtitle=textoidl('x [pc] / (L/10^{37} ergs s^{-1})^{1/2}') $&
       ,ytitle=textoidl('y [pc] / (L/10^{37} ergs s^{-1})^{1/2}')

xx0 = !p.clip[0] & xx1 = !p.clip[2]
yy0 = !p.clip[1] & yy1 = !p.clip[3]
color_bar,pos=posnorm([xx0,yy1+50,xx1,yy1+70],nx=nx,ny=ny),/norm,/up,lim=[minv_tau,maxv_tau] $&
         ,bartitle='Tau (Log Scale)',titlegap=30./ny,/white

stop
;Flux = prj_P_tot / (4.d0*!pi*d*d)

;ss = size(Flux)
;window,0, xs=ss[1], ys=ss[2] + 100
;tvscl,alog10(Flux)
;color_bar,lim = [min(alog10(Flux)),max(alog10(Flux))],/up,pos=[30,449,796-30,465],bartitle='log(Flux) @ Radio',titlegap=30

;window,1
;plot,Flux,/ylog,psym=3,xtitle='pixels',ytitle='Flux [erg/s/cm^2]',/xst,yst=2

;arcs_pix = dx/d*180./!pi*60.*60.  ; arcsecond per pixel w/ distance of 'd'
;;arcs_pix = dx/d*180./!pi*60.  ; arcminute per pixel w/ distance of 'd'

;Flux_mjy = 1.d26 * Flux  ; dimension of mJy (Jy = 10^-23 cgs)
;
;Flux_arcs2 = Flux_mjy / arcs_pix / arcs_pix

; Brightness Temperature
;T_b = 1.36d0*(!unit.c/nu)^2.d0 * Flux_arcs2
;T_b = 1.54d0*(!unit.c/nu)^2.d0 * Flux_arcs2

T_b = !unit.c^2.d0/(2.d0*nu*nu*!unit.k) * sufb
ss = size(T_b)

if (ps eq 1) then begin
mkeps,'radio_brT.eps',xsize=20.,ysize=20.*ny/nx
endif else begin
loadct,3,/sil
!p.background=255 & !p.color=0
window,1,xs=nx,ys=ny
endelse

minv_T = -2. 
maxv_T = 0.
nlev_T = 256
levs_T = fltarr(nlev_T+4)
levs_T[2:257] = findgen(nlev_T)/float(nlev_T)*(maxv_T - minv_T) + minv_T
levs_T[1] = minv_T - 1000. & levs_T[0] = levs_T[1] - 1000.
levs_T[258] = maxv_T + 1000. & levs_T[259] = levs_T[258] + 1000.

;minT = min(T_b) &  maxT = max(T_b)
;tv,bytscl(alog10(T_b),min=alog10(minT),max=alog10(maxT))
contour,alog10(T_b),x,y,/fill,/xst,/yst,/iso,levels=levs_T,xtitle='x [cm]', ytitle='y [cm]'

lab_bar3 = 'log(Brightness Temperature) [K] @ Radio'
color_bar,pos=posnorm([90,470,775,490],nx=nx,ny=ny),/norm,/up,lim=[minv_T,maxv_T],bartitle=lab_bar3,titlegap=30./ny,/white

if (ps eq 1) then epsfree

stop
end
