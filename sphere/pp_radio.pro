pro pp_radio,n,sample=sample $
   ,bh=bh,xc0=xc0,yc0=yc0,zc0=zc0 $
   ,ps=ps,out=out,bar=bar
device,decomposed=0

bh=1
sample=3

;--------------------------------------------------------------------------------------------
; parameters
pw = 2.5d0
gammin = 1.d0
;nu = 100.d9
nu = 1.d9
wnu = 2.d0*!pi*nu
d   = 1.d4 *!unit.pc
;--------------------------------------------------------------------------------------------

;n=200   ; sph_1e36 earlier case
if not keyword_set(out) then out='pp_radio_1e36_'+string(n,format='(I4.4)')+'.eps'

fname = 'JetSet_hdf5_plt_cnt_'+string(n,format='(I4.4)')
read_amr,fname,var='dens',parameters=params,tree=tree
time = params.time

jetOn = 9.76e4
time = time - jetOn

print,'time = ', time

; find position of BB
; temporary calculation to circular orbit. (in case of ellipse, need more tasks)
m1 = 4.d34
m2 = 2.d34
mtot = m1+m2
peri = 3.d12
th0 = !pi
th = th0 + sqrt(!unit.g * mtot / peri^3.d0) * time

th = th0
; center of mass, position of BB
cm2x = peri*m1/mtot*cos(th)
cm2y = peri*m1/mtot*sin(th)

; center of mass, position of Star
cm2xSt = peri*m2/mtot*cos(th+!pi)
cm2ySt = peri*m2/mtot*sin(th+!pi)

if keyword_set(bh) then begin
   xc0 = cm2x
   yc0 = cm2y
endif else begin
   if not keyword_set(xc0) then xc0 = 0.
   if not keyword_set(yc0) then yc0 = 0.
endelse
zc0 = 0.

if not keyword_set(sample) then sample=0
pre = dload(n,var='pres',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
jet = dload(n,var='jet',xc=xc0,yc=yc0,zc=zc0,sample=sample)

time = time - jetOn

sd = size(pre)

; assume equipartition (B^2 / 8 pi = 3 pres)
B = sqrt(3.d0*pre * 8.d0*!pi)

;jet[where(jet ge 1.e-6)] = 1.
E_tot = pre * jet    ;/ denxyz * !unit.mh

CC = E_tot/!unit.me/!unit.c/!unit.c *(-2.d0+pw)/gammin^(2.d0-pw) 

; Rybicki p180
sina = 1.d0 ; sin(a) = 1 ; maximum value
P_tot = sqrt(3.d0)*!unit.e^3.d0*CC*B*sina/(2.d0*!pi*!unit.me*!unit.c*!unit.c*(pw+1.d0))$
       *gamma(pw/4. + 19./12.)*gamma(pw/4. - 1./12.) $
       *(!unit.me*!unit.c*wnu/(3.d0*!unit.e*B*sina))^((1.d0-pw)/2.d0)  

dy = y[2]-y[1]

prj_P_tot = total(P_tot,2)*dy^3.d0

sufb = prj_P_tot / 4.d0 /!pi / dy / dy

maxv_sufb = 1.e-3
minv_sufb = 1.e-5

minv_sufb2 = minv_sufb * 2.35044e15  ; convert to [mJy/arcsec^2]
maxv_sufb2 = maxv_sufb * 2.35044e15  ; convert to [mJy/arcsec^2]

xs = sd[1] & ys = sd[2] & zs = sd[3]

loadct,3,/sil
x0 = 60 & y0 = 60
x1 = 10  & y1 = 10
if keyword_set(bar) then begin
   xs2 = xs + x0+x1+10 
   zs2 = zs + y0+y1+70 
endif else begin
   xs2 = xs + x0+x1
   zs2 = zs + y0+y1
endelse

if keyword_set(ps) then begin
   psxs = 20. & psys=20.*float(zs2)/xs2
   mkeps,out,xs=psxs,ys=psys
endif else window,1,xs=xs2,ys=zs2

plot,x/!unit.au,z/!unit.au,/xst,/yst,/nodata,position=posnorm([x0,y0,x0+xs,y0+zs],nx=xs2,ny=zs2),/norm $
    ,xtitle='x [AU]', ytitle='z [AU]',/iso
if keyword_set(bar) then $
   color_bar, pos=posnorm([x0,y0+zs+10,x0+xs,y0+zs+40],nx=xs2,ny=zs2),/norm,/up,lim=[minv_sufb2,maxv_sufb2] $
   ,bartitle=textoidl('Surface Brightness [mJy arcsec^{-2}]'),/log,/minor,titlegap=0.04

if keyword_set(ps) then $
   tv,bytscl(alog10(sufb),max=alog10(maxv_sufb),min=alog10(minv_sufb)),float(x0)/xs2 *psxs,float(y0)/zs2*psys,/centimeter $
   ,xsize=float(xs)/xs2*psxs, ysize=float(zs)/zs2*psys $
  else $ 
   tv,bytscl(alog10(sufb),max=alog10(maxv_sufb),min=alog10(minv_sufb)),float(x0)/xs2,float(y0)/zs2,/norm 

plot,x/!unit.au,z/!unit.au,/xst,/yst,/nodata,position=posnorm([x0,y0,x0+xs,y0+zs],nx=xs2,ny=zs2),/norm $
    ,xtitle='x [AU]', ytitle='z [AU]',/iso,/noerase
;legend,'X-Z',/top,/left,textcolor=255,box=0
legend,string(time/60./60.,format='(f6.2)')+' hrs',/right,/top,box=0,textcolor=255

;plots,ring(cm2xSt/!unit.au,0.,1.4e12/!unit.au),/data,line=2,color=0
if keyword_set(ps) then epsfree
stop
end
