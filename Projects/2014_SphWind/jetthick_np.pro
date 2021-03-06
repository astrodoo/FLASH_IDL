pro jetthick_np

;(-2.e12,0) (-2.6e12,2.395e12) (-3.8e12,4.791e12) (-5.e12,6.798e12)  for 1e36
 
xc0 = -2.e12
yc0 = 0. 
zc0 = 1.e12

sample=4
;n=474
;n=519
;n=850
;n=404
n=494
;n=900
;n=426
;n=453

fname='jetthick_'+strtrim(n,2)
;fname='jetthick_'+strtrim(n,2)+'_smp0.1'
crit=0.1
out=fname+'_crit'+string(crit,format='(f3.1)')
;dens = dload(n,var='dens',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;pres = dload(n,var='pres',xc=xc0,yc=yc0,zc=zc0,sample=sample)
jet  = dload(n,var='jet',xc=xc0,yc=yc0,zc=zc0,sample=sample)
velz = dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time,lref=lref)

LX = 1.9342e13+1.9158e13

z0ind_tmp = where(z ge 0) & z0ind = z0ind_tmp[0]

sz = size(velz,/dimension)

z2 = z[z0ind:*]
thickx = fltarr(sz[2]-z0ind)
thicky = fltarr(sz[2]-z0ind)
;x2 = fltarr(sz[2]-z0ind)
xx1 = fltarr(sz[2]-z0ind)
xx2 = fltarr(sz[2]-z0ind)
yy1 = fltarr(sz[2]-z0ind)
yy2 = fltarr(sz[2]-z0ind)

denj = fltarr(sz[2]-z0ind)
vj   = fltarr(sz[2]-z0ind)
vzj   = fltarr(sz[2]-z0ind)

xx = x#replicate(1.,512)
yy = transpose(y#replicate(1.,512))

;velz2 = velz > crit*3.e9
velz2 = velz*jet > crit*3.e9
velz2[where(velz2 eq crit*3.e9)]=0

lsmp2 = (x[2]-x[1])^2. 
np = fltarr(sz[2]-z0ind)
dx2_tmp = fltarr(sz[2]-z0ind)
maxl_tmp = fltarr(sz[2]-z0ind)
nsmp_tmp = fltarr(sz[2]-z0ind)

for i=z0ind,sz[2]-1 do begin
    velzxy = reform(velz2[*,*,i])
    lrefxy = reform(lref[*,*,i])

    jetind = where(velzxy ne 0,count)
    if (count ne 0) then begin
       xxx = xx[jetind]
       yyy = yy[jetind]

       thickx[i-z0ind] = max(xxx)-min(xxx)
       xx1[i-z0ind] = max(xxx)
       xx2[i-z0ind] = min(xxx)
       thicky[i-z0ind] = max(yyy)-min(yyy) 
       yy1[i-z0ind] = max(yyy)
       yy2[i-z0ind] = min(yyy)

       lrefj = lrefxy[jetind]
       nsmp = n_elements(jetind)

       dx = Lx / (2.^(lrefj+2))
       dx2ave = total(dx*dx) / nsmp
       np[i-z0ind] = nsmp*lsmp2 / dx2ave
  
       dx2_tmp[i-z0ind] = dx2ave
       maxl_tmp[i-z0ind] = max(lrefj)
       nsmp_tmp[i-z0ind] = nsmp

    endif else begin
       thickx[i-z0ind] = !values.f_nan
       xx1[i-z0ind] = !values.f_nan 
       xx2[i-z0ind] = !values.f_nan
       thicky[i-z0ind] = !values.f_nan
       yy1[i-z0ind] = !values.f_nan
       yy2[i-z0ind] = !values.f_nan

       np[i-z0ind] = !values.f_nan
       dx2_tmp[i-z0ind] = !values.f_nan
       maxl_tmp[i-z0ind] = !values.f_nan
       nsmp_tmp[i-z0ind] = !values.f_nan
    endelse
endfor

window,1
plot,z2,np
;h0 = 2.* 1.25e10

window,0
plot,z2,thicky,/xst,xtitle='z [cm]',ytitle='thickness of the jet'

zshift = 2.5e12
;zshift = 2.e11
;zshift=0.
;h0=3.28980e+10 ; 2.5d10 
;h0=8.e11 ;5.2e+10 ; 2.5d10 
;h0=5.e11 ;5.2e+10 ; 2.5d10 
;h0=0.8e11 ;5.2e+10 ; 2.5d10 
h0 = 3.e11 ; 2.5d10 
;h0=2.5d10
;P0=1300 ;563.
P0=1300 ;563.
dw=2d-14
vw=2.5d8
gam=4./3.
ll=3.e12
zz=findgen(1000)/1000.*2.e13
;hh = h0*(P0/(dw*vw*vw*cos(th_an)^4.))^(1./2./gam)/cos(th_an)
hh = h0*(P0/(dw*vw*vw))^(1./2./gam) * (zz^2./ll^2.+1)^(1./gam)

;th2 = 15. * !dtor
;tanth2 = tan(th2)
;hh2 = h0*(P0/(dw*vw*vw))^(1./2./gam) * $ 
;      ( (ll^2.+zz^2.*tanth2^2.+2.*ll*zz*tanth2+zz^2.)/(ll^2.+ll*zz*tanth2) )^(1./gam)

oplot,zz+zshift,hh,line=2
;oplot,zz+zshift,hh2,line=2

;save,filename=out+'.sav',z2,thickx,thicky,xx1,xx2,yy1,yy2
stop
end

pro jetthick_comb,ps=ps

fname='jetthick_494'
restore, fname+'_crit0.05.sav'
z21 = z2 & thickx1=thickx & thicky1=thicky 
restore, fname+'_crit0.1.sav'
z22 = z2 & thickx2=thickx & thicky2=thicky
restore, fname+'_crit0.2.sav'
z23 = z2 & thickx3=thickx & thicky3=thicky
restore, fname+'_crit0.3.sav'
z24 = z2 & thickx4=thickx & thicky4=thicky

!p.charsize=1.5
loadct,39,/sil
if keyword_set(ps) then mkeps,fname+'_comb',xs=20.,ys=15. $
else begin
     !p.background=255 & !p.color=0
     window,0,xs=800,ys=600
endelse

plot,z21,thicky1,xtitle='z [cm]' ,ytitle='thickness of the jet [cm]',xmargin=[10,5],/nodata
oplot,z21,thicky1,color=0, line=2, thick=2
oplot,z22,thicky2,color=50, line=3, thick=2
oplot,z23,thicky3,color=150, line=4, thick=2
oplot,z24,thicky4,color=250, line=5, thick=2

;oplot,z21,thickx1,color=50
;oplot,z22,thickx2,color=50
;oplot,z23,thickx3,color=50

; draw analytic expectations
;h0=3.28980e+10 ; 2.5d10 
zshift = 2.e12
h0=3.e11 ;5.2e+10 ; 2.5d10 
;P0=1300 ;563.
P0=1250 ;563.
dw=2d-14
vw=2.5d8
gam=4./3.
ll=3.e12
zz=findgen(1000)/1000.*2.e13
;hh = h0*(P0/(dw*vw*vw*cos(th_an)^4.))^(1./2./gam)/cos(th_an)
hh = h0*(P0/(dw*vw*vw))^(1./2./gam) * (zz^2./ll^2.+1)^(1./gam)
oplot,zz+zshift,hh,line=0,thick=3.

;oplot,th_an/!dtor,thick_an,line=2,color=fsc_color('magenta')
;legend,['simulation (para)','simulation (perp)','analytic solution'],/right,/bottom,box=0,line=[0,0,2] $
;      ,color=[50,0,0] ,textcolor=[50,0,0],charsize=1.5

legend,'crit'+['0.05','0.1','0.2','0.3'],/left,/top,box=0, $
      color=[0,50,150,250] ,textcolor=[0,50,150,250],charsize=1.5, line=[2,3,4,5]

;oplot,[4.5e11,4.5e11],10.^!y.crange,line=1
;oplot,[2.e12,2.e12],10.^!y.crange,line=1
;xyouts,2.e11,7.e11,/data,'I',charsize=2.,charthick=3.
;xyouts,1.1e12,7.e11,/data,'II',charsize=2.,charthick=3.
;xyouts,2.4e12,7.e11,/data,'III',charsize=2.,charthick=3.

if keyword_set(ps) then epsfree
stop
end

pro jetthick_cont,ps=ps

fname='jetthick_494'
n=494
sample=4
restore, fname+'_crit0.05.sav'
z21 = z2 & y11=yy1 & y12=yy2
restore, fname+'_crit0.1.sav'
z22 = z2 & y21=yy1 & y22=yy2
restore, fname+'_crit0.2.sav'
z23 = z2 & y31=yy1 & y32=yy2
restore, fname+'_crit0.3.sav'
z24 = z2 & y41=yy1 & y42=yy2
 

dens = dload(n,var='dens',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
jet  = dload(n,var='jet',xc=xc0,yc=yc0,zc=zc0,sample=sample)
velz = dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)

loadct,0,/sil
if keyword_set(ps) then mkeps,fname+'_cont',charsize=1.5,xs=20.,ys=20. $
else begin
     !p.background=255 & !p.color=0
     window,0,xs=512,ys=512
endelse

veljet = velz*jet

veljet2 = veljet[0:255,*,*] & veljet22 = total(veljet2,1)
dens2 = dens[0:255,*,*] & dens22 = total(dens2,1)

tvcoord, alog(dens22),y,z,/scale

loadct,39,/sil
plots,y11,z21,/data,color=0
plots,y12,z21,/data,color=0
plots,y21,z22,/data,color=50
plots,y22,z22,/data,color=50
plots,y31,z23,/data,color=150
plots,y32,z23,/data,color=150
plots,y41,z24,/data,color=250
plots,y42,z24,/data,color=250


if keyword_set(ps) then epsfree
stop
end

pro withzoom, ps=ps

;fname='jetthick_494'
fname='jetthick_900'
restore,file=fname+'_crit0.1.sav'
z2_4 = z2 & thicky_4 = thicky & yy1_4 = yy1 & yy2_4 = yy2
restore,file=fname+'_smp0.1_crit0.1.sav'
z2_0 = z2 & thicky_0 = thicky & yy1_0 = yy1 & yy2_0 = yy2

if not keyword_set(ps) then begin 
!p.background=255 & !p.color=0
window,0,xs=800,ys=600
endif else mkeps,fname+'_withzoom_1',xs=20.,ys=20.*6./8.
!p.multi=[0,1,2,0,1]
!x.charsize=1.5 & !y.charsize=1.5

!x.margin=[14,2] & !y.margin=[5,0]

sysp0 = !p & sysx0 = !x & sysy0 = !y
;plot,z2_4,thicky_4,/iso,/xst,xtitle='z [cm]', ytitle='thickness of jet [cm]' ,/nodata
plot,z2_4,thicky_4,/iso,/xst,xtitle='z [cm]', ytitle='thickness of jet [cm]' ,/nodata,xra=[0.,1.3e13]

x1 = min(z2_0) & x2 = max(z2_0)
y1 = min(thicky_0) & y2 = max(thicky_0)
y2 = 5.e11
polyfill, [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],/data,color=trp_color(fsc_color('yellow'),alpha=0.3)
plots, [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],/data,color=fsc_color('yellow'),thick=2

zshift = 2.5e12
;h0 = 3.e11 ; for 494 1e37
h0 = 5.e11 ; for 900 1e36
P0=1300 ;563.
dw=2d-14
vw=2.5d8
gam=4./3.
ll=3.e12
zz=findgen(1000)/1000.*2.e13
hh = h0*(P0/(dw*vw*vw))^(1./2./gam) * (zz^2./ll^2.+1)^(1./gam)
oplot,zz+zshift,hh,line=2

equiz = 1.2e11
plot,z2_0,thicky_0,/iso,/xst,xtitle='z [cm]', ytitle='thickness of jet [cm]'
oplot,[equiz,equiz],!y.crange,line=2

!p=sysp0 & !x=sysx0 & !y=sysy0
;plot,z2_4,thicky_4,/iso,/xst,xtitle='z [cm]', ytitle='thickness of jet [cm]' , /noerase
plot,z2_4,thicky_4,/iso,/xst,xtitle='z [cm]', ytitle='thickness of jet [cm]' , /noerase,xra=[0.,1.3e13]

!p.multi=0

if keyword_set(ps) then epsfree

if not keyword_set(ps) then begin 
!p.background=255 & !p.color=0
window,1,xs=800,ys=600
endif else mkeps,fname+'_withzoom_2',xs=20.,ys=20.*6./8.
!p.multi=[0,2,1,0,1]

!x.margin=[17,2] & !y.margin=[5,5]

sysp1 = !p & sysx1 = !x & sysy1 = !y
plot,yy1_4,z2_4,xtitle='y [cm]',ytitle='z [cm]',/iso,xr=[-4.e12,4.e12], xticks=2
oplot,yy2_4,z2_4

x1 = -5.e11 & x2 = 5.e11
y1 = min(z2_0) & y2 = max(z2_0)
polyfill, [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],/data,color=trp_color(fsc_color('yellow'),alpha=0.3)
plots, [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],/data,color=fsc_color('yellow'),thick=2

plot,yy1_0,z2_0,xtitle='y [cm]',ytitle='z [cm]',/iso,xr=[x1,x2], xticks=2,/xst,/yst
oplot,yy2_0,z2_0
oplot,!x.crange,[equiz,equiz],line=2

!p=sysp1 & !x=sysx1 & !y=sysy1
plot,yy1_4,z2_4,xtitle='y [cm]',ytitle='z [cm]',/iso,xr=[-4.e12,4.e12], xticks=2 ,/noerase
oplot,yy2_4,z2_4

if keyword_set(ps) then epsfree
!p.multi=0
stop
end
