pro jetthick2

;(-2.e12,0) (-2.6e12,2.395e12) (-3.8e12,4.791e12) (-5.e12,6.798e12)  for 1e36
;(-2.e12,0)  for 1e36
 
xc0 = -2.e12 
yc0 = 0. 
zc0 = 4.791e12 

sample=0
;n=474
;n=519
;n=443
;n=542
;n=395
n=404
fname='jetthick_404_smp0_3'
out=fname+'_fwhm'

dens = dload(n,var='dens',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
pres = dload(n,var='pres',xc=xc0,yc=yc0,zc=zc0,sample=sample)
jet  = dload(n,var='jet',xc=xc0,yc=yc0,zc=zc0,sample=sample)
velz = dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)

;stop
momf = dens*jet*velz
;momf = jet*velz

z0ind_tmp = where(z ge 0) & z0ind = z0ind_tmp[0]

sz = size(velz,/dimension)

z2 = z[z0ind:*]
thickx = fltarr(sz[2]-z0ind)
thicky = fltarr(sz[2]-z0ind)
jetd = fltarr(sz[2]-z0ind)
jetp = fltarr(sz[2]-z0ind)

for i=z0ind,sz[2]-1 do begin
    momfxy = reform(momf[*,*,i])

    momfx = total(momfxy,2)
    momfy = total(momfxy,1)

    momfx = (momfx-min(momfx))/(max(momfx)-min(momfx)) 
    momfy = (momfy-min(momfy))/(max(momfy)-min(momfy))

    jetx = x[where(momfx ge 0.5)]
    jety = y[where(momfy ge 0.5)]

    thickx[i-z0ind] = jetx[n_elements(jetx)-1] - jetx[0]
    thicky[i-z0ind] = jety[n_elements(jety)-1] - jety[0]

    jetd[i-z0ind] = max(dens[*,*,i])
    jetp[i-z0ind] = max(pres[*,*,i])
endfor

window,0
plot,z2,thicky,/xst

save,filename=out+'.sav',z2,thickx,thicky,jetd,jetp
stop
end

pro jetthick_comb,ps=ps

;fname='jetthick_443_smp0_'
fname='jetthick_404_smp0_'
restore, fname+'1_fwhm.sav'
z21 = z2 & thickx1=thickx & thicky1=thicky & jetd1=jetd & jetp1=jetp
restore, fname+'2_fwhm.sav'
z22 = z2 & thickx2=thickx & thicky2=thicky & jetd2=jetd & jetp2=jetp
restore, fname+'3_fwhm.sav'
z23 = z2 & thickx3=thickx & thicky3=thicky
 
loadct,39,/sil
if keyword_set(ps) then mkeps,'jetthick_comb_404_fwhm',charsize=1.5,xs=20.,ys=15. $
else begin
     !p.background=255 & !p.color=0
     window,0,xs=800,ys=600
endelse

plot,z21,thicky1,xr=[0.,5.e12],yr=[1.e10,1.e12],/ylog,xtitle='z [cm]' $
    ,ytitle='thickness of the jet [cm]'
oplot,z22,thicky2
oplot,z23,thicky3

oplot,z21,thickx1,color=50
oplot,z22,thickx2,color=50
oplot,z23,thickx3,color=50

; draw analytic expectations
;h0=3.28980e+10 ; 2.5d10 
;h0=5.2e+10 ; 2.5d10 
h0=25.e+10 ; 2.5d10 
;P0=1300 ;563.
P0=1000 ;563.
dw=2d-14
vw=2.5d8
gam=4./3.
ll=3.e12
zz=findgen(1000)/1000.*1.e13
;hh = h0*(P0/(dw*vw*vw*cos(th_an)^4.))^(1./2./gam)/cos(th_an)
hh = h0*(P0/(dw*vw*vw))^(1./2./gam) * (zz^2./ll^2.+1)^(1./gam)
oplot,zz,hh,line=2

;oplot,th_an/!dtor,thick_an,line=2,color=fsc_color('magenta')
legend,['simulation (para)','simulation (perp)','analytic solution'],/right,/bottom,box=0,line=[0,0,2] $
      ,color=[50,0,0] ,textcolor=[50,0,0],charsize=1.5

;oplot,[4.5e11,4.5e11],10.^!y.crange,line=1
;oplot,[2.e12,2.e12],10.^!y.crange,line=1
;xyouts,2.e11,7.e11,/data,'I',charsize=2.,charthick=3.
;xyouts,1.1e12,7.e11,/data,'II',charsize=2.,charthick=3.
;xyouts,2.4e12,7.e11,/data,'III',charsize=2.,charthick=3.

if keyword_set(ps) then epsfree
stop
end
