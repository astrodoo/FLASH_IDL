pro pp_jetthick,ps=ps

dir='/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_high/'
fname='jetthick_yz_850'
opt='crit7'
restore,filename=dir+fname+'_smp0_1_'+opt+'.sav'
z2_smp0_1_xy = z2 & thick_smp0_1_xy = thickxy
restore,filename=dir+fname+'_smp0_2_'+opt+'.sav'
z2_smp0_2_xy = z2 & thick_smp0_2_xy = thickxy
restore,filename=dir+fname+'_smp0_3_'+opt+'.sav'
z2_smp0_3_xy = z2 & thick_smp0_3_xy = thickxy

loadct,39
!p.background=255 & !p.color=0
if not keyword_set(ps) then window,1 $
   else mkeps,'pp_jetthick_1e36_850.eps',xs=20.,ys=20.*6./8.,charsize=1.5
plot,z2_smp0_1_xy,thick_smp0_1_xy,/xst,xtitle=textoidl('z [cm]'),ytitle=textoidl('thickness of the jet [cm]') $
    ,xr=[0.,1.e13],/yst, yr=[1.e10,1.e12],/ylog
oplot,z2_smp0_2_xy,thick_smp0_2_xy
oplot,z2_smp0_3_xy,thick_smp0_3_xy

; draw analytic expectations
h0=3.28980e+10 ; 2.5d10 
P0=1300 ;563.
dw=2d-14
vw=2.5d8
gam=4./3.
ll=3.e12
zz=findgen(1000)/1000.*1.e13
;hh = h0*(P0/(dw*vw*vw*cos(th_an)^4.))^(1./2./gam)/cos(th_an)
hh = h0*(P0/(dw*vw*vw))^(1./2./gam) * (zz^2./ll^2.+1)^(1./gam)
oplot,zz,hh,line=2

;oplot,th_an/!dtor,thick_an,line=2,color=fsc_color('magenta')
legend,['simulation','analytic solution'],/left,/top,box=0,line=[0,2],color=0 ,textcolor=0,charsize=1
if keyword_set(ps) then epsfree
stop
end
