pro jetthick3,opt=opt,ps=ps

if not keyword_set(opt) then opt='crit7'
bhx = -2.e12 
stx = 1.e12  
sep = 3.e12
avec = [bhx-stx,0]/sep
h0 = 2.*1.25e10

fname='jetthick_850'
;restore,filename=fname+'_smp4_'+opt+'.sav'
;z2_smp4 = z2 & thick_smp4 = thick & x2_smp4=x2
;ll = sqrt((x2_smp4-stx)^2. + (z2_smp4)^2.)
;bvec_smp4 = {x:(x2_smp4-stx)/ll, z:(z2_smp4)/ll}
;th_smp4 = acos(bvec_smp4.x*avec[0]+bvec_smp4.z*avec[1])
;
;restore,filename=fname+'_smp0_1_'+opt+'.sav'
;z2_smp0_1 = z2 & thick_smp0_1 = thick & x2_smp0_1=x2
;ll = sqrt((x2_smp0_1-stx)^2. + (z2_smp0_1)^2.)
;bvec_smp0_1 = {x:(x2_smp0_1-stx)/ll, z:(z2_smp0_1)/ll}
;th_smp0_1 = acos(bvec_smp0_1.x*avec[0]+bvec_smp0_1.z*avec[1])
;
;restore,filename=fname+'_smp0_2_'+opt+'.sav'
;z2_smp0_2 = z2 & thick_smp0_2= thick & x2_smp0_2=x2
;ll = sqrt((x2_smp0_2-stx)^2. + (z2_smp0_2)^2.)
;bvec_smp0_2 = {x:(x2_smp0_2-stx)/ll, z:(z2_smp0_2)/ll}
;th_smp0_2 = acos(bvec_smp0_2.x*avec[0]+bvec_smp0_2.z*avec[1])
;
;restore,filename=fname+'_smp0_3_'+opt+'.sav'
;z2_smp0_3 = z2 & thick_smp0_3 = thick & x2_smp0_3=x2
;ll = sqrt((x2_smp0_3-stx)^2. + (z2_smp0_3)^2.)
;bvec_smp0_3 = {x:(x2_smp0_3-stx)/ll, z:(z2_smp0_3)/ll}
;th_smp0_3 = acos(bvec_smp0_3.x*avec[0]+bvec_smp0_3.z*avec[1])
;
;restore,filename=fname+'_smp0_4_'+opt+'.sav'
;z2_smp0_4 = z2 & thick_smp0_4 = thick & x2_smp0_4=x2
;ll = sqrt((x2_smp0_4-stx)^2. + (z2_smp0_4)^2.)
;bvec_smp0_4 = {x:(x2_smp0_4-stx)/ll, z:(z2_smp0_4)/ll}
;th_smp0_4 = acos(bvec_smp0_4.x*avec[0]+bvec_smp0_4.z*avec[1])
;

fname='jetthick_yz_850'
restore,filename=fname+'_smp0_1_'+opt+'.sav'
z2_smp0_1_xy = z2 & thick_smp0_1_xy = thickxy
restore,filename=fname+'_smp0_2_'+opt+'.sav'
z2_smp0_2_xy = z2 & thick_smp0_2_xy = thickxy
restore,filename=fname+'_smp0_3_'+opt+'.sav'
z2_smp0_3_xy = z2 & thick_smp0_3_xy = thickxy

;opt='crit8'
;fname='jetthick_yz_850'
;restore,filename=fname+'_smp4_'+opt+'.sav'
;z2_smp4_xy = z2 & thick_smp4_xy = thickxy & x2_smp4_xy=x2
;denj_smp4_xy = denj & vj_smp4_xy = vj & vzj_smp4_xy = vzj
;ll = sqrt((x2_smp4_xy-stx)^2. + (z2_smp4_xy)^2.)
;bvec_smp4_xy = {x:(x2_smp4_xy-stx)/ll, z:(z2_smp4_xy)/ll}
;th_smp4_xy = acos(bvec_smp4_xy.x*avec[0]+bvec_smp4_xy.z*avec[1])
;
;restore,filename=fname+'_smp0_1_'+opt+'.sav'
;z2_smp0_1_xy = z2 & thick_smp0_1_xy = thickxy & x2_smp0_1_xy=x2
;denj_smp0_1_xy = denj & vj_smp0_1_xy = vj & vzj_smp0_1_xy = vzj
;ll = sqrt((x2_smp0_1_xy-stx)^2. + (z2_smp0_1_xy)^2.)
;bvec_smp0_1_xy = {x:(x2_smp0_1_xy-stx)/ll, z:(z2_smp0_1_xy)/ll}
;th_smp0_1_xy = acos(bvec_smp0_1_xy.x*avec[0]+bvec_smp0_1_xy.z*avec[1])
;
;restore,filename=fname+'_smp0_2_'+opt+'.sav'
;z2_smp0_2_xy = z2 & thick_smp0_2_xy= thickxy & x2_smp0_2_xy=x2
;denj_smp0_2_xy = denj & vj_smp0_2_xy = vj & vzj_smp0_2_xy = vzj
;ll = sqrt((x2_smp0_2_xy-stx)^2. + (z2_smp0_2_xy)^2.)
;bvec_smp0_2_xy = {x:(x2_smp0_2_xy-stx)/ll, z:(z2_smp0_2_xy)/ll}
;th_smp0_2_xy = acos(bvec_smp0_2_xy.x*avec[0]+bvec_smp0_2_xy.z*avec[1])
;
;restore,filename=fname+'_smp0_3_'+opt+'.sav'
;z2_smp0_3_xy = z2 & thick_smp0_3_xy = thickxy & x2_smp0_3_xy=x2
;denj_smp0_3_xy = denj & vj_smp0_3_xy = vj & vzj_smp0_3_xy = vzj
;ll = sqrt((x2_smp0_3_xy-stx)^2. + (z2_smp0_3_xy)^2.)
;bvec_smp0_3_xy = {x:(x2_smp0_3_xy-stx)/ll, z:(z2_smp0_3_xy)/ll}
;th_smp0_3_xy = acos(bvec_smp0_3_xy.x*avec[0]+bvec_smp0_3_xy.z*avec[1])

;restore,filename=fname+'_smp0_4_'+opt+'.sav'
;z2_smp0_4_xy = z2 & thick_smp0_4_xy = thickxy & x2_smp0_4_xy=x2
;denj_smp0_4_xy = denj & vj_smp0_4_xy = vj & vzj_smp0_4_xy = vzj
;ll = sqrt((x2_smp0_4_xy-stx)^2. + (z2_smp0_4_xy)^2.)
;bvec_smp0_4_xy = {x:(x2_smp0_4_xy-stx)/ll, z:(z2_smp0_4_xy)/ll}
;th_smp0_4_xy = acos(bvec_smp0_4_xy.x*avec[0]+bvec_smp0_4_xy.z*avec[1])

;restore,filename='moflux3_anal.sav'
;th_an = th & thick_an = h

loadct,39
!p.background=255 & !p.color=0
if not keyword_set(ps) then window,1 $
   else mkeps,'jetthick2_800_2'+opt+'.eps',xs=20.,ys=20.*6./8.
plot,z2_smp0_1_xy,thick_smp0_1_xy,/xst,xtitle=textoidl('z'),ytitle=textoidl('jet thick / 2r_{0}') $
    ,xr=[0.,1.e13],/yst,/nodata, yr=[1.e10,1.e12],/ylog
;oplot,z2_smp0_1,thick_smp0_1
;oplot,z2_smp0_2,thick_smp0_2
;oplot,th_smp0_3/!dtor,thick_smp0_3
;oplot,th_smp0_4/!dtor,thick_smp0_4

oplot,z2_smp0_1_xy,thick_smp0_1_xy
oplot,z2_smp0_2_xy,thick_smp0_2_xy
oplot,z2_smp0_3_xy,thick_smp0_3_xy

;oplot,z2_smp0_3_xy7,thick_smp0_3_xy7,color=100
;oplot,z2_smp0_3_xy7,thick_smp0_3_xy7,color=100
;oplot,z2_smp0_3_xy7,thick_smp0_3_xy7,color=100
;oplot,th_smp0_4_xy/!dtor,thick_smp0_4_xy,color=50
;
;plot,th_smp0_1/!dtor,thick_smp0_1,/xst,xtitle=textoidl('\theta [deg]'),ytitle=textoidl('jet thick / 2r_{0}') $
;    ,xr=[0.,80.],/yst,/nodata, yr=[0.,6.e11]
;oplot,th_smp0_1/!dtor,thick_smp0_1
;oplot,th_smp0_2/!dtor,thick_smp0_2
;oplot,th_smp0_3/!dtor,thick_smp0_3
;oplot,th_smp0_4/!dtor,thick_smp0_4

;oplot,th_smp0_1_xy/!dtor,thick_smp0_1_xy,color=50
;oplot,th_smp0_2_xy/!dtor,thick_smp0_2_xy,color=50
;oplot,th_smp0_3_xy/!dtor,thick_smp0_3_xy,color=50
;oplot,th_smp0_4_xy/!dtor,thick_smp0_4_xy,color=50
;
h0=3.28980e+10 ; 2.5d10 
P0=1300 ;563.
dw=2d-14
vw=2.5d8
gam=4./3.
ll=3.e12
zz=findgen(1000)/1000.*1.e13
;hh = h0*(P0/(dw*vw*vw*cos(th_an)^4.))^(1./2./gam)/cos(th_an)
hh = h0*(P0/(dw*vw*vw))^(1./2./gam) * (zz^2./ll^2.+1)^(1./gam)
oplot,zz,hh,color=250

;oplot,th_an/!dtor,thick_an,line=2,color=fsc_color('magenta')
legend,['x-z plane','x-y plane (perp. to wind)','analytic'],/left,/top,box=0,line=[0,0,2],color=[0,59,fsc_color('magenta')] $
      ,textcolor=[0,50,fsc_color('magenta')]
stop
legend,opt,/left,/center

if keyword_set(ps) then epsfree

momj_smp4 = denj_smp4_xy*vzj_smp4_xy^2.*thick_smp4_xy*thick_smp4
momj_smp0_1 = denj_smp0_1_xy*vzj_smp0_1_xy^2.*thick_smp0_1_xy*thick_smp0_1
momj_smp0_2 = denj_smp0_2_xy*vzj_smp0_2_xy^2.*thick_smp0_2_xy*thick_smp0_2
momj_smp0_3 = denj_smp0_3_xy*vzj_smp0_3_xy^2.*thick_smp0_3_xy*thick_smp0_3
momj_smp0_4 = denj_smp0_4_xy*vzj_smp0_4_xy^2.*thick_smp0_4_xy*thick_smp0_4
;momj_smp4_ = denj_smp4_xy*vj_smp4_xy^2.*thick_smp4_xy*thick_smp4
;momj_smp0_1_ = denj_smp0_1_xy*vj_smp0_1_xy^2.*thick_smp0_1_xy*thick_smp0_1
;momj_smp0_2_ = denj_smp0_2_xy*vj_smp0_2_xy^2.*thick_smp0_2_xy*thick_smp0_2
;momj_smp0_3_ = denj_smp0_3_xy*vj_smp0_3_xy^2.*thick_smp0_3_xy*thick_smp0_3
;momj_smp0_4_ = denj_smp0_4_xy*vj_smp0_4_xy^2.*thick_smp0_4_xy*thick_smp0_4

;window,3
;plot,th_smp4_xy/!dtor,momj_smp4, yr=[0.,9e26],/nodata,/yst,xtitle=textoidl('\theta [degree]') $
;    ,ytitle='momentum flux of jet'
;oplot,th_smp0_1_xy/!dtor,momj_smp0_1,color=50
;oplot,th_smp0_2_xy/!dtor,momj_smp0_2,color=50
;oplot,th_smp0_3_xy/!dtor,momj_smp0_3,color=50
;oplot,th_smp0_4_xy/!dtor,momj_smp0_4,color=50
;oplot,!x.crange,[3.3332790e+26,3.3332790e+26],line=2,color=0
;;oplot,th_smp0_1_xy/!dtor,momj_smp0_1_,color=250
;;oplot,th_smp0_2_xy/!dtor,momj_smp0_2_,color=250
;;oplot,th_smp0_3_xy/!dtor,momj_smp0_3_,color=250
;;oplot,th_smp0_4_xy/!dtor,momj_smp0_4_,color=250
;legend,['analytic','from data'],line=[2,0],color=[0,50],textcolor=[0,50],/right,/top,box=0
;
;window,4
;plot,th_smp4_xy/!dtor,momj_smp4, yr=[0.,6],/nodata,/yst,xtitle=textoidl('\theta [degree]') $
;    ,ytitle='thick_xy / thick_xz'
;oplot,th_smp0_1_xy/!dtor,thick_smp0_1/thick_smp0_1_xy
;oplot,th_smp0_2_xy/!dtor,thick_smp0_2/thick_smp0_2_xy
;oplot,th_smp0_3_xy/!dtor,thick_smp0_3/thick_smp0_3_xy
;oplot,th_smp0_4_xy/!dtor,thick_smp0_4/thick_smp0_4_xy
;
stop
end
