pro show_inc4_45,mkdata=mkdata

dir='/home/astrodoo/Work/Data/2015_BowPWN/out_PWN3d_600kms_45deg_highRes/'

;restore,dir+'cut2d.sav'


nfiles = [50,70,80,100]
if keyword_set(mkdata) then begin
   restore,file=dir+'cut2d.sav'
   dens = ds[nfiles].d
   time = ds[nfiles].t
   save,file='cut2d_inc4_45.sav',dens,x,y,time
endif else restore,file='cut2d_inc4_45.sav'

xcutind1 = (where(x ge -7.e18))[0] & xcutind2 = (where(x ge 1.3e19))[0]
ycutind1 = (where(y ge -5.e18))[0] & ycutind2 = (where(y ge 5.e18))[0]

dens2 = dens[xcutind1:xcutind2,ycutind1:ycutind2,*]
xx = x[xcutind1:xcutind2]
yy = y[ycutind1:ycutind2]

szd = size(dens2,/dimension)

loadct,0,/sil

pltx0 = 320. & plty0 = 170.
winxs = pltx0+szd[0]*2.+350. & winys= plty0+szd[1]*2.+30.

maxd = 5.e-24 & mind = 1.e-29

mkeps,'show_inc4_45', xs= 30., ys=30.*winys/winxs

tenc = 6.86367e+10
tvcoord,bytscl(alog10(dens2[*,*,0]),max=alog10(maxd),min=alog10(mind)),xx,yy,/axes $
    ,pos=[pltx0/winxs,(plty0+szd[1])/winys],/norm,imgsize=float(szd[0])/winxs,/black $
    ,ytickinterval=3.e18
legend,string((time[0]-tenc)/60./60./24./365.,format='(f7.2)')+ ' yr',box=0,pos=[-7e18,4.3e18],charsize=1.3
tvlct,r,g,b,/get
oplot,!x.crange, [0.,0.], line=2, color=fsc_color('cyan'),thick=3
;arrow,-3.e17,0.,-6.e17,0.,/data,thick=10,hsize=500.,color=0
inc = 45.*!dtor
xdiscont = -4.1181999e+18 + 6.e7*time[0]
oplot,!x.crange, tan(inc)*(!x.crange-xdiscont), line=3, color=fsc_color('crimson'),thick=3
arrow,-2.5e18,0.,-5.e18,0.,/data,thick=10,hsize=500.,color=0
legend,textoidl('\varphi=45^{\circ}'),pos=[-5e18,-3e18],/data,box=0
tvlct,r,g,b

tvcoord,bytscl(alog10(dens2[*,*,1]),max=alog10(maxd),min=alog10(mind)),xx,yy,/axes $
    ,pos=[(pltx0+szd[0])/winxs,(plty0+szd[1])/winys],/norm,imgsize=float(szd[0])/winxs,/black, ytickformat='(a1)' $
    ,ytickinterval=3.e18
legend,string((time[1]-tenc)/60./60./24./365.,format='(f7.2)')+ ' yr',box=0,pos=[-7e18,4.3e18],charsize=1.3
tvlct,r,g,b,/get
oplot,!x.crange, [0.,0.], line=2, color=fsc_color('cyan'),thick=3
;arrow,-3.e17,0.,-6.e17,0.,/data,thick=10,hsize=500.,color=0
xdiscont = -4.1181999e+18 + 6.e7*time[1]
oplot,!x.crange, tan(inc)*(!x.crange-xdiscont), line=3, color=fsc_color('crimson'),thick=3
arrow,-2.5e18,0.,-5.e18,0.,/data,thick=10,hsize=500.,color=0
tvlct,r,g,b

tvcoord,bytscl(alog10(dens2[*,*,2]),max=alog10(maxd),min=alog10(mind)),xx,yy,/axes $
    ,pos=[pltx0/winxs,plty0/winys],/norm,imgsize=float(szd[0])/winxs,/black $
    ,ytickinterval=3.e18
legend,string((time[2]-tenc)/60./60./24./365.,format='(f7.2)')+ ' yr',box=0,pos=[-7e18,4.3e18],charsize=1.3
tvlct,r,g,b,/get
oplot,!x.crange, [0.,0.], line=2, color=fsc_color('cyan'),thick=3
;arrow,-3.e17,0.,-6.e17,0.,/data,thick=10,hsize=500.,color=0
xdiscont = -4.1181999e+18 + 6.e7*time[2]
oplot,!x.crange, tan(inc)*(!x.crange-xdiscont), line=3, color=fsc_color('crimson'),thick=3
arrow,-2.5e18,0.,-5.e18,0.,/data,thick=10,hsize=500.,color=0
tvlct,r,g,b

tvcoord,bytscl(alog10(dens2[*,*,3]),max=alog10(maxd),min=alog10(mind)),xx,yy,/axes $
    ,pos=[(pltx0+szd[0])/winxs,plty0/winys],/norm,imgsize=float(szd[0])/winxs,/black, ytickformat='(a1)' $
    ,ytickinterval=3.e18
legend,string((time[3]-tenc)/60./60./24./365.,format='(f7.2)')+ ' yr',box=0,pos=[-7e18,4.3e18],charsize=1.3

xyouts,(pltx0+szd[0]-60)/winxs,(plty0-120)/winys,/norm,'x [cm]'
xyouts,0.03,0.5,/norm,'y [cm]',orientation=90

oplot,!x.crange, [0.,0.], line=2, color=fsc_color('cyan'),thick=3

xdiscont = -4.1181999e+18 + 6.e7*time[3]
oplot,!x.crange, tan(inc)*(!x.crange-xdiscont), line=3, color=fsc_color('crimson'),thick=3
arrow,-2.5e18,0.,-5.e18,0.,/data,thick=10,hsize=500.,color=0

arrow,3.e18,-2.5e18,2.e18,-1.2e18,/data,thick=3,hsize=300.,color=fsc_color('blue')
xyouts,3.e18,-3.5e18,/data,'K',color=fsc_color('blue')

;plots,0.,0.,/data,psym=7,color=255,symsize=2

loadct,0,/sil

color_bar,lim=[mind,maxd],/log,/right,bartitle=textoidl('density [g cm^{-3}]'), titlegap=0.09, pos=posnorm([pltx0+2.*szd[0]+20.,plty0,pltx0+2.*szd[0]+80.,plty0+2*szd[1]],nx=winxs,ny=winys),/norm

epsfree

stop
end
