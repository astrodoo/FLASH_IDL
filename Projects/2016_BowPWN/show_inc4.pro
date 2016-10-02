pro show_inc4,mkdata=mkdata

dir='/home/astrodoo/Work/Data/2015_BowPWN/PWN3d_05inc/'

nfiles = [22,81,227,370]
if keyword_set(mkdata) then begin
   restore,file=dir+'cut2d.sav'
   dens = ds[nfiles].d
   time = ds[nfiles].t
   save,file='cut2d_inc4.sav',dens,x,y,time
endif else restore,file='cut2d_inc4.sav'

xx = x & yy = y

restore,dir+'analBow_d24.sav'
x_24 = x & w1_24 = w1

restore,dir+'analBow_d23.sav'
x_23 = x & w1_23 = w1

szd = size(dens,/dimension)

loadct,0,/sil


pltx0 = 320. & plty0 = 170.
winxs = pltx0+szd[0]*2.+350. & winys= plty0+szd[1]*2.+30.

maxd = 9.e-23 & mind = 2.e-29

mkeps,'show_inc4', xs= 30., ys=30.*winys/winxs

tenc = 2.e10
tvcoord,bytscl(alog10(dens[*,*,0]),max=alog10(maxd),min=alog10(mind)),xx,yy,/axes,pos=[pltx0/winxs,(plty0+szd[1])/winys],/norm,imgsize=float(szd[0])/winxs,/black
legend,string((time[0]-tenc)/60./60./24./365.,format='(f7.2)')+ ' yr',box=0,pos=[-1e18,9.e17],charsize=1.3
tvlct,r,g,b,/get
oplot,!x.crange, [0.,0.], line=2, color=fsc_color('cyan'),thick=3
inc = 5.*!dtor
xdiscont = -2.e18 + 1.e8*time[0]
oplot,!x.crange, tan(inc)*(!x.crange-xdiscont), line=3, color=fsc_color('crimson'),thick=3

arrow,-3.e17,0.,-6.e17,0.,/data,thick=5,hsize=400.,color=0
legend,textoidl('\varphi=5^{\circ}'),pos=[2.e18,8e17],/data,box=0
tvlct,r,g,b


tvcoord,bytscl(alog10(dens[*,*,1]),max=alog10(maxd),min=alog10(mind)),xx,yy,/axes,pos=[(pltx0+szd[0])/winxs,(plty0+szd[1])/winys],/norm,imgsize=float(szd[0])/winxs,/black, ytickformat='(a1)'
legend,string((time[1]-tenc)/60./60./24./365.,format='(f7.2)')+ ' yr',box=0,pos=[-1e18,9.e17],charsize=1.3
tvlct,r,g,b,/get
oplot,!x.crange, [0.,0.], line=2, color=fsc_color('cyan'),thick=3
xdiscont = -2.e18 + 1.e8*time[1]
oplot,!x.crange, tan(inc)*(!x.crange-xdiscont), line=3, color=fsc_color('crimson'),thick=3
arrow,-3.e17,0.,-6.e17,0.,/data,thick=5,hsize=400.,color=0
tvlct,r,g,b

tvcoord,bytscl(alog10(dens[*,*,2]),max=alog10(maxd),min=alog10(mind)),xx,yy,/axes,pos=[pltx0/winxs,plty0/winys],/norm,imgsize=float(szd[0])/winxs,/black
legend,string((time[2]-tenc)/60./60./24./365.,format='(f7.2)')+ ' yr',box=0,pos=[-1e18,9.e17],charsize=1.3
tvlct,r,g,b,/get
oplot,!x.crange, [0.,0.], line=2, color=fsc_color('cyan'),thick=3
xdiscont = -2.e18 + 1.e8*time[2]
oplot,!x.crange, tan(inc)*(!x.crange-xdiscont), line=3, color=fsc_color('crimson'),thick=3
arrow,-3.e17,0.,-6.e17,0.,/data,thick=5,hsize=400.,color=0
tvlct,r,g,b

tvcoord,bytscl(alog10(dens[*,*,3]),max=alog10(maxd),min=alog10(mind)),xx,yy,/axes,pos=[(pltx0+szd[0])/winxs,plty0/winys],/norm,imgsize=float(szd[0])/winxs,/black, ytickformat='(a1)'
legend,string((time[3]-tenc)/60./60./24./365.,format='(f7.2)')+ ' yr',box=0,pos=[-1e18,9.e17],charsize=1.3

xyouts,(pltx0+szd[0]-60)/winxs,(plty0-120)/winys,/norm,'x [cm]'
xyouts,0.03,0.47,/norm,'y [cm]',orientation=90

oplot,!x.crange, [0.,0.], line=2, color=fsc_color('cyan'),thick=3
xdiscont = -2.e18 + 1.e8*time[3]
oplot,!x.crange, tan(inc)*(!x.crange-xdiscont), line=3, color=fsc_color('crimson'),thick=3

arrow,-3.e17,0.,-6.e17,0.,/data,thick=5,hsize=400.,color=0


xoff1 = 2.e16
oplot,x_24-xoff1,w1_24 *1.1, color=fsc_color('magenta'), thick=5
oplot,x_23,-w1_23, color=fsc_color('magenta'),thick=5

plots,0.,0.,/data,psym=7,color=255,symsize=2

loadct,0,/sil

color_bar,lim=[mind,maxd],/log,/right,bartitle=textoidl('density [g cm^{-3}]'), titlegap=0.09, pos=posnorm([pltx0+2.*szd[0]+20.,plty0,pltx0+2.*szd[0]+80.,plty0+2*szd[1]],nx=winxs,ny=winys),/norm

epsfree

stop
end
