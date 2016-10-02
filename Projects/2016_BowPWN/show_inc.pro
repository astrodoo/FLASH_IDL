pro show_inc,mkdata=mkdata

dir='/home/jianiye/Work/Data/2015_BowPWN/PWN3d_05inc/'

nfile = 370
if keyword_set(mkdata) then begin
   restore,file=dir+'cut2d.sav'

   x1_ind = 80 & x2_ind = 600
   y1_ind = 52 & y2_ind = 460

;   d = ds[nfile].d[x1_ind:x2_ind,y1_ind:y2_ind]
   d = ds[nfile].d
;   x = x[x1_ind:x2_ind]
;   y = y[y1_ind:y2_ind]

;   save,file='cut2d_'+strtrim(nfile,2)+'_zoom.sav',d,x,y,(time=ds[nfile].t)
   save,file='cut2d_'+strtrim(nfile,2)+'.sav',d,x,y,(time=ds[nfile].t)
;endif else restore,file='cut2d_'+strtrim(nfile,2)+'_zoom.sav'
endif else restore,file='cut2d_'+strtrim(nfile,2)+'.sav'

xx = x & yy = y

restore,dir+'analBow_d24.sav'
x_24 = x & w1_24 = w1

restore,dir+'analBow_d23.sav'
x_23 = x & w1_23 = w1

szd = size(d,/dimension)

loadct,0,/sil


pltx0 = 220. & plty0 = 120.
winxs = pltx0+szd[0]+250. & winys= plty0+szd[1]+30.
;window,xs=szd[0]+150,ys=szd[1]+150

mkeps,'show_inc', xs= 25., ys=25.*winys/winxs

tvcoord,alog10(d),xx,yy,/scale,/axes,pos=[pltx0/winxs,plty0/winys],/norm,imgsize=float(szd[0])/winxs,/black,xtitle='x [cm]', ytitle='y [cm]'

oplot,!x.crange, [0.,0.], line=2, color=fsc_color('cyan'),thick=3

inc = 5.*!dtor
xdiscont = -2.e18 + 1.e8*time
oplot,!x.crange, tan(inc)*(!x.crange-xdiscont), line=3, color=fsc_color('crimson'),thick=3

arrow,-3.e17,0.,-6.e17,0.,/data,thick=10,hsize=500.,color=0


xoff1 = 2.e16
oplot,x_24-xoff1,w1_24 *1.1, color=fsc_color('magenta'), thick=10
oplot,x_23,-w1_23, color=fsc_color('magenta'),thick=10

plots,0.,0.,/data,psym=7,color=255,symsize=2

loadct,0,/sil

color_bar,lim=[min(d),max(d)],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),bargap=0.01, titlegap=0.1

epsfree

stop
end
