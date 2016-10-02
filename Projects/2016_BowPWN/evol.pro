pro evol

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_600kms_2/'
fname1 = 'PWN2d_hdf5_plt_cnt_0005'
fname2 = 'PWN2d_hdf5_plt_cnt_0023'
fname3 = 'PWN2d_hdf5_plt_cnt_0050'
sample=2

xra=[-2.e17,1.5e18] & yra=[0.,1.5e18]
d1 = loaddata(dir+fname1,'dens',sample=sample,xrange=yra, yrange=xra, xCoord=y, yCoord=x,time=t1)
d1 = transpose(d1)
szd1 = size(d1,/dimension)
d1_2 = fltarr(szd1[0],szd1[1]*2)
d1_2[*,0:szd1[1]-1] = reverse(d1,2)
d1_2[*,szd1[1]:szd1[1]*2-1] = d1

y_2 = fltarr(szd1[1]*2)
y_2[0:szd1[1]-1] = -reverse(y)
y_2[szd1[1]:szd1[1]*2-1] = y

d2 = loaddata(dir+fname2,'dens',sample=sample,xrange=yra, yrange=xra,time=t2)
d2 = transpose(d2)
d2_2 = fltarr(szd1[0],szd1[1]*2)
d2_2[*,0:szd1[1]-1] = reverse(d2,2)
d2_2[*,szd1[1]:szd1[1]*2-1] = d2

xra3=[-2.e17,2.8e18] 
d3 = loaddata(dir+fname3,'dens',sample=sample,xrange=yra, yrange=xra3, yCoord=x3,time=t3)
d3 = transpose(d3)
szd3 = size(d3,/dimension)
d3_2 = fltarr(szd3[0],szd3[1]*2)
d3_2[*,0:szd3[1]-1] = reverse(d3,2)
d3_2[*,szd3[1]:szd3[1]*2-1] = d3


dir2 = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_600kms/l10_2/'
fname4 = 'PWN2d_hdf5_plt_cnt_0145'
fname5 = 'PWN2d_hdf5_plt_cnt_0205'
fname6 = 'PWN2d_hdf5_plt_cnt_0291'

xra4=[-5.e17,3.183e18] & yra4=[0.,4e18]
sample=3
d4 = loaddata(dir2+fname4,'dens',sample=sample,xrange=yra4, yrange=xra4, xCoord=y4, yCoord=x4,time=t4)
d4 = transpose(d4)
szd4 = size(d4,/dimension)
d4_2 = fltarr(szd4[0],szd4[1]*2)
d4_2[*,0:szd4[1]-1] = reverse(d4,2)
d4_2[*,szd4[1]:szd4[1]*2-1] = d4

y4_2 = fltarr(szd4[1]*2)
y4_2[0:szd4[1]-1] = -reverse(y4)
y4_2[szd4[1]:szd4[1]*2-1] = y4

d5 = loaddata(dir2+fname5,'dens',sample=sample,xrange=yra4, yrange=xra4,time=t5)
d5 = transpose(d5)
d5_2 = fltarr(szd4[0],szd4[1]*2)
d5_2[*,0:szd4[1]-1] = reverse(d5,2)
d5_2[*,szd4[1]:szd4[1]*2-1] = d5

xra6=[-5.e17,6.e18] 
d6 = loaddata(dir2+fname6,'dens',sample=sample,xrange=yra4, yrange=xra6, yCoord=x6,time=t6)
d6 = transpose(d6)
szd6 = size(d6,/dimension)
d6_2 = fltarr(szd6[0],szd6[1]*2)
d6_2[*,0:szd6[1]-1] = reverse(d6,2)
d6_2[*,szd6[1]:szd6[1]*2-1] = d6


;draw
pltx0=200. & plty0=100.
winxs=pltx0+szd1[0]*2+szd3[0]+20. & winys=plty0+szd1[1]*2.+120.
psxs=20. & psys=20.*winys/winxs
mkeps,'evol_init',xs=psxs, ys=psys
tvcoord,alog(d1_2),x,y_2,/scale,pos=[pltx0/winxs,plty0/winys],psx=szd1[0]/winxs,/norm,/axes,/black,xtickinterval=1.e18,ytitle='y [cm]'
legend,textoidl('t=')+string(fix(t1/60./60./24./365.),format='(I4)')+' yr',/left,/top,charsize=0.8,box=0
color_bar,lim=[min(d1_2),max(d1_2)],/up,/log,charsize=0.8
tvcoord,alog(d2_2),x,y_2,/scale,pos=[(pltx0+szd1[0])/winxs,plty0/winys],psx=szd1[0]/winxs,/norm,/axes,/black,xtickinterval=1.e18,ytickformat='(a1)',xtitle='x [cm]'
legend,textoidl('t=')+string(fix(t2/60./60./24./365.),format='(I4)')+' yr',/left,/top,charsize=0.8,box=0
color_bar,lim=[min(d2_2),max(d2_2)],/up,/log,charsize=0.8
tvcoord,alog(d3_2),x3,y_2,/scale,pos=[(pltx0+szd1[0]*2)/winxs,plty0/winys],psx=szd3[0]/winxs,/norm,/axes,/black,xtickinterval=1.e18,ytickformat='(a1)'
legend,textoidl('t=')+string(fix(t3/60./60./24./365.),format='(I4)')+' yr',/left,/top,charsize=0.8,box=0
color_bar,lim=[min(d3_2),max(d3_2)],/up,/log,charsize=0.8

xyouts,0.5,0.91,/norm,textoidl('\rho [g cm^{-3}]'),charsize=1.

epsfree


td=6.03e10+(5.e17/6.e7) & vp = 6.d7
pltx0=200. & plty0=100.
winxs=pltx0+szd4[0]*2+szd6[0]+20. & winys=plty0+szd4[1]*2.+120.
psxs=20. & psys=20.*winys/winxs
mkeps,'evol_disc',xs=psxs, ys=psys
tvcoord,alog(d4_2),x4,y4_2,/scale,pos=[pltx0/winxs,plty0/winys],psx=szd4[0]/winxs,/norm,/axes,/black,xtickinterval=2.e18,ytitle='y [cm]'
oplot,[(t4-td)*vp,(t4-td)*vp],!y.crange,line=2
legend,textoidl('t-t_{d}=')+string(fix((t4-td)/60./60./24./365.),format='(I4)')+' yr',/left,/top,charsize=0.8,box=0
color_bar,lim=[min(d4_2),max(d4_2)],/up,/log,charsize=0.8;,xtickinterval=2.
tvcoord,alog(d5_2),x4,y4_2,/scale,pos=[(pltx0+szd4[0])/winxs,plty0/winys],psx=szd4[0]/winxs,/norm,/axes,/black,xtickinterval=2.e18,ytickformat='(a1)',xtitle='x [cm]'
oplot,[(t5-td)*vp,(t5-td)*vp],!y.crange,line=2
legend,textoidl('t-t_{d}=')+string(fix((t5-td)/60./60./24./365.),format='(I4)')+' yr',/left,/top,charsize=0.8,box=0
color_bar,lim=[min(d5_2),max(d5_2)],/up,/log,charsize=0.8;,xtickinterval=2.
tvcoord,alog(d6_2),x6,y4_2,/scale,pos=[(pltx0+szd4[0]*2)/winxs,plty0/winys],psx=szd6[0]/winxs,/norm,/axes,/black,xtickinterval=2.e18,ytickformat='(a1)'
oplot,[(t6-td)*vp,(t6-td)*vp],!y.crange,line=2
legend,textoidl('t-t_{d}=')+string(fix((t6-td)/60./60./24./365.),format='(I4)')+' yr',/left,/top,charsize=0.8,box=0
color_bar,lim=[min(d6_2),max(d6_2)],/up,/log,charsize=0.8

xyouts,0.5,0.93,/norm,textoidl('\rho [g cm^{-3}]'),charsize=1.

epsfree







stop
end
