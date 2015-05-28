pro guitar
dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar2/'
fname = 'PWN2d_hdf5_plt_cnt_0117'

d=transpose(loaddata(dir+fname,'dens',sample=2,xCoord=y,yCoord=x,time=t,xra=[0.,7.e17],yra=[-1.e17,2.5e18]))

dsz = size(d,/dimension)
d2 = fltarr(dsz[0],2*dsz[1])
y2 = fltarr(2*dsz[1])

d2[*,dsz[1]:2*dsz[1]-1] = d
d2[*,0:dsz[1]-1] = reverse(d,2)

y2 = [-reverse(y),y]

;draw
pltx0=400. & plty0=200.
pltxs=dsz[0] & pltys=2*dsz[1]
winxs=pltx0+pltxs+400 & winys=plty0+pltys+350


dd = fltarr(dsz[0])
disc1 = 1.11002e+18 & disc2 = 1.43202e+18
d00 = 2.76e-26 & d01 = 1.23e-24 & d02 = 1.23e-25
dd[where(x lt disc1)] = d02
dd[where((x ge disc1) and (x lt disc2))] = d01
dd[where(x gt disc2)] = d00

;for analytic bubble
Esp = 1.d33
rb1 = (125./154./!pi)^0.2*(Esp/d00)^0.2*t^0.6 
rb2 = (125./154./!pi)^0.2*(Esp/d02)^0.2*(t-4.3e+09)^0.6 



mkeps,'guitar',xs=20.,ys=20.*winys/winxs

tvcoord,alog10(d2),x,y2,/scale,position=[pltx0/winxs,plty0/winys],psx=float(dsz[0])/winxs,/axes,/black $
       ,xtitle='x [cm]', ytitle='y [cm]',xtickinterval=1.e18
oplot,[disc1,disc1],!y.crange,line=2
oplot,[disc2,disc2],!y.crange,line=2

plots, 0.,0.,psym=7,color=255,/data,symsize=1.5
plots,ring(t*1.5d8,0.,rb1),/data
plots,ring((t-4.3e+09)*1.5d8,0.,rb2),/data




color_bar,/right,lim=[min(d),max(d)],/log,bartitle=textoidl('density [g cm^{-3}]'),bargap=0.01,titlegap=0.12

plot,x,dd,pos=posnorm([pltx0,plty0+pltys,pltx0+pltxs,plty0+pltys+300],nx=winxs,ny=winys) $
    ,/norm,/noerase,/xst,xtickformat='(a1)',/ylog,ytitle=textoidl('\rho_{0}')
oplot,[disc1,disc1],10.^!y.crange,line=2
oplot,[disc2,disc2],10.^!y.crange,line=2

epsfree

stop
end
