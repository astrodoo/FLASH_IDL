forward_function Rst

pro guitar_head2

spE = 1e33
d0 = 1.67d-24
d1 = d0*0.7
vs = 1.5d8
vp = 1.d10
Rst1 = Rst(d0=d0)
xoff1 = 8.e14

; read data
xra=[0.,3.e16]
yra=[-5.e15,7e16]
sample=1
dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_guitar_head/dl_0.1/'
fname='PWN2d_hdf5_plt_cnt_0500'

d = loaddata(dir+fname,'dens',sample=sample,xCoord=yy,yCoord=xx,xra=xra,yra=yra,time=tt)
d = transpose(d)

dir2 = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_guitar_head/dh_0.1/'
d2 = loaddata(dir2+fname,'dens',sample=sample,xCoord=yy,yCoord=xx,xra=xra,yra=yra)
d2 = transpose(d2)
d02 = d0/0.7
Rst2 = Rst(d0=d02)
xoff2 = 5.e14

dsz = size(d,/dimension)

xgr1 = -3.75d16 + vs*tt
xgr2 = -7.06d16 + vs*tt


;restore,file='analBow_0180.dat'
;x1 = x & w11 = w1 & w21 = w2
;restore,file='analBow_0180_d0.7.dat'
;x2 = x & w12 = w1 & w22 = w2

;analytic lines
dir='../data/'
restore,file=dir+'analBow2_0500.sav'
xan1 = x & w1an1 = w1
restore,file=dir+'analBow2_0500_0.7.sav'
xan07 = x & w1an07 = w1
restore,file=dir+'analBow2_0500_0.5.sav'
xan05 = x & w1an05 = w1
restore,file=dir+'analBow2_0500_0.3.sav'
xan03 = x & w1an03 = w1
restore,file=dir+'analBow2_0500_0.1.sav'
xan01 = x & w1an01 = w1

restore,file=dir+'analBow2_0500_3.sav'
xan3 = x & w1an3 = w1
restore,file=dir+'analBow2_0500_5.sav'
xan5 = x & w1an5 = w1
restore,file=dir+'analBow2_0500_7.sav'
xan7 = x & w1an7 = w1
restore,file=dir+'analBow2_0500_10.sav'
xan10 = x & w1an10 = w1

; drawing
x0 = 450. & y0 = 200.
yplt = 300.

xs = x0 + dsz[0]+420
ys = y0 + dsz[1]*2 + 100+2*yplt
psxs=20. & psys=psxs*ys/xs

loadct,0,/sil
mkeps,'guitar_head2',xs=psxs,ys=psys

maxd = max(d) & mind = min(d)
tvcoord,bytscl(alog10(d),max=alog10(maxd),min=alog10(mind)),xx,yy,psx=float(dsz[0])/xs $
       ,pos=[x0/xs,(y0+dsz[1]+yplt)/ys],/axes,/black,ytickinterval=1.e16
oplot,[xgr1,xgr1],!y.crange,line=1
oplot,[xgr2,xgr2],!y.crange,line=1
anacl = 'magenta'
loadct,39,/sil
oplot,xan1[where(xan1 le 6.e16)]-Rst(d0=d0)-xoff1,w1an1[where(xan1 le 6.e16)],color=50
oplot,xan07[where(xan07 le 5.e16)]-Rst(d0=d0*0.7)-xoff1,w1an07[where(xan07 le 5.e16)],color=100
oplot,xan05[where(xan05 le 4.e16)]-Rst(d0=d0*0.5)-xoff1,w1an05[where(xan05 le 4.e16)],color=150
oplot,xan03[where(xan03 le 3.e16)]-Rst(d0=d0*0.3)-xoff1,w1an03[where(xan03 le 3.e16)],color=200
oplot,xan01[where(xan01 le 2.e16)]-Rst(d0=d0*0.1)-xoff1,w1an01[where(xan01 le 2.e16)],color=250
arrow,3.e16,2.e16,1.5e16,2.e16,/data
xyouts,1.8e16,2.2e16,/data,textoidl('\nabla\rho_{0} < 0')

legend,[' 0.1 ',' 0.3 ',' 0.5 ',' 0.7 ','     ']+textoidl('\rho_{0,init}'),textcolor=[250,200,150,100,50],box=0,/left,/top,charsize=1.
;oplot,x1-Rst1-xoff1,w11
;oplot,x1-Rst1-xoff1,w21
;xyouts,1.e15,8.5e15,/data,textoidl('\rho_{1}')

dgrad1 = smooth(d[*,dsz[1]-1],200)
plot,xx,dgrad1,position=posnorm([x0,y0+dsz[1]*2+yplt,x0+dsz[0],y0+dsz[1]*2+yplt*2],nx=xs,ny=ys),/norm,/noerase,/xst,yst=2 $
    ,xtickformat='(a1)',ytickinterval=1.e-24,ytitle=textoidl('\rho_{0}')
oplot,[xgr1,xgr1],!y.crange,line=1
oplot,[xgr2,xgr2],!y.crange,line=1

loadct,0,/sil
maxd2 = max(d2) & mind2 = min(d2)
tvcoord,bytscl(alog10(reverse(d2,2)),max=alog10(maxd2),min=alog10(mind2)),xx,-yy,psx=float(dsz[0])/xs $
       ,pos=[x0/xs,(y0+yplt)/ys],/axes,/black,ytickinterval=1.e16,xtickformat='(a1)'
oplot,[xgr1,xgr1],!y.crange,line=1
oplot,[xgr2,xgr2],!y.crange,line=1
loadct,39,/sil
oplot,xan1[where(xan1 le 6.e16)]-Rst(d0=d0)-xoff1,-w1an1[where(xan1 le 6.e16)],color=50
oplot,xan3[where(xan3 le 5.e16)]-Rst(d0=d0*3)-xoff1,-w1an3[where(xan3 le 5.e16)],color=100
oplot,xan5[where(xan5 le 4.e16)]-Rst(d0=d0*5)-xoff1,-w1an5[where(xan5 le 4.e16)],color=150
oplot,xan7[where(xan7 le 3.e16)]-Rst(d0=d0*7)-xoff1,-w1an7[where(xan7 le 3.e16)],color=200
oplot,xan10[where(xan10 le 2.e16)]-Rst(d0=d0*10)-xoff1,-w1an10[where(xan10 le 2.e16)],color=250
arrow,3.e16,-2.e16,1.5e16,-2.e16,/data
xyouts,1.8e16,-2.4e16,/data,textoidl('\nabla\rho_{0} > 0')

legend,['  10 ','   7 ','   5 ','   3 ','     ']+textoidl('\rho_{0,init}'),textcolor=[250,200,150,100,50],box=0,/left,/bottom,charsize=1.
;oplot,x2-Rst2-xoff2,-w12
;oplot,x2-Rst2-xoff2,-w22
;xyouts,1.e15,-8.5e15,/data,textoidl('\rho_{2}=0.7\rho_{1}')

dgrad2 = smooth(d2[*,dsz[1]-1],200)
plot,xx,dgrad2,position=posnorm([x0,y0,x0+dsz[0],y0+yplt],nx=xs,ny=ys),/norm,/noerase,/xst,yst=2 $
    ,ytickinterval=1.e-23,ytitle=textoidl('\rho_{0}'),xtitle='x [cm]'
oplot,[xgr1,xgr1],!y.crange,line=1
oplot,[xgr2,xgr2],!y.crange,line=1


xyouts,0.04,0.45,/norm,'y [cm]',orientation=90

loadct,0,/sil

color_bar,/right,position=posnorm([x0+dsz[0]+20,y0+dsz[1]+yplt,x0+dsz[0]+70,y0+dsz[1]*2+yplt],nx=xs,ny=ys),/norm,lim=[mind,maxd],/log ;, bartitle=textoidl('\rho [g cm^{-3}]'),titlegap=0.12
color_bar,/right,position=posnorm([x0+dsz[0]+20,y0+yplt,x0+dsz[0]+70,y0+dsz[1]+yplt],nx=xs,ny=ys),/norm,lim=[mind2,maxd2],/log

xyouts,0.98,0.45,/normal,textoidl('\rho [g cm^{-3}]'),orientation=90.



epsfree

stop
end

function Rst,spE=spE,d0=d0,vs=vs,vp=vp
if not keyword_set(spE) then spE = 1e33
if not keyword_set(d0) then d0 = 1.67d-24
if not keyword_set(vs) then vs = 1.5d8
if not keyword_set(vp) then vp = 1.d10

return, sqrt(spE / (4.d0*!dpi*d0*vs^2.*vp))
end


pro comp_bow

; read data
xra=[0.,3.e16]
yra=[-5.e15,7e16]
sample=1

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_guitar_head/dl_0.1/'
fname='PWN2d_hdf5_plt_cnt_0500'
d1 = loaddata(dir+fname,'dens',sample=sample,xCoord=y1,yCoord=x1,xra=xra,yra=yra,time=tt)
d1 = transpose(d1)
szd1 = size(d1,/dimension)

bow1 = fltarr(szd1[0])
smd1 = smooth(d1,10)
for i=0,szd1[0]-1 do bow1[i] = y1[(where(smd1[i,*] eq max(smd1[i,*])))[0]]
Rst1 = Rst(d0=1.67d-24*0.7)
xoff1 = Rst1 + 3.e15
x12_ind = where(x1 ge -xoff1)
X12 = x1[x12_ind]
bow1 = bow1[x12_ind]

patch_ind = where((x1 ge 5.e15) and (x1 le 1.36e16))

smd1_patch = smooth(d1[patch_ind,*],100)
for j=0, n_elements(patch_ind)-1 do bow1[j+patch_ind[0]] = y1[(where(smd1_patch[j,*] eq max(smd1_patch[j,*])))[0]]


;window,0
;plot,x1,bow1,/iso
;oplot,[x1[patch_ind[0]],x1[patch_ind[0]]],!y.crange,line=2
;oplot,[x1[patch_ind[n_elements(patch_ind)-1]],x1[patch_ind[n_elements(patch_ind)-1]]],!y.crange,line=2
;stop

dir2 = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_guitar_head/dh_0.1/'
d2 = loaddata(dir2+fname,'dens',sample=sample,xCoord=y2,yCoord=x2,xra=xra,yra=yra)
d2 = transpose(d2)
szd2 = size(d2,/dimension)
bow2 = fltarr(szd2[0])
smd2 = smooth(d2,10)
for i=0,szd2[0]-1 do bow2[i] = y2[(where(smd2[i,*] eq max(smd2[i,*])))[0]]
Rst2 = Rst(d0=1.67d-24/0.7)
x22_ind = where(x2 ge -Rst2)
X22 = x2[x22_ind]
bow2 = bow2[x22_ind]

dir3 = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_guitar_head/nochange_d/'
d3 = loaddata(dir3+fname,'dens',sample=sample,xCoord=y3,yCoord=x3,xra=xra,yra=yra)
d3 = transpose(d3)
szd3 = size(d3,/dimension)
bow3 = fltarr(szd3[0])
smd3 = smooth(d3,10)
for i=0,szd3[0]-1 do bow3[i] = y3[(where(smd3[i,*] eq max(smd3[i,*])))[0]]
Rst3 = Rst(d0=1.67d-24)
xoff3 = Rst3 + 3.e15
x32_ind = where(x3 ge -xoff3)
X32 = x3[x32_ind]
bow3 = bow3[x32_ind]

loadct,39,/sil
;window,xs=szd1[0],ys=szd1[1]
;tvcoord,alog(d1),x1,y1,/scale
;plots,x1,bow1,color=fsc_color('magenta')
mkeps,'guitar_head_comp',xs=20.,ys=20./2.2
plot,x12,bow1,/iso,xra=[-5.e15,7.e16],yra=[0.,3.e16],/xst,/yst,/nodata,ytitle='y [cm]',xtitle='x [cm]',xmargin=[12.,0.9],ytickinterval=1.e16
oplot,x12,bow1,color=250
oplot,x22,bow2,color=50
oplot,x32,bow3,color=0

plots,0.,0.,/data,psym=dsym(16,/fill)

legend,[textoidl('\nabla\rho_{0}<0'),textoidl('\nabla\rho_{0}>0'),textoidl('\rho_{0}=\rho_{0,init}')],color=[250,50,0],textcolor=[250,50,0],/left,/top,box=0, line=0

epsfree



stop
end
