pro jet_snap,mkdata=mkdata

;file = '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L38_M30_lg/JetSet_hdf5_plt_cnt_2761'
file = '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L38_M30_lg/JetSet_hdf5_plt_cnt_2682'
;file = '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L38_M30_lw/JetSet_hdf5_plt_cnt_2682'

sample=5
;sample=4
xcut = -2.e12
ycut = 0.
;id = '1e38lw'
id = '1e38'

if keyword_set(mkdata) then begin
dyz = reform(loaddata(file,'dens',sample=sample,xrange=[xcut,xcut],xCoord=xx,yCoord=y,zCoord=z,time=time,lref=lref))
lyz = reform(lref)
pyz = reform(loaddata(file,'pres',sample=sample,xrange=[xcut,xcut]))
jyz = reform(loaddata(file,'jet',sample=sample,xrange=[xcut,xcut]))
v1yz = reform(loaddata(file,'velx',sample=sample,xrange=[xcut,xcut]))
v2yz = reform(loaddata(file,'vely',sample=sample,xrange=[xcut,xcut]))
v3yz = reform(loaddata(file,'velz',sample=sample,xrange=[xcut,xcut]))

dxz = reform(loaddata(file,'dens',sample=sample,yrange=[ycut,ycut],xCoord=x,yCoord=yy,zCoord=z,time=time,lref=lref))
lxz = reform(lref)
pxz = reform(loaddata(file,'pres',sample=sample,yrange=[ycut,ycut]))
jxz = reform(loaddata(file,'jet',sample=sample,yrange=[ycut,ycut]))
v1xz = reform(loaddata(file,'velx',sample=sample,yrange=[ycut,ycut]))
v2xz = reform(loaddata(file,'vely',sample=sample,yrange=[ycut,ycut]))
v3xz = reform(loaddata(file,'velz',sample=sample,yrange=[ycut,ycut]))


save,file='jet_snap_'+id+'.sav',dyz,lyz,pyz,jyz,v1yz,v2yz,v3yz, dxz,lxz,pxz,jxz,v1xz,v2xz,v3xz, x,y,z, time, id
endif else restore,file='jet_snap_'+id+'.sav'

zcut = 1.5e13
zcind = (where(z ge zcut))[0]
dyz = reform(dyz[*,0:zcind])
pyz = reform(pyz[*,0:zcind])
jyz = reform(jyz[*,0:zcind])
v1yz = reform(v1yz[*,0:zcind])
v2yz = reform(v2yz[*,0:zcind])
v3yz = reform(v3yz[*,0:zcind])
z = reform(z[0:zcind])

pltsz = size(dyz,/dimension)

vj = 3.d9
jvyz = jyz*abs(v3yz)/vj; > 0.95

ramp = dyz*(v1yz*v1yz + v2yz*v2yz)

pltx0=150 & plty0=100
loadct,0,/sil
winxs = pltx0+pltsz[0]*4. & winys = plty0+pltsz[1]+100.
swindow, xs=winxs, ys=winys
xyouts,winxs/2.+30., 50,/dev,'y [cm]'

mind=min(dyz) & maxd=max(dyz)
;minp=min(pyz) & maxp=max(pyz)
minp = 5. & maxp=1.e3
minjv=min(jvyz) & maxjv=max(jvyz)
tvcoord,bytscl(alog10(dyz),min=alog10(mind),max=alog10(maxd)),y,z,pos=[pltx0,plty0,pltx0+pltsz[0],plty0+pltsz[1]],/dev,/axes,ytitle='z [cm]'
legend,id,/left,/top,box=0
color_bar,lim=[mind,maxd],/log,/top,bartitle='density',titlegap=0.04
loadct,3
tvcoord,bytscl(jvyz,min=minjv,max=maxjv),y,z,pos=[pltx0+pltsz[0],plty0,pltx0+2*pltsz[0],plty0+pltsz[1]],/dev,/axes,ytickformat='(a1)'
color_bar,lim=[minjv,maxjv],/top,bartitle='jet*velz/jetv',titlegap=0.04
loadct,13
tvcoord,bytscl(alog10(pyz),min=alog10(minp),max=alog(maxp)),y,z,pos=[pltx0+2*pltsz[0],plty0,pltx0+3*pltsz[0],plty0+pltsz[1]],/dev,/axes,ytickformat='(a1)'
color_bar,lim=[minp,maxp],/log,/top,bartitle='pressure',titlegap=0.04
tvcoord,bytscl(alog10(ramp),min=alog10(minp),max=alog(maxp)),y,z,pos=[pltx0+3*pltsz[0],plty0,pltx0+4*pltsz[0],plty0+pltsz[1]],/dev,/axes,ytickformat='(a1)'
color_bar,lim=[minp,maxp],/log,/top,bartitle='ram pressure',titlegap=0.04

stop
end
