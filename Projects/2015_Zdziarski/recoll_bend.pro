pro recoll_bend
device,decomposed=0

;dir= '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L6E37_M30_lw_lg/'
;fname = 'JetSet_hdf5_plt_cnt_2963'
dir= '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L38_M30_lw_lg/'
fname = 'JetSet_hdf5_plt_cnt_3036'
dir= '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/sphere_1e37/'
fname = 'JetSet_hdf5_plt_cnt_0463'
sample=2

;jet = loaddata(dir+fname,'jet',sample=sample,xCoord=x,yCoord=y,zCoord=z,time=time)
;v3 = loaddata(dir+fname,'velz',sample=sample)

xrange = [-7.e12,1.e12] & yrange = [-6.e12,6.e12] & zrange = [0.,1.5e13]
jet = loaddata(dir+fname,'jet',sample=sample,xCoord=x,yCoord=y,zCoord=z,time=time,xra=xrange,yra=yrange, zra=zrange)
v3 = loaddata(dir+fname,'velz',sample=sample,xra=xrange,yra=yrange,zra=zrange)

zcut0 = 0.
zcind = (where(z ge zcut0))[0]

nz = n_elements(z)
z = reform(z[zcind:nz-1])
jet = reform(jet[*,*,zcind:nz-1])
v3 = reform(v3[*,*,zcind:nz-1])

sz = size(jet,/dimension)
velj = 3.d9
jv = jet*abs(v3)/velj
jcrit = 0.95
;jv2 = jv>jcrit

jyl = fltarr(sz[2]) & jyr = fltarr(sz[2]) & jx = fltarr(sz[2])
for i=0,sz[2]-1 do begin
    iflag = 1
    jythick = 0.
    for j=0,sz[0]-1 do begin
        jin_ind = where(jv[j,*,i] ge jcrit,count)
        if (count ne 0) then begin
           jyl_x = y[jin_ind[0]]
           jyr_x = y[jin_ind[n_elements(jin_ind)-1]]
           jythick_x = jyr_x-jyl_x
           if (jythick_x gt jythick) then begin
              jythick = jythick_x
              jx[i] = x[j]
              jyl[i] = jyl_x
              jyr[i] = jyr_x
              iflag = 0
           endif
        endif

        if (iflag) then begin
           jx[i] = -2.e12
           jyl[i]=0.
           jyr[i]=0.
        endif
    endfor
endfor

;save,file='recoll_bend_6e37.sav',z,jyl,jyr,jx,time
save,file='recoll_bend_1e38.sav',z,jyl,jyr,jx,time

loadct,0,/sil
window,0
plot,jyl,z,xtitle='y [cm]',ytitle='z [cm]',/iso,xrange=[-3.e12,3.e12],yrange=[0.,1.5e13],/xst,/yst
oplot,jyr,z

end

pro recoll_comp
; compare between data cut by -2.e12 and data trace the jet

;resfname = 'recoll_bend_6e37.sav'
resfname = 'recoll_bend_1e38.sav'
restore,file=resfname
id = strmid(resfname,12,4)
z1 = z & jyl1 = jyl & jyr1 = jyr & jx1=jx

;dir= '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L6E37_M30_lw_lg/'
;fname = 'JetSet_hdf5_plt_cnt_2963'
dir= '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L38_M30_lw_lg/'
fname = 'JetSet_hdf5_plt_cnt_3036'

d = reform(loaddata(dir+fname,'dens',sample=4,xCoord=x,yCoord=y,zCoord=z,yrange=[0.,0.]))

szd = size(d,/dimension)

loadct,0,/sil
mkeps,'recoll_bend_comp1_'+id,xs=20.,ys=20.
;window,0, xs=szd[0]+100, ys=szd[1]+100
;tvcoord,alog(d),x,z,/scale,/axes
contour,alog10(d),x,z,/cell_fill,nlevels=255,xra=[-5.e12,1.e12],yra=[0.,1.e13],/iso,xtitle='x [cm]',ytitle='z [cm]'
oplot,[-2.e12,-2.e12],!y.crange,line=2
oplot,jx1,z1,psym=1

;legend,'L=6E37',/right,/top,box=0,textcolor=0
legend,'L=1E38',/right,/top,box=0,textcolor=0

epsfree

mkeps,'recoll_bend_comp2_'+id,xs=20.,ys=30.
restore,file='recoll_6e37.sav'
z2 = z & jyl2 = jyl & jyr2 = jyr

;window,1,xs=500,ys=800
plot,jyl1,z1,xra=[-5.e12,5.e12],yra=[0.,1.5e13],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]'
oplot,jyr1,z1
oplot,jyl2,z2,line=2
oplot,jyr2,z2,line=2

legend,['trace jet','cut by x=-2.e12 cm'],line=[0,2],color=0,textcolor=0,/left,/bottom,box=0
;legend,'L=6E37',/right,/top,box=0,textcolor=0
legend,'L=1E38',/right,/top,box=0,textcolor=0
epsfree

stop

end

pro recoll_comp2
device,decomposed=0

restore,'recoll_bend_6e37.sav'
z1=z & jyl1=jyl & jyr1=jyr
restore,'recoll_bend_1e38.sav'
z2=z & jyl2=jyl & jyr2=jyr

ycut=1.5e13

ycutind1 = (where(z1 ge ycut))[0]
ycutind2 = (where(z2 ge ycut))[0]

jyl1_ycut = jyl1[ycutind1]
jyr1_ycut = jyr1[ycutind1]
jyl2_ycut = jyl2[ycutind2]
jyr2_ycut = jyr2[ycutind2]

linyl1=jyl1_ycut/ycut*z1[0:ycutind1]
linyr1=jyr1_ycut/ycut*z1[0:ycutind1]

linyl2=jyl2_ycut/ycut*z2[0:ycutind2]
linyr2=jyr2_ycut/ycut*z2[0:ycutind2]

loadct,39,/sil

winyx = 6./8.
pltx0 = 0.2 & plty0 = 0.12
pltdx = 0.38 & pltdy = pltdx* ycut/1.e13/winyx
mkeps,'recoll_bend_comp2',xs=20,ys=20.*winyx

plot,jyr1,z1,xra=[-5.e12,5.e12],yra=[0.,ycut],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]',position=[pltx0,plty0,pltx0+pltdx,plty0+pltdy],/norm,xtickinterval=4.e12
oplot,jyl1,z1
oplot,linyl1,z1[0:ycutind1],line=2,color=254
oplot,linyr1,z1[0:ycutind1],line=2,color=254
legend,'L=6e37',/right,/top,textcolor=0,box=0
legend,['jet','line to edge of jet'],color=[0,254],textcolor=[0,254],line=[0,2],box=0,position=[pltx0,0.98],/norm
plot,jyr2,z2,xra=[-5.e12,5.e12],yra=[0.,ycut],/xst,/yst,/iso,xtitle='y [cm]',position=[pltx0+pltdx,plty0,pltx0+2*pltdx,plty0+pltdy],/norm,/noerase,ytickformat='(a1)',xtickinterval=4.e12
oplot,jyl2,z2
oplot,linyl2,z2[0:ycutind2],line=2,color=254
oplot,linyr2,z2[0:ycutind2],line=2,color=254
legend,'L=1e38',/right,/top,textcolor=0,box=0

epsfree

mkeps,'recoll_bend_comp2_dev',xs=20,ys=20.*6./8.
dev1 = (linyl1*linyl1-jyl1[0:ycutind1]*jyl1[0:ycutind1]) / (linyl1*linyl1)
dev2 = (linyl2*linyl2-jyl2[0:ycutind2]*jyl2[0:ycutind2]) / (linyl2*linyl2)
plot,z1[0:ycutind1],dev1,xtitle='z [cm]',ytitle='deviation',/xst,/yst,xtickinterval=5.e12
oplot,z2[0:ycutind2],dev2,color=254
legend,['L=6e37','L=1e38'],/right,/top,box=0,textcolor=[0,254],color=[0,254],line=0
epsfree
stop
end
