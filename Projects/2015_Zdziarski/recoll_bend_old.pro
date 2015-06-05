pro recoll_bend
;id='3e37'
;id='1e38'
;id='3e37zoom'
id='1e38zoom'

;dir= '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L6E37_M30_lw_lg/'
;fname = 'JetSet_hdf5_plt_cnt_2963'
;dir= '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L38_M30_lw_lg/'
;fname = 'JetSet_hdf5_plt_cnt_3036'
;dir= '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/sphere_1e37/'
;fname = 'JetSet_hdf5_plt_cnt_0463'
if ((id eq '1e38') or (id eq '1e38zoom')) then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet38/'
   fname='JetSet_hdf5_plt_cnt_3032'
endif else if ((id eq '3e37') or (id eq '3e37zoom')) then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet_3E37/'
   fname='JetSet_hdf5_plt_cnt_3024'
endif

if ((id eq '3e37') or (id eq '1e38')) then begin
   sample=4
   xrange = [-7.e12,1.e12] & yrange = [-6.e12,6.e12] & zrange = [0.,1.e13]
endif else if ((id eq '3e37zoom') or (id eq '1e38zoom')) then begin
  sample=3
  xrange = [-3.5e12,5.e11] & yrange = [-3.e12,3.e12] & zrange = [0.,5.e12]
endif
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
;;jv2 = jv>jcrit

jyl = fltarr(sz[2]) & jyr = fltarr(sz[2]) & jx_org = fltarr(sz[2]) & jx = fltarr(sz[2])
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
              jx_org[i] = x[j]
;              jyl[i] = jyl_x
;              jyr[i] = jyr_x
              iflag = 0
           endif
        endif

        if (iflag) then begin
           jx_org[i] = -2.e12
;           jyl[i]=0.
;           jyr[i]=0.
        endif
    endfor
endfor

;for i=0,sz[2]-1 do begin
;   x_ind = where(jv[*,sz[1]/2,i] eq max(jv[*,sz[1]/2,i]),count)
;   if (count gt 1) then x_ind = x_ind[0]
;   jx[i] = x[x_ind]
;
;   jin_ind = where(jv[x_ind,*,i] ge jcrit,count2)
;   if (count2 ne 0) then begin
;      jyl[i] = y[jin_ind[0]]
;      jyr[i] = y[jin_ind[n_elements(jin_ind)-1]]
;   endif else begin
;      jyl[i] = 0.
;      jyr[i] = 0.
;   endelse
;endfor

fit = poly_fit(z,jx_org,2)
jxfit = fit[0] + fit[1]*z + fit[2]*z*z
plots,jxfit,z,psym=1, color=fsc_color('red')

for i=0,sz[2]-1 do begin
   x_ind = (where(x ge jxfit[i]))[0]
   jx[i] = x[x_ind]

   jin_ind = where(jv[x_ind,*,i] ge jcrit,count2)
   if (count2 ne 0) then begin
      jyl[i] = y[jin_ind[0]]
      jyr[i] = y[jin_ind[n_elements(jin_ind)-1]]
   endif else begin
      jyl[i] = 0.
      jyr[i] = 0.
   endelse
endfor

loadct,0,/sil
window,0,xs=sz[0],ys=sz[2]
tvcoord,reform(jv[*,sz[1]/2,*]),x,z,/scale
plots,jx_org,z,psym=1,color=fsc_color('yellow')
plots,jx,z,psym=1, color=fsc_color('red')

loadct,0,/sil
window,1
plot,jyl,z,xtitle='y [cm]',ytitle='z [cm]',/iso,xrange=[-3.e12,3.e12],yrange=[0.,1.e13],/xst,/yst
oplot,jyr,z

save,file='recoll_bend_'+id+'.sav',z,jyl,jyr,jx,time,jx_org

stop
end

pro recoll_comp
; compare between data cut by -2.e12 and data trace the jet

;id='3e37'
id='1e38'

resfname = 'recoll_bend_'+id+'.sav'
restore,file=resfname
id = strmid(resfname,12,4)
z1 = z & jyl1 = jyl & jyr1 = jyr & jx1=jx & jx_org1=jx_org

if (id eq '1e38') then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet38/'
   fname='JetSet_hdf5_plt_cnt_3032'
endif else if (id eq '3e37') then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet_3E37/'
   fname='JetSet_hdf5_plt_cnt_3024'
endif

sample=4
xrange = [-5.e12,1.e12] &  zrange = [0.,1.e13]
d = reform(loaddata(dir+fname,'dens',sample=sample,xCoord=x,yCoord=y,zCoord=z,xrange=xrange,yrange=[0.,0.],zrange=zrange))

szd = size(d,/dimension)

loadct,0,/sil
;mkeps,'recoll_bend_comp1_'+id,xs=20.,ys=20.
pltx0=100 & plty0=100
;winxs=pltx0+szd[0]+20 & winys=plty0+szd[1]+150
;window,0, xs=winxs, ys=winys
window,0,xs=500,ys=800
;maxd=1.e-12 & mind=1.e-18
maxd=max(d) & mind=min(d)
;tvcoord,bytscl(alog10(d),max=alog10(maxd),min=alog10(mind)),x,z,/scale,/axes,xtitle='x [cm]',ytitle='z [cm]',pos=[pltx0,plty0]
contour,alog10(d),x,z,/fill,nlevels=255,xra=[-5.e12,1.e12],yra=[0.,1.e13],/iso,xtitle='x [cm]',ytitle='z [cm]',/xst,/yst $
       ,xtickinterval=2.e12
oplot,[-2.e12,-2.e12],!y.crange,line=2,color=fsc_color('yellow')
oplot,jx_org1,z1,psym=1,color=fsc_color('yellow')
oplot,jx1,z1,color=fsc_color('red')

legend,'L='+id,/right,/top,box=0,textcolor=fsc_color('cyan')

loadct,0,/sil
color_bar,lim=[mind,maxd],/log,/up,bartitle='density',titlegap=0.05,bargap=0.02

;epsfree

stop

end

pro recoll_comp2
device,decomposed=0

restore,'recoll_bend_3e37.sav'
z1=z & jyl1=jyl & jyr1=jyr
restore,'recoll_bend_3e37zoom.sav'
z11=z & jyl11=jyl & jyr11=jyr
restore,'recoll_bend_1e38.sav'
z2=z & jyl2=jyl & jyr2=jyr
restore,'recoll_bend_1e38zoom.sav'
z22=z & jyl22=jyl & jyr22=jyr

;ycut=1.5e13
;ycut=1.e13
ycut=max(z1)

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

;winyx = 6./8.
;pltx0 = 0.2 & plty0 = 0.12
;pltdx = 0.38 & pltdy = pltdx* ycut/1.e13/winyx
;
pltx0 = 180. & plty0=100.
pltxs = 400. & pltys=ycut/1.e13 *pltxs
winxs = pltx0+2*pltxs+30. & winys = plty0+pltys+100.

mkeps,'recoll_bend_comp2_zoom',xs=20.,ys=20.*winys/winxs

;plot,jyr1,z1,xra=[-5.e12,5.e12],yra=[0.,ycut],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]',position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm,xtickinterval=4.e12
;oplot,jyl1,z1
plot,jyr11,z11,xra=[-1.e12,1.e12],yra=[0.,5.e12],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]',position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm,xtickinterval=4.e12
oplot,jyl11,z11
oplot,linyl1,z1[0:ycutind1],line=2,color=254
oplot,linyr1,z1[0:ycutind1],line=2,color=254
legend,'L=3e37',/right,/top,textcolor=0,box=0
legend,['jet','line to edge of jet'],color=[0,254],textcolor=[0,254],line=[0,2],box=0,position=[pltx0/winxs,0.98],/norm
;plot,jyr2,z2,xra=[-5.e12,5.e12],yra=[0.,ycut],/xst,/yst,/iso,xtitle='y [cm]',position=posnorm([pltx0+pltxs,plty0,pltx0+2*pltxs,plty0+pltys],nx=winxs,ny=winys),/norm,/noerase,ytickformat='(a1)',xtickinterval=4.e12
;oplot,jyl2,z2
plot,jyr22,z22,xra=[-1.e12,1.e12],yra=[0.,5.e12],/xst,/yst,/iso,xtitle='y [cm]',position=posnorm([pltx0+pltxs,plty0,pltx0+2*pltxs,plty0+pltys],nx=winxs,ny=winys),/norm,/noerase,ytickformat='(a1)',xtickinterval=4.e12
oplot,jyl22,z22
oplot,linyl2,z2[0:ycutind2],line=2,color=254
oplot,linyr2,z2[0:ycutind2],line=2,color=254
legend,'L=1e38',/right,/top,textcolor=0,box=0

epsfree

mkeps,'recoll_bend_comp2_dev_zoom',xs=20,ys=20.*6./8.
dev1 = (linyl1*linyl1-jyl1[0:ycutind1]*jyl1[0:ycutind1]) / (linyl1*linyl1)
dev2 = (linyl2*linyl2-jyl2[0:ycutind2]*jyl2[0:ycutind2]) / (linyl2*linyl2)
plot,z1[0:ycutind1],dev1,xtitle='z [cm]',ytitle='deviation',/xst,/yst,xtickinterval=5.e12,yra=[-1.,1.]
oplot,z2[0:ycutind2],dev2,color=254
legend,['L=3e37','L=1e38'],/right,/top,box=0,textcolor=[0,254],color=[0,254],line=0
epsfree
stop
end
