forward_function jet_anal
pro recoll_bend
id='3e37'
;id='1e38'
;id='3e37zoom'
;id='1e38zoom'

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

;den = loaddata(dir+fname,'dens',sample=sample,xra=xrange,yra=yrange,zra=zrange)
v3 = loaddata(dir+fname,'velz',sample=sample,xra=xrange,yra=yrange,zra=zrange)

zcut0 = 0.
zcind = (where(z ge zcut0))[0]

nz = n_elements(z)
z = reform(z[zcind:nz-1])
jet = reform(jet[*,*,zcind:nz-1])
v3 = reform(v3[*,*,zcind:nz-1])
;den = reform(den[*,*,zcind:nz-1])

sz = size(jet,/dimension)
velj = 3.d9
jv = jet*abs(v3)/velj

jcrit = 0.95
;;jv2 = jv>jcrit

jyl = fltarr(sz[2]) & jyr = fltarr(sz[2]) & jx = fltarr(sz[2]) ;& jxmax = fltarr(sz[2])

jythick_x = fltarr(sz[0])
for i=0,sz[2]-1 do begin
    jythick_x[*] = 0.
    for j=0,sz[0]-1 do begin
        jin_ind = where(jv[j,*,i] ge jcrit,count)
        if (count ne 0) then begin
           jyl_x = y[jin_ind[0]]
           jyr_x = y[jin_ind[n_elements(jin_ind)-1]]
           jythick_x[j] = jyr_x-jyl_x
        endif 
    endfor
; find jet x center by weighted y-width
    if (total(jythick_x) eq 0.) then jx_org=-2.e12 $
       else jx_org = total(x*jythick_x)/total(jythick_x)

    x_ind = (where(x ge jx_org))[0]
    jx[i] = x[x_ind]

    jin_ind = where(jv[x_ind,*,i] ge jcrit,count2)
    if (count2 ne 0) then begin
       jyl[i] = y[jin_ind[0]]
       jyr[i] = y[jin_ind[n_elements(jin_ind)-1]]
    endif else begin
       jyl[i] = 0.
       jyr[i] = 0.
    endelse
; jet x center by maximum jv
;    x_ind2 = (where(jv[*,sz[1]/2,i] eq max(jv[*,sz[1]/2,i])))[0]
;    jxmax[i] = x[x_ind2]
endfor

loadct,0,/sil
window,0,xs=sz[0],ys=sz[2]
tvcoord,reform(jv[*,sz[1]/2,*]),x,z,/scale
;tvcoord,alog10(reform(den[*,sz[1]/2,*])),x,z,/scale
;plots,jx_org,z,psym=1,color=fsc_color('yellow')
plots,jx,z,psym=1, color=fsc_color('red')
;plots,jxmax,z,psym=1, color=fsc_color('blue')


loadct,0,/sil
window,1
plot,jyl,z,xtitle='y [cm]',ytitle='z [cm]',/iso,xrange=[-3.e12,3.e12],yrange=[0.,1.e13],/xst,/yst
oplot,jyr,z

save,file='recoll_bend_'+id+'.sav',z,jyl,jyr,jx,time,jxmax

stop
end

pro recoll_comp
; compare between data cut by -2.e12 and data trace the jet

;id='3e37'
;id='1e38'

;resfname = 'recoll_bend_'+id+'.sav'
;restore,file=resfname
;id = strmid(resfname,12,4)
restore,file='recoll_bend_3e37.sav'
z1 = z & jyl1 = jyl & jyr1 = jyr & jx1=jx ;& jxmax1=jxmax
restore,file='recoll_bend_1e38.sav'
z2 = z & jyl2 = jyl & jyr2 = jyr & jx2=jx ;& jxmax1=jxmax

;if (id eq '1e38') then begin
dir2='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet38/'
fname2='JetSet_hdf5_plt_cnt_3032'
;endif else if (id eq '3e37') then begin
dir1='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet_3E37/'
fname1='JetSet_hdf5_plt_cnt_3024'
;endif

sample=4
xrange = [-5.e12,1.e12] &  zrange = [0.,1.e13]
d1 = reform(loaddata(dir1+fname1,'dens',sample=sample,xCoord=xx1,yCoord=yy1,zCoord=zz1,xrange=xrange,yrange=[0.,0.],zrange=zrange))
d2 = reform(loaddata(dir2+fname2,'dens',sample=sample,xCoord=xx2,yCoord=yy2,zCoord=zz2,xrange=xrange,yrange=[0.,0.],zrange=zrange))

szd = size(d1,/dimension)

jx_anal38 = jet_anal(Lj=1.d38,jyl=jyl1,jx=jx1,z=z1)
jx_anal37 = jet_anal(Lj=3.d37,jyl=jyl1,jx=jx1,z=z1)

loadct,0,/sil
pltx0=100. & plty0=60.
winxs=pltx0+szd[0]*2+30 & winys=plty0+szd[1]+90

mkeps,'recoll_bend_comp1',xs=20.,ys=20.*winys/winxs
;window,0, xs=winxs, ys=winys
;window,0,xs=500,ys=800
;maxd=1.e-12 & mind=1.e-18
maxd=max(d1) & mind=min(d1)
;tvcoord,bytscl(alog10(d),max=alog10(maxd),min=alog10(mind)),x,z,/scale,/axes,xtitle='x [cm]',ytitle='z [cm]',pos=[pltx0,plty0]
contour,alog10(d1),xx1,zz1,/fill,nlevels=255,xra=[-5.e12,1.e12],yra=[0.,1.e13],/iso,xtitle='x [cm]',ytitle='z [cm]',/xst,/yst $
       ,xtickinterval=2.e12,pos=posnorm([pltx0,plty0,pltx0+szd[0],plty0+szd[1]],nx=winxs,ny=winys),/norm
oplot,[-2.e12,-2.e12],!y.crange,line=2,color=fsc_color('yellow')
oplot,(2.*jx1-2.e12)/3.,z1,psym=1,color=fsc_color('yellow'),symsize=0.5
oplot,jx_anal37,z1,color=fsc_color('red'),thick=4.
legend,'L=3e37',/right,/top,box=0,textcolor=fsc_color('cyan')

loadct,0,/sil

contour,alog10(d2),xx2,zz2,/fill,nlevels=255,xra=[-5.e12,1.e12],yra=[0.,1.e13],/iso,xtitle='x [cm]',/xst,/yst $
       ,xtickinterval=2.e12,pos=posnorm([pltx0+szd[0],plty0,pltx0+szd[0]*2,plty0+szd[1]],nx=winxs,ny=winys),/norm,/noerase,ytickformat='(a1)'
oplot,[-2.e12,-2.e12],!y.crange,line=2,color=fsc_color('yellow')
oplot,(2.*jx2-2.e12)/3.,z1,psym=1,color=fsc_color('yellow'),symsize=0.5
oplot,jx_anal38,z2,color=fsc_color('red'),thick=4.
legend,'L=1e38',/right,/top,box=0,textcolor=fsc_color('cyan')

loadct,0,/sil
color_bar,lim=[mind,maxd],/log,/up,bartitle='density',titlegap=0.05,position=posnorm([pltx0,plty0+szd[1]+10.,pltx0+2*szd[0],plty0+szd[1]+30.],nx=winxs,ny=winys),/norm

epsfree

stop

end

pro recoll_comp2
device,decomposed=0

restore,'recoll_bend_3e37.sav'
z1=z & jyl1=jyl & jyr1=jyr
;restore,'recoll_bend_3e37zoom.sav'
;z11=z & jyl11=jyl & jyr11=jyr
restore,'recoll_bend_1e38.sav'
z2=z & jyl2=jyl & jyr2=jyr
;restore,'recoll_bend_1e38zoom.sav'
;z22=z & jyl22=jyl & jyr22=jyr

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

mkeps,'recoll_bend_comp2',xs=20.,ys=20.*winys/winxs

plot,jyr1,z1,xra=[-5.e12,5.e12],yra=[0.,ycut],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]',position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm,xtickinterval=4.e12
oplot,jyl1,z1
oplot,linyl1,z1[0:ycutind1],line=2,color=254
oplot,linyr1,z1[0:ycutind1],line=2,color=254
legend,'L=3e37',/right,/top,textcolor=0,box=0
legend,['identified jet','line to the end of jet'],color=[0,254],textcolor=[0,254],line=[0,2],box=0,position=[pltx0/winxs,0.98],/norm
plot,jyr2,z2,xra=[-5.e12,5.e12],yra=[0.,ycut],/xst,/yst,/iso,xtitle='y [cm]',position=posnorm([pltx0+pltxs,plty0,pltx0+2*pltxs,plty0+pltys],nx=winxs,ny=winys),/norm,/noerase,ytickformat='(a1)',xtickinterval=4.e12
oplot,jyl2,z2
oplot,linyl2,z2[0:ycutind2],line=2,color=254
oplot,linyr2,z2[0:ycutind2],line=2,color=254
legend,'L=1e38',/right,/top,textcolor=0,box=0

epsfree

mkeps,'recoll_bend_comp2_dev',xs=20,ys=20.*6./8.
dev1 = (linyl1*linyl1-jyl1[0:ycutind1]*jyl1[0:ycutind1]) / (linyl1*linyl1)
dev2 = (linyl2*linyl2-jyl2[0:ycutind2]*jyl2[0:ycutind2]) / (linyl2*linyl2)
plot,z1[0:ycutind1],dev1,xtitle='z [cm]',ytitle='deviation',/xst,/yst,yra=[-1.,1.]
oplot,z2[0:ycutind2],dev2,color=254
oplot,!x.crange,[0.,0.],line=2
legend,['L=3e37','L=1e38'],/right,/top,box=0,textcolor=[0,254],color=[0,254],line=0
epsfree
stop
end


function jet_anal,jyl=jyl,jx=jx,Lj=Lj,z=z

;vw0 = 2.14d8
vw0 = 2.5d8
;Mdw = 4.4d20
Mdw = 6.3d20

vj = 3.d9
;Lj = 1.d38

sep = 3.d12

;restore,'recoll_bend_1e38.sav'
z2=z & jyl2=jyl & jx2 = jx

zcut=max(z2)
zcutind2 = (where(z2 ge zcut))[0]

jyl2_zcut = jyl2[zcutind2]

;tanalp = abs(jyl2_zcut/zcut)
if (Lj eq 1.d38) then tanalp = 0.15 $
  else if (Lj eq 3.d37) then tanalp = abs(jyl2_zcut/zcut)
print,tanalp,atan(tanalp)/!dtor

x0=-2.d12
jx_anal = x0 - vw0*vj*Mdw/(4.*!pi*Lj) * tanalp * (z2-sep*atan(z2/sep))

;window,0
;plot,z2,x0-jx2,/iso
;oplot,z2,x0-jx_anal

return,jx_anal

end
