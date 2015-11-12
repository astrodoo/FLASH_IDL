forward_function jet_anal, Fth, PQ_Limits
pro recoll_bend,zoom=zoom
;id='1e37'
id='3e37'
;id='1e38'

if (id eq '1e38') then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet38/'
   fname='JetSet_hdf5_plt_cnt_3032'
endif else if (id eq '3e37') then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet_3E37/'
   fname='JetSet_hdf5_plt_cnt_3024'
endif else if (id eq '1e37') then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/M30lateral0.01/'
   fname='JetSet_hdf5_plt_cnt_0439'
endif

bhx = -2.e12

if ((id eq '3e37') or (id eq '1e38')) then begin
   sample=4
   dx=6.e12/2.
   dz=1.e13
   xrange = [bhx-dx,bhx+dx] & yrange = [-dx,dx] & zrange = [0.,dz]
endif else if (id eq '1e37') then begin
   sample=3
   dx=9.6e11/2.
   dz=1.6e12 
   xrange = [bhx-dx,bhx+dx] & yrange = [-dx,dx] & zrange = [0.,dz]
endif

if keyword_set(zoom) then begin
   dx=9.6e11/2.
   dz=1.6e12 
   xrange = [bhx-dx,bhx+dx] & yrange = [-dx,dx] & zrange = [0.,dz]
   sample=2
endif

jet = loaddata(dir+fname,'jet',sample=sample,xCoord=x,yCoord=y,zCoord=z,time=time,xra=xrange,yra=yrange, zra=zrange)

den = loaddata(dir+fname,'dens',sample=sample,xra=xrange,yra=yrange,zra=zrange)
v3 = loaddata(dir+fname,'velz',sample=sample,xra=xrange,yra=yrange,zra=zrange)

zcut0 = 0.
zcind = (where(z ge zcut0))[0]

ycut0 = 0.
ycind = (where(y ge ycut0))[0]

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
plot,jyl,z,xtitle='y [cm]',ytitle='z [cm]',/iso,xrange=[-dx,dx],yrange=[0.,dz],/xst,/yst
oplot,jyr,z

dxz = reform(den[*,ycind,*])
window,2
tvcoord,alog10(dxz),x,z,/scale,/axes

if keyword_set(zoom) then $
   save,file='recoll_bend_'+id+'_zoom.sav',z,x,dxz,jyl,jyr,jx,time $ ;,jxmax 
  else $
   save,file='recoll_bend_'+id+'.sav',z,x,dxz,jyl,jyr,jx,time ;,jxmax

stop
end

pro recoll_comp
; compare between data cut by -2.e12 and data trace the jet

restore,file='recoll_bend_3e37.sav'
z1 = z & x1=x & jyl1 = jyl & jyr1 = jyr & jx1=jx & dxz1=dxz ;& jxmax1=jxmax
restore,file='recoll_bend_1e38.sav'
z2 = z & x2=x & jyl2 = jyl & jyr2 = jyr & jx2=jx & dxz2=dxz ;& jxmax1=jxmax
restore,file='recoll_bend_1e37.sav'
z3 = z & x3=x & jyl3 = jyl & jyr3 = jyr & jx3=jx & dxz3=dxz ;& jxmax1=jxmax

szd = size(dxz1,/dimension)

jx_anal38 = jet_anal(Lj=1.d38,jyl=jyl1,jx=jx1,z=z1)
jx_anal37 = jet_anal(Lj=3.d37,jyl=jyl1,jx=jx1,z=z1)

loadct,0,/sil
pltx0=130. & plty0=70.
xgap=120.
winxs=pltx0+szd[0]*3+30+xgap & winys=plty0+szd[1]+100

mkeps,'recoll_bend_comp1',xs=30.,ys=30.*winys/winxs


bhx=-2.e12
dx=6.e12/2.
dz=1.e13
dx2=9.6e11/2.
dz2=1.6e12 
 
x3=x3-bhx
x1=x1-bhx
x2=x2-bhx

;modify
inlin = where((z3 ge 1.2e12) and (z3 le 1.4e12))
lin3 = linfit(z3[inlin],jx3[inlin])

jx3[where(z3 gt 1.4e12)] = lin3[1]*z3[where(z3 gt 1.4e12)]+lin3[0] 

maxd2=max(dxz3) & mind2=min(dxz3)
tvcoord,bytscl(alog10(dxz3),max=alog10(maxd2),min=alog10(mind2)),x3,z3,/scale,/axes,xtitle='x [cm]',ytitle='z [cm]' $
       ,pos=[pltx0/winxs,plty0/winys],/norm,psx=szd[0]/winxs,/black,xtickinterval=3.e11
;contour,alog10(dxz3),x3,z3,/fill,nlevels=255,xra=[bhx-dx2,bhx+dx2],yra=[0.,dz2],/iso,xtitle='x [cm]',ytitle='z [cm]',/xst,/yst $
;       ,xtickinterval=2.e12,pos=posnorm([pltx0,plty0,pltx0+szd[0],plty0+szd[1]],nx=winxs,ny=winys),/norm
oplot,[0.,0.],!y.crange,line=2,color=fsc_color('yellow')
oplot,jx3-bhx,z3,psym=1,color=fsc_color('yellow'),symsize=0.5
oplot,xjet_37-bhx,zjet_37,color=fsc_color('red'),thick=4.
;legend,'L=1e37',/right,/top,box=0,textcolor=fsc_color('cyan')
xyouts,(pltx0+szd[0]-70.)/winxs,(plty0+szd[1]-30.)/winys,/norm,'P1e37',color=fsc_color('cyan')
loadct,0,/sil
color_bar,lim=[mind2,maxd2],/log,/up,bartitle='density',titlegap=0.07,position=posnorm([pltx0,plty0+szd[1],pltx0+szd[0],plty0+szd[1]+20.],nx=winxs,ny=winys),/norm

loadct,0,/sil
maxd=max(dxz1) & mind=min(dxz1)
;tvcoord,bytscl(alog10(d),max=alog10(maxd),min=alog10(mind)),x,z,/scale,/axes,xtitle='x [cm]',ytitle='z [cm]',pos=[pltx0,plty0]
contour,alog10(dxz1),x1,z1,/fill,nlevels=255,/iso,xtitle='x [cm]',ytitle='z [cm]',/xst,/yst $
;contour,alog10(dxz1),x1,z1,/fill,nlevels=255,xra=[-dx,dx],yra=[0.,dz],/iso,xtitle='x [cm]',ytitle='z [cm]',/xst,/yst $
       ,xtickinterval=2.e12,pos=posnorm([pltx0+szd[0]+xgap,plty0,pltx0+2*szd[0]+xgap,plty0+szd[1]],nx=winxs,ny=winys),/norm,/noerase
oplot,[0.,0.],!y.crange,line=2,color=fsc_color('yellow')
oplot,(2.*jx1-2.e12)/3.-bhx,z1,psym=1,color=fsc_color('yellow'),symsize=0.5
oplot,jx_anal37-bhx,z1,color=fsc_color('red'),thick=4.
;legend,'L=3e37',/right,/top,box=0,textcolor=fsc_color('cyan')
xyouts,(pltx0+2*szd[0]+xgap-70.)/winxs,(plty0+szd[1]-30.)/winys,/norm,'P3e37',color=fsc_color('cyan')

loadct,0,/sil

contour,alog10(dxz2),x2,z2,/fill,nlevels=255,/iso,xtitle='x [cm]',/xst,/yst $
;contour,alog10(dxz2),x2,z2,/fill,nlevels=255,xra=[-dx,dx],yra=[0.,dz],/iso,xtitle='x [cm]',/xst,/yst $
       ,xtickinterval=2.e12,pos=posnorm([pltx0+2*szd[0]+xgap,plty0,pltx0+szd[0]*3+xgap,plty0+szd[1]],nx=winxs,ny=winys),/norm,/noerase,ytickformat='(a1)'
oplot,[0.,0.],!y.crange,line=2,color=fsc_color('yellow')
oplot,(2.*jx2-2.e12)/3.-bhx,z1,psym=1,color=fsc_color('yellow'),symsize=0.5
oplot,jx_anal38-bhx,z2,color=fsc_color('red'),thick=4.
;legend,'L=1e38',/right,/top,box=0,textcolor=fsc_color('cyan')
xyouts,(pltx0+3*szd[0]+xgap-70.)/winxs,(plty0+szd[1]-30.)/winys,/norm,'P1e38',color=fsc_color('cyan')

loadct,0,/sil
color_bar,lim=[mind,maxd],/log,/up,bartitle='density',titlegap=0.07,position=posnorm([pltx0+szd[0]+xgap,plty0+szd[1],pltx0+3*szd[0]+xgap,plty0+szd[1]+20.],nx=winxs,ny=winys),/norm

epsfree

stop

end

pro recoll_comp2
device,decomposed=0

;restore,'recoll_3e37.sav'
;z1=z & jyl1=jyl & jyr1=jyr
;restore,'recoll_1e38.sav'
;z2=z & jyl2=jyl & jyr2=jyr
restore,'recoll_1e37.sav'
z3=z & jyl3=jyl & jyr3=jyr

inlin = where((z3 ge 9e11) and (z3 le 1.1e12))
lin3 = linfit(z3[inlin], jyl3[inlin])
jyl3[where(z3 ge 1.1e12)] = lin3[1]*z3[where(z3 ge 1.1e12)] + lin3[0]


inlin2 = where(z3 lt 5.e11)
lin32 = linfit(z3[inlin2],jyl3[inlin2])
jyl3lin2 = lin32[1]*z3+lin32[0]

z3end = z3[n_elements(z3)-1]
jyl3end = jyl3[n_elements(z3)-1]
jyl3lin = jyl3end/z3end *z3 + z3[0]

bhx=-2.e12
dx=6.e12/2.
dz=1.e13
dx2=9.6e11/2.
dz2=1.6e12 
 
;z1=z1-z1[0] & jyl1=jyl1-jyl1[0] & jyr1=jyr1-jyr1[0]
;z2=z2-z2[0] & jyl2=jyl2-jyl2[0] & jyr2=jyr2-jyr2[0]
;loadct,39,/sil
;window,0
;
;plot,jyl3,z3,xra=[-dx2,dx2],yra=[0.,1.5e12],/iso,/xst,/yst,/nodata
;oplot,jyl3,z3
;oplot,-jyl3,z3
;oplot,jyl3lin,z3,line=2
;oplot,-jyl3lin,z3,line=2
;oplot,!x.crange,[7.e11,7.e11],line=1
;
;stop
;oplot,smooth(jyr1,20),z1
;oplot,smooth(jyl1,20),z1
;oplot,smooth(jyr2,20),z2
;oplot,smooth(jyl2,20),z2


restore,'recoll_bend_3e37.sav'
z1=z & jyl1=jyl & jyr1=jyr
restore,'recoll_bend_1e38.sav'
z2=z & jyl2=jyl & jyr2=jyr

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


;winyx = 6./8.
;pltx0 = 0.2 & plty0 = 0.12
;pltdx = 0.38 & pltdy = pltdx* ycut/1.e13/winyx
;
pltx0 = 180. & plty0=100.
pltxs = 400. & pltys=ycut/1.e13 *pltxs
xgap=40.
winxs = pltx0+3*pltxs+30.+xgap & winys = plty0+pltys+100.

;loadct,39,/sil
;mkeps,'recoll_bend_comp2',xs=30.,ys=30.*winys/winxs
;
;plot,jyl3,z3,xra=[-5.e11,5.e11],yra=[0.,1.5e12],/iso,/xst,/yst,/nodata,position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm $
;    ,xtickinterval=3.e11, ytickinterval=5.e11,xtitle='y [cm]', ytitle='z [cm]'
;oplot,jyl3,z3
;oplot,-jyl3,z3
;oplot,jyl3lin,z3,line=2,color=254
;oplot,-jyl3lin,z3,line=2,color=254
;oplot,jyl3lin2,z3,line=2,color=50
;oplot,-jyl3lin2,z3,line=2,color=50
;oplot,!x.crange,[7.e11,7.e11],line=1
;xyouts,(pltx0+10.)/winxs,(plty0+pltys-40)/winys,/norm,'P1e37'
;
;plot,jyr1,z1,xra=[-5.e12,5.e12],yra=[0.,ycut],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]',position=posnorm([pltx0+pltxs+xgap,plty0,pltx0+2*pltxs+xgap,plty0+pltys],nx=winxs,ny=winys),/norm,xtickinterval=4.e12,/noerase
;oplot,jyl1,z1
;oplot,linyl1,z1[0:ycutind1],line=2,color=254
;oplot,linyr1,z1[0:ycutind1],line=2,color=254
;;legend,'L=3e37',/right,/top,textcolor=0,box=0
;xyouts,(pltx0+pltxs+xgap+20.)/winxs,(plty0+pltys-40)/winys,/norm,'P3e37'
;
;legend,['identified jet','line to the end of jet'],color=[0,254],textcolor=[0,254],line=[0,2],box=0,position=[(pltx0-30)/winxs,0.98],/norm
;
;plot,jyr2,z2,xra=[-5.e12,5.e12],yra=[0.,ycut],/xst,/yst,/iso,xtitle='y [cm]',position=posnorm([pltx0+pltxs*2+xgap,plty0,pltx0+3*pltxs+xgap,plty0+pltys],nx=winxs,ny=winys),/norm,/noerase,ytickformat='(a1)',xtickinterval=4.e12
;oplot,jyl2,z2
;oplot,linyl2,z2[0:ycutind2],line=2,color=254
;oplot,linyr2,z2[0:ycutind2],line=2,color=254
;;legend,'L=1e38',/right,/top,textcolor=0,box=0
;xyouts,(pltx0+2*pltxs+xgap+20.)/winxs,(plty0+pltys-40)/winys,/norm,'P1e38'
;
;epsfree

mkeps,'recoll_bend_comp2_dev',xs=20,ys=20.*6./8.
dev1 = (linyl1*linyl1-jyl1[0:ycutind1]*jyl1[0:ycutind1]) / (linyl1*linyl1)
dev1[where(z1 le 5.e11)] = dev1[where(z1 le 5.e11)]*0.1

dev2 = (linyl2*linyl2-jyl2[0:ycutind2]*jyl2[0:ycutind2]) / (linyl2*linyl2)
dev2[where(z1 le 5.e11)] = dev2[where(z1 le 5.e11)]*0.1

plot,z1[0:ycutind1],smooth(dev1,30),xtitle='z [cm]',ytitle='fractional deviation',/xst,/yst,yra=[-1.,1.]
oplot,z2[0:ycutind2],smooth(dev2,30),color=254
oplot,!x.crange,[0.,0.],line=2
legend,['P3e37','P1e38'],/right,/top,box=0,textcolor=[0,254],color=[0,254],line=0
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


; calculating analytic jet trajectory for pressure equilibrium model
pro analjet37

restore,file='recoll_bend_1e37.sav'
z1 = z & x1=x & jyl1 = jyl & jyr1 = jyr & jx1=jx & dxz1=dxz ;& jxmax1=jxmax

jetpos,th0=0,Lj=1e37,xjet=xjet_37,yjet=zjet_37

window,0
tvcoord,alog10(dxz1),x,z,/scale,/axes
oplot,xjet_37,zjet_37,psym=1

save_append,file='recoll_bend_1e37.sav',xjet_37,zjet_37
stop
end

; getting analytical lines
function Fth,th1,th2
common Power, intpow
   return, (1.+tan(th1)^2.)*cos(th2)^intpow
end

function PQ_Limits,th1
common Param, th_start
;print, 'th_start=', th_start
   return, [th_start,th1]
end

pro jetpos,th0=th0,xjet=xjet,yjet=yjet,nth=nth,theta=theta, Lj=Lj, gam=gam
common Param, th_start
common Power, intpow

if keyword_set(gam) then intpow = 0.8 else intpow = 0.5
print,'Lj=  ',Lj,'      th0 = ', th0 
if keyword_set(gam) then print,'gamma = 5/3'

if not keyword_set(nth) then nth=10000

case Lj of
   1e35: h0=3.e10
   1e36: h0=5.5e10
   1e37: if keyword_set(gam) then h0=1.e11 else h0=1.2e11
   else: print, 'out of range in Energy'
endcase

;h0=5e10 ; L35
;h0=6.e10 ; L36
;h0=1.2e11 ; L37

;Lj = 1.d35
;Lj = 1.d36
;Lj = 1.d37

l = 3.d12
Mdot = 4.415d20
vw = 2.5d8
vj = 3.d9
Mj = 30.

th0 = -th0*!dtor
th_start=th0
Const = 0.25/!pi * Mdot*vw * h0 * vj * (1.+9./2./Mj^2) / Lj * cos(th0)^(-0.25)

th = findgen(nth)/float(nth)*(!pi/2.-th0) + th0
theta=th

int2val = fltarr(nth)
for i=0, nth-1 do begin
  AB_Limits = [th0, th[i]]
  int2val[i] = int_2d('Fth', AB_Limits, 'PQ_Limits', 96, /double)
endfor

;print,'int2val:', int2val[nth-1]

dx = Const * int2val

stx0 = 1.e12 & sty0 = 0.
jetx0= stx0-l*cos(abs(th0))  & jety0=l*sin(th0)

jetx = jetx0 - dx
jety = l*cos(th0)*tan(th)

jetx_rot =  cos(abs(th0))*(jetx-stx0) + sin(abs(th0))*(jety-sty0) + stx0
jety_rot = -sin(abs(th0))*(jetx-stx0) + cos(abs(th0))*(jety-sty0) + sty0
xjet=jetx_rot & yjet=jety_rot

;case Lj of
;   1e35: restore,file='/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e35/jettrace2_1e35_700_smp3.sav' 
;   1e36: restore,file='/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/_WRO_sphere_1e36_high/jettrace2_1e36_851_smp3.sav' 
;   1e37: restore,file='/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/_WRO_sphere_1e37_old/jettrace2_1e37_400_smp4.sav' 
;   else: print, 'out of range in Energy'
;endcase
;
;window,4
;plot,xjet,yjet,xra=[-6.e12,-1.e12],/xst,/iso,yra=[0.,8e12],/yst
;oplot,x2,z2
end



pro showpres

;id='1e37'
;id='3e37'
id='1e38'

if (id eq '1e38') then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet38/'
   fname='JetSet_hdf5_plt_cnt_3032'
endif else if (id eq '3e37') then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet_3E37/'
   fname='JetSet_hdf5_plt_cnt_3024'
endif else if (id eq '1e37') then begin
   dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/M30lateral0.01/'
   fname='JetSet_hdf5_plt_cnt_0439'
endif

bhx = -2.e12

if ((id eq '3e37') or (id eq '1e38')) then begin
   sample=2
   dx=6.e12/2.
   dz=1.e13
   xrange = [bhx-dx,bhx+dx] & yrange = [0.,0.] & zrange = [0.,dz]
endif else if (id eq '1e37') then begin
   sample=1
   dx=9.6e11/2.
   dz=1.6e12 
   xrange = [bhx-dx,bhx+dx] & yrange = [0.,0.] & zrange = [0.,dz]
endif

pres = reform(loaddata(dir+fname,'pres',sample=sample,xra=xrange,yra=yrange,zra=zrange,xCoord=x,yCoord=y,zCoord=z,time=time))

sz = size(pres,/dimension)

loadct,39,/sil
!p.charsize=3.
pltx0=200. & plty0=100
swindow,xs=pltx0+sz[0]+200,ys=plty0+sz[1]+50
tvcoord,alog10(pres),x-bhx,z,/scale,/axes,xtitle='x [cm]',ytitle='z [cm]',position=[pltx0,plty0],/dev
legend,'P'+id,/left,/top,textcolor=fsc_color('magenta'),box=0

loadct,39,/sil
color_bar,lim=[min(pres),max(pres)],/log,/right,bartitle='pressure'

stop
end


pro playpres

id='1e37'

dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/M30lateral0.01/'
fnames = file_search(dir+'*plt_cnt*')

nn = fix(strmid(fnames,3,4,/rev))

startn = 250

fnames_id = where(nn ge startn)  

fnames = fnames[fnames_id]
nfiles = n_elements(fnames)

;   fname='JetSet_hdf5_plt_cnt_0439'
;   fname='JetSet_hdf5_plt_cnt_0359'

bhx = -2.e12

if ((id eq '3e37') or (id eq '1e38')) then begin
   sample=2
   dx=6.e12/2.
   dz=1.e13
   xrange = [bhx-dx,bhx+dx] & yrange = [0.,0.] & zrange = [0.,dz]
endif else if (id eq '1e37') then begin
;   sample=3
   dx=9.6e11/2.
   dz=1.6e12 
   xrange = [bhx-dx,bhx+dx] & yrange = [0.,0.] & zrange = [0.,dz]
endif


spawn,'mkdir png_playpres'
outdir='png_playpres/'

for i=0, nfiles-1 do begin

print,i,' of ', nfiles

read_amr,fnames[i],/nodata,param=param,tree=tree
sample1 = max([max(tree.lrefine)-7,0])
sample2 = max([sample1-1,0])

preszoom = reform(loaddata(fnames[i],'pres',sample=sample2,xra=xrange,yra=yrange,zra=zrange,xCoord=xx,yCoord=yy,zCoord=zz,time=time))
pres = reform(loaddata(fnames[i],'pres',sample=sample1,yra=[0.,0.],xCoord=x,yCoord=y,zCoord=z,time=time))
sz = size(pres,/dimension)

zcut_ind = (where(z ge 0))[0] 
nz2 = (sz[1]-zcut_ind)*2
z2 = fltarr(nz2) & pres2 = fltarr(sz[0],nz2)
z2[0:nz2/2-1] = -reverse(z[zcut_ind:*]) & z2[nz2/2:*] = z[zcut_ind:*]
pres2[*,0:nz2/2-1] = reverse(pres[*,zcut_ind:*],2) & pres2[*,nz2/2:*] = pres[*,zcut_ind:*]


szz = size(preszoom,/dimension)

;minp = min(pres) & maxp = max(pres)
minp = 0.754456 & maxp = 234738.

loadct,39,/sil
!p.charsize=2.
pltx0=150. & plty0=100
window,xs=pltx0+sz[0]+120,ys=plty0+sz[1]+50

tvcoord,bytscl(alog10(pres2),min=alog10(minp),max=alog10(maxp)),x-bhx,z2,/axes,xtitle='x [cm]', ytitle='z [cm]',position=[pltx0,plty0],/dev
xyouts,pltx0-130,plty0+sz[1],/dev,'P'+id,color=255
tj = 29910.
xyouts,pltx0-130,plty0+sz[1]-20,/dev,'t='+string((time-tj)/60./60.,format='(f5.2)')+' hrs',color=255

loadct,39,/sil
color_bar,lim=[minp,maxp],/log,/right,bartitle='pressure'

if (time ge tj) then $
   tvcoord,bytscl(alog10(preszoom),min=alog10(minp),max=alog10(maxp)),xx-bhx,zz,/axes,position=[pltx0+sz[0]-szz[0],plty0+sz[1]-szz[1]],/dev,xtickinterval=4.e11

snapshot,outdir+'playpres_P'+id+'_'+string(i,format=('(I3.3)'))
endfor

stop
end

