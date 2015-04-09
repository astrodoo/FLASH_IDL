pro jettrace_rot, n, sample=sample

device,decomposed=0

xc0 = -2.e12 
yc0 = 0. 
zc0 = 0. 

;n=534 ;589
if (n_elements(sample) eq 0) then sample=2

;velz = dload(700,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(490,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(891,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(1000,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
jet = dload(n,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
velz = dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(520,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)

;if keyword_set(rotback) then begin
;   print, 'start rotating'
   m1 = 4.d34
   m2 = 2.d34
   mtot = m1+m2
   peri = 3.d12
   dth = sqrt(!unit.g * mtot / peri^3.d0) * time
  
   stx = peri*m2/mtot *cos(dth)
   sty = peri*m2/mtot *sin(dth)
;
;   dth = dth / !dtor
;   jet_rot = rot_3d(jet,rotaxis=3,degree=dth,/interp) 
;   jet=jet_rot
;   print,'complete rotating'
;endif


;y0ind_tmp = where(y ge 0) & y0ind = y0ind_tmp[0]
;velz2 = reform(velz[*,y0ind,*])
;jet2 = reform(jet[*,y0ind,*])

z0ind_tmp = where(z ge 0) & z0ind = z0ind_tmp[0]
;sz = size(jet2,/dimension)
jetpos = replicate({x:0.,y:0.,z:0.},n_elements(z)-z0ind)
;z2 = z[z0ind:*]
jetpos.z = z[z0ind:*]

critv=0.1
critj=0.9

for i=z0ind,n_elements(z)-1 do begin
;    maxv_xind = where(velz2[*,i] eq max(velz2[*,i]),count)
    maxj_ind = where( ((jet[*,*,i] ge critj) and (velz[*,*,i] ge 3.e9*critv)), count)
    if (count ge 2) then begin
        index = array_indices(jet,maxj_ind)
        dst2 = (index[0,*]-stx)^2. + (index[1,*]-sty)^2.
        shortdst_ind = where(dst2 eq min(dst2))

;        jetpos[i-z0ind].x = x[index[0,shortdst_ind]]
;        jetpos[i-z0ind].y = y[index[1,shortdst_ind]]
        maxj_ind = reform(maxj_ind[shortdst_ind])
    endif
    if ((count eq 0) or (maxj_ind eq 0)) then begin
       jetpos[i-z0ind].x = !values.f_nan
       jetpos[i-z0ind].y = !values.f_nan
    endif else begin
       index = array_indices(jet,maxj_ind)
       jetpos[i-z0ind].x = x[index[0]]
       jetpos[i-z0ind].y = y[index[1]]
    endelse
endfor

jetposxyz = fltarr(3,n_elements(jetpos))
jetposxyz[0,*] = jetpos.x
jetposxyz[1,*] = jetpos.y
jetposxyz[2,*] = jetpos.z
jetposcyl = cv_coord(from_rect=jetposxyz,/to_cylin)

!p.background=255 & !p.color=0
loadct,0,/sil
window,0
plot_3dbox,jetpos.x,jetpos.y,jetpos.z,/xy_plane,/xz_plane,/yz_plane,psym=1,xtitle='x',ytitle='y',ztitle='z'

loadct,39,/sil
window,1,xs=800,ys=300
plot,jetpos.x,jetpos.y,psym=1,/iso,xtitle='x',ytitle='y',/nodata
plots,jetpos.x,jetpos.y,psym=1,/data,color=bytscl(jetposxyz[2,*])
color_bar,lim=[min(jetposxyz[2,*]),max(jetposxyz[2,*])],/up,bartitle='z',titlegap=0.05

window,2,xs=500,ys=600
plot,jetposcyl[1,*],jetposcyl[2,*],/iso,psym=1,xtitle='r',ytitle='z'

save,filename='jettrace_1e36_rot_cyl_'+strtrim(n,2)+'.sav',jetposxyz,jetposcyl

stop
end

;---------------------------------------------------------------------------------------
pro jettrace2,ps=ps

fname1 = '/d/d3/yoon/outputs/out_mHD_Binary_sphere/sphere_1e35/jettrace_1e35_700.sav'
fname2 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_high/jettrace_1e36_high_891.sav'
fname3 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37/jettrace_1e37_490.sav'
fname4 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_rot/Stampede/jettrace_1e36_rot_570.sav'
fname5 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37_rot/Stampede/jettrace_1e37_rot_520.sav'
;fname6 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_rot/Stampede/jettrace_1e36_rot_cyl_570.sav'
fname6 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_rot/Brian/jettrace_1e36_rot_cyl_589.sav'
restore,filename=fname1 & x2_35 = x2 & z2_35 = z2
restore,filename=fname2 & x2_36 = x2 & z2_36 = z2
restore,filename=fname3 & x2_37 = x2 & z2_37 = z2
restore,filename=fname4 & x2_36_rot = x2 & z2_36_rot = z2
restore,filename=fname5 & x2_37_rot = x2 & z2_37_rot = z2
restore,filename=fname6 & r_36_rot = jetposcyl[1,*] & z_36_rot = jetposcyl[2,*]

loadct,39,/sil
if not keyword_set(ps) then begin
!p.background=255 & !p.color=0
window,1,xs=800,ys=700
endif else mkeps,'pp_jettrace',xs=20.,ys=17.
plot,x2_35,z2_35,xr=[-1.e13,2.e12],/xst,yr=[0.,1.e13],yst=2,psym=1,/iso,xtitle='x [cm]',ytitle='z [cm]',charsize=1.5
oplot,x2_36,z2_36,psym=1,color=50
oplot,x2_37,z2_37,psym=1,color=250
oplot,x2_36_rot,z2_36_rot,psym=dsym(8,/fill),color=50,symsize=1.
oplot,x2_37_rot,z2_37_rot,psym=dsym(8,/fill),color=250,symsize=1.


xc0=1.e12 & yc0=0.
xc1=!x.crange[0]
plotsym,3
plots,xc0,yc0,/data,psym=8,symsize=3

; fitted results
ind35 = where(x2_35 le -2.5e12)
lin35 = linfit(x2_35[ind35],z2_35[ind35])
x0=x2_35[ind35] & x0=reform(x0[0])
y0 = lin35[1]*x0+lin35[0]
oplot,[x0,!x.crange[0]],[lin35[1]*x0+lin35[0],lin35[0]+lin35[1]*!x.crange[0]]
oplot,[x0,xc1],[y0, -tan(16.69*!dtor)*(xc1-x0)+y0],line=2

ind36 = where((z2_36 ge 2.e12) and (z2_36 le 8.e12))
xx2 = x2_36[ind36] & xx2h = xx2[where(finite(xx2) eq 1)]
zz2 = z2_36[ind36] & zz2h = zz2[where(finite(xx2) eq 1)]
lin36 = linfit(xx2h,zz2h)
x0=reform(xx2h[0])
y0=lin36[0]+lin36[1]*x0
oplot,[x0,!x.crange[0]],[lin36[0]+lin36[1]*x0,lin36[0]+lin36[1]*!x.crange[0]],color=50
oplot,[x0,xc1],[y0,-tan(66.*!dtor)*(xc1-x0)+y0],color=50,line=2

ind36r = where((z2_36_rot ge 2.e12) and (z2_36_rot le 8.e12))
xx2 = x2_36_rot[ind36r] & xx2h = xx2[where(finite(xx2) eq 1)]
zz2 = z2_36_rot[ind36r] & zz2h = zz2[where(finite(xx2) eq 1)]
lin36r = linfit(xx2h,zz2h)
x0=reform(xx2h[0])
y0=lin36r[0]+lin36r[1]*x0
oplot,[x0,!x.crange[0]],[lin36r[0]+lin36r[1]*x0,lin36r[0]+lin36r[1]*!x.crange[0]],color=50

oplot,-r_36_rot,z_36_rot,psym=2,color=fsc_color('cyan')


ind37 = where(z2_37 ge 1.e12)
lin37 = linfit(x2_37[ind37],z2_37[ind37])
x0=x2_37[ind37] & x0=reform(x0[0])
y0=lin37[0]+lin37[1]*x0
oplot,[x0,!x.crange[0]],[lin37[0]+lin37[1]*x0,lin37[0]+lin37[1]*!x.crange[0]],color=250
oplot,[x0,xc1],[y0, -tan(78.68*!dtor)*(xc1-x0)+y0],color=250,line=2

ind37r = where((z2_37_rot ge 3.e12) and (z2_37_rot le 8.e12))
xx2 = x2_37_rot[ind37r] & xx2h = xx2[where(finite(xx2) eq 1)]
zz2 = z2_37_rot[ind37r] & zz2h = zz2[where(finite(xx2) eq 1)]
lin37r = linfit(xx2h,zz2h)
x0=reform(xx2h[0])
y0=lin37r[0]+lin37r[1]*x0
;oplot,[x0,!x.crange[0]],[lin37r[0]+lin37r[1]*x0,lin37r[0]+lin37r[1]*!x.crange[0]],color=250
oplot,[x0,(!y.crange[1]-lin37r[0])/lin37r[1]],[lin37r[0]+lin37r[1]*x0,!y.crange[1]],color=250


;oplot,[xc0,xc1],[yc0, -tan(16.69*!dtor)*(xc1-xc0)+yc0],line=2
;oplot,[xc0,xc1],[yc0, -tan(66.*!dtor)*(xc1-xc0)+yc0],color=50,line=2
;oplot,[xc0,xc1],[yc0, -tan(78.68*!dtor)*(xc1-xc0)+yc0],color=250,line=2

legend,textoidl('P_{jet}=')+['1e37','1e36','1e35','1e37_rot','1e36_rot'],psym=[1,1,1,dsym(8,/fill),dsym(8,/fill)],color=[250,50,0,250,50],textcolor=[250,50,0,250,50],/right,/top,charsize=1.5

if keyword_set(ps) then epsfree
stop
end

pro jettrace_anal,ps=ps

fname0 = '/d/d7/yoon/out_FLASH3.3_mhd/comp_sphere/jetpos_analxz.sav'
fname1 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_high/jettrace_1e36_high_891.sav'
;fname4 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_rot/Stampede/jettrace_1e36_rot_570.sav'
fname2 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_rot/Stampede/jettrace_1e36_rot_cyl_570.sav'
;fname3 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_rot/Brian/jettrace_1e36_rot_cyl_589.sav'
fname3 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_rot/Brian/jettrace_1e36_rot_cyl_534.sav'
fname4 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_rot/aci_vcorrect/jettrace_1e36_rot_cyl_534.sav'

restore,filename=fname0 & x0 = analx & z0 = analz                        ;analytic sol
restore,filename=fname1 & x1 = x2 & z1 = z2                              ;norot
restore,filename=fname2 & x2 = -jetposcyl[1,*]  & z2 = jetposcyl[2,*]    ;rot
restore,filename=fname3 & x3 = -jetposcyl[1,*]  & z3 = jetposcyl[2,*]    ;rot + corr
restore,filename=fname4 & x4 = -jetposcyl[1,*]  & z4 = jetposcyl[2,*]    ;rot + corr2


loadct,39,/sil

if not keyword_set(ps) then begin
!p.background=255 & !p.color=0
window,0,xs=450,ys=600
endif else mkeps,'jettrace_anal_1',xs=20.,ys=20.*600./450.

plot,x0,z0,xra=[-6.e12,-1.e12],/xst,/iso,yra=[0.,8e12],/yst,xticks=3,xtitle='x',ytitle='z',/nodata
oplot,x1,z1,color=50,psym=1
oplot,x2,z2,color=150,psym=1
oplot,x3,z3,color=250,psym=1
oplot,x4,z4,color=100,psym=1
oplot,x0,z0
legend,['norot','rot','rot+cor','rot+cor2'],color=[50,150,250,100],psym=1,/left,/bottom,box=0

if keyword_set(ps) then epsfree

;jetx0 = -2.d12
;bndang = (!pi/2. - atan(z0/(jetx0-x0))) /!dtor

bndang = (!pi/2 + atan(deriv(x0,z0)))/!dtor

intz1_bndang = interpol(bndang,z0,z1,/spline)
intz1_x0     = interpol(x0,z0,z1,/spline)
dev_x1 = abs(x1-intz1_x0)

intz2_bndang = interpol(bndang,z0,z2,/spline)
intz2_x0     = interpol(x0,z0,z2,/spline)
dev_x2 = abs(x2-intz2_x0)

intz3_bndang = interpol(bndang,z0,z3,/spline)
intz3_x0     = interpol(x0,z0,z3,/spline)
dev_x3 = abs(x3-intz3_x0)

intz4_bndang = interpol(bndang,z0,z4,/spline)
intz4_x0     = interpol(x0,z0,z4,/spline)
dev_x4 = abs(x4-intz4_x0)

if not keyword_set(ps) then $
 window,1,xs=800,ys=600 $
else mkeps,'jettrace_anal_2',xs=20.,ys=20.*600/800.

plot,intz1_bndang,dev_x1,psym=1,xtitle='bending angle [degree]' ,ytitle='deviation |x_anal - x_data|',/nodata,xmargin=[11.,3.]
oplot,intz1_bndang,dev_x1,psym=1,color=50
oplot,intz2_bndang,dev_x2,psym=1,color=150
oplot,intz3_bndang,dev_x3,psym=1,color=250
oplot,intz4_bndang,dev_x4,psym=1,color=100
legend,['norot','rot','rot+cor','rot+cor2'],color=[50,150,250,100],psym=1,/left,/top,box=0

if keyword_set(ps) then epsfree
stop
end
