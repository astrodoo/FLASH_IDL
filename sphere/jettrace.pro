pro jettrace, rotback=rotback

xc0 = -2.e12 
yc0 = 0. 
zc0 = 0. 

sample=3
;n=891
cd,'/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37_gam166' 
n=654
;velz = dload(700,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(490,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
jet = dload(n,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
velz = dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(1000,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(570,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(520,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(700,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
;jet = dload(349,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)

if keyword_set(rotback) then begin
   print, 'start rotating'
   m1 = 4.d34
   m2 = 2.d34
   mtot = m1+m2
   peri = 3.d12
   dth = sqrt(!unit.g * mtot / peri^3.d0) * time

   dth = dth / !dtor
   jet_rot = rot_3d(jet,rotaxis=3,degree=dth,/interp) 
   jet=jet_rot
   print,'complete rotating'
endif

y0ind_tmp = where(y ge 0) & y0ind = y0ind_tmp[0]
velz2 = reform(velz[*,y0ind,*])
jet2 = reform(jet[*,y0ind,*])

z0ind_tmp = where(z ge 0) & z0ind = z0ind_tmp[0]
sz = size(jet2,/dimension)
z2 = z[z0ind:*]
x2 = fltarr(sz[1]-z0ind)

critv=0.1
critj=0.7

for i=z0ind,sz[1]-1 do begin
;    maxv_xind = where(velz2[*,i] eq max(velz2[*,i]),count)
;    maxj_xind = where(jet2[*,i] eq max(jet2[*,i]),count)
    maxj_xind = where( (jet2[*,i] ge critj) and (velz2[*,i] ge 3.e9*critv) ,count)
    if (count ge 2) then maxj_xind = reform(maxj_xind[n_elements(maxj_xind)-1]) ;closest point to star 
    if ((count eq 0) or (maxj_xind eq 0)) then begin
       x2[i-z0ind] = !values.f_nan
    endif else begin
       x2[i-z0ind] = x[maxj_xind] 
    endelse
endfor

;x2[where(z2 le 0.5e12)]=!values.f_nan ; 1e37 case

window,0
plot,x2,z2,/iso

;save,filename='jettrace_1e35_700.sav',x2,z2
;save,filename='jettrace_1e36_800.sav',x2,z2
;save,filename='jettrace_1e36_700.sav',x2,z2
;save,filename='jettrace_1e36_high_851.sav',x2,z2
save,filename='jettrace_1e36_high_891.sav',x2,z2
;save,filename='jettrace_1e36_high_1000.sav',x2,z2
;save,filename='jettrace_1e37_490.sav',x2,z2
;save,filename='jettrace_1e37_349.sav',x2,z2
;save,filename='jettrace_1e36_rot_570.sav',x2,z2
;save,filename='jettrace_1e37_rot_520.sav',x2,z2

stop
end

;---------------------------------------------------------------------------------------
pro jettrace2,ps=ps

fname1 = '/d/d3/yoon/outputs/out_mHD_Binary_sphere/sphere_1e35/jettrace_1e35_700.sav'
fname2 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_high/jettrace_1e36_high_891.sav'
fname3 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37/jettrace_1e37_490.sav'
fname4 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_rot/Stampede/jettrace_1e36_rot_570.sav'
fname5 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37_rot/Stampede/jettrace_1e37_rot_520.sav'
fname6 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36/jettrace_1e36_700.sav'
fname7 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37/jettrace_1e37_349.sav'

restore,filename=fname1 & x2_35 = x2 & z2_35 = z2
restore,filename=fname2 & x2_36 = x2 & z2_36 = z2
restore,filename=fname3 & x2_37 = x2 & z2_37 = z2
restore,filename=fname4 & x2_36_rot = x2 & z2_36_rot = z2
restore,filename=fname5 & x2_37_rot = x2 & z2_37_rot = z2

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
