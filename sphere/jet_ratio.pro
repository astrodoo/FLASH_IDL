pro jet_ratio

;spherical wind
;sample=1
;n=870
;xc0=-2.e12 & yc0=0. & zc0=0.

;uniform wind
sample=1
n=800
xc0=0. & yc0=0. & zc0=0.

jcrit=(vcrit=0.7)

velz = dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample)
jet  = dload(n,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample)
ss = size(jet,/dimension)

y0ind_tmp = where(y ge 0) & y0ind = y0ind_tmp[0]
velz_xz = reform(velz[*,y0ind,*])
jet_xz  = reform(jet[*,y0ind,*])

z0ind_tmp = where(z ge 0) & z0ind = z0ind_tmp[0]

z2 = z[z0ind:*]
thickxy = fltarr(ss[2]-z0ind)
thickxz = fltarr(ss[2]-z0ind)
x2 = fltarr(ss[2]-z0ind)

for i=z0ind,ss[2]-1 do begin
; xz plane
;   thickxz_ind = where(velz_xz[*,i] ge vcrit*3.e9,cnt)
;   thickxz_ind = where(jet_xz[*,i] ge jcrit,cnt)
   thickxz_ind = where((velz_xz[*,i] ge vcrit*3.e9) and (jet_xz[*,i] ge jcrit),cnt)
   if (cnt eq 0) then thickxz[i-z0ind] = !values.f_nan $
     else if (thickxz_ind[0] eq 0) then thickxz[i-z0ind] = !values.f_nan $
     else thickxz[i-z0ind]=x[thickxz_ind[n_elements(thickxz_ind)-1]]-x[thickxz_ind[0]]
   
;xy plane
   maxv_xind = where(velz_xz[*,i] eq max(velz_xz[*,i]),cnt)
;   maxv_xind = where(jet_xz[*,i] eq max(jet_xz[*,i]),cnt)
   if (cnt ge 2) then maxv_xind = reform(maxv_xind[0])

;   thickxy_ind = where(velz[maxv_xind,*,i] ge crit*3.e9,cnt)
;   thickxy_ind = where(jet[maxv_xind,*,i] ge crit,cnt)
   thickxy_ind = where((velz[maxv_xind,*,i] ge vcrit*3.e9) and (jet[maxv_xind,*,i] ge jcrit),cnt)
   if ((cnt eq 0) or (maxv_xind eq 0)) then begin
      x2[i-z0ind]=!values.f_nan
      thickxy[i-z0ind]=!values.f_nan
   endif else begin 
      x2[i-z0ind]=x[maxv_xind]
      thickxy[i-z0ind] = y[thickxy_ind[n_elements(thickxy_ind)-1]]-y[thickxy_ind[0]]
   endelse
endfor
 
;thickxz[where(z2 le 2.e11)] = !values.f_nan
;thickxy[where(z2 le 2.e11)] = !values.f_nan

window,0,xs=800,ys=800
;mkeps,'jet_ratio_uni_800.eps',xs=20.,ys=20.
;mkeps,'jet_ratio_sph_800.eps',xs=20.,ys=20.
multiplot,[1,2],xgap=0.01
plot,z2,thickxz,ytitle='Jet thickness [cm]',xr=[0.,2.e12],/xst
oplot,z2,thickxy,line=2
legend,['thick_xz','thick_xy'],line=[0,2],box=0,/left,/top
multiplot
plot,z2,thickxz/thickxy,xtitle='z [cm]',ytitle='thick_xz / thick_xy',xr=[0.,2.e12],/xst
multiplot,/reset
;epsfree
stop
end
