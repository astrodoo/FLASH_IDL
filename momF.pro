pro momF

n=800
sample=0
xc0=-3.8e12 & yc0=0. & zc0=4.791e12

dens=dload(n,var='dens',sample=sample,xc=xc0,yc=yc0,zc=zc0,x,y,z,time)
velz=dload(n,var='velz',sample=sample,xc=xc0,yc=yc0,zc=zc0)
ss = size(dens,/dimension)

dx=x[2]-x[1]
dy=y[2]-y[1]

z0ind_tmp = where(z ge 0) & z0ind = z0ind_tmp[0]

crit=0.7

z2 = z[z0ind:ss[2]-1]
area  = fltarr(ss[2]-z0ind)
ave_v = fltarr(ss[2]-z0ind)
ave_d = fltarr(ss[2]-z0ind)
ave_m = fltarr(ss[2]-z0ind)

for i=z0ind,ss[2]-1 do begin
   dens_cut = reform(dens[*,*,i])
   velz_cut = reform(velz[*,*,i])
   in_ind = where(velz_cut ge crit*3.e9,cnt)
   
   if (cnt eq 0) then begin
      area[i-z0ind]  = !values.f_nan   
      ave_v[i-z0ind] = !values.f_nan   
      ave_d[i-z0ind] = !values.f_nan   
      ave_m[i-z0ind] = !values.f_nan   
   endif else begin
      area[i-z0ind]  = cnt * dx*dy
      ave_v[i-z0ind] = total(velz_cut[in_ind])/cnt
      ave_d[i-z0ind] = total(dens_cut[in_ind])/cnt
      ave_m[i-z0ind] = ave_d[i-z0ind]*ave_v[i-z0ind]^2.*area[i-z0ind]
   endelse
endfor

loadct,39,/sil
!p.background=255 & !p.color=0
window,0
;mkeps,'momF_0.eps'
;plot,z2,area,xr=[0.,6.e12],/xst,xtitle='z [cm]',ytitle=textoidl('area [cm^{2}]'),xmargin=[33,4]
plot,z2,area,/xst,xtitle='z [cm]',ytitle=textoidl('area [cm^{2}]'),xmargin=[33,4]
axis,-1.25e12,yaxis=0,yr=[1.e9,4.e9],/yst,/save,color=50,ytitle='velz [cm/s]'
oplot,z2,ave_v,color=50
axis,-2.5e12,yaxis=0,yr=[1.e-17,5.e-15],/yst,/save,color=250,ytitle=textoidl('\rho [g/cm^{3}]')
oplot,z2,ave_d,color=250
;epsfree

window,1
;mkeps,'momF2_0.eps'
;plot,z2,ave_m,xtitle='z [cm]',ytitle='momentum flux (z-direction)',xr=[0.,6.e12],/xst,xmargin=[11,4]
plot,z2,ave_m,xtitle='z [cm]',ytitle='momentum flux (z-direction)',/xst,xmargin=[11,4]
;epsfree

save,filename='momF_0_3.sav',z2,area,ave_v,ave_d,ave_m
stop
end

pro momF_comb

restore,filename='momF_0_1.sav'
z2_1=z2 & area_1=area & ave_v1 = ave_v & ave_d1 = ave_d & ave_m1 = ave_m
restore,filename='momF_0_2.sav'
z2_2=z2 & area_2=area & ave_v2 = ave_v & ave_d2 = ave_d & ave_m2 = ave_m
restore,filename='momF_0_3.sav'
z2_3=z2 & area_3=area & ave_v3 = ave_v & ave_d3 = ave_d & ave_m3 = ave_m

loadct,39,/sil
;!p.background=255 & !p.color=0
;window,0
mkeps,'momF_0.eps',xs=20.,ys=20.*0.7
plot,z2_1,area_1,xr=[0.,5.5e12],/xst,xtitle='z [cm]',ytitle=textoidl('area [cm^{2}]'),xmargin=[33,2],yr=[0.,1.e23]
oplot,z2_2,area_2
oplot,z2_3,area_3
axis,-1.1e12,yaxis=0,yr=[1.e9,4.e9],/yst,/save,color=50,ytitle='velz [cm/s]'
oplot,z2_1,ave_v1,color=50
oplot,z2_2,ave_v2,color=50
oplot,z2_3,ave_v3,color=50
axis,-2.2e12,yaxis=0,yr=[1.e-17,8.e-14],/yst,/save,color=250,ytitle=textoidl('\rho [g/cm^{3}]'),/ylog
oplot,z2_1,ave_d1,color=250
oplot,z2_2,ave_d2,color=250
oplot,z2_3,ave_d3,color=250
epsfree

;window,1
mkeps,'momF2_0.eps',xs=20.,ys=20.*0.8
plot,z2_1,ave_m1,xtitle='z [cm]',ytitle='momentum flux (z-direction)',xr=[0.,5.5e12],/xst,xmargin=[11,2]
oplot,z2_2,ave_m2
oplot,z2_3,ave_m3
epsfree

stop
end
