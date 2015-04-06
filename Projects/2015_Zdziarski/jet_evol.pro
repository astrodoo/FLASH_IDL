pro jet_evol,mk_rddata=mk_rddata,mkdata=mkdata


zcut0 = 0. & zcut1 = 1.5e13
if keyword_set(mk_rddata) then begin
   restore,file='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/L38_M30_lg/clip2d.sav'
   ndata = size(data.dyz,/dimension)

   zc0_ind = (where(z ge zcut0))[0]
   zc1_ind = (where(z ge zcut1))[0]
   z = reform(z[zc0_ind:zc1_ind])
   nz = n_elements(z)

   data_rd = replicate({dyz:fltarr(ndata[0],nz),jyz:fltarr(ndata[0],nz),v3yz:fltarr(ndata[0],nz),time:0.},ndata[2])

   data_rd.dyz = reform(data.dyz[*,zc0_ind:zc1_ind])
   data_rd.jyz = reform(data.jyz[*,zc0_ind:zc1_ind])
   data_rd.v3yz = reform(data.v3yz[*,zc0_ind:zc1_ind])
   data_rd.time = data.time

;   dyz = reform(dyz[*,zc0_ind:zc1_ind])
;   jyz = reform(jyz[*,zc0_ind:zc1_ind])
;   v3yz = reform(v3yz[*,zc0_ind:zc1_ind])

   save,file='clip2d_reduced.sav', (data=data_rd),y,z,xcut
;endif else restore,file='tmp_evol.sav'
endif else restore,file='clip2d_reduced.sav'

print,'complete to read the clip2d_reduced.sav data'

;identifying jet

if keyword_set(mkdata) then begin
ndata = size(data,/dimension)

nz2 = n_elements(z)
velj = 3.d9
jcrit = 0.95

jyl = fltarr(nz2,ndata[0]) & jyr = fltarr(nz2,ndata[0])
for i=0,ndata[0]-1 do begin

   jv = data[i].jyz*abs(data[i].v3yz)/velj

   for j=0,nz2-1 do begin
       jyind = where(jv[*,j] ge jcrit,count)
       if (count ne 0) then begin
          jyl[j,i] = y[jyind[0]]
          jyr[j,i] = y[jyind[n_elements(jyind)-1]]
       endif else begin
          jyl[j,i] = 0.
          jyr[j,i] = 0.
       endelse
   endfor
endfor

save,file='jet_evol.sav',jyl,jyr,y,z,(time=data.time)
endif else restore,file='jet_evol.sav'

print,'complete to read the jet_evol.sav data'

outdir = 'png_jet_evol'
spawn,'mkdir '+outdir
!p.background=255 & !p.color=0.
pltx0=130 &plty0=80
pltsz = size(data.dyz,/dimension)
winxs = pltx0+pltsz[0]+400 & winys = plty0+pltsz[1]+100
mind = 1.e-16 & maxd = 5.e-12

ndata = n_elements(time)
for i=0,ndata-1 do begin

print,i,' of  ',ndata 
 
loadct,0,/sil
window,0,xs=winxs,ys=winys,/pixmap

tvcoord,bytscl(alog10(data[i].dyz),min=alog10(mind),max=alog10(maxd)),y,z,pos=[pltx0,plty0,pltx0+pltsz[0],plty0+pltsz[1]] $
       ,/axes,xtitle='y [cm]',ytitle='z [cm]',/black
color_bar,lim=[mind,maxd],/log,/top,bartitle='density',titlegap=0.06

jyl2 = reform(jyl[*,i]) & jyr2 = reform(jyr[*,i])
not_j = where((jyl2 eq 0) and (jyr2 eq 0) and (z ge 3.e11),count) 

if (count ne 0) then begin
   jyl2[not_j] = !values.f_nan 
   jyr2[not_j] = !values.f_nan
endif

plot,jyl2,z,xtitle='y [cm]',/iso,xrange=[-4.e12,4.e12],yrange=[zcut0,zcut1],/xst,/yst $
    ,position=[pltx0+pltsz[0],plty0,pltx0+2*pltsz[0],plty0+pltsz[1]],/dev,/noerase,ytickformat='(a1)'$
    ,xtickinterval=2.e12,title='identified jet'
oplot,jyr2,z

xyouts,winxs-200,winys-30,/dev,'time:  '+string(time[i]/60./60.,format='(f5.2)')+' hrs'

snapshot,outdir+'/jet_evol_'+string(i,format='(I3.3)')
endfor
stop
end

pro jet_evol_comb

restore,file='jet_evol.sav'

loadct,0,/sil
;!p.background=255 & !p.color=0
mkeps,'jet_evol_comb',xs=15.,ys=20.
;window,0
plot,y,z,/nodata,xra=[-4.e12,4.e12],yra=[0.,1.5e13],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z[cm]',xmargin=[12,1] $
    ,title='jet evolution'

ndata = n_elements(time)

loadct,13,/sil
for i=0,ndata-1 do begin

  jyl2 = reform(jyl[*,i]) & jyr2 = reform(jyr[*,i])
  not_j = where((jyl2 eq 0) and (jyr2 eq 0) and (z ge 3.e11),count) 

  if (count ne 0) then begin
     if (not_j[0] ge 50) then begin
        jyl2[not_j[0]-50:*] = !values.f_nan 
        jyr2[not_j[0]-50:*] = !values.f_nan
     endif else if (not_j[0] ge 5) then begin
        jyl2[not_j[0]-5:*] = !values.f_nan
        jyr2[not_j[0]-5:*] = !values.f_nan
     endif
  endif

  oplot,jyl2,z,color=float(i)/float(ndata)*255
  oplot,jyr2,z,color=float(i)/float(ndata)*255
endfor

color_bar,lim=([time[0],time[ndata-1]]-time[0])/60./60.,/right,bartitle='time (hrs)',titlegap=0.07

epsfree
stop
end
