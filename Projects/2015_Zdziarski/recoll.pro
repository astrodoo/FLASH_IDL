pro recoll,mkdata=mkdata
device,decomposed=0

restore,file='jet_snap_1e38.sav'

zcut0 = 0.
zcind = (where(z ge zcut0))[0]

nz = n_elements(z)
z = reform(z[zcind:nz-1])
jyz = reform(jyz[*,zcind:nz-1])
v3yz = reform(v3yz[*,zcind:nz-1])

nz2 = n_elements(z)
velj = 3.d9
jv = jyz*abs(v3yz)/velj
jcrit = 0.95

jyl = fltarr(nz2) & jyr = fltarr(nz2)
for i=0,nz2-1 do begin
    jyind = where(jv[*,i] ge jcrit,count)
    if (count ne 0) then begin
       jyl[i] = y[jyind[0]]
       jyr[i] = y[jyind[n_elements(jyind)-1]]
    endif else begin
       jyl[i] = 0.
       jyr[i] = 0.
    endelse
endfor

save,file='recoll_1e38.sav',z,jyl,jyr,time

loadct,0,/sil
window,0
plot,jyl,z,xtitle='y [cm]',ytitle='z [cm]',/iso,xrange=[-3.e12,3.e12],yrange=[0.,1.5e13],/xst,/yst
oplot,jyr,z

end

pro jet_resol
restore,file='jet_snap_1e38.sav'
zcut0 = 0.
zcind = (where(z ge zcut0))[0]
nz = n_elements(z)
z = reform(z[zcind:nz-1])
lyz = reform(lyz[*,zcind:nz-1])
nz2 = n_elements(z)

restore,file='recoll_1e38.sav'

j_thick = jyr-jyl

j_lref = fltarr(nz2)
j_ncell = fltarr(nz2)
boxsz = 1.2e13
for i=0,nz2-1 do begin
    ylind = (where(y ge jyl[i]))[0]
    yrind = (where(y ge jyr[i]))[0]
    j_lref[i] = mean(lyz[ylind:yrind,i])
    j_ncell[i] = j_thick[i] / (boxsz / 2.^(j_lref[i]+2.))
endfor

drefz = [3.e11,6.e11,1.e12,2.e12,4.e12,1.e13]

mkeps,'jet_resol',xs=25.,ys=25.*6./8.
multiplot,[1,3]
plot,z,j_thick,/iso,xra=[0.,1.5e13],/xst,ytitle='jet thick',ytickinterval=1.e12
for i=0,5 do oplot,[drefz[i],drefz[i]],!y.crange,line=2
multiplot
plot,z,j_lref,xra=[0.,1.5e13],/xst,ytitle='refine in jet',yra=[4,12],yst=2
for i=0,5 do oplot,[drefz[i],drefz[i]],!y.crange,line=2
multiplot
plot,z,j_ncell,xra=[0.,1.5e13],/xst,xtitle='z [cm]',ytitle='ncell in jet'
for i=0,5 do oplot,[drefz[i],drefz[i]],!y.crange,line=2
multiplot,/reset

epsfree

stop
end
