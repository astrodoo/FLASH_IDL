pro jettrace,n,mkdata=mkdata,sample=sample,xc0=xc0,yc0=yc0,zc0=zc0

if not keyword_set(xc0) then xc0 = -2.e12 
if not keyword_set(yc0) then yc0 = 0.
if not keyword_set(zc0) then zc0 = 5.e11
if (n_elements(sample) eq 0) then sample=2

fname='jettrace_1e36_'+strtrim(n,2)
;fname='jettrace2_1e35_'+strtrim(n,2)+'_smp'+strtrim(sample,2)
;fname='jettrace2_1e36_60deg_'+strtrim(n,2)+'_smp'+strtrim(sample,2)

if keyword_set(mkdata) then begin
   d  = dload(n,var='dens',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
   d  = reform(d[*,256,*])
   j  = dload(n,var='jet',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
   j  = reform(j[*,256,*])
   vz = dload(n,var='velz',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=sample,time)
   vz = reform(vz[*,256,*])

   save,file=fname+'.sav',d,j,vz,x,z,time
endif else restore,file=fname+'.sav' 

loadct,0,/sil
window,0,xs=1024, ys=1024
tvscl,alog(d),0
tvscl,j,1
tvscl,vz,2

sz = n_elements(z)
z0ind_tmp = where(z ge 0) & z0ind = z0ind_tmp[0]
z2 = z[z0ind:sz-1]
x2 = fltarr(sz-z0ind)

jj = vz*j
for i=z0ind,sz-1 do begin
    maxv_xind = where(jj[*,i] eq max(jj[*,i]),count)
;    maxv_xind = where(vz[*,i] eq max(vz[*,i]),count)
    if (count ge 2) then maxv_xind = reform(maxv_xind[n_elements(maxv_xind)-1]) ;closest point to star 
    x2[i-z0ind] = x[maxv_xind] 
    if (j[maxv_xind,i] le 1.e-2) then x2[i-z0ind] = !values.f_nan
endfor

window,3
plot,x2,z2,/iso

save_append,file=fname+'.sav',x2,z2

stop
end

;---------------------------------------------------------------------------------------
pro jettrace_draw, ps=ps, mkdata=mkdata

fout='jettrace_comb.sav'

if keyword_set(mkdata) then begin

; E35 data
fn_e35_1 = '/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e35/jettrace2_1e35_700_smp1.sav'
fn_e35_2 = '/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e35/jettrace2_1e35_700_smp2.sav'
fn_e35_3 = '/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e35/jettrace2_1e35_700_smp3.sav'

restore,filename=fn_e35_1 & x2_35_1 = x2 & z2_35_1 = z2
restore,filename=fn_e35_2 & x2_35_2 = x2[1:*] & z2_35_2 = z2[1:*]
restore,filename=fn_e35_3 & x2_35_3 = x2 & z2_35_3 = z2

x2_35_dummy1 = x2_35_2[where(x2_35_2 lt min(x2_35_1))]
x2_35_dummy2 = x2_35_3[where(x2_35_3 lt min(x2_35_2))]
x2_35 = [x2_35_1,x2_35_dummy1,x2_35_dummy2]
z2_35_dummy1 = z2_35_2[where(x2_35_2 lt min(x2_35_1))]
z2_35_dummy2 = z2_35_3[where(x2_35_3 lt min(x2_35_2))]
z2_35 = [z2_35_1,z2_35_dummy1,z2_35_dummy2]

x2_35[where(z2_35 gt 3.e12)] = !values.f_nan
z2_35[where(z2_35 gt 3.e12)] = !values.f_nan
x2_35 = x2_35[where(finite(x2_35))]
z2_35 = z2_35[where(finite(x2_35))]


; E36 data
;fn_e36_1 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/_WRO_sphere_1e36_high/jettrace2_1e36_851_smp3.sav'
;fn_e36_2 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/_WRO_sphere_1e36_high/jettrace2_1e36_851_smp4.sav'
fn_e36_1 = 'jettrace_1e36_459_smp3.sav'

restore,filename=fn_e36_1 & x2_36 = x2 & z2_36 = z2
;restore,filename=fn_e36_2 & x2_36_2 = x2 & z2_36_2 = z2

;z2_36_dummy = z2_36_2[where(z2_36_2 gt max(z2_36_1))]
;z2_36 = [z2_36_1,z2_36_dummy]
;x2_36_dummy = x2_36_2[where(z2_36_2 gt max(z2_36_1))]
;x2_36 = [x2_36_1,x2_36_dummy]

; eliminate bad pixels
x2_36[where((z2_36 ge 6.e12))] = !values.f_nan
z2_36[where((z2_36 ge 6.e12))] = !values.f_nan
x2_36 = x2_36[where(finite(x2_36))]
z2_36 = z2_36[where(finite(x2_36))]


; E36 rot & acc
fn_e36_rot = 'jettrace_1e36_rot_cyl_72.sav'
;fn_e36_rot = '/d/d7/yoon/out_FLASH3.3_mhd/brian/out_sphere_1e36_rot/jettrace_1e36_rot_cyl_589.sav'
;fn_e36_acc = '/d/d7/yoon/out_FLASH3.3_mhd/aci/out_sphere_1e36_accel/jettrace_1e36_rot_cyl_260.sav'

restore,filename=fn_e36_rot & x2_36_rot = -jetposcyl[1,*] & z2_36_rot = jetposcyl[2,*]
;restore,filename=fn_e36_acc & x2_36_acc = -jetposcyl[1,*] & z2_36_acc = jetposcyl[2,*]
;zcrit_acc = 3.7e12
;x2_36_acc = x2_36_acc[where(z2_36_acc le zcrit_acc)]
;z2_36_acc = z2_36_acc[where(z2_36_acc le zcrit_acc)]

; eliminate bad pixels
x2_36_rot_tmp = x2_36_rot & z2_36_rot_tmp = z2_36_rot

x2_36_rot_tmp[where((z2_36_rot ge 6.e12))] = !values.f_nan
z2_36_rot_tmp[where((z2_36_rot ge 6.e12))] = !values.f_nan

x2_36_rot_tmp[where((x2_36_rot le -3.5e12) and (z2_36_rot le 4.e12))] = !values.f_nan
x2_36_rot_tmp[where((x2_36_rot le -2.8e12) and (z2_36_rot le 3.e12))] = !values.f_nan
;z2_36_rot[where((x2_36_rot le -1.e12) and (z2_36_rot ge 4.e12))] = !values.f_nan


x2_36_rot = x2_36_rot[where(finite(x2_36_rot_tmp))]
z2_36_rot = z2_36_rot[where(finite(x2_36_rot_tmp))]


; E37 data 
fn_e37_1   = 'jettrace_1e37_380_smp3.sav'
fn_e37_2   = 'jettrace_1e37_380_smp4.sav'

restore,filename=fn_e37_1 & x2_37_1 = x2 & z2_37_1 = z2
restore,filename=fn_e37_2 & x2_37_2 = x2 & z2_37_2 = z2

z2_37_dummy = z2_37_2[where(z2_37_2 gt max(z2_37_1))]
z2_37 = [z2_37_1,z2_37_dummy]
x2_37_dummy = x2_37_2[where(z2_37_2 gt max(z2_37_1))]
x2_37 = [x2_37_1,x2_37_dummy]

x2_37 = x2_37[where(finite(x2_37))]
z2_37 = z2_37[where(finite(x2_37))]


; E37 gam166 data
fn_e37_g5_1   = 'jettrace_1e37_gam166_354_smp3.sav'
fn_e37_g5_2   = 'jettrace_1e37_gam166_354_smp4.sav'

restore,filename=fn_e37_g5_1 & x2_37_g5_1 = x2 & z2_37_g5_1 = z2
restore,filename=fn_e37_g5_2 & x2_37_g5_2 = x2 & z2_37_g5_2 = z2

z2_37_g5_dummy = z2_37_g5_2[where(z2_37_g5_2 gt max(z2_37_g5_1))]
z2_37_g5 = [z2_37_g5_1,z2_37_g5_dummy]
x2_37_g5_dummy = x2_37_g5_2[where(z2_37_g5_2 gt max(z2_37_g5_1))]
x2_37_g5 = [x2_37_g5_1,x2_37_g5_dummy]

x2_37_g5 = x2_37_g5[where(finite(x2_37_g5))]
z2_37_g5 = z2_37_g5[where(finite(x2_37_g5))]


jetpos,th0=0,Lj=1e35,xjet=xjet_35,yjet=zjet_35
jetpos,th0=0,Lj=1e36,xjet=xjet_36,yjet=zjet_36
jetpos,th0=0,Lj=1e37,xjet=xjet_37,yjet=zjet_37
jetpos,th0=0,Lj=1e37,xjet=xjet_37_g5,yjet=zjet_37_g5,/gam

save,file=fout,x2_35,z2_35,x2_36,z2_36,x2_37,z2_37,xjet_35,zjet_35,xjet_36,zjet_36,xjet_37,zjet_37 $
               ,x2_37_g5, z2_37_g5, xjet_37_g5, zjet_37_g5, x2_36_rot,z2_36_rot ;$
;               ,x2_36_acc,z2_36_acc

endif else restore, file=fout

;--------------------------------------------------------------------------------
; drawing
;

;restore, file='/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37_gam166/c2var_xyjet.sav'
;xjet_37_g5=xjet_tp & zjet_37_g5=yjet_tp

loadct,39,/sil
if not keyword_Set(ps) then begin 
  !p.background=255 & !p.color=0
  window,1,xs=600,ys=600 
endif else mkeps,'jettrace_comb',xs=18.,ys=20.

plot,x2_35,z2_35,xra=[-5.e12,1.e12],yra=[0.,8.e12],/nodata,/iso,xticks=3,xtitle='x [cm]', ytitle='z [cm]' $
    ,/xst,/yst,xtickinterval=5.e12
;plot,x2_35,z2_35,xra=[-7.e12,3.e12],yra=[0.,1.e13],/nodata,/iso,xticks=3,xtitle='x [cm]', ytitle='z [cm]' $
;    ,xmargin=[13.,0.1],/xst,/yst,xtickinterval=5.e12

oplot,x2_35,z2_35,psym=dsym(6),symsize=0.5
oplot,xjet_35,zjet_35
oplot,x2_36,z2_36,psym=dsym(5),color=50,symsize=0.5
oplot,xjet_36,zjet_36,color=50
oplot,x2_37,z2_37,psym=dsym(8),color=254,symsize=0.5
oplot,xjet_37,zjet_37,color=254

oplot,x2_37_g5,z2_37_g5,psym=dsym(8),color=200,symsize=0.5
oplot,xjet_37_g5,zjet_37_g5,color=200

oplot,x2_36_rot,z2_36_rot,psym=dsym(4),color=100,symsize=0.5
;oplot,x2_36_acc,z2_36_acc,psym=dsym(7),color=150,symsize=0.7, thick=0.5

plots,1.e12,0.,psym=dsym(16), symsize=3

xc0=1.e12 & yc0=0.
xc1=!x.crange[0]
oplot,[xc0,xc1],[yc0, -tan(18*!dtor)*(xc1-xc0)+yc0],line=2
oplot,[xc0,xc1],[yc0, -tan(62.*!dtor)*(xc1-xc0)+yc0],color=50,line=2
oplot,[xc0,xc1],[yc0, -tan(82*!dtor)*(xc1-xc0)+yc0],color=250,line=2

;legend,textoidl('L_{j}=')+[textoidl('10^{35}'),textoidl('10^{36}'),textoidl('10^{37}')]+textoidl(' erg s^{-1}'),color=[0,50,254],textcolor=[0,50,254] $
;      ,box=0,/right,/top
legend,'SphWind_'+['E35','E36','E37','E37_gam166','E36_rot'],color=[0,50,254,200,100],textcolor=[0,50,254,200,100] $
      ,box=0,/right,/top,psym=['dsym(6)','dsym(5)','dsym(8)','dsym(8)','dsym(4)'],charsize=1.2

if keyword_Set(ps) then epsfree
stop
end







;==========================================================================
; getting analytical lines
forward_function Fth, PQ_Limits
function Fth,th1,th2
common Power, intpow
   return, (1.+tan(th1)^2.)*cos(th2)^intpow
end

function PQ_Limits,th1
common Param, th_start
;print, 'th_start=', th_start
   return, [th_start,th1]
end

pro jetpos,th0=th0,xjet=xjet,yjet=yjet,nth=nth,theta=theta, Lj=Lj, gam53=gam53
common Param, th_start
common Power, intpow

if keyword_set(gam53) then begin
   print,'gamma = 5/3'  
   intpow = 0.8 
endif else begin
   print,'gamma = 4/3'
   intpow = 0.5
endelse

if (n_elements(th0) eq 0) then th0 = 0.

print,'Lj=  ',Lj,'      th0 = ', th0 

if not keyword_set(nth) then nth=10000

case Lj of
   1e35: h0=3.e10
   1e36: h0=6.5e10
   1e37: if keyword_set(gam53) then h0=6.5e10 else h0=6.5e10
   else: print, 'out of range in Energy'
endcase

print,'h0 = ',h0

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

save,file='jettrace_anal.sav',xjet,yjet

case Lj of
   1e35: restore,file='/d/d2/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e35/jettrace2_1e35_700_smp3.sav' 
;   1e36: restore,file='/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/_WRO_sphere_1e36_high/jettrace2_1e36_851_smp3.sav' 
   1e36: restore,file='/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36/jettrace_1e36_417.sav'
   1e37: restore,file='/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/_WRO_sphere_1e37_old/jettrace2_1e37_400_smp4.sav' 
   else: print, 'out of range in Energy'
endcase

window,4
plot,xjet,yjet,xra=[-6.e12,-1.e12],/xst,/iso,yra=[0.,8e12],/yst
oplot,x2,z2
end





pro errorchk, ps=ps

restore, file='jettrace_comb.sav'

; for off-axis jets data
;restore, file='/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_30deg/jettrace2_1e36_30deg_530_smp3.sav'
;restore, file='/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_60deg/jettrace2_1e36_60deg_460_smp3.sav'

stx = 1.d12

;d35 = sqrt( (xjet_35-stx)^2. + zjet_35*zjet_35 )
;th_35 = acos(abs(xjet_35-stx)/d35)
;d36 = sqrt( (xjet_36-stx)^2. + zjet_36*zjet_36 )
;th_36 = acos(abs(xjet_36-stx)/d36)
;d37 = sqrt( (xjet_37-stx)^2. + zjet_37*zjet_37 )
;th_37 = acos(abs(xjet_37-stx)/d37)

nflag = 9990

slop35 = deriv(xjet_35,zjet_35)
ang_35 = 90.-atan(-slop35)/!dtor 
ang_35 = ang_35[0:nflag] & xjet_35 = xjet_35[0:nflag] & zjet_35 = zjet_35[0:nflag]
slop36 = deriv(xjet_36,zjet_36)
ang_36 = 90.-atan(-slop36)/!dtor 
ang_36 = ang_36[0:nflag] & xjet_36 = xjet_36[0:nflag] & zjet_36 = zjet_36[0:nflag]
slop37 = deriv(xjet_37,zjet_37)
ang_37 = 90.-atan(-slop37)/!dtor 
ang_37 = ang_37[0:nflag] & xjet_37 = xjet_37[0:nflag] & zjet_37 = zjet_37[0:nflag]

;slop30deg = deriv(xjet_30deg,zjet_30deg)
;ang_30deg = -60.-atan(-slop30deg)/!dtor  
;ang_30deg = ang_30deg[0:nflag] & xjet_30deg = xjet_30deg[0:nflag] & zjet_30deg = zjet_30deg[0:nflag]
;ang_30deg[where(ang_30deg lt 0)] += 180.
;ang_30deg[0]=0.

xjet_36_rot = xjet_36 & zjet_36_rot=zjet_36
;xjet_36_acc = xjet_36 & zjet_36_acc=zjet_36

;xjet_60deg = xjet_36_60deg & zjet_60deg = zjet_36_60deg
;x2_60deg = x2_36_60deg & z2_60deg = z2_36_60deg
;x2_60deg = x2_60deg(where(finite(x2_60deg))) & z2_60deg = z2_60deg(where(finite(x2_60deg)))

;slop60deg = deriv(xjet_60deg,zjet_60deg)
;ang_60deg = -30.-atan(-slop60deg)/!dtor 
;ang_60deg = ang_60deg[0:nflag] & xjet_60deg = xjet_60deg[0:nflag] & zjet_60deg = zjet_60deg[0:nflag]
;ang_60deg[where(ang_60deg lt 0)] += 180.
;ang_60deg[0]=0.

intpx35 = interpol(x2_35,z2_35,zjet_35)
xjet_35[where(intpx35 le -1.1e13)] = !values.f_nan
zjet_35[where(intpx35 le -1.1e13)] = !values.f_nan
ang_35[where(intpx35 le -1.1e13)] = !values.f_nan
intpx35[where(intpx35 le -1.1e13)] = !values.f_nan
xjet_35 = xjet_35[where(finite(xjet_35))]
zjet_35 = zjet_35[where(finite(zjet_35))]
intpx35 = intpx35[where(finite(intpx35))]
ang_35 = ang_35[where(finite(ang_35))]
dev35 = abs(xjet_35 - intpx35)

intpx36 = interpol(x2_36,z2_36,zjet_36)
xjet_36[where(intpx36 le -1.5e13)] = !values.f_nan
zjet_36[where(intpx36 le -1.5e13)] = !values.f_nan
intpx36[where(intpx36 le -1.5e13)] = !values.f_nan
xjet_36 = xjet_36[where(finite(xjet_36))]
zjet_36 = zjet_36[where(finite(zjet_36))]
intpx36 = intpx36[where(finite(intpx36))]
dev36 = abs(xjet_36 - intpx36)

intpx36_rot = interpol(x2_36_rot,z2_36_rot,zjet_36_rot)
;xjet_36_rot[where(intpx36_rot le -1.5e13)] = !values.f_nan
;zjet_36_rot[where(intpx36_rot le -1.5e13)] = !values.f_nan
;intpx36_rot[where(intpx36_rot le -1.5e13)] = !values.f_nan
xjet_36_rot = xjet_36_rot[where(finite(xjet_36_rot))]
zjet_36_rot = zjet_36_rot[where(finite(zjet_36_rot))]
intpx36_rot = intpx36_rot[where(finite(intpx36_rot))]
dev36_rot = abs(xjet_36_rot - intpx36_rot)

;intpx36_acc = interpol(x2_36_acc,z2_36_acc,zjet_36_acc)
;;xjet_36_acc[where(intpx36_acc le -1.5e13)] = !values.f_nan
;;zjet_36_acc[where(intpx36_acc le -1.5e13)] = !values.f_nan
;;intpx36_acc[where(intpx36_acc le -1.5e13)] = !values.f_nan
;xjet_36_acc = xjet_36_acc[where(finite(xjet_36_acc))]
;zjet_36_acc = zjet_36_acc[where(finite(zjet_36_acc))]
;intpx36_acc = intpx36_acc[where(finite(intpx36_acc))]
;dev36_acc = abs(xjet_36_acc - intpx36_acc)

intpx37 = interpol(x2_37,z2_37,zjet_37)
xjet_37[where(zjet_37 ge 1.8e13)] = !values.f_nan
intpx37[where(zjet_37 ge 1.8e13)] = !values.f_nan
zjet_37[where(zjet_37 ge 1.8e13)] = !values.f_nan
xjet_37 = xjet_37[where(finite(xjet_37))]
zjet_37 = zjet_37[where(finite(zjet_37))]
intpx37 = intpx37[where(finite(intpx37))]
dev37 = abs(xjet_37 - intpx37)

;intpx30deg = interpol(x2_30deg,z2_30deg,zjet_30deg)
;dev30deg   = abs(xjet_30deg - intpx30deg)

;intpx60deg = interpol(x2_60deg,z2_60deg,zjet_60deg)
;dev60deg   = abs(xjet_60deg - intpx60deg)

sep=3.e12
loadct,39,/sil
if not keyword_set(ps) then begin
   !p.background=255 & !p.color=0
   window,1
endif else mkeps,'jettrace2_errorchk',xs=20.,ys=20.*6./8.

plot,ang_35, dev35/abs(xjet_35),xtitle='jet bending angle [degree]',ytitle='fractional deviation',xrange=[0.,70.],yrange=[-0.05,1.0],/xst,/yst,/nodata
oplot,ang_35, dev35/abs(xjet_35), line=2
oplot,ang_36, dev36/abs(xjet_36),color=50, line=3
oplot,ang_37, dev37/abs(xjet_37),color=250, line=4

oplot,ang_36, dev36_rot/abs(xjet_36_rot),color=100,line=5
;oplot,ang_36, dev36_acc/abs(xjet_36_acc),color=150,line=1

;dev30deg[where((ang_30deg ge 40) and (dev30deg le 6.e12))] = 6.e12
;oplot,ang_30deg, dev30deg/abs(xjet_30deg),color=30
;dev60deg[where((ang_60deg ge 60) and (dev60deg le 3.5e12))] = 3.5e12
;oplot,ang_60deg, dev60deg/abs(xjet_60deg),color=210

oplot,[20.,20.],!y.crange,line=0
legend,'SphWind_'+['E35','E36','E37','E36_rot','E36_acc'],color=[0,50,254,100,150],textcolor=[0,50,254,100,150] $
      ,box=0,/right,/bottom, charsize=1., line=[2,3,4,5,1]
;legend,'SphWind_'+['E35','E36','E37','E36_rot','E36_acc','E36_30deg','E36_60deg'],color=[0,50,254,100,150,30,210],textcolor=[0,50,254,100,150,30,210] $
;      ,box=0,/left,/top, charsize=1.

if keyword_Set(ps) then epsfree

stop
end
