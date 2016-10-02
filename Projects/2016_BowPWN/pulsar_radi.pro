pro pulsar_radi,savepng=savepng
;/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/600kms/v7
;/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/600kms/v8
;/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/600kms/v9
;/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/600kms/v10

title='v10'
restore,file='save_pwn_'+title+'.sav'
nt = n_elements(time)
nr = n_elements(r)

radi = fltarr(nt)

loadct,0,/sil
!p.background=255 & !p.color=0
if keyword_set(savepng) then $
   window,0,xs=800,ys=900,/pixmap $
 else window,0,xs=800,ys=900

x0=100 & x1=770 & y0=50
ny = 275 

if keyword_set(savepng) then spawn,'mkdir png_radi'

for i=0,nt-1 do begin
   print,i, ' of ', nt
;radius finding
   dd = reform(d[*,i])
;   ddr = deriv(r,dd)
;   ind_tmp = where(ddr ge 0) & ind = ind_tmp[0]
;   dtmp = dd[ind:*] & rtmp = r[ind:*]

;   radi_tmpInd = where(dtmp eq max(dtmp)) & radi_Ind = radi_tmpInd[0]
   radi_tmpInd = where(dd eq 1.67e-24) & radi_Ind = radi_tmpInd[0]

;   radi[i] = rtmp[radi_Ind]
   if (radi_ind eq -1) then break else radi[i] = r[radi_Ind]

;   plot,r,dd,/ylog,position=[x0,y0+2*ny,x1,y0+3*ny],/dev,xtickformat='(a1)',ytitle='density',/xst,yra=[5.e-26,1.e-16],/yst
;   plot,r,dd,/ylog,position=[x0,y0+2*ny,x1,y0+3*ny],/dev,xtickformat='(a1)',ytitle='density',/xst,yra=[1.e-27,1.e-23],/yst
   plot,r,dd,/ylog,position=[x0,y0+2*ny,x1,y0+3*ny],/dev,xtickformat='(a1)',ytitle='density',/xst,yra=[1.e-30,1.e-23],/yst
   oplot,[radi[i],radi[i]],10.^!y.crange,line=2
   legend,/right,/top,string(time[i]/60./60./24./365.,format='(F9.2)')+' yr',color=0,box=0
   plot,r,p[*,i],/ylog,position=[x0,y0+ny,x1,y0+2*ny],/dev,xtickformat='(a1)',/noerase,ytitle='pressure',/xst,yra=[1.e-16,1.e-8],/yst
;   plot,r,p[*,i],/ylog,position=[x0,y0+ny,x1,y0+2*ny],/dev,xtickformat='(a1)',/noerase,ytitle='pressure',/xst,yra=[1.e-16,1.e-7],/yst
   oplot,[radi[i],radi[i]],10.^!y.crange,line=2
;   plot,r,v[*,i],position=[x0,y0,x1,y0+ny],/dev,/noerase,ytitle='Vr',xtitle='r [cm]',/xst,yra=[-5.e5,1.2e7],/yst
;   plot,r,v[*,i],position=[x0,y0,x1,y0+ny],/dev,/noerase,ytitle='Vr',xtitle='r [cm]',/xst,yra=[-5.e5,1.2e8],/yst
;   plot,r,v[*,i],position=[x0,y0,x1,y0+ny],/dev,/noerase,ytitle='Vr',xtitle='r [cm]',/xst,yra=[-5.e5,1.2e9],/yst
   plot,r,v[*,i],position=[x0,y0,x1,y0+ny],/dev,/noerase,ytitle='Vr',xtitle='r [cm]',/xst,yra=[-5.e5,1.2e10],/yst
   oplot,[radi[i],radi[i]],!y.crange,line=2

if keyword_set(savepng) then draw,'png_radi/radi_'+string(i,format='(I4.4)')
endfor

radi = radi[where(radi ne 0)] 
time = time[where(radi ne 0)] 

save,file='pulsar_radi_'+title+'.sav',time,radi

window,1,xs=800,ys=600
plot,time,radi,/xlog,/ylog
stop
end

pro pulsar_radi_comb

restore,file='/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/v7/pulsar_radi_v7.sav'
t_v7 = time & r_v7 = radi
restore,file='/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/v8/pulsar_radi_v8.sav'
t_v8 = time & r_v8 = radi
restore,file='/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/v9/pulsar_radi_v9.sav'
t_v9 = time & r_v9 = radi
restore,file='/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/v10/pulsar_radi_v10.sav'
t_v10 = time & r_v10 = radi

loadct,39,/sil
!p.background=255 & !p.color=0
window,1,xs=800,ys=600
plot,t_v7,r_v7,/xlog,/ylog,xtitle='time [s]',ytitle='radius [cm]',xra=[1.e7,1.e11],/xst,yra=[3.e15,2.e18]
oplot,t_v7,10^(3./5.*alog10(t_v7/t_v7[1])+alog10(r_v7[1])),line=2

oplot,t_v8,r_v8,color=50
oplot,t_v8,10^(3./5.*alog10(t_v8/t_v8[1])+alog10(r_v8[1])),line=2,color=50
oplot,t_v9,r_v9,color=150
oplot,t_v9,10^(3./5.*alog10(t_v9/t_v9[1])+alog10(r_v9[1])),line=2,color=150
oplot,t_v10,r_v10,color=250
oplot,t_v10,10^(3./5.*alog10(t_v10/t_v10[1])+alog10(r_v10[1])),line=2,color=250

legend,['v10','v9','v8','v7'],line=0,color=[250,150,50,0],textcolor=[250,150,50,0],box=0,/left,/top
stop
end
