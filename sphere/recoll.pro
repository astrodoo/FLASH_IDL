pro recoll,n,mkdata=mkdata

fname = 'JetSet_hdf5_plt_cnt_'+string(n,format='(I4.4)')
print,'reading file: ',fname
outname = 'recoll_M30_'+string(n,format='(I4.4)')

;xrange = [-2.5e12,-1.5e12]
;yrange = [-8.e11,8.e11]
;zrange = [-1.e11,2.e12]

xrange = [-2.5e12,-1.5e12]
yrange = [-5.e11,5.e11]
zrange = [-1.e11,2.e12]
;sample=3 ; for M10
sample=2 ; for M30

if keyword_set(mkdata) then begin
   d = loaddata(fname,'dens',xra=xrange,yra=yrange $
       ,zra=zrange,time=time,sample=sample,xcoords=x,ycoords=y,zcoords=z)
   vz = loaddata(fname,'velz',xra=xrange,yra=yrange $
       ,zra=zrange,sample=sample)
   j = loaddata(fname,'jet',xra=xrange,yra=yrange $
       ,zra=zrange,sample=sample)
   save, file=outname+'.sav', d,vz,j,x,y,z,time
endif else restore, file=outname+'.sav'      ;mkdata

xjind_tmp = where(x ge -2.e12) & xjind = xjind_tmp[0]
yjind_tmp = where(y ge 0.) & yjind = yjind_tmp[0]
zjind_tmp = where(z ge 0) & zjind = zjind_tmp[0]

sz = size(d,/dimension)
loadct,0,/sil
window,0,xs=sz[0]*2+sz[1],ys=sz[2]
tvscl,alog(reform(d[*,sz[1]/2,*])),0
xyouts,10,sz[2]-20,/dev,strmid(outname,7,3)
xyouts,10,sz[2]-50,/dev,'density in X-Z'

loadct,13,/sil
;nl = 30
;lcolor = findgen(nl)/float(nl)*256.
;zind = zjind+1
;for i=0,nl-1 do begin
;   zind = [zind, zjind+(i+1)*20]
;   plots,[0,511],[zind[i],zind[i]],/dev,thick=1,color=lcolor[i]
;endfor

; measuring the thickness
jj = vz*j
jcrit = 2.8e9
jj2 = jj > jcrit

jxl = fltarr(sz[2]-zjind) & jxr = fltarr(sz[2]-zjind)
jyl = fltarr(sz[2]-zjind) & jyr = fltarr(sz[2]-zjind)
zj  = z[zjind:*]
for i=zjind,sz[2]-1 do begin
   jx = reform(jj[*,yjind,i])  
   jxind = where(jx ge jcrit,count)
   if (count ne 0) then begin
      jxl[i-zjind] = x[jxind[0]] 
      jxr[i-zjind] = x[jxind[n_elements(jxind)-1]]
   endif else begin
      jxl[i-zjind] = 0.
      jxr[i-zjind] = 0.
   endelse
   jy = reform(jj[xjind,*,i])  
   jyind = where(jy ge jcrit,count)
   if (count ne 0) then begin
      jyl[i-zjind] = y[jyind[0]] 
      jyr[i-zjind] = y[jyind[n_elements(jyind)-1]]
   endif else begin
      jyl[i-zjind] = 0.
      jyr[i-zjind] = 0.
   endelse
endfor

save,file=outname+'_jet.sav', zj,jxl,jxr,jyl,jyr, time

tvlct,r,g,b,/get
tvlct,255,255,255,0
tvlct,0,0,0,1
tvcoord,reform(jj2[*,sz[1]/2,*]),x,z,pos=[sz[0],0],/scale
oplot,jxl,zj,color=1,thick=3
oplot,jxr,zj,color=1,thick=3
xyouts,sz[0]+10,sz[2]-50,/dev,'j*vz>crit in X-Z',color=254
tvcoord,reform(jj2[xjind,*,*]),y,z,pos=[sz[0]*2,0],/scale
oplot,jyl,zj,color=1,thick=3
oplot,jyr,zj,color=1,thick=3
xyouts,sz[0]*2+10,sz[2]-50,/dev,'j*vz>crit in Y-Z',color=254
tvlct,r,g,b

;mkeps,outname,xs=20.,ys=20.
;contour, reform(jj[*,*,zind[0]]),x,y,/iso,xra=[-2.1e12,-1.9e12],yra=[-1.e11,1.e11],/xst,/yst,/nodata,levels=jcrit, xtitle='x [cm]', ytitle='y [cm]'
;oplot,!x.crange,[0.,0.],line=1
;oplot,[-2.e12,-2.e12],!y.crange,line=1
;for i=0,nl-1 do contour, reform(jj[*,*,zind[i]]),x,y,levels=jcrit,/overplot,color=lcolor[i]
;epsfree

;stop
;=======
jyz = reform(jj[xjind,*,zjind:*])
yy = y
zz = z[zjind:*]
maxv = max(jyz)
minv = min(jyz)

nlev = 255
levs = findgen(nlev)/float(nlev)*(maxv-minv) + minv
mkeps, outname+'_colorbar',xs=20.,ys=20.*6./8.
contour,jyz,yy,zz,/iso,xra=[-5.e11,5.e11],yra=[0.,1.2e12],/xst,/yst,/fill,levels=levs,xtitle='x [cm]',ytitle='y [cm]'
color_bar,/right,lim=[minv,maxv],bartitle='jj'

epsfree
stop
end

pro comp_jet
dir = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/re-coll/M10lateral0.01/'
restore,file=dir+'recoll_M10_0500_jet.sav'
zj10 = zj[where((jyl ne 0) and (jyr ne 0))] & jyl10 = jyl[where((jyl ne 0) and (jyr ne 0))] & jyr10 = jyr[where((jyl ne 0) and (jyr ne 0))]

;jyl10[where(jyl10 eq 0)] = !values.f_nan
;jyr10[where(jyr10 eq 0)] = !values.f_nan

dir = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/re-coll/M30lateral0.01/'
restore,file=dir+'recoll_M30_0440_jet.sav'
zj30 = zj[where((jyl ne 0) and (jyr ne 0))] & jyl30 = jyl[where((jyl ne 0) and (jyr ne 0))] & jyr30 = jyr[where((jyl ne 0) and (jyr ne 0))]

;jyl30[where(jyl30 eq 0)] = !values.f_nan
;jyr30[where(jyr30 eq 0)] = !values.f_nan

loadct,39,/sil
recz = 7.e11
mkeps,'recoll_compjet',xs=20.,ys=20.*6./8.
plot,jyl10,zj10,xr=[-8.e11,8.e11],yr=[0.,1.3e12],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]'
;oplot,jyr10,zj10
oplot,-jyl10,zj10
;recz10 = 3.1e11
oplot,!x.crange,[recz,recz],line=1
und10_ind = where(zj10 le recz)
fit10 = linfit(jyl10[und10_ind],zj10[und10_ind])
oplot,[0.,-4.e11],[0.,fit10[1]*(-4.e11)],line=2
xyouts,2.e11,4.e11,/data, textoidl('\alpha = ') + string(90.-atan(abs(fit10[1]))/!dtor,format='(f5.2)') +' deg'

oplot,jyl30,zj30,color=50
oplot,-jyl30,zj30,color=50
;oplot,jyr30,zj30,color=50
;recz = 6.8e11
oplot,!x.crange,[recz,recz],line=1
und30_ind = where(zj30 le recz)
fit30 = linfit(jyl30[und30_ind],zj30[und30_ind])
oplot,[0.,-2.e11],[0.,fit30[1]*(-2.e11)],line=2,color=50
xyouts,2.e11,3.2e11,/data, textoidl('\alpha = ') + string(90.-atan(abs(fit30[1]))/!dtor,format='(f5.2)') +' deg', color=50
legend,['M10','M30'],line=0,color=[0,50],textcolor=[0,50],/right,/bottom,box=0
epsfree

;jthk10 = jyr10-jyl10
;jthk30 = jyr30-jyl30

jthk10 = 2*abs(jyl10)
jthk30 = 2*abs(jyl30)

mkeps,'recoll_compjet_thick',xs=20.,ys=20.*6./8.
;cutz10 = 4.1e11
cutz10 = 1.3e12
plot,zj10[where(zj10 le cutz10)],jthk10[where(zj10 le cutz10)], xra=[0.,1.e12],yra=[0.,5.e11], /xst,/yst,xtitle='z [cm]', ytitle='jet thickness [cm]',xtickinterval=4.e11
cutz30 = 1.3e12
oplot,zj30[where(zj30 le cutz30)],jthk30[where(zj30 le cutz30)],color=50
oplot,[recz,recz],!y.crange,line=1
legend,['M10','M30'],line=0,color=[0,50],textcolor=[0,50],/left,/top,box=0
epsfree

stop
end

pro comp_jet_time
;dir = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/re-coll/M10lateral0.01/'
dir = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/re-coll/M30lateral0.01/'

;fs = 'recoll_M10_*_jet.sav'
fs = 'recoll_M30_*_jet.sav'

fnames = file_search(dir+fs)

Mstr = strmid(fs,7,3)

nfile = n_elements(fnames)

max_nzj = 0
for i=0,nfile-1 do begin
   restore,fnames[i]
   nzj = n_elements(zj)
   if (nzj gt max_nzj) then max_nzj = nzj
endfor

jet = replicate({zj: fltarr(max_nzj), jyl: fltarr(max_nzj), jyr: fltarr(max_nzj), time: time},nfile)

for i=0,nfile-1 do begin
   restore,file=fnames[i]
   nzj = n_elements(zj)
   jet[i].zj[0:nzj-1] = zj
   jet[i].jyl[0:nzj-1] = jyl
   jet[i].jyr[0:nzj-1] = jyr
   jet[i].time = time
endfor

loadct,39,/sil

col = findgen(nfile) / float(nfile-1) *254

;window,0
for i =0,nfile-1 do begin
   jet[i].zj[where((jet[i].jyl eq 0) or (jet[i].jyr eq 0))]=!values.f_nan
   jet[i].jyl[where((jet[i].jyl eq 0) or (jet[i].jyr eq 0))]=!values.f_nan
   jet[i].jyr[where((jet[i].jyl eq 0) or (jet[i].jyr eq 0))]=!values.f_nan
endfor
mkeps,'recoll_compjet_time_'+Mstr,xs=20.,ys=20.*7./8.
plot,jet[0].jyl,jet[0].zj,xr=[-8.e11,8.e11],yr=[0.,1.5e12],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]', color=col[0]
oplot,jet[0].jyr,jet[0].zj, color= col[0]
for i=0,nfile-1 do begin
   oplot,jet[i].jyl,jet[i].zj, color=col[i]
   oplot,jet[i].jyr,jet[i].zj, color=col[i]
endfor

timestr = strtrim(jet[*].time,2)
legend,'time = ' + timestr,/right,/bottom,box=0, textcolor=col, color=col

legend,Mstr,/left,/top,box=0
epsfree
;recz = 6.8e11
;mkeps,'recoll_compjet',xs=20.,ys=20.
;plot,jyl10,zj10,xr=[-5.e11,5.e11],yr=[0.,1.e12],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]'
;oplot,jyr10,zj10
;recz10 = 3.1e11
;oplot,!x.crange,[recz,recz],line=1
;und10_ind = where(zj10 le recz)
;fit10 = linfit(jyr10[und10_ind],zj10[und10_ind])
;oplot,[0.,4.e11],[0.,fit10[1]*4.e11],line=2
;xyouts,1.5e11,3.5e11,/data, textoidl('\alpha = ') + string(90.-atan(fit10[1])/!dtor,format='(f5.2)') +' deg'

;oplot,jyl30,zj30,color=50
;oplot,jyr30,zj30,color=50
;recz = 6.8e11
;oplot,!x.crange,[recz,recz],line=1
;und30_ind = where(zj30 le recz)
;fit30 = linfit(jyr30[und30_ind],zj30[und30_ind])
;oplot,[0.,2.e11],[0.,fit30[1]*2.e11],line=2,color=50
;xyouts,1.5e11,4.e11,/data, textoidl('\alpha = ') + string(90.-atan(fit30[1])/!dtor,format='(f5.2)') +' deg', color=50
;legend,['M10','M30'],line=0,color=[0,50],textcolor=[0,50],/left,/bottom,box=0
;;epsfree
;
;jthk10 = jyr10-jyl10
;jthk30 = jyr30-jyl30
;
;mkeps,'recoll_compjet_thick',xs=20.,ys=20.*6./8.
;;cutz10 = 4.1e11
;cutz10 = 1.e12
;plot,zj10[where(zj10 le cutz10)],jthk10[where(zj10 le cutz10)], xra=[0.,1.e12],yra=[0.,5.e11], /xst,/yst,xtitle='z [cm]', ytitle='jet thickness [cm]',xtickinterval=4.e11
;cutz30 = 1.e12
;oplot,zj30[where(zj30 le cutz30)],jthk30[where(zj30 le cutz30)],color=50
;oplot,[recz,recz],!y.crange,line=1
;legend,['M10','M30'],line=0,color=[0,50],textcolor=[0,50],/left,/top,box=0
;epsfree

stop
end
