pro recoll,n


fname = 'JetSet_hdf5_plt_cnt_'+string(n,format='(I4.4)')
print,'reading file: ',fname
outname = 'recoll_M10_'+string(n,format='(I4.4)')

xrange = [-2.5e12,-1.5e12]
yrange = [-5.e11,5.e11]
zrange = [-1.e11,1.5e12]
sample=2

mkdata=0
if (mkdata) then begin
   d = loaddata(fname,'dens',xra=xrange,yra=yrange $
       ,zra=zrange,time=time,sample=sample,xcoords=x,ycoords=y,zcoords=z)
   vz = loaddata(fname,'velz',xra=xrange,yra=yrange $
       ,zra=zrange,sample=sample)
   j = loaddata(fname,'jet',xra=xrange,yra=yrange $
       ,zra=zrange,sample=sample)
   save, file=outname+'.sav', d,vz,j,x,y,z,time
endif else restore, file=outname+'.sav'      ;mkdata

xjind = where(x ge -2.e12) & xjind = reform(xjind[0])
yjind = where(y ge 0.) & yjind = reform(yjind[0])
zjind = where(z ge 0) & zjind = reform(zjind[0])

sz = size(d,/dimension)
loadct,0,/sil
window,0,xs=sz[0]*2+sz[1],ys=sz[2]
tvscl,alog(reform(d[*,sz[1]/2,*])),0
xyouts,10,sz[2]-20,/dev,strmid(outname,7,3)
xyouts,10,sz[2]-50,/dev,'density in X-Z'

loadct,39,/sil
nl = 10 ;30
lcolor = findgen(nl)/float(nl)*256.
zind = zjind+1
for i=0,nl-1 do begin
   zind = [zind, zjind+(i+1)*20]
   plots,[0,511],[zind[i],zind[i]],/dev,thick=1,color=lcolor[i]
endfor

; measuring the thickness
jj = vz*j
jcrit = 2.8e9
jj2 = jj > jcrit

jxl = fltarr(sz[2]-zind[0]) & jxr = fltarr(sz[2]-zind[0])
jyl = fltarr(sz[2]-zind[0]) & jyr = fltarr(sz[2]-zind[0])
zj  = z[zind[0]:*]
for i=zind[0],sz[2]-1 do begin
   jx = reform(jj[*,yjind,i])  
   jxind = where(jx ge jcrit,count)
   if (count ne 0) then begin
      jxl[i-zind[0]] = x[jxind[0]] 
      jxr[i-zind[0]] = x[jxind[n_elements(jxind)-1]]
   endif else begin
      jxl[i-zind[0]] = 0.
      jxr[i-zind[0]] = 0.
   endelse
   jy = reform(jj[xjind,*,i])  
   jyind = where(jy ge jcrit,count)
   if (count ne 0) then begin
      jyl[i-zind[0]] = y[jyind[0]] 
      jyr[i-zind[0]] = y[jyind[n_elements(jyind)-1]]
   endif else begin
      jyl[i-zind[0]] = 0.
      jyr[i-zind[0]] = 0.
   endelse
endfor

save,file=outname+'_jet.sav', zj,jxl,jxr,jyl,jyr

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

mkeps,outname,xs=20.,ys=20.
contour, reform(jj[*,*,zind[0]]),x,y,/iso,xra=[-2.1e12,-1.9e12],yra=[-1.e11,1.e11],/xst,/yst,/nodata,levels=jcrit, xtitle='x [cm]', ytitle='y [cm]'
oplot,!x.crange,[0.,0.],line=1
oplot,[-2.e12,-2.e12],!y.crange,line=1
for i=0,nl-1 do contour, reform(jj[*,*,zind[i]]),x,y,levels=jcrit,/overplot,color=lcolor[i]
epsfree

stop
end

pro comp_jet
restore,file='recoll_M10_0319_jet.sav'
zj10 = zj & jyl10 = jyl & jyr10 = jyr

jyl10[where(jyl10 eq 0)] = !values.f_nan
jyr10[where(jyr10 eq 0)] = !values.f_nan

restore,file='recoll_M30_0365_jet.sav'
zj30 = zj & jyl30 = jyl & jyr30 = jyr

jyl30[where(jyl30 eq 0)] = !values.f_nan
jyr30[where(jyr30 eq 0)] = !values.f_nan

loadct,39,/sil
mkeps,'recoll_compjet',xs=20.,ys=20.
plot,jyl10,zj10,xr=[-5.e11,5.e11],yr=[0.,1.e12],/xst,/yst,/iso,xtitle='y [cm]',ytitle='z [cm]'
oplot,jyr10,zj10
oplot,jyl30,zj30,color=50
oplot,jyr30,zj30,color=50
legend,['M10','M30'],line=0,color=[0,50],textcolor=[0,50],/left,/top,box=0
epsfree

stop
end
