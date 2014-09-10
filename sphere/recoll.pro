pro recoll,n


fname = 'JetSet_hdf5_plt_cnt_'+string(n,format='(I4.4)')
print,'reading file: ',fname
outname = 'recoll_'+string(n,format='(I4.4)')

xrange = [-2.5e12,-1.5e12]
yrange = [-5.e11,5.e11]
zrange = [-1.e11,1.e12]
sample=2

mkdata=0
if (mkdata) then begin
   d = loaddata(fname,'dens',xra=xrange,yra=yrange $
       ,zra=zrange,time=time,sample=sample,xcoords=x,ycoords=y,zcoords=z)
   vz = loaddata(fname,'velz',xra=xrange,yra=yrange $
       ,zra=zrange,sample=sample)
   j = loaddata(fname,'jet',xra=xrange,yra=yrange $
       ,zra=zrange,sample=sample)
   save, file=outname, d,vz,j,x,y,z,time
endif else restore, file=outname      ;mkdata

zjind = where(z ge 0) & zjind = reform(zjind[0])

sz = size(d,/dimension)
window,0,xs=sz[0],ys=sz[2]
tvscl,alog(reform(d[*,sz[1]/2,*]))

zind = zjind+1
for i=0,20 do begin
   zind = [zind, zjind+(i+1)*20]
   plots,[0,511],[zind[i],zind[i]],/dev,thick=1
endfor

jj = vz*j
window,1
contour, reform(jj[*,*,zind[0]]),x,y,/iso,xra=[-2.2e12,-1.8e12],yra=[-2.e11,2.e11],/xst,/yst,/nodata,levels=3.e8
for i=0,20 do contour, reform(jj[*,*,zind[i]]),x,y,levels=3.e8,/overplot


stop
end
