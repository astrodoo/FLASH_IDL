pro timechk,sample=sample

xcut=-2.e12
ycut=0.

time0 = systime(1)

d=loaddata('JetSet_hdf5_plt_cnt_0439','dens',sample=sample,xra=[-1.5e13,1.e12],yra=[-1.e13,1.e13],zra=[-3.e11,1.9e13],xcoord=x,ycoord=y,zcoord=z)
xind = (where(x ge xcut))[0]
yind = (where(y ge ycut))[0]

dyz = reform(d[xind,*,*])
dxz = reform(d[*,yind,*])

time1 = systime(1)

dyz=reform(loaddata('JetSet_hdf5_plt_cnt_0439','dens',sample=sample,xra=[xcut,xcut],yra=[-1.e13,1.e13],zra=[-3.e11,1.9e13],xcoord=x,ycoord=y,zcoord=z))
dxz=reform(loaddata('JetSet_hdf5_plt_cnt_0439','dens',sample=sample,xra=[-1.5e13,1.e12],yra=[ycut,ycut],zra=[-3.e11,1.9e13],xcoord=x,ycoord=y,zcoord=z))
y = reform(y)

time2 = systime(2)

print,'read 3d: ',time1-time0
print,'read 2d: ',time2-time1
stop
end
