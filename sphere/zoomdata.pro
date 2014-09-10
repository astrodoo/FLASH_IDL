pro zoomdata,fname,xrange=xrange,yrange=yrange,zrange=zrange,sample=sample

dx = 4.e11
xrange = [-2.e12-dx,-2.e12+dx]
yrange = [-dx,dx]
zrange = [-dx,dx]

if not keyword_set(sample) then sample=0

dens=loaddata(fname,'dens',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,time=time,xCoord=x,yCoord=y,zCoord=z)
pres=loaddata(fname,'pres',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)
velz=loaddata(fname,'velz',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)
jet =loaddata(fname,'jet',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)

save,filename='zoomdata_'+strmid(fname,3,4,/reverse)+'.dat', time,x,y,z,dens,pres,velz,jet

stop
end
