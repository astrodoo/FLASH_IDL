pro recoll_snap,mkdata=mkdata

fl = '/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/M30lateral0.01/JetSet_hdf5_plt_cnt_0439'

dx = 9.6e11/2.
bhx = -2.e12
xrange = [bhx-dx,bhx+dx]
yrange = [-5.e11,5.e11]
zrange = [0.,1.6e12]
sample=3 

dxz = reform(loaddata(fl,'dens',xra=xrange,yra=yrange,zra=zrange,time=time,sample=sample,xcoords=x,ycoords=y,zcoords=z)
x = x -bhx



stop
end
