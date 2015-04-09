pro mkddxz,n,sample=sample, xrange=xrange, yrange=yrange, zrange=zrange, xc0=xc0, yc0=yc0, zc0=zc0, jeton=jeton
device,decomposed=0

if not keyword_set(sample) then sample=3

if not keyword_set(xrange) then xrange=[-1.4e13,3.5e12]
if not keyword_set(yrange) then yrange=[-5.e12,5.e12]
if not keyword_set(zrange) then zrange=[-1.2e13,1.2e13]
if not keyword_set(xc0) then xc0 = -2.e12
if not keyword_set(yc0) then yc0 = 0.
if not keyword_set(zc0) then zc0 = 0.
if not keyword_set(jeton) then jeton = 9.76e4

fout = 'ddxz_'+string(n,format='(I4.4)')+ '.sav'

fname = 'JetSet_hdf5_plt_cnt_'+string(n,format='(I4.4)')

d = loaddata(fname,'dens',xrange=xrange,yrange=yrange $
   ,zrange=zrange,time=time,sample=sample,xcoords=x,ycoords=y,zcoords=z)

time = time - jetOn

sd = size(d,/dimension)

indy_tmp = where(y ge yc0) & indy = indy_tmp[0]
indz_tmp = where(z ge zc0) & indz = indz_tmp[0]

ddxz = reform(d[*,indy,*])

ddxz_dummy = reform(ddxz[*,indz:*])
z_dummy = z[indz:*]

ddxz_cp = fltarr(sd[0],2*(sd[2]-indz))
z_cp = fltarr(2*(sd[2]-indz))

ddxz_cp[*,sd[2]-indz:*] = ddxz_dummy
ddxz_cp[*,0:sd[2]-indz-1] = reverse(ddxz_dummy,2)
z_cp[sd[2]-indz:*] = z_dummy
z_cp[0:sd[2]-indz-1] = -reverse(z_dummy)
sd2 = size(ddxz_cp,/dimension)

window,xs=sd2[0],ys=sd2[1]
tvscl,alog(ddxz_cp)

ddxz = ddxz_cp
z    = z_cp
save,file=fout,ddxz, x, z, time 
stop
end
