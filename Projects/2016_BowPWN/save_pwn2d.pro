pro save_pwn2d,step=step,startn=startn, endn=endn

if not keyword_set(step) then step=1
if not keyword_set(startn) then startn=0

fnames=file_search('PWN2d_hdf5_plt_cnt_*')
nfiles = n_elements(fnames)

if not keyword_set(endn) then endn=nfiles-1

tmp = loaddata(fnames[0],'dens',xCoord=x,yCoord=y,time=tt,xra=[0.,1.e18],yra=[-1.e18,1.e18])
sz  = size(tmp,/dimension)

ndata = (endn-startn+1)/step
time = fltarr(ndata)
d    = fltarr(sz[0],ndata)
p    = fltarr(sz[0],ndata)
v    = fltarr(sz[0],ndata)

y0ind_tmp = where(y ge 0) & y0ind = y0ind_tmp[0]

k=0
for i=startn,endn,step do begin
   print,'*************************************** ',i,' of ',endn
   dd = loaddata(fnames[i],'dens',xCoord=x,yCoord=y,time=tt,xra=[0.,1.e18],yra=[-1.e18,1.e18])
   pp = loaddata(fnames[i],'pres',xra=[0.,1.e18],yra=[-1.e18,1.e18])
   vv = loaddata(fnames[i],'velx',xra=[0.,1.e18],yra=[-1.e18,1.e18])

   d[*,k] = dd[*,y0ind]
   p[*,k] = pp[*,y0ind]
   v[*,k] = vv[*,y0ind]

   time[k] = tt

   k=k+1
endfor

r=x
save,file='save_pwn2d_v9.sav',r,d,p,v,time
stop
end
