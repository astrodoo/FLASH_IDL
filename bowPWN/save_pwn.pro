pro save_pwn,step=step,startn=startn, endn=endn
;/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/600kms/v7
;/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/600kms/v8
;/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/600kms/v9
;/d1/yoon/outputs_FLASH3.3_mhd/out_ds_PWN1d/600kms/v10

if not keyword_set(step) then step=1
if not keyword_set(startn) then startn=0

fnames=file_search('PWN1d_hdf5_plt_cnt_*')
nfiles = n_elements(fnames)

if not keyword_set(endn) then endn=nfiles-1

tmp = loaddata(fnames[0],'dens',xCoord=r,time=tt)
nn  = n_elements(tmp)

ndata = (endn-startn+1)/step
time = fltarr(ndata)
d    = fltarr(nn,ndata)
p    = fltarr(nn,ndata)
v    = fltarr(nn,ndata)

k=0
for i=startn,endn,step do begin
   print,'*************************************** ',i,' of ',endn
   dd = loaddata(fnames[i],'dens',xCoord=x,time=tt)
   pp = loaddata(fnames[i],'pres')
   vv = loaddata(fnames[i],'velx')

   d[*,k] = dd
   p[*,k] = pp
   v[*,k] = vv

   time[k] = tt

   k=k+1
endfor

save,file='save_pwn_v10.sav',r,d,p,v,time
stop
end
