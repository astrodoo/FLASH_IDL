pro relbub

fname = 'relbub_rel0.99c.sav'
restore,file=fname
print,'done to restore the data'

;d0 = 1.67e-24
d0 = 1.67e-4    ;for relativistic case 

sz = size(d,/dimension)
shIn  = fltarr(sz[2])
shOut = fltarr(sz[2])

xcind = (where(x ge 0))[0]
for i=0, sz[2]-1 do begin
   shIn_ind = (where(d[xcind,*,i] ge d0))[0]
   shIn[i] = y[shIn_ind]
   shOut_ind = (where(d[xcind,*,i] eq d0))[0]
   shOut[i] = y[shOut_ind]
endfor

out = strmid(fname,0,strlen(fname)-4)+'_bubble.sav'
window,0
plot,time,shIn
oplot,time,shOut

save,file=out,shIn,shOut,time

stop
end

pro bubcomb

restore,file='data/relbub_norelv9_bubble.sav'
nrv9_shIn = shIn & nrv9_shOut = shOut & nrv9_t = time
restore,file='data/relbub_norelv10_bubble.sav'
nrv10_shIn = shIn & nrv10_shOut = shOut & nrv10_t = time
restore,file='data/relbub_norel0.6c_bubble.sav'
nr6c_shIn = shIn & nr6c_shOut = shout & nr6c_t = time
restore,file='data/relbub_norel0.9c_bubble.sav'
nr9c_shIn = shIn & nr9c_shOut = shOut & nr9c_t = time
restore,file='data/relbub_norel0.99c_bubble.sav'
nr99c_shIn = shIn & nr99c_shOut = shOut & nr99c_t = time

L = 1.d13
restore,file='data/relbub_rel0.33c_bubble.sav'
r33c_shIn = shIn*L & r33c_shOut = shout*L & r33c_t = time * L/!unit.c
restore,file='data/relbub_rel0.6c_bubble.sav'
r6c_shIn = shIn*L & r6c_shOut = shout*L & r6c_t = time * L/!unit.c
restore,file='data/relbub_rel0.9c_bubble.sav'
r9c_shIn = shIn*L & r9c_shOut = shout*L & r9c_t = time * L/!unit.c
restore,file='data/relbub_rel0.99c_bubble.sav'
r99c_shIn = shIn*L & r99c_shOut = shout*L & r99c_t = time * L/!unit.c


ana_t = nrv9_t
Edot = 6.66d35
d0 = 1.67d-24
;ana_sh = (25./14./!pi)^0.2*(Edot/d0)^0.2*ana_t^0.6
ana_sh = (125./154./!pi)^0.2*(Edot/d0)^0.2*ana_t^0.6

loadct,39,/sil
mkeps,'relbub_comb',xs=20,ys=20.*6./8.


plot,nrv9_t,nrv9_shIn,xra=[1.e9,5.e9],yra=[1.e17,7.e17],/xst,/yst,/nodata,xtitle='time [s]',ytitle='radius of bubble [cm]',/xlog,/ylog

;oplot,nrv9_t,nrv9_shIn
oplot,nrv10_t,nrv10_shIn
;oplot,nrv10_t,nrv10_shOut
oplot,nr6c_t,nr6c_shIn,color=50
;oplot,nr6c_t,nr6c_shOut
oplot,nr9c_t,nr9c_shIn,color=150
;oplot,nr9c_t,nr9c_shOut
oplot,nr99c_t,nr99c_shIn,color=210
;oplot,nr99c_t,nr99c_shOut

r033c=1.05
oplot,r33c_t,r33c_shIn*r033c,line=2
;oplot,r33c_t,r33c_shOut*r033c
r06c=1.02
oplot,r6c_t,r6c_shIn*r06c,line=2,color=50
;oplot,r6c_t,r6c_shOut*r06c
r09c=0.9
oplot,r9c_t,r9c_shIn*r09c,line=2,color=150
;oplot,r9c_t,r9c_shOut*r09c
r099c=0.8
oplot,r99c_t,r99c_shIn*r099c,line=2,color=210
;oplot,r99c_t,r99c_shOut*r099c

oplot,ana_t,ana_sh,color=254,thick=10

strvw = textoidl('v_{wind}=')
legend,strvw+['0.33c','0.6c','0.9c','0.99c'],color=[0,50,150,210],textcolor=[0,50,150,210],/left,/top,box=0
legend,['analytic solution','non-relativistic algorithm','relativistic algoritm'],color=[254,0,0] $
      ,textcolor=[254,0,0],/right,/bottom,box=0,line=[0,0,2],thick=[10,2,2]
epsfree
stop
end









pro mkdata

id = 'rel0.99c'
;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/tst_pulsar/v0.6c/'
dir = '/d/d7/yoon/out_FLASH2.5_rhd/out_PWNe/tst_rhd_PWN2d/tst_pulsar/v0.99c/denmod1st/'

fname = file_search(dir+'*plt_cnt*')
nfile= n_elements(fname)
sample=2

dd=loaddata(fname[0],'dens',sample=sample,xCoord=y,yCoord=x,time=t)
dd=transpose(dd)

sz = size(dd,/dimension)
d = fltarr(sz[0],sz[1],nfile)
time = fltarr(nfile)
d[*,*,0] = dd
time[0] = t

for i=1,nfile-1 do begin
   print, i,'  of  ', nfile
   dd=loaddata(fname[i],'dens',sample=sample,time=t) 
   dd = transpose(dd)
   d[*,*,i]=dd
   time[i]=t
endfor

save,file='relbub_'+id+'.sav',d,time,x,y

stop
end
