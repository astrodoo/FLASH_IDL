pro jettst,n,var=var,dx=dx
device,decomposed=0

fname = 'JetSet_hdf5_plt_cnt_'+string(n,format='(I4.4)')

if not keyword_set(var) then var='dens'
if not keyword_set(dx) then dx=2.e11

d = loaddata(fname,var,xc=0,yc=0,zc=0,xra=[-dx,dx],yra=[-dx,dx],zra=[-dx,dx],time=time,sample=0)

sd = size(d)
dd = reform(d[*,sd[2]/2,*])
if (sd[1] gt 600) then begin 
     xs = sd[1] & ys=sd[3] 
     dd2 = dd
endif else begin
     xs = 600 & ys = 600
     dd2 = congrid(dd, 600,600) 
endelse

window,0,xs=xs,ys=ys
tvscl,alog(dd2)
stop
end
