pro line1d,n, nowindow=nowindow,d0=d0,v0=v0,init=init,nostop=nostop

if (n_elements(n) eq 0) then n=0

fname = 'Line1D_hdf5_plt_cnt_'+string(n,format='(I4.4)')

d = loaddata(fname,'dens',xCoord=x,time=t)
v = loaddata(fname,'velx')

if not keyword_set(nowindow) then begin
   d0 = d
   v0 = v
endif 

if keyword_set(init) then begin
   d0 = loaddata('Line1D_hdf5_plt_cnt_0000','dens',xCoord=x)
   v0 = loaddata('Line1D_hdf5_plt_cnt_0000','velx')
endif

xplt = 700 & yplt = 350
x0   = 100 & y0   = 100
loadct,39,/sil
!p.background=255 & !p.color=0
if not keyword_set(nowindow) then $
   window,0,xs=x0+xplt+110,ys=y0+2.*yplt+30

plot,x,d,/ylog,pos=[x0,y0+yplt,x0+xplt,y0+2*yplt],/dev,/xst,yst=8,xtickformat='(a1)',ytitle='density'
oplot,x,d0,line=2
legend,strtrim(n,0),/right,/top,box=0

plot,x,v0,pos=[x0,y0+yplt,x0+xplt,y0+2*yplt],/dev,yst=4,color=50,/noerase,xst=5,line=2,yra=[0.,max([v,v0])]
oplot,x,v,color=50
axis,yaxis=1,color=50,ytitle='velocity [cm/s]'
plot,x,4.*!pi*d*v*x*x,pos=[x0,y0,x0+xplt,y0+yplt],/dev,/xst,/noerase,xtitle='r [cm]',ytitle='Mass loss rate'
oplot,x,4.*!pi*d0*v0*x*x,line=2

if not keyword_set(nostop) then stop
end

pro play_line1d,skip=skip

if not keyword_set(skip) then skip=1
fnames = file_search('Line1D_hdf5_plt_cnt*')

n = strmid(fnames,3,4,/reverse)

nfiles= n_elements(fnames)

nowindow = 0
d0=0. & v0=0.
for i=0,nfiles-1,skip do begin
   line1d,i,/nostop,nowindow=nowindow,d0=d0,v0=v0
   nowindow=1
endfor
end

pro monitor_line1d, fname=fname, update=update

if not keyword_set(update) then update=1.

nfiles_old = 0
nowindow=0

true = 1
while (true) do begin
fnames = file_search('Line1D_hdf5_plt_cnt*')
nfiles = n_elements(fnames)

if (nfiles ne nfiles_old) then begin
  i = fix(strmid(fnames[nfiles-1],3,4,/rev))
  line1d,i,/nostop,/init,nowindow=nowindow
  nowindow=1
endif
nfiles_old = nfiles

print,nfiles
wait,update
endwhile

stop
end

