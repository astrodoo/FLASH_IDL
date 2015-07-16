pro play_lmxb
dir = '/d/d2/yoon/out_Flash_jet/out_LMXB/out_jet_xrbs_wind_100/den2dcut/'
restore, dir+'den2d_all.sav'

szd = size(dens,/dimension)

xc = 6.e19
yc = 9.874e20
mind = min(dens[*,*,szd[2]-1]) & maxd = max(dens[*,*,szd[2]-1])

pltx0=100. & plty0=80.
winxs=pltx0+szd[0]+100. & winys=plty0+szd[1]+20.
loadct,3,/sil

spawn,'mkdir '+dir+'png_data'

!p.background=255 & !p.color=0
for i=0,szd[2]-1 do begin
   print,i, ' of ', szd[2] 
   window,0,xs= winxs,ys=winys;,/pixmap
   tvcoord,bytscl(alog10(dens[*,*,i]),min=alog10(mind),max=alog10(maxd)),x,y,/axes,pos=[pltx0,plty0],/dev,xtitle='x [cm]', ytitle='y [cm]',xtickinterval=2.e20, ytickinterval=2.e20,/black
   plots,xc,yc,psym=7,color=255
   legend,'t = '+string(time[i]/60./60./24./365./1.e3,format='(I4)')+' kyr',/right,/top,box=0,textcolor=0
   color_bar,lim=[mind,maxd],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),titlegap=0.12
   snapshot,dir+'png_data/lmxb_'+string(i,format='(I4.4)')
endfor

stop
end



pro den2dcut

dir = '/d/d2/yoon/out_Flash_jet/out_LMXB/out_jet_xrbs_wind_100/'
fname = 'jet_xrbs_hdf5_plt_cnt_'

spawn,'mkdir '+dir+'den2dcut'

fnames = file_search(dir+fname+'*')
nfiles = n_elements(fnames)
id = strmid(fnames,3,4,/rev)

xc = 6.e19
yc = 9.874e20
zc = 9.874e20

for i=0, nfiles-1 do begin
    print,i,' of ',nfiles
    den=loaddata(fnames[i],'dens',xCoord=x,yCoord=y,time=t,sample=2,xra=[1.e18,7.e20],yra=[yc-3.e20,yc+3.e20],zra=[zc,zc])
    save,file=dir+'den2dcut/den2d_'+id[i]+'.sav',den,x,y,t
endfor

stop
end

pro merge2d
dir = '/d/d2/yoon/out_Flash_jet/out_LMXB/out_jet_xrbs_wind_100/den2dcut/'
fname = 'den2d_'
fnames = file_search(dir+fname+'*')
nfiles= n_elements(fnames)

restore,file=fnames[0]
szd = size(den,/dimension)

dens = fltarr(szd[0],szd[1],nfiles)
time = fltarr(nfiles)
for i=0, nfiles-1 do begin
    print,i,' of ',nfiles
    restore,file=fnames[i]
    dens[*,*,i] = den[*,*]
    time[i] = t
endfor

save,file=dir+'den2d_all.sav',dens,time,x,y

stop
end
