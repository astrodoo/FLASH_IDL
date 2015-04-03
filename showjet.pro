pro showjet,n,var=var,sample=sample,block=block,xc0=xc0,yc0=yc0,zc0=zc0,chk=chk,nodload=nodload  $
           ,xrange=xrange,yrange=yrange,zrange=zrange, ct=ct, nolog=nolog
device,decomposed=0

if keyword_set(chk) then begin
   fname = 'JetSet_hdf5_chk_'+string(n,format='(I4.4)') 
   print, 'read check point file: '+fname   
endif else begin 
   fname = 'JetSet_hdf5_plt_cnt_'+string(n,format='(I4.4)')
   chk=0
endelse

if not keyword_set(var) then var='dens'
if not keyword_set(xc0) then xc0 = -2.e12
if not keyword_set(yc0) then yc0 = 0.
if not keyword_set(zc0) then zc0 = 0.
if not keyword_set(sample) then sample=0
if not keyword_set(ct) then ct=0
if not keyword_set(nolog) then nolog=0

if not keyword_set(xrange) then str_xra = '' $
  else str_xra = 'xra=['+strtrim(xrange[0],2)+','+strtrim(xrange[1],1)+']'
if not keyword_set(yrange) then str_yra = '' $
  else str_yra = 'yra=['+strtrim(yrange[0],2)+','+strtrim(yrange[1],1)+']'
if not keyword_set(zrange) then str_zra = '' $
  else str_zra = 'zra=['+strtrim(zrange[0],2)+','+strtrim(zrange[1],1)+']'

if keyword_set(xrange) then begin
   if keyword_set(yrange) then yrange = ','+yrange
   if keyword_set(zrange) then zrange = ','+zrange
endif else if (keyword_seT(yrange) and keyword_set(zrange)) then zrange = ','+zrange 

str_xyzra = str_xra+str_yra+str_zra
if (strmid(str_xyzra,0,1,/reverse) eq ']') then str_xyzra = str_xyzra+','

if (where(strmatch(['velx','vely','velz','jet'],var) eq 1) ne -1) then nolog=1 

read_amr,fname,var='dens',parameters=params,tree=tree ,/nodata
time = params.time

print,'time = ', time

if (keyword_set(nodload) or keyword_set(xrange) or keyword_set(yrange) or keyword_set(zrange)) then $
   strexe = execute("data = loaddata(fname,'"+var+"',"+str_xyzra+"sample=sample,lref=lref,xCoord=x,yCoord=y,zCoord=z,time=time)") $
  else data = dload(n,var=var,xc=xc0,yc=yc0,zc=zc0,x,y,z,time,sample=sample,lref=lref,chk=chk)


sz = size(data,/dimension)

pltx0=100 & plty0=80
winxs = pltx0+sz[0]+sz[1]+30 & winys = plty0+sz[2]+sz[1]+30

scrsz = get_screen_size()

loadct,ct,/sil
if (((winxs-100) ge scrsz[0]) or ((winys-100) ge scrsz[1])) then  $
   swindow, xs=winxs, ys=winys $
  else window, xs=winxs, ys=winys


xind = (where(x ge xc0))[0] & yind=(where(y ge yc0))[0] & zind=(where(z ge zc0))[0]
dxy = reform(data[*,*,zind])
dyz = reform(data[xind,*,*])
dxz = reform(data[*,yind,*])

; draw slice cuts
maxd = max(data) & mind = min(data)
if not keyword_set(nolog) then begin
   tvcoord, bytscl(alog10(dxz),min=alog10(mind),max=alog10(maxd)),x,z,pos=[pltx0,plty0,pltx0+sz[0],plty0+sz[2]],/dev,/axes,xtitle='x [cm]',ytitle='z [cm]'
   if keyword_set(block) then showblock,tree=tree,param=params,yc=yc0 
   tvcoord, bytscl(alog10(dxy),min=alog10(mind),max=alog10(maxd)),x,y,pos=[pltx0,plty0+sz[2],pltx0+sz[0],plty0+sz[1]+sz[2]],/dev,/axes,ytitle='y [cm]',xtickformat='(a1)'
   if keyword_set(block) then showblock,tree=tree,param=params,zc=zc0 
   tvcoord, bytscl(alog10(dyz),min=alog10(mind),max=alog10(maxd)),y,z,pos=[pltx0+sz[0],plty0,pltx0+sz[0]+sz[1],plty0+sz[2]],/dev,/axes,xtitle='y [cm]',ytickformat='(a1)'
   if keyword_set(block) then showblock,tree=tree,param=params,xc=xc0 
   color_bar,lim=[mind,maxd],/log,/up,bartitle=var,pos=[pltx0+sz[0]+20, plty0+sz[2], pltx0+sz[0]+sz[1], plty0+sz[2]+20],titlegap=30.
endif else begin
   tvcoord, bytscl(dxz,min=mind,max=maxd),x,z,pos=[pltx0,plty0,pltx0+sz[0],plty0+sz[2]],/dev,/axes,xtitle='x [cm]',ytitle='z [cm]'
   if keyword_set(block) then showblock,tree=tree,param=params,yc=yc0 
   tvcoord, bytscl(dxy,min=mind,max=maxd),x,y,pos=[pltx0,plty0+sz[2],pltx0+sz[0],plty0+sz[1]+sz[2]],/dev,/axes,ytitle='y [cm]',xtickformat='(a1)'
   if keyword_set(block) then showblock,tree=tree,param=params,zc=zc0 
   tvcoord, bytscl(dyz,min=mind,max=maxd),y,z,pos=[pltx0+sz[0],plty0,pltx0+sz[0]+sz[1],plty0+sz[2]],/dev,/axes,xtitle='y [cm]',ytickformat='(a1)'
   if keyword_set(block) then showblock,tree=tree,param=params,xc=xc0 
   color_bar,lim=[mind,maxd],/up,bartitle=var,pos=[pltx0+sz[0]+20, plty0+sz[2], pltx0+sz[0]+sz[1], plty0+sz[2]+20],titlegap=30.
endelse

loadct,0,/sil
xyouts, pltx0+sz[0]+10, plty0+sz[2]+sz[1]-20,/dev,fname
xyouts, pltx0+sz[0]+10, plty0+sz[2]+sz[1]-40,/dev,'time = '+string(time/60./60.,format='(f5.2)')+'hrs'


stop
end
