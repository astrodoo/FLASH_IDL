pro showPWN2d,n,var=var,sample=sample,block=block,xrange=xrange,yrange=yrange,ct=ct,nolog=nolog,chk=chk
device,decomposed=0

if keyword_set(chk) then begin
   fname = 'PWN2d_hdf5_chk_'+string(n,format='(I4.4)') 
   print, 'read check point file: '+fname   
endif else begin 
   fname = 'PWN2d_hdf5_plt_cnt_'+string(n,format='(I4.4)')
   chk=0
endelse

if not keyword_set(var) then var='dens'
if not keyword_set(sample) then sample=0
if not keyword_set(ct) then ct=0
if not keyword_set(nolog) then nolog=0

if not keyword_set(xrange) then str_xra = '' $
  else str_xra = 'yra=['+strtrim(xrange[0],2)+','+strtrim(xrange[1],1)+']'
if not keyword_set(yrange) then str_yra = '' $
  else str_yra = 'xra=['+strtrim(yrange[0],2)+','+strtrim(yrange[1],1)+']'

if (keyword_set(xrange) and keyword_set(yrange)) then str_yra = ','+str_yra

str_xyra = str_xra+str_yra
if (strmid(str_xyra,0,1,/reverse) eq ']') then str_xyra = str_xyra+','

if (where(strmatch(['velx','vely','velz','jet'],var) eq 1) ne -1) then nolog=1

read_amr,fname,var='dens',parameters=params,tree=tree ,/nodata
time = params.time

print,'time = ', time

strexe = execute("data = loaddata(fname,'"+var+"',"+str_xyra+"sample=sample,lref=lref,xCoord=y,yCoord=x,time=time)")

data=transpose(data)
tree.bndbox=reverse(tree.bndbox,2)

sz = size(data,/dimension)

pltx0=100. & plty0=60.
winxs=pltx0+sz[0]+20 & winys=plty0+sz[1]+100

scrsz = get_screen_size()

loadct,ct,/sil
if (((winxs+100) ge scrsz[0]) or ((winys+100) ge scrsz[1])) then  $
   swindow, xs=winxs, ys=winys $
  else window, xs=winxs, ys=winys

; draw slice cuts
;maxd = max(data) & mind = min(data)
maxd = 5.d-25 & mind=1.d-26
if not keyword_set(nolog) then begin
   tvcoord, bytscl(alog10(data),min=alog10(mind),max=alog10(maxd)),x,y,pos=[pltx0,plty0,pltx0+sz[0],plty0+sz[1]],/dev,/axes,xtitle='x [cm]',ytitle='y [cm]'
   if keyword_set(block) then showblock,tree=tree,param=params
   color_bar,lim=[mind,maxd],/log,/up,bartitle=var,pos=[pltx0, plty0+sz[1]+10, pltx0+sz[0], plty0+sz[1]+30],titlegap=0.04
endif else begin
   tvcoord, bytscl(data,min=mind,max=maxd),x,y,pos=[pltx0,plty0,pltx0+sz[0],plty0+sz[1]],/dev,/axes,xtitle='x [cm]',ytitle='y [cm]'
   if keyword_set(block) then showblock,tree=tree,param=params
   color_bar,lim=[mind,maxd],/up,bartitle=var,pos=[pltx0, plty0+sz[1]+10, pltx0+sz[0], plty0+sz[1]+30],titlegap=0.04
endelse

stop
end



