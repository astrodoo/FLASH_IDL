pro asywind, n=n, lout=l, ld=ld, lvr=lvr, lion=lion, time=time

fname = file_search('*plt_cnt_*')

nfile = n_elements(fname)
skip=10

;for i=0,nfile-1,skip
i=n
read_amr,fname[i],tree=tree,/nodata
maxlref = max(tree.lrefine)
n = fix(strmid(fname[i],3,4,/rev))

d=dload(n,var='dens',sample=maxlref-7,xc=0,yc=0,zc=0,x,y,z,time)
if arg_present(lvr) then begin
  vx=dload(n,var='velx',sample=maxlref-7,xc=0,yc=0,zc=0)
  vy=dload(n,var='vely',sample=maxlref-7,xc=0,yc=0,zc=0)
endif

dsz = size(d,/dimension)
print,dsz[2]/2
dxy = d[*,*,dsz[2]/2]

if arg_present(lvr) then begin
  vx_xy = vx[*,*,dsz[2]/2]
  vy_xy = vy[*,*,dsz[2]/2]
endif
;!p.background=255. & !p.color=0
;window,0,xs=1500,ys=800
loadct,0,/sil
window,xs=512,ys=512
tvcoord,alog(dxy),x,y,/scale

; find position of BB
; temporary calculation to circular orbit. (in case of ellipse, need more tasks)
m1 = 4.d34
m2 = 2.d34
mtot = m1+m2
peri = 3.d12
th0 = !pi
th = th0 + sqrt(!unit.g * mtot / peri^3.d0) * time

r = x[dsz[0]-1] - peri*m2/mtot

; read lines: two ends of the line
; center of mass, position of Star
cm2xSt = peri*m2/mtot*cos(th+!pi)
cm2ySt = peri*m2/mtot*sin(th+!pi)
; center of mass, position of BB
cm2xj = peri*m1/mtot*cos(th)
cm2yj = peri*m1/mtot*sin(th)

if arg_present(lvr) then begin
   x2dSt = (x-cm2xSt)#replicate(1.,dsz[1])
   y2dSt = replicate(1.,dsz[0])#(y-cm2ySt)
   rSt   = sqrt(x2dSt*x2dSt + y2dSt*y2dSt)
   vrxy  = vx_xy*x2dSt/rSt + vy_xy*y2dSt/rSt
endif

Lx   = 1.d36
x2dj = (x-cm2xj)#replicate(1.,dsz[1])
y2dj = replicate(1.,dsz[0])#(y-cm2yj)
rj   = sqrt(x2dj*x2dj + y2dj*y2dj)
ionxy  = Lx*!unit.mh/dxy/rj/rj

pos0 = fltarr(2) & pos1 = fltarr(2)
pos0[0] = r * cos(th) + cm2xSt
pos0[1] = r * sin(th) + cm2ySt
pos1[0] = r * cos(th+!pi) + cm2xSt
pos1[1] = r * sin(th+!pi) + cm2ySt

plots,ring(cm2xSt,cm2ySt,r),/data,color=255
tvlct,rr,gg,bb,/get
plots,[pos0[0],pos1[0]],[pos0[1],pos1[1]],/data,color=fsc_color('yellow')
tvlct,rr,gg,bb
timestr = string(sqrt(!unit.g * mtot / peri^3.d0)/2./!pi*time,format='(f4.2)')+' Period'
legend, timestr,/right,/top,box=0

lineread,dxy,x=x,y=y,pos0,pos1,/data,lout=l,value=ld,/no_window, npoint=500
if arg_present(lvr) then $
lineread,vrxy,x=x,y=y,pos0,pos1,/data,lout=l,value=lvr,/no_window, npoint=500
lineread,ionxy,x=x,y=y,pos0,pos1,/data,lout=l,value=lion,/no_window, npoint=500

end

pro batch, fileindex

;fileindex = [0,50,100]
fileindex = [0,100,200,300]
sfile = n_elements(fileindex)

npoint=500
data = replicate({img:fltarr(3,512,512), l:fltarr(npoint), ld:fltarr(npoint), lvr:fltarr(npoint), lion:fltarr(npoint), time:0.}, sfile)
for i=0, sfile-1 do begin
   asywind, n=fileindex[i], lout=l, ld=ld, lvr=lvr, lion=lion, time=time
   data[i].img = tvrd(true=1)
   data[i].l = l
   data[i].ld = ld
   data[i].lvr = lvr
   data[i].lion = lion
   data[i].time = time
endfor

save,file='asywind_data.sav', data
stop
end

pro makeplot,ps=ps,out=out

restore,file='asywind_data.sav'

if not keyword_set(out) then out=''
sd = size(data,/dimension)

loadct,39,/sil
if not keyword_set(ps) then begin
  !p.background=255 & !p.color=0.
  window,0,xs=700,ys=900
endif else mkeps,'asywind_plot'+out,xs=20.,ys=20.
  
width = 0.29 & x0=0.15 & y0=0.1 & x1=0.95
; density
plot, data[0].l-data[0].l[n_elements(data[0].l)-1]/2. ,data[0].ld,/ylog, ytitle='density',pos=[x0,y0+width*2,x1,y0+width*3],/norm,xtickformat="(a1)"
oplot,data[1].l-data[1].l[n_elements(data[1].l)-1]/2. ,data[1].ld,color=50
oplot,data[2].l-data[2].l[n_elements(data[2].l)-1]/2. ,data[2].ld,color=250
;oplot,data[2].l-data[2].l[n_elements(data[2].l)-1]/2. ,data[2].ld,color=150
;oplot,data[3].l-data[3].l[n_elements(data[3].l)-1]/2. ,data[3].ld,color=250

oplot,[-1.4e12,-1.4e12],10.^!y.crange,line=1
oplot,[1.4e12,1.4e12],10.^!y.crange,line=1
oplot,[-3e12,-3e12],10.^!y.crange,line=1

m1 = 4.d34
m2 = 2.d34
mtot = m1+m2
peri = 3.d12
timestr = string(sqrt(!unit.g * mtot / peri^3.d0)/2./!pi*data.time,format='(f4.2)')+' Period'

legend, timestr[0:2],/left,/top,box=0,color=[0,50,250],textcolor=[0.,50.,250]
;legend, timestr[0:3],/left,/top,box=0,color=[0,50,150,250],textcolor=[0.,50.,150,250]

; radial velocity
plot, data[0].l-data[0].l[n_elements(data[0].l)-1]/2. ,data[0].lvr,ytitle='radial velocity',pos=[x0,y0+width,x1,y0+width*2],/norm,/noerase, xtickformat="(a1)", yra=[0.,3.e8]
oplot,data[1].l-data[1].l[n_elements(data[1].l)-1]/2. ,data[1].lvr,color=50
oplot,data[2].l-data[2].l[n_elements(data[2].l)-1]/2. ,data[2].lvr,color=250
;oplot,data[2].l-data[2].l[n_elements(data[2].l)-1]/2. ,data[2].lvr,color=150
;oplot,data[3].l-data[3].l[n_elements(data[3].l)-1]/2. ,data[3].lvr,color=250

oplot,[-1.4e12,-1.4e12],!y.crange,line=1
oplot,[1.4e12,1.4e12],!y.crange,line=1
oplot,[-3e12,-3e12],!y.crange,line=1

; mass loss rate
plot, data[0].l-data[0].l[n_elements(data[0].l)-1]/2. , 4.*!pi*(data[0].l-data[0].l[n_elements(data[0].l)-1]/2.)^2.*data[0].lvr*data[0].ld,xtitle='r [cm]' $
    , ytitle='mass loss rate',pos=[x0,y0,x1,y0+width],/norm,/noerase,yra=[0.,4.5e20]
oplot,data[1].l-data[1].l[n_elements(data[1].l)-1]/2. , 4.*!pi*(data[1].l-data[1].l[n_elements(data[1].l)-1]/2.)^2.*data[1].lvr*data[1].ld,color=50
oplot,data[2].l-data[2].l[n_elements(data[2].l)-1]/2. , 4.*!pi*(data[2].l-data[2].l[n_elements(data[2].l)-1]/2.)^2.*data[2].lvr*data[2].ld,color=250
;oplot,data[2].l-data[2].l[n_elements(data[2].l)-1]/2. , 4.*!pi*(data[2].l-data[2].l[n_elements(data[2].l)-1]/2.)^2.*data[2].lvr*data[2].ld ,color=150
;oplot,data[3].l-data[3].l[n_elements(data[3].l)-1]/2. , 4.*!pi*(data[3].l-data[3].l[n_elements(data[3].l)-1]/2.)^2.*data[3].lvr*data[3].ld ,color=250

oplot,[-1.4e12,-1.4e12],!y.crange,line=1
oplot,[1.4e12,1.4e12],!y.crange,line=1
oplot,[-3e12,-3e12],!y.crange,line=1
if keyword_set(ps) then epsfree

simg = size(data.img,/dimension)
for i=0,sd[0]-1 do begin
   window,i+1,xs=simg[1],ys=simg[2]
   tv,data[i].img,true=1
   if keyword_set(ps) then draw,'asywind_'+out+''+string(i,format='(I2.2)')
endfor

ionparam=0
if (ionparam) then begin
if keyword_set(ps) then mkeps,'asywind_plot2'+out,xs=20.,ys=20.*6./8. $
  else window, sd[0]+1
  
plot, data[0].l-data[0].l[n_elements(data[0].l)-1]/2. ,data[0].lion,xtitle='r [cm]',ytitle='ion paramter',/ylog
oplot,data[1].l-data[1].l[n_elements(data[1].l)-1]/2. ,data[1].lion,color=50
;oplot,data[2].l-data[2].l[n_elements(data[2].l)-1]/2. ,data[2].lion,color=250
oplot,data[2].l-data[2].l[n_elements(data[2].l)-1]/2. ,data[2].lion,color=150
oplot,data[3].l-data[3].l[n_elements(data[3].l)-1]/2. ,data[3].lion,color=250

oplot,!x.crange,[100.,100.],line=1
oplot,[-1.4e12,-1.4e12],10.^!y.crange,line=1
oplot,[1.4e12,1.4e12],10.^!y.crange,line=1
oplot,[-3e12,-3e12],10.^!y.crange,line=1

if keyword_set(ps) then epsfree
endif ; ionparam

stop
end
