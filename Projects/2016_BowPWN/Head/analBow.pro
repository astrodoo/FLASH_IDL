forward_function dw_dx_funct_1 ;, dw_dx_funct_2

pro analBow, time=time, nx=nx, Ls=Ls, d0=d0, vs=vs, vp=vp, P0=P0,  $         ; input            
                  x, w1, dw1dx, info=info, out=out,bubp=bubp                 ; output

common dw_Params, L_s, d_0, P_0, v_s, tt, aa, x00, R_bub, v_bub, bubblep

if not keyword_set(nx) then nx=1000
if not keyword_set(time) then time=2.5e10
if not keyword_set(Ls) then Ls = 1.d36          ; source luminosity
if not keyword_set(d0) then d0 = 1.67d-24          ; ambient density
if not keyword_set(vs) then vs = 1.d8              ; source velocity 
if not keyword_set(vp) then vp = 1.d10             ; pulsar wind velocity
if not keyword_set(P0) then P0 = 3.d-12            ; ambient pressure
if keyword_set(bubp) then bubblep=1 else bubblep=0

; stagnation point
Rst = sqrt(Ls / (4.d0*!dpi*d0*vs^2.d0*vp))

gam  = 1.33333d0
a    = d0^(1.d0-gam)*vp^(2.d0*gam)*vs^(2.d0-2.d0*gam) ; constant 'a'
print, 'a = ',a

L_s = Ls
d_0 = d0
P_0 = P0
v_s = vs
tt  = time
aa  = a

info = {Ls:Ls,d0:d0,P0:P0,vs:vs,vp:vp,time:time}

bubble, time=tt, vs=v_s, rho0=d_0, Ls=L_s, R_b=R_b, v_b=v_b

R_bub = R_b
v_bub = v_b

if not keyword_set(nx) then nx=100
x = (dindgen(nx)/double(nx) + 0.001d0) * (vs*time)
x00 = x[0]

w1 = dblarr(nx)
;w2 = dblarr(nx)
dw1dx = dblarr(nx)
;dw2dx = dblarr(nx)

for i=0, nx-1 do begin
   x0 = x[0]
   w10 = x[0]
   ddeabm, 'dw_dx_funct_1', x0, w10, x[i]
   w1[i] = w10
endfor

for i=0,nx-1 do dw1dx[i]=dw_dx_funct_1(x[i],w1[i])

x = x - Rst

;plot,x,w1,/iso

if keyword_set(out) then save,file=out,x,(w=w1),info
end


;=====================================================================================
pro bubble, time=time, vs=vs, rho0=rho0, Ls=Ls,  $    ; input
            R_b=R_b, v_b=v_b, pos=pos                 ; output

pos = vs*time

;C1 = 25.d0/14.d0/!dpi
C1 = 125.d0/154.d0/!dpi

R_b = C1^(1.d0/5.d0) * (Ls/rho0)^(1.d0/5.d0) * time^(3.d0/5.d0)
v_b = C1^(1.d0/5.d0) * (Ls/rho0)^(1.d0/5.d0) * 3.d0/5.d0*time^(-2.d0/5.d0)
end

;=====================================================================================
function dw_dx_funct_1, x, w

common dw_Params, L_s, d_0, P_0, v_s, tt, aa, x00, R_bub, v_bub, bubblep

;M0 = ((d_0*v_s^2.)/(5./3.*P_0))^(1./2.)
M0 = ((d_0*v_s^2.)/(4./3.*P_0))^(1./2.)

c = L_s/(5.d0/2.d0*aa^(3.d0/5.d0)*(d_0*v_s^2.d0)^(2.d0/5.d0))
;C2 = (5.d0/2.d0)^(-1.d0/4.d0) * (c*aa^(3.d0/5.d0)/sqrt(2.d0*aa^(3.d0/5.d0)))^(1.d0/2.d0) * (d_0*v_s^2.d0)^(-1.d0/10.d0)  ;1.d0
C2 = 2^(-1.d0/4.d0) * sqrt(L_s) * aa^(-3.d0/16.d0) * (d_0*v_s^2.d0)^(-3.d0/16.d0) * 4.^(-1./4.)

;bubble, time=tt, vs=v_s, rho0=d_0, Ls=L_s, R_b=R_b, v_b=v_b

M_b = v_bub/(5.d0/3.d0*P_0/d_0)^(1.d0/2.d0)
P_min = P_0*(5.d0/4.d0*M_b^2.d0-1.d0/4.d0)

;R_i = (1.d0 - (2.d0/3.d0*M_b^2.d0+2.d0)/(8.d0/3.d0*M_b^2.d0))^(1.d0/3.d0)*R_bub
;P_min2 = 2.d0/5.d0*L_s*tt/(4.d0/3.d0*!dpi*R_i^3.d0)

;p_min=p_min * v_0*time/(v_0*time-R_b)

P_min1=P_min

v_i=v_s-v_bub
;M_i=((d_0*v_i^2.d0)/(5.d0/3.d0*P_0))^(1.d0/2.d0)
M_i=((d_0*v_i^2.d0)/(4.d0/3.d0*P_0))^(1.d0/2.d0)

beta_min=v_bub/(v_s-v_bub) 
;theta_min=(M_i^2.d0*beta_min^2.d0-1.d0)/(4.d0/3.d0*m_i^2.d0)/beta_min
theta_min=(M_i^2.d0*beta_min^2.d0-1.d0)/(7.d0/6.d0*m_i^2.d0)/beta_min

;M1=M_i
M1=M0

;return,max([3.d0*(C2^(10.d0/3.d0) - P_0*w^(10.d0/3.))/(M1*P_0*w^(10.d0/3.d0)*(20.d0*C2^(10.d0/3.d0)/P_0/w^(10.d0/3.d0)+5.d0)^(1.d0/2.d0)),theta_min])
if (bubblep) then $
   return,max([3.d0*(C2^(8.d0/3.d0) - P_0*w^(8.d0/3.d0))/(M1*P_0*w^(8.d0/3.d0)*(14.d0*C2^(8.d0/3.d0)/P_0/w^(8.d0/3.d0)+2.d0)^(1.d0/2.d0)),theta_min])  $
  else $
   return,[3.d0*(C2^(8.d0/3.d0) - P_0*w^(8.d0/3.d0))/(M1*P_0*w^(8.d0/3.d0)*(14.d0*C2^(8.d0/3.d0)/P_0/w^(8.d0/3.d0)+2.d0)^(1.d0/2.d0))]

END

pro multibow

nx=1000
time=2.5e10
Ls = 1.d36          ; source luminosity
;d0 = 1.67d-24          ; ambient density
d0 = 1.67d-23          ; ambient density
vs = 1.d8              ; source velocity 
vp = 1.d10             ; pulsar wind velocity
P0 = 3.d-12            ; ambient pressure

;time=1.e11
;Ls = 4.3d33          ; source luminosity
;;d0 = 1.67d-24          ; ambient density
;d0 = 1.67d-23          ; ambient density
;vs = 6.d6              ; source velocity 
;vp = 1.d10             ; pulsar wind velocity
;P0 = 3.d-12            ; ambient pressure

d1 = d0/20.
nd = 100
;d1 = d0/200.
;nd = 1000
dd = findgen(nd)/float(nd)*(d0-d1)+d1

analBow, time=time, nx=nx, Ls=Ls, d0=dd[0], vs=vs, vp=vp, P0=P0,  $         ; input            
         x, w1, dw1dx, info=info, /nosav                                  ; output

nx = n_elements(x)
ds = replicate( {x:fltarr(nx),w1:fltarr(nx),dw1dx:fltarr(nx), d0:0.}, nd)
ds[0].x = x
ds[0].w1=w1
ds[0].dw1dx=dw1dx
ds[0].d0=dd[0]

for i=1,nd-1 do begin
   print,i,' of ',nd
   analBow, time=time, nx=nx, Ls=Ls, d0=dd[i], vs=vs, vp=vp, P0=P0, x, w1, dw1dx, /nosav                                  ; output
   ds[i].x = x
   ds[i].w1= w1
   ds[i].dw1dx=dw1dx
   ds[i].d0=dd[i]
endfor

save,file='multibow.sav',ds,info

stop
end

pro draw_multibow

restore,file='multibow_model2.sav'
;restore,file='multibow.sav'

;ds = ds[5:500]

nds = (size(ds,/dimension))[0]

loadct,39,/sil

xra=[-5.e17,2.e18]
yra=[0.,1.5e18]
xtickint = 1.e18
ytickint = 1.e18

;xra=[-5.e17,2.5e18]
;yra=[0.,2.5e18]
;xtickint = 1.e18
;ytickint = 1.e18

;xra=[-1e17,2.5e17]
;yra=[0.,2.5e17]
;xtickint = 2.e17
;ytickint = 1.e17

d0_24 = 1.67e-24

pltx0=250. & plty0=120.
;pltxs=800. & pltys=pltxs* 2.5e18/ 3.e18
pltxs=800. & pltys=pltxs* yra[1]/(xra[1]-xra[0])
winxs=pltx0+pltxs+200. & winys=plty0+pltys+20.
mkeps, 'draw_multibow_model2',xs=20.,ys=20.* winys/winxs
plot,ds[0].x,ds[0].w1,/iso,xtitle='x [cm]',ytitle='y [cm]' $
    ,position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm $
    ,xtickinterval=xtickint,ytickinterval=ytickint,xra=xra,yra=yra,/xst,/yst

for i=1,nds-1 do oplot,ds[i].x,ds[i].w1 ,color=fix(255.*i/float(nds))

color_bar,lim=[min(ds.d0),max(ds.d0)]/d0_24,/right,bartitle=textoidl('\rho_{0}/\rho_{0,24}'),titlegap=0.09

epsfree

;stop
xcrit1 = 1.5e18
xcrit2 = max(ds.x)
;xcrit1 = 2.e17
;xcrit2 = 2.5e17

openang = fltarr(nds)
w1xcrit = fltarr(nds)
for i=0,nds-1 do begin
   xx = double( ds[i].x[ where((ds[i].x ge xcrit1) and (ds[i].x le xcrit2)) ] ) 
   ww1 = double( ds[i].w1[ where((ds[i].x ge xcrit1) and (ds[i].x le xcrit2)) ] )
   linf = linfit(xx,ww1)
   openang[i] = linf[1]
   w1xcrit[i] = ww1[0]
endfor

pltx0=110. & plty0=100.
pltxs=800. & pltys=600.
;winxs=pltx0+pltxs+120. & winys=plty0+pltys+20.
winxs=pltx0+pltxs+20. & winys=plty0+pltys+20.

loadct,39,/sil
mkeps,'draw_multibow2_model2',xs=20.,ys=20.*winys/winxs
plot,ds.d0/d0_24,openang/!dtor,/yst,/xst,xtitle=textoidl('\rho_{0}/\rho_{0,24}'),ytitle=textoidl('\alpha [degree]') $
    ,position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm, thick=5., yra=[5.,25.]
;plot,ds.d0/d0_24,openang/!dtor,/yst,/xst,xtitle=textoidl('\rho_{0}/\rho_{0,24}'),ytitle=textoidl('\alpha [degree]') $
;    ,position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm, thick=5., yra=[3.,9.]

;plot,ds.d0/d0_24,openang/!dtor,yst=9,/xst,xtitle=textoidl('\rho_{0}/\rho_{0,24}'),ytitle=textoidl('\alpha [degree]') $
;    ,position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm, thick=5.

;plot,ds.d0/d0_24,w1xcrit/w1xcrit[nds-1],xst=5,yst=5,/nodata $
;    ,position=posnorm([pltx0,plty0,pltx0+pltxs,plty0+pltys],nx=winxs,ny=winys),/norm, /noerase
;oplot,ds.d0/d0_24,w1xcrit/w1xcrit[nds-1],thick=5.,color=50.
    
;axis,yaxis=1,yst=1,color=50.,ytitle=textoidl('w_{\rho_{1}}/w_{\rho_{0}}')
epsfree
stop
end
