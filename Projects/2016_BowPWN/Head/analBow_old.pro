forward_function dw_dx_funct_1, dw_dx_funct_2

pro analBow, ps=ps

; parameters

Ls = 1d33            ; source luminosity
rho0 = 1.67d-24 *0.7     ; ambient density
vs = 1.5d8             ; source velocity 
vp = 1.d10            ; pulsar wind velocity
P0 = 3.d-12           ; ambient pressure


;time = 8.d10
gam  = 1.6666d0

a    = rho0^(1.d0-gam)*vp^(2.d0*gam)*vs^(2.d0-2.d0*gam)
print, 'a = ',a
;a    = 6.24d38

fname = 'PWN2d_hdf5_plt_cnt_0180_d0.7'
xra=[0.,1.e16]
yra=[-2.e15,2.e16]
sample=0
read_amr,fname,Param=Param,/nodata
time=Param.time

mkdata = 1

if mkdata then begin
   call_dw_dx, x, w1, w2, dw1dx, dw2dx, $
              time=time, Ls=Ls, rho0=rho0, vs=vs, P0=P0, a=a, nx=1000

   bubble, time=time, Ls=Ls, rho0=rho0, vs=vs, $
           R_b=R_b, pos=pos
   
   save,file='analBow_0180_d0.7.dat',x,w1,w2,R_b,pos
;   save,file='analBow.dat',x,w1,R_b,pos
endif else restore, file='analBow_0180_d0.7.dat'

;d = loaddata(fname,'dens',sample=1,xra=[0.,1.e18],yra=[-1.e17,3.e18],xCoord=xx,yCoord=yy)
d = loaddata(fname,'dens',sample=sample,xCoord=xx,yCoord=yy,xra=xra,yra=yra)
sz = size(d,/dimension)


;xoffset = 2.d16
xoffset = 4.e14
; stagnation point
Rst = sqrt(Ls / (4.d0*!dpi*rho0*vs^2.d0*vp))

tvlct, r,g,b,/get
if keyword_Set(ps) then mkeps, 'pp_analBow', xs=20., ys=20.*(sz[0]+100)/(sz[1]+300) $
   else window,0,xs=sz[1]+300,ys=sz[0]+100
tvcoord,transpose(alog(d)),yy,xx,/scale,/on, xtitle='x [cm]', ytitle='y [cm]',imgsize=0.63, position=[0.19,0.24]

oplot,x-Rst-xoffset,w1, color=fsc_color('magenta'), line=2
oplot,x-Rst-xoffset,w2, color=fsc_color('magenta'), line=2
rings = ring(pos, 0, R_b)
xring = reform(rings[0,*]) & yring = reform(rings[1,*])
oplot,xring,yring,color=0

tvlct,r,g,b
color_bar,lim=[min(d),max(d)],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),titlegap=0.12,charsize=1.5

if keyword_set(ps) then epsfree
stop
end

;=====================================================================================
pro call_dw_dx, time=time, nx=nx, Ls=Ls, rho0=rho0, vs=vs, P0=P0, a=a,  $         ; input            
                x, w1, w2, dw1dx, dw2dx                                          ; output

common dw_Params, L_s, d_0, P_0, v_s, tt, aa, x00, R_bub, v_bub

L_s = Ls
d_0 = rho0
P_0 = P0
v_s = vs
tt  = time
aa  = a

bubble, time=tt, vs=v_s, rho0=d_0, Ls=L_s, R_b=R_b, v_b=v_b

R_bub = R_b
v_bub = v_b

if not keyword_set(nx) then nx=100
x = (dindgen(nx)/double(nx) + 0.001d0) * (vs*time)
x00 = x[0]

w1 = dblarr(nx)
w2 = dblarr(nx)
dw1dx = dblarr(nx)
dw2dx = dblarr(nx)

for i=0, nx-1 do begin
   x0 = x[0]
   w10 = x[0]
   ddeabm, 'dw_dx_funct_1', x0, w10, x[i]
   w1[i] = w10
endfor

for i=0,nx-1 do dw1dx[i]=dw_dx_funct_1(x[i],w1[i])

for i=0,nx-1 do dw2dx[i]=dw_dx_funct_2(x[i])
;t0=replicate(0.001d0*1.d19,1001)
print,'start calculating w2 by qromo'
x0=replicate(x[0],nx)
for i=1,nx-1 do x0[i]=x[i-1]
;w2_test=qromo('dw_dx_funct_2',x0,x,/double,eps=1.e-2,jmax=14)
w2_test=qromo('dw_dx_funct_2',x0,x,eps=1.e-2,jmax=14)
w2[0]=w1[0]

for i=1,nx-1 do w2[i]=w2[i-1]+w2_test[i]

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
;a=a, time=time ;, con=con,c1=c1,p0=p0,m1=m1,c2=c2,R_b=R_b,P_min1=P_min1,beta_min=beta_min,theta_min=theta_min,v_b=v_b,R_i=R_i,p_min2=p_min2


common dw_Params, L_s, d_0, P_0, v_s, tt, aa, x00, R_bub, v_bub

M0 = ((d_0*v_s^2.)/(5./3.*P_0))^(1./2.)

c = L_s/(5.d0/2.d0*aa^(3.d0/5.d0)*(d_0*v_s^2.d0)^(2.d0/5.d0))
C2 = (5.d0/2.d0)^(-1.d0/4.d0) * (c*aa^(3.d0/5.d0)/sqrt(2.d0*aa^(3.d0/5.d0)))^(1.d0/2.d0) * (d_0*v_s^2.d0)^(-1.d0/10.d0)  ;1.d0

;bubble, time=tt, vs=v_s, rho0=d_0, Ls=L_s, R_b=R_b, v_b=v_b

M_b = v_bub/(5.d0/3.d0*P_0/d_0)^(1.d0/2.d0)
P_min = P_0*(5.d0/4.d0*M_b^2.d0-1.d0/4.d0)

R_i = (1.d0 - (2.d0/3.d0*M_b^2.d0+2.d0)/(8.d0/3.d0*M_b^2.d0))^(1.d0/3.d0)*R_bub
P_min2 = 2.d0/5.d0*L_s*tt/(4.d0/3.d0*!dpi*R_i^3.d0)

;p_min=p_min * v_0*time/(v_0*time-R_b)

P_min1=P_min

v_i=v_s-v_bub
M_i=((d_0*v_i^2.d0)/(5.d0/3.d0*P_0))^(1.d0/2.d0)

beta_min=v_bub/(v_s-v_bub) 
theta_min=(M_i^2.d0*beta_min^2.d0-1.d0)/(4.d0/3.d0*m_i^2.d0)/beta_min

M1=M_i

return,max([3.d0*(C2^(10.d0/3.d0) - P_0*w^(10.d0/3.))/(M1*P_0*w^(10.d0/3.d0)*(20.d0*C2^(10.d0/3.d0)/P_0/w^(10.d0/3.d0)+5.d0)^(1.d0/2.d0)),theta_min])

END


function dw_dx_funct_2, x

common dw_Params, L_s, d_0, P_0, v_s, tt, aa, x00, R_bub, v_bub

M0 = ((d_0*v_s^2.d0)/(5.d0/3.d0*P_0))^(1.d0/2.d0)

;m_i=((rho_0*v_i^2.)/(5./3.*p0))^(1./2.)

;time=3.27602d+13 ; for file=1660
;bubble, time=tt, vs=v_s, rho0=d_0, Ls=L_s, R_b=R_b, v_b=v_b

M_b = v_bub/(5.d0/3.d0*P_0/d_0)^(1.d0/2.d0)
P_min = P_0*(5.d0/4.d0*M_b^2.d0 - 1.d0/4.d0)

v_i=v_s-v_bub
M_i=((d_0*v_i^2.d0)/(5.d0/3.d0*P_0))^(1.d0/2.d0)

M1=M_i

; this is dw2/dx in terms of dw1/dx at w
; want it in terms of x I think.

;x0 = v_s*tt/double(nxx)*1.d-3
x0 = x00
w10 = x0

ddeabm,'dw_dx_funct_1',x0,w10,x

df=dw_dx_funct_1(x,w10)

return,(8.d0/9.d0*df^2.d0+1.d0/M1^2.d0+8.d0/9.d0*df*(df^2.d0+9.d0/4.d0/M1^2.d0)^(1.d0/2.d0))^(1.d0/2.d0)

END

