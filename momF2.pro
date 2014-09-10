pro momF2,n,sample=sample,xc0=xc0,zc0=zc0, z2=z2,momFzj=momFzj,momFzt=momFzt,time=time

if (n_elements(n) eq 0) then n=800
if not keyword_set(sample) then sample=0
if not keyword_set(xc0) then xc0=-2.e12
if not keyword_set(zc0) then zc0=0.
yc0=0.  

dens=dload(n,var='dens',sample=sample,xc=xc0,yc=yc0,zc=zc0,x,y,z,time)
velz=dload(n,var='velz',sample=sample,xc=xc0,yc=yc0,zc=zc0)
jet =dload(n,var='jet',sample=sample,xc=xc0,yc=yc0,zc=zc0)
ss = size(dens,/dimension)

dx=x[2]-x[1]
dy=y[2]-y[1]
dA = dx*dy

z0ind_tmp = where(z ge 0) & z0ind = z0ind_tmp[0]

z2     = z[z0ind:ss[2]-1]
momFzj = fltarr(ss[2]-z0ind)
momFzt = fltarr(ss[2]-z0ind)

for i=z0ind,ss[2]-1 do begin
   momFzj[i-z0ind] = total(dens[*,*,i]*velz[*,*,i]^2.*jet[*,*,i])*dA
   momFzt[i-z0ind] = total(dens[*,*,i]*velz[*,*,i]^2.)*dA
endfor

;!p.background=255 & !p.color=0
;window,0
;plot,z2,momFzt,xtitle='z [cm]',ytitle=textoidl('Momentum Flux_{z}'),yr=[1.e24,1.e29],/yst,/ylog
;oplot,z2,momFzj,line=2
;legend,['total','jet'],line=[0,2],color=0,textcolor=0,box=0,/right,/center
end

pro momF2_comb
xc0=-2.e12 & zc0=0.
xc1=-2.6e12 & zc1=2.395e12

dir='momF2_data'
spawn,'mkdir '+dir

for n=0,900,4 do begin

momF2,n,sample=4,xc0=xc0,zc0=zc0,z2=z2_p4,momFzj=momFzj_p4,momFzt=momFzt_p4,time=time
momF2,n,sample=0,xc0=xc0,zc0=zc0,z2=z2_p01,momFzj=momFzj_p01,momFzt=momFzt_p01
momF2,n,sample=0,xc0=xc1,zc0=zc1,z2=z2_p02,momFzj=momFzj_p02,momFzt=momFzt_p02

!p.background=255 & !p.color=0
loadct,39,/sil
window,0,xs=800,ys=600,/pixmap
;mkeps,'momF2_comb.eps'
plot,z2_p4,momFzt_p4,xtitle='z [cm]',ytitle=textoidl('Momentum Flux_{z}') $
    ,xr=[0.,1.e13],yr=[1.e24,1.e29],/yst,/ylog,line=2
;oplot,z2_p01,momFzt_p01,color=50
;oplot,z2_p02,momFzt_p02,color=50

oplot,z2_p4,momFzj_p4
oplot,z2_p01,momFzj_p01,color=50
oplot,z2_p02,momFzj_p02,color=50
legend,['total','jet'],line=[2,0],color=0,textcolor=0,box=0,/right,/bottom
xyouts,!p.clip[2]-150,!p.clip[3]+5,/dev,'time='+string((time-9.76e4)/60./60.,format='(f5.1)')+' hr'
;legend,['smp4','smp0'],color=[0,50],textcolor=[0,50],box=0,/right,/top
;epsfree

fname='momF2_comb_'+string(n,format='(I3.3)')
draw,dir+'/'+fname
save,filename=dir+'/'+fname + '.sav' $
    ,z2_smp4,momFzj_smp4,momFzt_smp4,z2_smp01,momFzj_smp01,momFzt_smp01 $
    ,z2_smp02,momFzj_smp02,momFzt_smp02

endfor

stop
end
