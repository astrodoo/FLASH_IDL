pro tmpwind2

restore,filename='therm1d_15deg_1000.sav'
;restore,filename='therm1d_00deg_1000'
x00=x & d00=den & p00=pres & vx00=vx
restore,filename='therm1d_15deg_d1e14_1000.sav'
;restore,filename='therm1d_05deg_1000'
x05=x & d05=den & p05=pres & vx05=vx
restore,filename='therm1d_15deg_d5e15_1000.sav'
;restore,filename='therm1d_15deg_1000'
x15=x & d15=den & p15=pres & vx15=vx

rs = 1.4e12
dy = 0.29
loadct,39,/sil
!p.background=255 & !p.color=0
window,0,xs=800,ys=1000
plot,x00/rs,d00 ,pos=[0.15,0.1+2.*dy,0.97,0.1+3.*dy],/norm,xtickformat='(a1)',xst=2 $
    ,ytitle='density',/ylog,yst=2
oplot,x05/rs,d05,color=50
oplot,x15/rs,d15,color=240
legend,textoidl('\rho = ')+['2e-14','1e-14','5e-15'],color=[0,50,240],textcolor=[0,50,240],box=0,/right,/top
;legend,'ang = '+[' 0 deg',' 5 deg','15 deg'],color=[0,50,240],textcolor=[0,50,240],box=0,/right,/top
plot,x/rs,p00,pos=[0.15,0.1+dy,0.97,0.1+2.*dy],/norm,xtickformat='(a1)',/noerase,xst=2 $
    ,ytitle='pressure',yst=2
oplot,x05/rs,p05,color=50
oplot,x15/rs,p15,color=240

plot,x/rs,vx00/1.e5,pos=[0.15,0.1,0.97,0.1+dy],/norm,xtitle='r/Rs',/noerase,xst=2 $ 
    ,ytitle='velx [km/s]',yst=2,yr=[0.,1000.]
oplot,x05/rs,vx05/1.e5,color=50
oplot,x15/rs,vx15/1.e5,color=240

stop
end
