forward_function bubrad

pro show45deg

dir='/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN3d/out_PWN3d_600kms_45deg_highRes/'
file='PWN3d_hdf5_plt_cnt_0071'

;d=reform(loaddata(dir+file,'dens',sample=1,xcoord=x,ycoord=y,xra=[-1.e18,3.e18],yra=[-1.5e18,1.5e18],zra=[0.,0.],time=t))
;save,file=dir+'tmp2d_0071.sav', d,x,y,t

restore,file=dir+'tmp2d_0071.sav'

t0 = 6.86367e10
t = t-t0

bubrad = (bubrad(t))[0]
bubcen = (bubrad(t))[1]

dsm=smooth(d,10)
dsz=size(d,/dimension)

;drawing
pltx0=150. & plty0=80.
winxs=pltx0+dsz[0]+150 & winys=plty0+dsz[1]+20
winxs_cm=25. & winys_cm=winxs_cm*winys/winxs

mkeps,'show45deg',xs=winxs_cm,ys=winys_cm
;!p.background=255 & !p.color=0 & !p.charsize=2.
;window,xs=winxs,ys=winys
lthick=5

tvcoord,alog10(d),x,y,/scale,/axes,/black,xtickinterval=1.e18,pos=[pltx0/winxs,plty0/winys],/norm $
       ,imgsize=dsz[0]/winxs,xtitle='x [cm]', ytitle='y [cm]'

pulsarRad = 1.24e16
;plots,ring(0.,0., pulsarRad),color=255,thick=1

plots,ring(bubcen,0,bubrad,th=[45.,225.]*!dtor),color=fsc_color('magenta'),thick=lthick
contour,alog10(dsm),x,y,level=[-24.5,-24.45,-24.4,-24.1,-23.9,-23.76,-23.75,-23.7],/overplot

; analytic bow shock
restore,file='analBow3_45deg.sav'
bowx = x & bowy = w1
xoffset = 4.e14
Ls = 6.66d35            ; source luminosity
rho0 = 1.67d-24         ; ambient density
vs = 6.d7             ; source velocity 
vp = 1.d10            ; pulsar wind velocity
P0 = 3.d-12           ; ambient pressure
Rst = sqrt(Ls / (4.d0*!dpi*rho0*vs^2.d0*vp))

oplot,x-Rst-xoffset,w1,color=fsc_color('magenta'),thick=lthick
oplot,x-Rst-xoffset,-w1,color=fsc_color('magenta'),thick=lthick

plots,0.,0.,/data,psym=7,color=255,symsize=2

;xyouts,-8.e17,1.3e18,/data, 't= '+string(t/60./60./24./365.,format='(f6.2)')+' yr'


;plots,[0.,0.],psym=7,symsize=2,color=255,/data
arrow,5.e16,-8.e17,5.e16,-3.5e17,thick=3.,/data
loadct,0,/sil
color_bar,lim=[min(d),max(d)],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),bargap=0.01


epsfree

stop
end

function bubrad, time, ambden=ambden

Esp = 6.66d35
vw  = 6.d7
vp  = 1.d10

if keyword_set(ambden) then d1=ambden else d1 = 1.67d-25
rad = (125./154./!pi)^0.2*(Esp/d1)^0.2*time^0.6 
loc = time*vw
;bub1r = ring(l1,0.,r1)
;print,'location 1st bubble: ',l1

return, [rad,loc]
end

pro play45deg

dir='/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN3d/out_PWN3d_600kms_45deg_highRes/'
files = file_search(dir+'PWN3d_hdf5_plt_cnt_????')

outdir='png_play45deg'
spawn,'mkdir '+outdir

startn=0
endn=100

; analytic bow shock
restore,file='analBow3_45deg.sav'
bowx = x & bowy = w1
xoffset = 4.e14
Ls = 6.66d35            ; source luminosity
rho0 = 1.67d-24         ; ambient density
vs = 6.d7             ; source velocity 
vp = 1.d10            ; pulsar wind velocity
P0 = 3.d-12           ; ambient pressure
Rst = sqrt(Ls / (4.d0*!dpi*rho0*vs^2.d0*vp))


maxd = 6.7e-24 & mind = 9.2e-30
; loop
for i=startn, endn do begin

print, i,' of ',endn

d=reform(loaddata(files[i],'dens',sample=1,xcoord=x,ycoord=y,xra=[-1.e18,3.e18],yra=[-1.5e18,1.5e18],zra=[0.,0.],time=t))

t0 = 6.86367e10
tb = t-t0

if (tb ge 0.) then begin
   bubrad = (bubrad(tb))[0]
   bubcen = (bubrad(tb))[1]
endif else begin
   bubrad = 0. & bubcen = 0.
endelse

bubrad2 = (bubrad(t,ambden=1.67d-24))[0]
bubcen2 = (bubrad(t,ambden=1.67d-24))[1]

;dsm=smooth(d,10)
dsz=size(d,/dimension)

;drawing
pltx0=150. & plty0=80.
winxs=pltx0+dsz[0]+150 & winys=plty0+dsz[1]+20
;winxs_cm=25. & winys_cm=winxs_cm*winys/winxs

!p.background=255 & !p.color=0 & !p.charsize=2.
window,xs=winxs,ys=winys,/pixmap
lthick=3

tvcoord,bytscl(alog10(d),min=alog10(mind),max=alog10(maxd)),x,y,/axes,/black,xtickinterval=1.e18,pos=[pltx0,plty0],/dev $
       ,xtitle='x [cm]', ytitle='y [cm]'

plots,ring(bubcen,0,bubrad,th=[45.,225.]*!dtor),color=fsc_color('magenta'),thick=lthick
;contour,alog10(dsm),x,y,level=[-24.5,-24.45,-24.4,-24.1,-23.9,-23.76,-23.75,-23.7],/overplot

bowxx = bowx-Rst-xoffset
bowx_in = where(bowxx lt bubcen2-bubrad2)

if ((n_elements(bowx_in) gt 1) and (bubcen2 ne 0)) then begin
   oplot,bowxx[bowx_in],bowy[bowx_in],color=fsc_color('magenta'),thick=lthick
   oplot,bowxx[bowx_in],-bowy[bowx_in],color=fsc_color('magenta'),thick=lthick
endif

plots,0.,0.,/data,psym=7,color=255,symsize=2

xyouts,-8.e17,1.3e18,/data, 't= '+string(tb/60./60./24./365.,format='(f8.2)')+' yr'


;plots,[0.,0.],psym=7,symsize=2,color=255,/data
;arrow,5.e16,-8.e17,5.e16,-3.5e17,thick=3.,/data
loadct,0,/sil
color_bar,lim=[mind,maxd],/log,/right,bartitle=textoidl('density [g cm^{-3}]'),bargap=0.01

snapshot,outdir+'/play45deg_'+string(i,format='(I4.4)')
endfor ; end loop


stop
end
