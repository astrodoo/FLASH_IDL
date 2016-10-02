pro viewangle,angle=angle,mkdata=mkdata

if not keyword_set(angle) then angle=60.

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_60inc/'
fname = 'phot_3d.sav'
restore,file=dir+fname

if keyword_set(mkdata) then begin
phot3d_rot = rot_3d(phot3d,rotaxis=2,degree=90.-angle,/interp)

prj_rot = total(phot3d_rot,3)

;save,file=dir+'surfb3d_rot_'+strtrim(fix(angle),2)+'.sav',surfb3d_rot,x,y,z,angle,prj_rot
save,file=dir+'phot3d_rot_'+strtrim(fix(angle),2)+'.sav',x,y,z,angle,prj_rot
endif else restore,file=dir+'phot3d_rot_'+strtrim(fix(angle),2)+'.sav'

sz = size(phot3d,/dimension)

prj = total(phot3d,3)
loadct,1,/sil
winxs=sz[0] & winys=sz[1]*3
window,0,xs=winxs, ys=winys
tvscl,prj,0
tvscl,prj_rot,1
tvscl,alog((prj_rot)>1.e-18),2

stop
end

pro drawview

;dir='/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_elmfist_0.5/'
dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar/'

;restore,file=dir+'phot_3d.sav'
;prj = total(phot3d,3)
;save,file=dir+'phot3d_rot_90.sav',x,y,z,prj
;stop

restore,file=dir+'phot3d_rot_90.sav'
x90 = x & y90 = y

dir2 = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_60inc/'
restore,file=dir2+'phot3d_rot_60.sav'
prj_60 = prj_rot
x60 = x & y60 = y

;dir3 = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_30inc/'
;restore,file=dir3+'surfb3d_rot_30.sav'
;prj_30 = prj_rot
;x30 = x & y30 = y

print,'complete to read the data'
sz = size(prj,/dimension)

xl = min(x90) & xr = max(x90)
xoff60 = mean(x60)-mean(x60)*cos(30.*!dtor)
x60 = x60-xoff60
ind60x = where((x60 ge xl) and (x60 le xr))
ind60x = [ind60x[0]-1,ind60x]
x60_2 = x60[ind60x]
prj60_2 = prj_60[ind60x,*]

;xoff30 = mean(x30)-mean(x30)*cos(60.*!dtor)
;x30 = x30-xoff30
;ind30x = where((x30 ge xl) and (x30 le xr))
;ind30x = [ind30x[0]-1,ind30x]
;x30_2 = x30[ind30x]
;prj30_2 = prj_30[ind30x,*]


; discontinuity positions
vs = 1.5e8 
time=3.43e10
disx1 = vs*time - 3.3e18
disx2 = vs*time - 3.913e18
disx3 = vs*time - 4.43e18

pltx0 = 100. & plty0 = 50.
winxs=sz[0]+pltx0+200 & winys=plty0+2*sz[1]+20
loadct,13,/sil
window,xs=winxs,ys=winys
;minv=1.e-15 & maxv=max(prj)
minv=1.e-3 & maxv=max(prj)
tvcoord, bytscl(alog10(prj),min=alog10(minv),max=alog10(maxv)),x90,y90,pos=[pltx0,plty0+sz[1]],/dev,xtickformat='(a1)',/axes,ytitle='y [cm]'
loadct,0,/sil
oplot,[0.,0.],!y.crange,line=2
oplot,[disx1,disx1],!y.crange,line=2
oplot,[disx2,disx2],!y.crange,line=2
oplot,[disx3,disx3],!y.crange,line=2
legend,'view angle= 90 deg',/right,/top,box=0
loadct,13,/sil
tvcoord, bytscl(alog10(prj60_2),min=alog10(minv),max=alog10(maxv)),x90,y90,pos=[pltx0,plty0],/dev,/axes,ytitle='y [cm]',xtitle='x [cm]'
loadct,0,/sil
oplot,[0.,0.],!y.crange,line=2
oplot,[disx1,disx1],!y.crange,line=2
oplot,[disx2,disx2],!y.crange,line=2
oplot,[disx3,disx3],!y.crange,line=2
legend,'view angle= 60 deg',/right,/top,box=0
loadct,13,/sil
color_bar,lim=[minv,maxv]*!unit.h*!unit.c/656.3d-9,/log,pos=[pltx0+sz[0]+10,plty0,pltx0+sz[0]+30,plty0+sz[1]*2],/right,bartitle=textoidl('Surface Brightness [ergs cm^{-2} s^{-1} sr^{-1}]'),charsize=2.,titlegap=0.14

stop
end

pro drawview_ps

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar/'

restore,file=dir+'phot3d_rot_90.sav'
x90 = x & y90 = y

dir2 = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_60inc/'
restore,file=dir2+'phot3d_rot_60.sav'
prj_60 = prj_rot
x60 = x & y60 = y

print,'complete to read the data'
sz = size(prj,/dimension)

xl = min(x90) & xr = max(x90)
xoff60 = mean(x60)-mean(x60)*cos(30.*!dtor)
x60 = x60-xoff60
ind60x = where((x60 ge xl) and (x60 le xr))
ind60x = [ind60x[0]-1,ind60x]
x60_2 = x60[ind60x]
prj60_2 = prj_60[ind60x,*]

; discontinuity positions
vs = 1.5e8 
time=3.43e10
disx1 = vs*time - 3.3e18
disx2 = vs*time - 3.913e18
disx3 = vs*time - 4.43e18

pltx0 = 130. & plty0 = 70.
winxs=sz[0]+pltx0+180 & winys=plty0+2*sz[1]+20
mkeps,'viewangle',xs=30.,ys=30.*winys/winxs

loadct,13,/sil
;window,xs=winxs,ys=winys
;minv=1.e-15 & maxv=max(prj)
minv=1.e-3 & maxv=max(prj)
tvcoord, bytscl(alog10(prj),min=alog10(minv),max=alog10(maxv)),x90,y90,pos=[pltx0/winxs,(plty0+sz[1])/winys],/norm,psx=sz[0]/winxs ;,xtickformat='(a1)',/axes,ytitle='y [cm]'
loadct,0,/sil
plot,x90,y90,xra=[min(x90),max(x90)],yra=[min(y90),max(y90)],/xst,/yst,pos=posnorm([pltx0,plty0+sz[1],pltx0+sz[0],plty0+2*sz[1]],nx=winxs,ny=winys),/norm,/noerase,/nodata,color=0,ytitle='y [cm]'
plot,x90,y90,xra=[min(x90),max(x90)],yra=[min(y90),max(y90)],/xst,/yst,pos=posnorm([pltx0,plty0+sz[1],pltx0+sz[0],plty0+2*sz[1]],nx=winxs,ny=winys),/norm,/noerase,/nodata,color=255 $
    ,xtickformat='(a1)',ytickformat='(a1)'
oplot,[0.,0.],!y.crange,line=2,color=255
oplot,[disx1,disx1],!y.crange,line=2,color=255
oplot,[disx2,disx2],!y.crange,line=2,color=255
oplot,[disx3,disx3],!y.crange,line=2,color=255
legend,textoidl('view angle=90^{\circ}'),/left,/top,box=0,textcolor=255
loadct,13,/sil
tvcoord, bytscl(alog10(prj60_2),min=alog10(minv),max=alog10(maxv)),x90,y90,pos=[pltx0/winxs,plty0/winys],/norm,psx=sz[0]/winxs ;,/axes,ytitle='y [cm]',xtitle='x [cm]'
loadct,0,/sil
plot,x90,y90,xra=[min(x90),max(x90)],yra=[min(y90),max(y90)],/xst,/yst,pos=posnorm([pltx0,plty0,pltx0+sz[0],plty0+sz[1]],nx=winxs,ny=winys),/norm,/noerase,/nodata,color=0,ytitle='y [cm]',xtitle='x [cm]'
plot,x90,y90,xra=[min(x90),max(x90)],yra=[min(y90),max(y90)],/xst,/yst,pos=posnorm([pltx0,plty0,pltx0+sz[0],plty0+sz[1]],nx=winxs,ny=winys),/norm,/noerase,/nodata,color=255 $
    ,xtickformat='(a1)',ytickformat='(a1)'
oplot,[0.,0.],!y.crange,line=2,color=255
oplot,[disx1,disx1],!y.crange,line=2,color=255
oplot,[disx2,disx2],!y.crange,line=2,color=255
oplot,[disx3,disx3],!y.crange,line=2,color=255
legend,textoidl('view angle=60^{\circ}'),/left,/top,box=0,textcolor=255
loadct,13,/sil

srtoasc = (180./!pi)^2.*60.*60.
color_bar,lim=[minv,maxv]/srtoasc*1.e2,/log,pos=posnorm([pltx0+sz[0]+10,plty0,pltx0+sz[0]+30,plty0+sz[1]*2],nx=winxs,ny=winys),/norm,/right $
         ,bartitle=textoidl('Surface Brightness [photons cm^{-2} s^{-1} arcsec^{-2}]'),charsize=2.,titlegap=0.11,/minor

epsfree

stop
end
