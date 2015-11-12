pro viewangle,angle=angle,mkdata=mkdata

if not keyword_set(angle) then angle=60.

;dir='/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_elmfist/'
;dir='/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_elmfist_0.5/'
dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_0.5/'
;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar2/'
fname = 'surfb_3d.sav'
restore,file=dir+fname

if keyword_set(mkdata) then begin
surfb3d_rot = rot_3d(surfb3d,rotaxis=2,degree=90.-angle,/interp)

prj_rot = total(surfb3d_rot,3)

;save,file=dir+'surfb3d_rot_'+strtrim(fix(angle),2)+'.sav',surfb3d_rot,x,y,z,angle,prj_rot
save,file=dir+'surfb3d_rot_'+strtrim(fix(angle),2)+'.sav',x,y,z,angle,prj_rot
endif else restore,file=dir+'surfb3d_rot_'+strtrim(fix(angle),2)+'.sav'

sz = size(surfb3d,/dimension)

prj = total(surfb3d,3)
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
dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_0.5/'
restore,file=dir+'surfb_3d.sav'
prj = total(surfb3d,3)

restore,file=dir+'surfb3d_rot_30.sav'
prj_30 = prj_rot
restore,file=dir+'surfb3d_rot_45.sav'
prj_45 = prj_rot
restore,file=dir+'surfb3d_rot_60.sav'
prj_60 = prj_rot

print,'complete to read the data'
sz = size(surfb3d,/dimension)

pltx0 = 120. & plty0 = 80.
winxs=sz[0]+pltx0+200 & winys=plty0+sz[1]*4+20
loadct,39,/sil
swindow,xs=winxs,ys=winys
minv=min(prj) & maxv=max(prj)
tvcoord, bytscl(prj,min=minv,max=maxv),x,y,pos=[pltx0,plty0+3*sz[1]],/dev,xtickformat='(a1)',/axes,ytitle='y [cm]'
legend,'view angle= 90 deg',/right,/top,box=0
tvcoord, bytscl(prj_60,min=minv,max=maxv),x,y,pos=[pltx0,plty0+2*sz[1]],/dev,xtickformat='(a1)',/axes,ytitle='y [cm]'
legend,'view angle= 60 deg',/right,/top,box=0
tvcoord, bytscl(prj_45,min=minv,max=maxv),x,y,pos=[pltx0,plty0+sz[1]],/dev,xtickformat='(a1)',/axes, ytitle='y [cm]'
legend,'view angle= 45 deg',/right,/top,box=0
tvcoord, bytscl(prj_30,min=minv,max=maxv),x,y,pos=[pltx0,plty0],/dev,xtitle='x [cm]',/axes, ytitle='y [cm]'
legend,'view angle= 30 deg',/right,/top,box=0
color_bar,lim=[minv,maxv],pos=[pltx0+sz[0]+10,plty0,pltx0+sz[0]+30,plty0+sz[1]*4],/right,bartitle='projected surface brightness',charsize=2.,titlegap=0.14

stop
end

pro drawview_ps

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_0.5/'
restore,file=dir+'surfb_3d.sav'
minprj = 1.e-18
prj = total(surfb3d,3)>minprj

restore,file=dir+'surfb3d_rot_60.sav'
prj_60 = prj_rot>minprj

print,'complete to read the data'
sz = size(surfb3d,/dimension)

pltx0 = 220. & plty0 = 120.
winxs=sz[0]+pltx0+200 & winys=plty0+sz[1]*2+20
;loadct,39,/sil
loadct,3,/sil
;swindow,xs=winxs,ys=winys
mkeps,'guitar_halpha',xs=20.,ys=20.*winys/winxs
minv=min(prj) & maxv=max(prj)
tvcoord, bytscl(alog10(prj),min=alog10(minv),max=alog10(maxv)),x,y,pos=[pltx0/winxs,(plty0+sz[1])/winys],/norm,psx=sz[0]/winxs,xtickformat='(a1)',/axes,ytitle='y [cm]',color=0,ytickinterval=4.e17
plot,x,y,xra=[x[0],x[n_elements(x)-1]],yra=[y[0],y[n_elements(y)-1]],pos=posnorm([pltx0,plty0+sz[1],pltx0+sz[0],plty0+sz[1]*2],nx=winxs,ny=winys) $
    ,/nodata,xtickformat='(a1)',ytickformat='(a1)',/xst,/yst,color=255,/noerase
legend,'view angle= 90'+textoidl('^{\circ}'),/left,/top,box=0,textcolor=255,charsize=1.
tvcoord, bytscl(alog10(prj_60),min=alog10(minv),max=alog10(maxv)),x,y,pos=[pltx0/winxs,plty0/winys],/norm,psx=sz[0]/winxs,xtickformat='(a1)',/axes,ytitle='y [cm]',color=0,xtitle='x [cm]',ytickinterval=4.e17
plot,x,y,xra=[x[0],x[n_elements(x)-1]],yra=[y[0],y[n_elements(y)-1]],pos=posnorm([pltx0,plty0,pltx0+sz[0],plty0+sz[1]],nx=winxs,ny=winys) $
    ,/nodata,xtickformat='(a1)',ytickformat='(a1)',/xst,/yst,color=255,/noerase
legend,'view angle= 60'+textoidl('^{\circ}'),/left,/top,box=0,textcolor=255,charsize=1.
color_bar,lim=[minv,maxv],/log,pos=posnorm([pltx0+sz[0]+10,plty0,pltx0+sz[0]+30,plty0+sz[1]*2],nx=winxs,ny=winys),/norm,/right,bartitle='projected surface brightness',titlegap=0.12

epsfree

stop
end
