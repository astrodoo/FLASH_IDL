pro viewangle,angle=angle,mkdata=mkdata

if not keyword_set(angle) then angle=60.

dir='/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_highdenTest/'
fname = 'PWN2d_hdf5_plt_cnt_0117_3d.sav'
restore,file=dir+fname

if keyword_set(mkdata) then begin
rotd = rot_3d(den3d_xyz,rotaxis=2,degree=90.-angle,/interp)

den3d_rot = rotd

save,file=dir+'den3d_rot_'+strtrim(fix(angle),2)+'.sav',den3d_rot,x,y,z,angle
endif else restore,file=dir+'den3d_rot_'+strtrim(fix(angle),2)+'.sav'

sz = size(den3d_rot,/dimension)

loadct,1,/sil
window,0,xs=sz[0],ys=sz[1]*2
tvscl,alog(total(double(den3d_xyz)*double(den3d_xyz),3)),0
tvscl,alog(total(double(den3d_rot)*double(den3d_rot),3)),1
legend,'Emission Measure (view=90'+textoidl('^{\circ})'),pos=[10,sz[1]*2-10],/dev,textcolor=fsc_color('magenta'),box=0
legend,'Emission Measure (view='+strtrim(fix(angle),2)+textoidl('^{\circ})'),pos=[10,sz[1]-10],/dev,textcolor=fsc_color('magenta'),box=0

stop
end

pro drawall

dir='/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_highdenTest/'
fname = 'PWN2d_hdf5_plt_cnt_0117_3d.sav'
restore,file=dir+fname

restore,file=dir+'den3d_rot_60.sav'
den3d_rot60 = den3d_rot & ang60 = angle
restore,file=dir+'den3d_rot_30.sav'
den3d_rot30 = den3d_rot & ang30 = angle

sz = size(den3d_rot,/dimension)

em_90 = total(double(den3d_xyz)*double(den3d_xyz),3)

em_60 = total(double(den3d_rot60)*double(den3d_rot60),3)
em_60_dummy = em_60
em_60_dummy[where(em_60 eq min(em_60))] = max(em_60)
em_60[where(em_60 eq min(em_60))] = min(em_60_dummy)   ; set the minimum to the second minimum

em_30 = total(double(den3d_rot30)*double(den3d_rot30),3) 
em_30_dummy = em_30
em_30_dummy[where(em_30 eq min(em_30))] = max(em_30)
em_30[where(em_30 eq min(em_30))] = min(em_30_dummy)   ; set the minimum to the second minimum 

loadct,1,/sil
window,0,xs=sz[0],ys=sz[1]*3
tvscl,alog(em_90),0
tvscl,alog(em_60),1
tvscl,alog(em_30),2
legend,'Emission Measure (view=90'+textoidl('^{\circ})'),pos=[10,sz[1]*3-10],/dev,textcolor=fsc_color('magenta'),box=0
legend,'Emission Measure (view='+strtrim(fix(ang60),2)+textoidl('^{\circ})'),pos=[10,sz[1]*2-10],/dev,textcolor=fsc_color('magenta'),box=0
legend,'Emission Measure (view='+strtrim(fix(ang30),2)+textoidl('^{\circ})'),pos=[10,sz[1]-10],/dev,textcolor=fsc_color('magenta'),box=0

stop
end



