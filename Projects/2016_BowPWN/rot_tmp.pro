pro rot_tmp

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_30inc/'
fname = 'surfb_3d.sav'
restore,file=dir+fname

sz = size(surfb3d,/dimension)

rot3d = fltarr(sz[0],sz[1],sz[2])
prj = fltarr(sz[0],sz[1],61)

prj[*,*,0] = total(surfb3d,3)

pngdir = 'png_rot_tmp'
spawn,'mkdir '+pngdir
i=90
loadct,39,/sil
window,xs=sz[0],ys=sz[1],/pixmap
tvscl,reform(alog(prj[*,*,0]))
xyouts,sz[0]-100,sz[1]-30,/dev,'inc='+string(i,format='(I2.2)')+' deg'
snapshot,pngdir+'/'+'png_'+string(90-i,format='(I2.2)')

for i=89,30,-1 do begin 
   print, i
   rot3d = rot_3d(surfb3d,rotaxis=2,degree=90.-i,/interp)
   prj[*,*,90-i]=total(rot3d,3)
   tvscl,reform(alog(prj[*,*,90-i]))
   xyouts,sz[0]-100,sz[1]-30,/dev,'inc='+string(i,format='(I2.2)')+' deg'
   snapshot,pngdir+'/'+'png_'+string(90-i,format='(I2.2)')
endfor

save,file='rot_tmp.sav',prj

stop
end

