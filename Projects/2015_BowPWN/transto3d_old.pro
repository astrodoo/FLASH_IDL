pro transto3d,n,var=var,sample=sample, xrange=xrange, yrange=yrange, hdf5=hdf5, nth=nth

; For guitar nebula ============================================================================
n=337 ;117
xrange=[-1.e17,2.5e18]
yrange=[0.,7.e17]
sample=4 ;3
; ==============================================================================================
if keyword_set(chk) then begin
   fname = 'PWN2d_hdf5_chk_'+string(n,format='(I4.4)') 
   print, 'read check point file: '+fname   
endif else begin 
   fname = 'PWN2d_hdf5_plt_cnt_'+string(n,format='(I4.4)')
   chk=0
endelse
if not keyword_set(sample) then sample=4
if not keyword_set(var) then var='dens'

if not keyword_set(xrange) then str_xra = '' $
  else str_xra = 'yra=['+strtrim(xrange[0],2)+','+strtrim(xrange[1],1)+']'
if not keyword_set(yrange) then str_yra = '' $
  else str_yra = 'xra=['+strtrim(yrange[0],2)+','+strtrim(yrange[1],1)+']'

if (keyword_set(xrange) and keyword_set(yrange)) then str_yra = ','+str_yra

str_xyra = str_xra+str_yra
if (strmid(str_xyra,0,1,/reverse) eq ']') then str_xyra = str_xyra+','

strexe = execute("data = loaddata(fname,'"+var+"',"+str_xyra+"sample=sample,xCoord=yy,yCoord=x,time=time)")

data=transpose(data)
sz = size(data,/dimension)

if not keyword_set(nth) then nth = sz[1]
th = findgen(nth)/float(nth)*2.*!pi

;den3d_xyz = fltarr(sz[0],sz[1],sz[1])
den3d_xyz = fltarr(sz[0],sz[1]*2,sz[1]*2)
for i=0,sz[0]-1 do begin
   print, i,'     of',sz[0]
   dcutx_rth = reform(data[i,*]) # replicate(1.,nth)   ; [r, th]
   tv_polar,dcutx_rth,yy,th,xout=y,yout=z,imgout=dcutx_yz,/extrapol,/no_roff,/no_window

   den3d_xyz[i,*,*] = dcutx_yz > min(data[i,*])
endfor

save,file=fname+'_3d.sav',den3d_xyz,x,y,z


if keyword_set(hdf5) then begin
; write hdf5 file
print,'start to write hdf5'
den3d_xyz = transpose(den3d_xyz)
sz_den3d = size(den3d_xyz,/dimension)
fileID = h5f_create(fname+'_smp4_3dxyz')
datatypeID = H5T_IDL_CREATE(den3d_xyz)
dataspaceID = H5S_CREATE_SIMPLE([sz_den3d[0],sz_den3d[1],sz_den3d[2]])
datasetID = H5D_CREATE(fileID, 'Density', datatypeID, dataspaceID)
H5D_WRITE, datasetID, den3d_xyz

;uniform_coord = [[min(x),max(x)],[min(y),max(y)],[min(z),max(z)]]
;uniform_coord = transpose(uniform_coord)
;sz_coord = size(uniform_coord,/dimension)
;datatypeID2 = H5T_IDL_CREATE(uniform_coord)
;dataspaceID2 = H5S_CREATE_SIMPLE([sz_coord[0],sz_coord[1]])
;datasetID2 = H5D_CREATE(fileID, 'Coordinate', datatypeID2, dataspaceID2)
;H5D_WRITE, datasetID2, uniform_coord

H5F_CLOSE, fileID
endif

stop
end
