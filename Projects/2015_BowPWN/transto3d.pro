pro transto3d, hdf5=hdf5, nth=nth

;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar2/' 
;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_elmfist_0.5/'
;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_0.5/'
;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar/'
dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_60inc/'
;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_30inc/'
restore,dir+'phot.sav'
fname = 'phot'
yy=y

sz = size(phot,/dimension)

if not keyword_set(nth) then nth = sz[1]
th = findgen(nth)/float(nth)*2.*!pi

phot3d = fltarr(sz[0],sz[1]*2,sz[1]*2)
for i=0,sz[0]-1 do begin
   print, i,'     of',sz[0]
   dcutx_rth = reform(phot[i,*]) # replicate(1.,nth)   ; [r, th]
   tv_polar,dcutx_rth,yy,th,xout=y,yout=z,imgout=dcutx_yz,/extrapol,/no_roff,/no_window

   phot3d[i,*,*] = dcutx_yz > min(phot[i,*])
endfor

save,file=dir+fname+'_3d.sav',phot3d,x,y,z

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
