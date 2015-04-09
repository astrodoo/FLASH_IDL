pro transto3d,fname,sample=sample, xrange=xrange, yrange=yrange

if not keyword_set(sample) then sample=4

;xra = [-2.17e17,2.51e18] & yra=[0.,8.e17]
d=loaddata(fname,'dens',sample=sample,xCoord=r,yCoord=z)   ; [r, z]
;d=loaddata(fname,'dens',sample=sample,xCoord=r,yCoord=z,xra=yra,yra=xra)   ; [r, z]

sd = size(d,/dimension)

nth = sd[0]
th = findgen(nth)/float(nth)*2.*!pi

;tmpd3d = replicate({dd:d},nth)

;i=400
stop
den3d_xyz = fltarr(sd[1],sd[0]*2,sd[0]*2)
for i=0,sd[1]-1 do begin
   print, i,'     of',sd[1]
   dcutz = reform(d[*,i]) # replicate(1.,nth)   ; [r, th]
   tv_polar,dcutz,r,th,xout=x,yout=y,imgout=dcutz_xy,/extrapol,/no_roff,/no_window

   den3d_xyz[i,*,*] = dcutz_xy > min(dcutz)
endfor
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

stop
end
