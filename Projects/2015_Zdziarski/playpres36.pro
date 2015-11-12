pro playpres36

id='1e36'

dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/sphere_1e36/'
fnames = file_search(dir+'JetSet_hdf5_plt_cnt*')

;nn = fix(strmid(fnames,3,4,/rev))
;startn = 250
;fnames_id = where(nn ge startn)  

;fnames = fnames[fnames_id]
nfiles = n_elements(fnames)

bhx = -2.e12

sample=1
xrange=[-5.e12,1e12]
yrange= [0.,0.]
zrange= [0.,7.e12]

spawn,'mkdir png_playpres36'
outdir='png_playpres36/'

for i=0, nfiles-1 do begin

;i=257

print,i,' of ', nfiles

;read_amr,fnames[i],/nodata,param=param,tree=tree
;sample1 = max([max(tree.lrefine)-7,0])
;sample2 = max([sample1-1,0])

preszoom = reform(loaddata(fnames[i],'pres',sample=sample,xra=xrange,yra=yrange,zra=zrange,xCoord=xx,yCoord=yy,zCoord=zz,time=time))
;pres = reform(loaddata(fnames[i],'pres',sample=sample1,yra=[0.,0.],xCoord=x,yCoord=y,zCoord=z,time=time))
;sz = size(pres,/dimension)

;zcut_ind = (where(z ge 0))[0] 
;nz2 = (sz[1]-zcut_ind)*2
;z2 = fltarr(nz2) & pres2 = fltarr(sz[0],nz2)
;z2[0:nz2/2-1] = -reverse(z[zcut_ind:*]) & z2[nz2/2:*] = z[zcut_ind:*]
;pres2[*,0:nz2/2-1] = reverse(pres[*,zcut_ind:*],2) & pres2[*,nz2/2:*] = pres[*,zcut_ind:*]


szz = size(preszoom,/dimension)

;minp = min(pres) & maxp = max(pres)
minp = 1.12050   & maxp= 1655.45 

loadct,39,/sil
!p.charsize=2.
pltx0=150. & plty0=100
window,xs=pltx0+szz[0]+120,ys=plty0+szz[1]+50

;tvcoord,bytscl(alog10(pres2),min=alog10(minp),max=alog10(maxp)),x-bhx,z2,/axes,xtitle='x [cm]', ytitle='z [cm]',position=[pltx0,plty0],/dev
tvcoord,bytscl(alog10(preszoom),min=alog10(minp),max=alog10(maxp)),xx-bhx,zz,/axes,position=[pltx0,plty0],/dev
xyouts,pltx0-130,plty0+szz[1],/dev,'P'+id,color=255
tj = 9.76e4
xyouts,pltx0-130,plty0+szz[1]-20,/dev,'t='+string((time-tj)/60./60.,format='(f5.2)')+' hrs',color=255

loadct,39,/sil
color_bar,lim=[minp,maxp],/log,/right,bartitle='pressure'

if (time ge tj) then $

snapshot,outdir+'playpres_P'+id+'_'+string(i,format=('(I3.3)'))

endfor

stop
end

