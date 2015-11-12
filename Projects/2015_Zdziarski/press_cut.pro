pro press_cut

id='3e37'

dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/re-coll/boundaryJet_3E37/'
fname='JetSet_hdf5_plt_cnt_3024'

bhx = -2.e12

sample=2
dx=6.e12/2.
dz=1.e13
xrange = [bhx-dx,bhx+dx] & yrange = [-dx,dx] & zrange = [0.,dz]

if keyword_set(zoom) then begin
   dx=9.6e11/2.
   dz=1.6e12 
   xrange = [bhx-dx,bhx+dx] & yrange = [-dx,dx] & zrange = [0.,dz]
   sample=2
endif

pres = loaddata(dir+fname,'pres',sample=sample,xCoord=x,yCoord=y,zCoord=z,time=time,xra=xrange,yra=yrange, zra=zrange)
psz = size(pres,/dimension)

maxp = 4.e3 & minp = 1.e-2

loadct,39,/sil
swindow,xs=psz[0]+300,ys=psz[1]+200
tvcoord,bytscl(reform(alog10(pres[*,*,200])),max=alog10(maxp),min=alog10(minp)),x,y,/axes
color_bar,lim=[minp,maxp],/log,/right,bartitle='pres'
    
stop
end
