pro press_cut

id='3e37'

dir='/home/jianiye/Work/Data/2015_Zdziarski/boundaryJet_3E37/'
fname='JetSet_hdf5_plt_cnt_3024'

bhx = -2.e12

sample=3
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

outdir= 'png_press_cut'
spawn,'mkdir '+outdir
zend=8.e12
zend_ind = (where(z ge zend))[0]
nzcut = 200
mult = zend_ind/nzcut
window,xs=psz[0]+230,ys=psz[1]+150
for i=0,nzcut-1 do begin
    print,i,' of ',nzcut
    pres_cut = reform(pres[*,*,i*mult])
    tvcoord,bytscl(reform(alog10(pres[*,*,i*mult])),max=alog10(maxp),min=alog10(minp)),x,y,/axes,xtitle='x [cm]',ytitle='y [cm]'
    legend,'z='+string(z[i*mult],format='(e10.2)')+' cm',/right,/top,textcolor=0,box=0
    color_bar,lim=[minp,maxp],/log,/right,bartitle='pressure'
    snapshot,outdir+'/'+'png_pxycut_'+string(i,format='(I3.3)')
endfor
    
stop
end
