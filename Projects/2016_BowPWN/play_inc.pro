pro play_inc
device,decomposed=0

dir = '/home/jianiye/Work/Data/2015_BowPWN/PWN3d_05inc/'
fname = 'cut2d.sav'

restore,file=dir+fname

nds = n_elements(ds)
dsz = size(ds.d,/dimension)

outdir = 'png_PWN3d_05inc'
spawn, 'mkdir '+outdir

maxd = 7.5e-23 & mind = 3.e-29

t_enc = 2.e18/1.e8

loadct,0,/sil
!p.charsize=2.
!p.background=255
!p.color=0
px0=150 & py0=70
window,xs=dsz[0]+px0+180,ys=dsz[1]+py0+20,/pixmap

for i=0,nds-1 do begin
   print, i, ' of ', nds

   tvcoord,bytscl(alog10(ds[i].d),max=alog10(maxd),min=alog10(mind)),x,y,/axes,/black,xtitle='x [cm]',ytitle='y [cm]',pos=[px0,py0]
   arrow,-3.e17,0.,-6.e17,0.,/data,thick=3

   legend,string((ds[i].t-t_enc)/!unit.year,format='(f7.2)')+ ' yr',/left,/top,box=0,textcolor=fsc_color('yellow')

   loadct,0,/sil
   color_bar,lim=[mind,maxd],/right,bartitle=textoidl('density [g cm^{-3}]'),titlegap=0.09,bthick=0.02
   snapshot,outdir+'/'+'PWN3d_05inc_'+string(i,format='(I3.3)')
endfor

stop
end
