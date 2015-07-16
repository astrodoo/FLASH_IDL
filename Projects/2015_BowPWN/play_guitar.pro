pro play_guitar

device,decomposed=0

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_0.5/'

;fname = 'PWN2d_hdf5_plt_cnt_0337'
fname = 'PWN2d_hdf5_plt_cnt'
fnames = file_search(dir+fname+'*')
nfiles= n_elements(fnames)

xrange=[-1.e17,2.5e18]
yrange=[0.,7.e17]
sample=4 ;3

dens = loaddata(fnames[nfiles-1],'dens',sample=sample,xra=yrange,yra=xrange,xCoord=y,yCoord=x,time=t)
dens = transpose(dens) 

sz = size(dens,/dimension)
maxd = max(dens)
mind = min(dens)

pltx0=150. & plty0=80.
winxs=sz[0]+pltx0+130 & winys=sz[1]*2+plty0+220

dens2 = fltarr(sz[0],2*sz[1])
;dens2[*,0:sz[1]-1] = reverse(dens,2)
;dens2[*,sz[1]:sz[1]*2-1] = dens

y2 = fltarr(sz[1]*2)
y2[0:sz[1]-1] = -y
y2[sz[1]:sz[1]*2-1] = y


;making png
spawn,'mkdir '+dir+'png_data'

!p.background=255 & !p.color=0
!p.charsize=2.

for i= 0, nfiles-1 do begin

   print,i, ' of ', nfiles
   dens = loaddata(fnames[i],'dens',sample=sample,xra=yrange,yra=xrange,xCoord=y,yCoord=x,time=t)
   dens = transpose(dens) 
   dens2[*,0:sz[1]-1] = reverse(dens,2)
   dens2[*,sz[1]:sz[1]*2-1] = dens

   denprof,x,t,ambient=amb,ch1x=ch1x,ch2x=ch2x,ch3x=ch3x,ch4x=ch4x

   loadct,3,/sil
   window,0,xs=winxs, ys=winys ,/pixmap

   tvcoord,bytscl(alog10(dens2),min=alog10(mind),max=alog10(maxd)),x,y2,/axes,pos=[pltx0,plty0],/dev,xtitle='x [cm]',ytitle='y [cm]',color=0
   legend,'t= '+string(fix(t/60./60./24./365.),format='(I4)')+' yr',/left,/top,textcolor=0,box=0
   oplot,[ch1x,ch1x],!y.crange,line=1
   oplot,[ch2x,ch2x],!y.crange,line=1
   oplot,[ch4x,ch4x],!y.crange,line=1
   color_bar,lim=[mind,maxd],/log,/right,bartitle='density',charsize=1.5,bargap=0.01,titlegap=0.07

   ;plot,x,smooth(dens[*,178],10),/xst,pos=[pltx0,plty0+sz[1]*2,pltx0+sz[0],plty0+sz[1]*2+200],/dev,/noerase $
   ;    , xtickformat='(a1)',ytitle=textoidl('\rho_{0}'),yst=2,color=0.
   plot,x,amb,/xst,pos=[pltx0,plty0+sz[1]*2,pltx0+sz[0],plty0+sz[1]*2+200],/dev,/noerase $
       , xtickformat='(a1)',ytitle=textoidl('\rho_{0}'),color=0.,yra=[0.,1.8e-25],/yst
   oplot,[ch1x,ch1x],!y.crange,line=1
   oplot,[ch2x,ch2x],!y.crange,line=1
   oplot,[ch4x,ch4x],!y.crange,line=1
  
   denp = (amb[where(x ge 0.)])[0]
   plots,0,denp,psym=dsym(8,/fill),symsize=2,color=trp_color(fsc_color('blue'),alpha=0.3,/silent,/whitebg)
   plots,0,denp,psym=dsym(8),symsize=2,color=fsc_color('blue')

   snapshot,dir+'png_data/guitar_'+string(i,format='(I4.4)')
endfor

stop
end



pro denprof,x,t,ambient=dd,ch1x=ch1x,ch2x=ch2x,ch3x=ch3x,ch4x=ch4x
  nx = n_elements(x)
  dd = fltarr(nx)

  vp = 1.5e8

  ch1x = vp*t-3.3e18
  ch2x = vp*t-3.8e18
  ch3x = vp*t-3.95e18
  ch4x = vp*t-4.e18

  x1_ind = where(x ge ch1x,ct1)
  x2_ind = where((x ge ch2x) and (x lt ch1x),ct2)
  x3_ind = where((x ge ch3x) and (x lt ch2x),ct3)
  x4_ind = where((x ge ch4x) and (x lt ch3x),ct4)
  x5_ind = where(x lt ch4x,ct5)

  if (ct1 ne 0) then dd[x1_ind] = 1.38e-25
  if (ct2 ne 0) then dd[x2_ind] = 2.76e-26
  if (ct3 ne 0) then dd[x3_ind] = (x[x3_ind]-ch2x)/(ch3x-ch2x)*(1.38e-26-2.76e-26) + 2.76e-26
  if (ct4 ne 0) then dd[x4_ind] = 1.38e-26
  if (ct5 ne 0) then dd[x5_ind] = 1.38e-25
end

