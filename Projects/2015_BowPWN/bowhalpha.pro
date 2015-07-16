pro bowhalpha

; read the data of ionization state in collisitional ionization equilibrium(CIE) from Mappings3
readcol,'H0001.cie',id,temph,h0,h1,format='(i,f,d,d)',/sil
h1fr_data = h1/(h0 + h1)

readcol,'He0001.cie',idhe,temphe,he0,he1,he2,format='(i,f,d,d,d)',/sil


; Osterbrock table 4.2 (pp 73)
;emit_ne_np_data = [3.72d-25, 2.2d-25, 1.24d-25, 6.62d-26]
emit_T = [2500.,5000.,10000.,20000.]
rec_hb_data = [9.07d-14,5.37d-14,3.03d-14,1.62d-14]
jha_jhb_data = [3.3,3.05,2.87,2.76]

;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar2/'
;fname = 'PWN2d_hdf5_plt_cnt_0117'
;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_elmfist_0.5/'
dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_0.5/'

fname = 'PWN2d_hdf5_plt_cnt_0337'

xrange=[-1.e17,2.5e18]
yrange=[0.,7.e17]
sample=4 ;3

dens = loaddata(dir+fname,'dens',sample=sample,xra=yrange,yra=xrange,xCoord=y,yCoord=x,time=time,lref=lref)
pres = loaddata(dir+fname,'pres',sample=sample,xra=yrange,yra=xrange)
dens = transpose(dens) & pres = transpose(pres) & lref= transpose(lref)

sz = size(dens,/dimension)

temp = !unit.mh/!unit.k * pres/dens; *0.5 ; assuming mu=1

h1fr = interpol(h1fr_data,temph,temp)

;emit_ne_np = interpol(emit_ne_np_data,emit_T,temp)
;min_emit_ne_np = 1.d-27
;emit_ne_np = emit_ne_np > min_emit_ne_np

jha_jhb = interpol(jha_jhb_data,emit_T,temp)
min_jha_jhb = 0.1
jha_jhb = jha_jhb > min_jha_jhb

rec_hb = interpol(rec_hb_data,emit_T,temp)
min_rec_hb = 1.d-20
rec_hb = rec_hb > min_rec_hb

nel = h1fr * dens/!unit.mh      ; nel = np

;emit_ha = emit_ne_np * jha_jhb * nel * nel
dx = double(x[2])-double(x[1])
;emit_tot = emit_ha * dx^3.d0
;surfb = emit_tot / (4.d0*!pi*dx*dx)

Ehb = !unit.h*!unit.c/486.1d-9
EM = nel*nel*dx^3.d0
L_ha = EM * rec_hb * Ehb * jha_jhb

surfb = L_ha/(4.d0*!pi*dx*dx)

; substract background
backg = surfb[*,sz[1]-1]#replicate(1,sz[1])
lrefmin = min(lref[0,*])
yref_ind_dummy = where(lref[0,*] gt lrefmin) & yref_ind = yref_ind_dummy[n_elements(yref_ind_dummy)-1]
backg[*,0:yref_ind] = surfb[*,yref_ind]#replicate(1,yref_ind+1)

disx1 = 1.5e8*time - 3.945e18
disx2 = 1.5e8*time - 3.3e18

; fine tunning
disx1_ind = (where(x ge disx1))[0]
disx2_ind = (where(x ge disx2))[0]
backg_1 = surfb[disx1_ind/2,sz[1]-1]
backg_2 = surfb[(disx1_ind+disx2_ind)/2,sz[1]-1]
backg_3 = surfb[(disx2_ind + sz[0])/2,sz[1]-1]
backrow = fltarr(sz[0])
backrow[0:disx1_ind-1] =backg_1
backrow[disx1_ind:disx2_ind-1] =backg_2
backrow[disx2_ind:sz[0]-1] =backg_3
lref_dummy = lref[0,*]
lref_dummy[where(lref_dummy eq min(lref_dummy))] = max(lref_dummy)
lref2min = min(lref_dummy)
yref2_ind_dummy = where(lref[0,*] gt lref2min) & yref2_ind = yref2_ind_dummy[n_elements(yref2_ind_dummy)-1]
backg[*,0:yref2_ind] = backrow#replicate(1,yref2_ind+1)

surfb2 = surfb-backg
surfb2 = surfb2>0.

;dx = 7.e16
;disx1 = 1.43201e+18
;disx2 = 1.11001e+18
;inwall_x = where((x le disx1+dx) and (x ge disx2-dx))

;surfb[inwall_x,*] = surfb[inwall_x,*]*1.e-20

loadct,39,/sil
window,xs=sz[0],ys=sz[1]*2
tvcoord,surfb,x,y,pos=[0,sz[1]],/dev,/scale
tvcoord,surfb2,x,y,pos=[0,0],/dev,/scale
oplot,[disx1,disx1],!y.crange,line=2
oplot,[disx2,disx2],!y.crange,line=2


save,file=dir+'surfb.sav',x,y,(surfb=surfb2)
;loadct,0,/sil
;swindow,xs=sz[0]+150,ys=sz[1]*4,x_scr=1400,y_scr=1000
;tvscl,alog(dens),0
;color_bar,lim=[min(dens),max(dens)],/log,pos=[sz[0],sz[1]*3,sz[0]+20,sz[1]*4],/right,bartitle='dens',titlegap=0.05
;loadct,1,/sil
;tvscl,alog(pres),1
;color_bar,lim=[min(pres),max(pres)],/log,pos=[sz[0],sz[1]*2,sz[0]+20,sz[1]*3],/right,bartitle='pres',titlegap=0.05
;loadct,3,/sil
;tvscl,alog(temp),2
;color_bar,lim=[min(temp),max(temp)],/log,pos=[sz[0],sz[1],sz[0]+20,sz[1]*2],/right,bartitle='temperature',titlegap=0.05
;loadct,7,/sil
;tvscl,surfb,3
;color_bar,lim=[min(surfb),max(surfb)],pos=[sz[0],0,sz[0]+20,sz[1]],/right,bartitle='surface brightness',titlegap=0.06

stop
end


pro denshow

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_0.5/'

fname = 'PWN2d_hdf5_plt_cnt_0337'

xrange=[-1.e17,2.5e18]
yrange=[0.,7.e17]
sample=4 ;3

dens = loaddata(dir+fname,'dens',sample=sample,xra=yrange,yra=xrange,xCoord=y,yCoord=x,time=time,lref=lref)
dens = transpose(dens) 

sz = size(dens,/dimension)

dens2 = fltarr(sz[0],2*sz[1])
dens2[*,0:sz[1]-1] = reverse(dens,2)
dens2[*,sz[1]:sz[1]*2-1] = dens

y2 = fltarr(sz[1]*2)
y2[0:sz[1]-1] = -y
y2[sz[1]:sz[1]*2-1] = y

pltx0=200. & plty0=100.
winxs=sz[0]+pltx0+200 & winys=sz[1]*2+plty0+250
mkeps,'guitar_den',xs=20.,ys=20.*winys/winxs
;window,0,xs=winxs, ys=winys

tvcoord,alog10(dens2),x,y2,/scale,/axes,pos=[pltx0/winxs,plty0/winys],/norm,psx=sz[0]/winxs,xtitle='x [cm]',ytitle='y [cm]',color=0
color_bar,lim=[min(dens2),max(dens2)],/log,/right,bartitle='density',charsize=1.5,bargap=0.01,titlegap=0.11

plot,x,smooth(dens[*,178],10),/xst,pos=posnorm([pltx0,plty0+sz[1]*2,pltx0+sz[0],plty0+sz[1]*2+200],nx=winxs,ny=winys),/norm,/noerase $
    , xtickformat='(a1)',ytitle=textoidl('\rho_{0}'),yst=2

epsfree

stop
end
