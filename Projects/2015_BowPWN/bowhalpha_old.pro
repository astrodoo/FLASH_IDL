pro bowhalpha

; read the data of ionization state in collisitional ionization equilibrium(CIE) from Mappings3
readcol,'H0001.cie',id,temph,h0,h1,format='(i,f,d,d)',/sil
h1fr_data = h1/(h0 + h1)

readcol,'He0001.cie',idhe,temphe,he0,he1,he2,format='(i,f,d,d,d)',/sil


; Osterbrock table 4.2 (pp 73)
emit_ne_np_data = [3.72d-25, 2.2d-25, 1.24d-25, 6.62d-26]
emit_T = [2500.,5000.,10000.,20000.]
jha_jhb_data = [3.3,3.05,2.87,2.76]

;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar2/'
;fname = 'PWN2d_hdf5_plt_cnt_0117'
dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_elmfist_0.5/'
fname = 'PWN2d_hdf5_plt_cnt_0337'



xrange=[-1.e17,2.5e18]
yrange=[0.,7.e17]
sample=3

dens = loaddata(dir+fname,'dens',sample=sample,xra=yrange,yra=xrange,xCoord=y,yCoord=x,time=time)
pres = loaddata(dir+fname,'pres',sample=sample,xra=yrange,yra=xrange)
dens = transpose(dens) & pres = transpose(pres)

sz = size(dens,/dimension)

temp = !unit.mh/!unit.k * pres/dens; *0.5 ; assuming mu=1

h1fr = interpol(h1fr_data,temph,temp)

emit_ne_np = interpol(emit_ne_np_data,emit_T,temp)
min_emit_ne_np = 1.d-27
emit_ne_np = emit_ne_np > min_emit_ne_np

jha_jhb = interpol(jha_jhb_data,emit_T,temp)
min_jha_jhb = 0.1
jha_jhb = jha_jhb > min_jha_jhb

nel = h1fr * dens/!unit.mh      ; nel = np

emit_ha = emit_ne_np * jha_jhb * nel * nel
dx = double(x[2])-double(x[1])
emit_tot = emit_ha * dx^3.d0
surfb = emit_tot / (4.d0*!pi*dx*dx)

loadct,0,/sil
swindow,xs=sz[0]+150,ys=sz[1]*4,x_scr=1400,y_scr=1000
tvscl,alog(dens),0
color_bar,lim=[min(dens),max(dens)],/log,pos=[sz[0],sz[1]*3,sz[0]+20,sz[1]*4],/right,bartitle='dens',titlegap=0.05
loadct,1,/sil
tvscl,alog(pres),1
color_bar,lim=[min(pres),max(pres)],/log,pos=[sz[0],sz[1]*2,sz[0]+20,sz[1]*3],/right,bartitle='pres',titlegap=0.05
loadct,3,/sil
tvscl,alog(temp),2
color_bar,lim=[min(temp),max(temp)],/log,pos=[sz[0],sz[1],sz[0]+20,sz[1]*2],/right,bartitle='temperature',titlegap=0.05
loadct,7,/sil
tvscl,surfb,3
color_bar,lim=[min(surfb),max(surfb)],pos=[sz[0],0,sz[0]+20,sz[1]],/right,bartitle='surface brightness',titlegap=0.06

stop
end
