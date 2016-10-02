pro bowhalpha

; read the data of ionization state in collisitional ionization equilibrium(CIE) from Mappings3
readcol,'H0001.cie',id,temph,h0,h1,format='(i,f,d,d)',/sil
h1fr_data = h1/(h0 + h1)

;readcol,'He0001.cie',idhe,temphe,he0,he1,he2,format='(i,f,d,d,d)',/sil


; Osterbrock table 4.2 (pp 73)
;emit_ne_np_data = [3.72d-25, 2.2d-25, 1.24d-25, 6.62d-26]
;emit_T = [2500.,5000.,10000.,20000.]
;rec_hb_data = [9.07d-14,5.37d-14,3.03d-14,1.62d-14]
;jha_jhb_data = [3.3,3.05,2.87,2.76]

incid = 90

if (incid eq 90) then begin
   dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar/' 
   inc = 90.*!dtor
endif else if (incid eq 60) then begin
   dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_60inc/'
   inc = 60.*!dtor
endif
;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_30inc/'

fname = 'PWN2d_hdf5_plt_cnt_0343'

;inc = 60. *!dtor
;inc = 90. *!dtor
xrange=[-1.e17,2.5e18]/sin(inc)
yrange=[0.,7.e17]
sample=5 
;sample=4

dens = loaddata(dir+fname,'dens',sample=sample,xra=yrange,yra=xrange,xCoord=y,yCoord=x,time=time,lref=lref)
pres = loaddata(dir+fname,'pres',sample=sample,xra=yrange,yra=xrange)
dens = transpose(dens) & pres = transpose(pres) & lref= transpose(lref)

sz = size(dens,/dimension)

temp = !unit.mh/!unit.k * pres/dens *0.1 

h1fr = interpol(h1fr_data,temph,temp)

;h1fr = replicate(1.,sz[0],sz[1])

rec_h = 2.59d-13*(temp/1.d4)^(-0.845)   ; [cm^3 s^-1]

nel = h1fr * dens/!unit.mh      ; nel = np

dx = double(x[2])-double(x[1])

EM = nel*nel*dx^3.d0
;dst = 1.9d3*!unit.pc
;phot = EM * rec_h / (4.d0*!pi* dst*dst)
phot = EM * rec_h / (4.d0*!pi* dx*dx)

vs = 1.5e8 /sin(inc)
disx1 = vs*time - 3.3e18/sin(inc)
disx2 = vs*time - 3.913e18/sin(inc)
disx3 = vs*time - 4.43e18/sin(inc)

disx1_ind = (where(x ge disx1))[0]
disx2_ind = (where(x ge disx2))[0]
disx3_ind = (where(x ge disx3))[0]

ldisx3 = reform(lref[disx3_ind,*])
minl = min(ldisx3)
y1_ind = where(ldisx3 gt minl)   & y1_ind = y1_ind[n_elements(y1_ind)-1]
y2_ind = where(ldisx3 gt minl+1) & y2_ind = y2_ind[n_elements(y2_ind)-1]
y3_ind = where(ldisx3 gt minl+2) & y3_ind = y3_ind[n_elements(y3_ind)-1]

ldisx2 = reform(lref[disx2_ind,*])
minl2 = min(ldisx2)
y4_ind = where(ldisx2 gt minl2)   & y4_ind = y4_ind[n_elements(y4_ind)-1]
y5_ind = where(ldisx2 gt minl2+1) & y5_ind = y5_ind[n_elements(y5_ind)-1]
y6_ind = where(ldisx2 gt minl2+2) & y6_ind = y6_ind[n_elements(y6_ind)-1]

ldisx1 = reform(lref[disx1_ind,*])
minl1 = min(ldisx1)
y7_ind = where(ldisx1 gt minl1)   & y7_ind = y7_ind[n_elements(y7_ind)-1]
y8_ind = where(ldisx1 gt minl1+1) & y8_ind = y8_ind[n_elements(y8_ind)-1]
y9_ind = where(ldisx1 gt minl1+2) & y9_ind = y9_ind[n_elements(y9_ind)-1]

loadct,39,/sil
window,0,xs=sz[0],ys=sz[1]*2
tvcoord,alog(dens),x,y,/scale,pos=[0,sz[1]]
oplot,[disx1,disx1],!y.crange,line=2
oplot,[disx2,disx2],!y.crange,line=2
oplot,[disx3,disx3],!y.crange,line=2
oplot,!x.crange,[y[y1_ind],y[y1_ind]],line=2
oplot,!x.crange,[y[y2_ind],y[y2_ind]],line=2
;oplot,!x.crange,[y[y3_ind],y[y3_ind]],line=2
tvcoord,lref,x,y,/scale,pos=[0,0]
oplot,[disx1,disx1],!y.crange,line=2
oplot,[disx2,disx2],!y.crange,line=2
oplot,[disx3,disx3],!y.crange,line=2
oplot,!x.crange,[y[y1_ind],y[y1_ind]],line=2
oplot,!x.crange,[y[y2_ind],y[y2_ind]],line=2
;oplot,!x.crange,[y[y3_ind],y[y3_ind]],line=2

;===============================  substract background ============================
; for 90 inc
if (incid eq 90) then begin
   backg = phot[*,sz[1]-1]#replicate(1,sz[1])
   backg[*,0:y1_ind] = phot[*,y1_ind]#replicate(1,y1_ind+1)
   backg[0:disx2_ind,0:y2_ind] = phot[0:disx2_ind,y2_ind]#replicate(1,y2_ind+1)
   backg[0:disx3_ind+20,0:y3_ind] = phot[0:disx3_ind+20,y3_ind]#replicate(1,y3_ind+1)
endif else if (incid eq 60) then begin
; for 60 inc
   backg = phot[*,sz[1]-1]#replicate(1,sz[1])
   backg[*,0:y1_ind] = phot[*,y1_ind]#replicate(1,y1_ind+1)
   backg[0:disx2_ind,0:y5_ind] = phot[0:disx2_ind,y5_ind]#replicate(1,y5_ind+1)
   backg[0:disx3_ind+20,0:y2_ind] = phot[0:disx3_ind+20,y2_ind]#replicate(1,y2_ind+1)
endif

; for 30 inc
;backg = surfb[*,sz[1]-1]#replicate(1,sz[1])
;backg[*,0:y7_ind] = surfb[*,y7_ind]#replicate(1,y7_ind+1)
;;backg[0:disx2_ind,0:y5_ind] = surfb[0:disx2_ind,y5_ind]#replicate(1,y5_ind+1)
;;backg[0:disx3_ind+20,0:y2_ind] = surfb[0:disx3_ind+20,y2_ind]#replicate(1,y2_ind+1)

phot2 = phot-backg
phot2 = phot2>0.
phot2_dum = phot2 & phot2_dum[where(phot2 eq 0)] = max(phot2)
phot2[where(phot2 eq 0)] = min(phot2_dum)

;============================= fine-tuning ===================================== 
; for 90 inc
if (incid eq 90) then begin
   phot2[191:200,35:61] = min(phot2_dum)
   phot2[103:194,31:102] = min(phot2_dum)
   phot2[192:220,65:170] = min(phot2_dum)
   phot2[231:351,95:147] = min(phot2_dum)
   phot2[530:555,109:127] = min(phot2_dum)
   phot2[514:536,120:137] = min(phot2_dum)
   phot2[554:621,63:145] = min(phot2_dum)
   phot2[613:sz[0]-1,69:143] = min(phot2_dum)
   phot2[*,133:sz[1]-1] = min(phot2_dum)
endif else if (incid eq 60) then begin
; for 60 inc
   phot2[227:240,49:73]= min(phot2_dum)
   phot2[216:233,35:57]= min(phot2_dum)
   phot2[0:221,33:sz[1]-1]= min(phot2_dum)
   phot2[227:247,57:143]= min(phot2_dum)
   phot2[236:256,65:127]= min(phot2_dum)
   phot2[320:435,98:126]= min(phot2_dum)
   phot2[0:20,*]= min(phot2_dum)
   phot2[626:sz[0]-1,73:sz[1]-1]= min(phot2_dum)
   phot2[*,126:sz[1]-1]= min(phot2_dum)
endif

; for 30 inc
;phot2[322:427,63:137] = min(phot2_dum)
;phot2[375:403,44:78] = min(phot2_dum)
;phot2[*,122:sz[1]-1] = min(phot2_dum)
;phot2[104:136,20:144] = min(phot2_dum)
;phot2[1027:1090,94:134] = min(phot2_dum)
;phot2[990:1012,111:129] = min(phot2_dum)
;phot2[730:851,89:125] = min(phot2_dum)
;phot2[554:683,93:129] = min(phot2_dum)
;phot2[1052:sz[0]-1,82:sz[1]-1] = min(phot2_dum)


window,2,xs=sz[0],ys=sz[1]*3
tvcoord,alog(phot),x,y,/scale,pos=[0,2*sz[1]]
oplot,[disx1,disx1],!y.crange,line=2,color=255
oplot,[disx2,disx2],!y.crange,line=2,color=255
oplot,[disx3,disx3],!y.crange,line=2,color=255

loadct,39,/sil
tvcoord,alog(backg),x,y,/scale,pos=[0,sz[1]]
oplot,[disx1,disx1],!y.crange,line=2,color=255
oplot,[disx2,disx2],!y.crange,line=2,color=255
oplot,[disx3,disx3],!y.crange,line=2,color=255

loadct,39,/sil
tvcoord,alog(phot2),x,y,/scale,pos=[0,0]
oplot,[disx1,disx1],!y.crange,line=2,color=255
oplot,[disx2,disx2],!y.crange,line=2,color=255
oplot,[disx3,disx3],!y.crange,line=2,color=255
oplot,!x.crange,[y[y1_ind],y[y1_ind]],line=2
oplot,!x.crange,[y[y2_ind],y[y2_ind]],line=2
;oplot,!x.crange,[y[y3_ind],y[y3_ind]],line=2


save,file=dir+'phot.sav',x,y,(phot=phot2)

px0 = 100. & py0= 80.
winx = px0+sz[0]+30 & winy = py0+sz[1]+200
window,1,xs=winx, ys=winy
tvcoord,alog10(phot2),x,y,/scale,pos=[px0,py0],/dev,/axes
color_bar,lim=[min(phot2),max(phot2)],/log,/up

stop
end


pro denshow

;dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar_0.5/'

dir = '/d/d2/yoon/out_FLASH3.3_mhd/out_PWNe/out_PWN2d/out_PWN2d_guitar/' 

fname = 'PWN2d_hdf5_plt_cnt_0343'

xrange=[-1.e17,2.5e18]
yrange=[0.,7.e17]
sample=5 

dens = loaddata(dir+fname,'dens',sample=sample,xra=yrange,yra=xrange,xCoord=y,yCoord=x,time=time,lref=lref)
dens = transpose(dens) 

sz = size(dens,/dimension)

dens2 = fltarr(sz[0],2*sz[1])
dens2[*,0:sz[1]-1] = reverse(dens,2)
dens2[*,sz[1]:sz[1]*2-1] = dens

y2 = fltarr(sz[1]*2)
y2[0:sz[1]-1] = -y
y2[sz[1]:sz[1]*2-1] = y

loadct,0,/sil
pltx0=200. & plty0=100.
winxs=sz[0]+pltx0+200 & winys=sz[1]*2+plty0+220
mkeps,'guitar_den',xs=20.,ys=20.*winys/winxs

tvcoord,alog10(dens2),x,y2,/scale,/axes,pos=[pltx0/winxs,plty0/winys],/norm,psx=sz[0]/winxs,xtitle='x [cm]',ytitle='y [cm]',color=0
color_bar,lim=[min(dens2),max(dens2)],/log,/right,bartitle='density',charsize=1.5,bargap=0.01,titlegap=0.11

;plot,x,smooth(dens[*,178],10),/xst,pos=posnorm([pltx0,plty0+sz[1]*2,pltx0+sz[0],plty0+sz[1]*2+200],nx=winxs,ny=winys),/norm,/noerase $
;    , xtickformat='(a1)',ytitle=textoidl('\rho_{0}'),yst=2,/ylog

vs = 1.5e8 
disx1 = vs*time - 3.3e18
disx2 = vs*time - 3.913e18
disx3 = vs*time - 4.43e18

d0=1.24e-25
d1=2.47e-26
d2=1.24e-26
delt=1.3e17
d3=1.24e-25

backd = fltarr(sz[0])
backd[where(x ge disx1)]=d0
backd[where((x ge disx2) and (x lt disx1))]=d1
backd[where((x ge disx3) and (x lt disx2))]=(d2-d1)/delt*((disx2-x[where((x ge disx3) and (x lt disx2))])<delt) + d1
backd[where(x lt disx3)]=d0

plot,x,backd,/xst,pos=posnorm([pltx0,plty0+sz[1]*2,pltx0+sz[0],plty0+sz[1]*2+200],nx=winxs,ny=winys),/norm,/noerase $
    , xtickformat='(a1)',ytitle=textoidl('\rho_{0}'),/ylog,yra=[8.e-27,3.e-25],/yst

epsfree

stop
end
