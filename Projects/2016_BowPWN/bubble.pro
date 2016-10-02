pro bubble,step=step

if not keyword_set(step) then step=5

fnames = file_search('PWN2d_hdf5_plt_cnt_*')
nfile = n_elements(fnames)-1

xra=[0.,4.e18]
yra=[-5.e17,1.e19]

d = loaddata(fnames[0],'dens',time=t,xCoord=x,yCoord=y,sample=3)
;d = loaddata(fnames[0],'dens',time=t,xCoord=x,yCoord=y,xra=xra,yra=yra,sample=3)
sz = size(d,/dimension)

nout = ceil(float(nfile)/step)
dens = fltarr(sz[0],sz[1],nout)
time = fltarr(nout)

k=0
for i=0,nfile-1,step do begin
    print, fnames[i], ' of ', fnames[nfile-1], ' step:',step
;    d = loaddata(fnames[i],'dens',time=t,xra=xra,yra=yra,sample=3)
    d = loaddata(fnames[i],'dens',time=t,sample=3)

    dens[*,*,k] = d
    time[k] = t
    k = k+1
endfor

save,file='bubble_comb.dat',dens,time,x,y

stop
end

pro bubradi,ps=ps

den0 = 1.67e-24

;dirPWN2d = '/d/d7/yoon/out_FLASH3.3_mhd/out_PWN2d/out_PWN2d_600kms/'
dirPWN2d = ''
restore,file=dirPWN2d+'bubble_comb.dat'
d2=dens & t2=time & x2=x & y2=y

sz = size(d2,/dimension)
shout = fltarr(sz[2])
shin  = fltarr(sz[2])

;vv = 3.e7
;vv = 6.e7
vv = 9.e7

for i=0,sz[2]-1 do begin
   bubc = t2[i]*vv
   bubcInd_tmp = where(y2 ge bubc) & bubcInd = bubcInd_tmp[0]

   if (bubcInd eq -1) then bubcInd = sz[2]-1
   shinInd_tmp = where(d2[*,bubcInd,i] ge den0) & shinInd = shinInd_tmp[0]
   shin[i] = x2[shinInd]

   shoutInd_tmp = where(d2[*,bubcInd,i] eq den0) & shoutInd = shoutInd_tmp[0]
   shout[i] = x2[shoutInd]
endfor

Edot = 6.66d35
ana1 = (125./154./!pi)^0.2*(Edot/den0)^0.2*t2^0.6
ana2 = (25./14./!pi)^0.2*(Edot/den0)^0.2*t2^0.6

loadct,39,/sil
if not keyword_set(ps) then window,0 $
  else mkeps,'bubradi_PWN2d',xs=20.,ys=20.*6./8.

;plot,t2,shin,/xlog,/ylog,xtitle='time [s]',ytitle='radius [cm]',xra=[4.e8,1.5e11],/xst,yra=[6.e16,6.e18],/yst
plot,t2,shin,/xlog,/ylog,xtitle='time [s]',ytitle='radius [cm]',/yst,/xst
oplot,t2,shout
oplot,t2,ana1,line=2,color=50
oplot,t2,ana2,line=2,color=254

legend,[textoidl('25/14\pi'),textoidl('125/154\pi')],/left,/top,box=0,textcolor=[254,50],color=[254,50]

if keyword_set(ps) then epsfree

save,file='bubradi.sav',t2,shin,shout
;dir3d = '/d/d7/yoon/out_Flash_jet/out_guitar/'
;restore,file=dir3d+'bubble_comb_guitar.dat'
;dxy3=denxy & dxz3=denxz & x3=x & y3=y & z3=z & t3=t

stop
end

pro bubcomb,ps=ps

restore, file='/d/d7/yoon/out_FLASH3.3_mhd/out_PWN2d/out_PWN2d_300kms_l10/bubradi.sav'
t300=t2 & shin300=shin & shout300=shout 
restore, file='/d/d7/yoon/out_FLASH3.3_mhd/out_PWN2d/out_PWN2d_600kms/bubradi.sav'
t600=t2 & shin600=shin & shout600=shout 
restore, file='/d/d7/yoon/out_FLASH3.3_mhd/out_PWN2d/out_PWN2d_900kms/bubradi.sav'
t900=t2 & shin900=shin & shout900=shout 

ncut = n_elements(t900)-10
t900 = reform(t900[0:ncut]) & shin900=reform(shin900[0:ncut]) & shout900=reform(shout900[0:ncut])

den0 = 1.67e-24  & Edot = 6.66d35
;tana = t600
tmax = 2.e11 & tmin=5.e8
nt=300
tana = findgen(nt)/float(nt)*(tmax-tmin) + tmin
ana1 = (125./154./!pi)^0.2*(Edot/den0)^0.2*tana^0.6

loadct,39,/sil

if not keyword_set(ps) then begin
!p.background=255 & !p.color=0
window
endif else mkeps,'bubble_comb',xs=20.,ys=20.*6./8.

plot,tana,ana1,color=0,/xlog,/ylog, xra=[5.e8,2.e11], yra=[5.e16,5.e18],/xst,/yst, xtitle='time [s]', ytitle='radius [cm]',/nodata
oplot,t300,shout300,color=50
oplot,t300,shin300,color=50
oplot,t600,shout600,color=150
oplot,t600,shin600,color=150
oplot,t900,shout900,color=250
oplot,t900,shin900,color=250

polyfill,[t300,reverse(t300)],[shin300,reverse(shout300)],color=50,/line_fill,orientation=90
polyfill,[t600,reverse(t600)],[shin600,reverse(shout600)],color=150,/line_fill,orientation=45
polyfill,[t900,reverse(t900)],[shin900,reverse(shout900)],color=250,/line_fill

oplot,tana,ana1,color=0,thick=5

vtt = textoidl('v_{*} = ') + ['300','600','900'] + textoidl(' km s^{-1}')
legend, ['analytic model',vtt], color=[0,50,150,250],textcolor=[0,50,150,250],/left,/top,box=0,thick=[5,1,1,1],line=0

if keyword_set(ps) then epsfree
stop
end
