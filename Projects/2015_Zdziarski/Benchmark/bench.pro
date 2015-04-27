pro bench

restore,'scaling_aci320.sav'
ncpu_aci320 = nproc
cputime_aci320 = cputime - cputime[0]
nblk_aci320 = nblk
step_aci320 = step - step[0]
phytime_aci320 = phytime - phytime[0]

restore,'scaling_sta320.sav'
ncpu_sta320 = nproc
cputime_sta320 = cputime - cputime[0]
nblk_sta320 = nblk
step_sta320 = step - step[0]
phytime_sta320 = phytime -phytime[0]

restore,'scaling_sta640.sav'
ncpu_sta640 = nproc
cputime_sta640 = cputime - cputime[0]
nblk_sta640 = nblk
step_sta640 = step - step[0]
phytime_sta640 = phytime - phytime[0]

restore,'scaling_sta960.sav'
ncpu_sta960 = nproc
cputime_sta960 = cputime -cputime[0]
nblk_sta960 = nblk
step_sta960 = step - step[0]
phytime_sta960 = phytime - phytime[0]

print,'ncpu aci320 = ', ncpu_aci320
print,'ncpu sta320 = ', ncpu_sta320
print,'ncpu sta640 = ', ncpu_sta640
print,'ncpu sta960 = ', ncpu_sta960

loadct,39,/sil
;!p.background=255 & !p.color=0
;window,0

mkeps,'scaling1',xs=20,ys=20.*6./8.

plot,step_aci320, cputime_aci320, xtitle='nstep - nstep0', ytitle='cpu time [sec]'
oplot,step_sta320, cputime_sta320,color=50
oplot,step_sta640, cputime_sta640,color=200
oplot,step_sta960, cputime_sta960,color=250
legend,['aci320','sta320','sta640','sta960'],line=0,color=[0,50,200,250],textcolor=[0,50,200,250],/left,/top,box=0

epsfree

xcrit = 100.
xind_aci320 = (where(phytime_aci320 ge xcrit))[0]
fit_aci320 = linfit(phytime_aci320[xind_aci320:*],cputime_aci320[xind_aci320:*])

xind_sta320 = (where(phytime_sta320 ge xcrit))[0]
fit_sta320 = linfit(phytime_sta320[xind_sta320:*],cputime_sta320[xind_sta320:*])

xind_sta640 = (where(phytime_sta640 ge xcrit))[0]
fit_sta640 = linfit(phytime_sta640[xind_sta640:*],cputime_sta640[xind_aci320:*])

xind_sta960 = (where(phytime_sta960 ge xcrit))[0]
fit_sta960 = linfit(phytime_sta960[xind_sta960:*],cputime_sta960[xind_sta960:*])

print, '** fitting result [cputime / physical time] **'
print, 'aci320 : ', fit_aci320[1] 
print, 'sta320 : ', fit_sta320[1] 
print, 'sta640 : ', fit_sta640[1] 
print, 'sta960 : ', fit_sta960[1] 

mkeps,'scaling2',xs=20.,ys=20.*6./8.
plot,phytime_aci320, cputime_aci320, xtitle='physical time in code unit [sec]', ytitle='cpu time [sec]'
;oplot,[100.,200.],[fit_aci320[1]*100+fit_aci320[0], fit_aci320[1]*200+fit_aci320[0]],color=254
oplot,phytime_sta320, cputime_sta320,color=50
oplot,phytime_sta640, cputime_sta640,color=200
oplot,phytime_sta960, cputime_sta960,color=250
;oplot,!x.crange,[IO_cputime[1]-cputime[0],IO_cputime[1]-cputime[0]], line=2
legend,['aci320','sta320','sta640','sta960'],line=0,color=[0,50,200,250],textcolor=[0,50,200,250],/left,/top,box=0

epsfree

mkeps,'scaling3',xs=20.,ys=20.*6./8.
plot,phytime_aci320, cputime_aci320 * ncpu_aci320, xtitle='physical time in code unit [sec]', ytitle='cpu time * ncpu'
oplot,phytime_sta320, cputime_sta320 * ncpu_sta320,color=50
oplot,phytime_sta640, cputime_sta640 * ncpu_sta640,color=200
oplot,phytime_sta960, cputime_sta960 * ncpu_sta960,color=250

;oplot,!x.crange,[IO_cputime[1]-cputime[0],IO_cputime[1]-cputime[0]], line=2
legend,['aci320','sta320','sta640','sta960'],line=0,color=[0,50,200,250],textcolor=[0,50,200,250],/left,/top,box=0

epsfree







stop
end
