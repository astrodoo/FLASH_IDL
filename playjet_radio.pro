pro playjet_radio,fname,var=var,startn=startn,endn=endn,step=step,png=png
device, decomposed=0

if not keyword_set(fname) then fname='JetSet_hdf5_plt_cnt_'
if not keyword_set(var) then var='dens'
if not keyword_set(step) then step=1

fnames = file_search(fname+'*')
nfiles = n_elements(fnames)

if keyword_set(startn) then begin
   str_start = string(startn,format='(I4.4)')
   start_j_tmp = where(strpos(fnames,str_start) ne -1)
   start_j = start_j_tmp[0]
   if start_j eq -1 then begin
      print,'out of range of startnumber'
      stop
   endif
endif else  begin
   start_j = 0 & startn = 0
endelse

if keyword_set(endn) then begin
   str_end = string(endn,format='(I4.4)')
   end_j_tmp = where(strpos(fnames,str_end) ne -1)
   end_j = end_j_tmp[0]
   if end_j eq -1 then begin
      print,'out of range of startnumber'
      stop
   endif
endif else end_j = nfiles-1
nfiles = end_j - start_j + 1

; Making PNG files
if keyword_set(png) then begin
   spawn,'mkdir png_radio'
   opt=',/pixmap' 
endif else opt=''

;maxd = 2.2d-13 & mind = 1.1d-15
;maxd = 2.76d-13 & mind = 2.9d-16
;maxd = 2.76d-13 & mind = 1.d-16
maxd = 8.e-14 & mind = 3.d-18

m1 = 4.d34
m2 = 2.d34
mtot = m1+m2
peri = 3.d12
th0 = !pi
;th = th0 + sqrt(!unit.g * mtot / peri^3.d0) * time
th = th0
; center of mass, position of BB
cm2x = peri*m1/mtot*cos(th)
cm2y = peri*m1/mtot*sin(th)

xc0 = cm2x        ; position -> jet centered
yc0 = cm2y
zc0 = 0.

jeton = 9.76e4
;--------------------------------------------------------------------------------------------
; parameters
pw = 2.5d0
gammin = 1.d0
;nu = 100.d9
nu = 1.d9
wnu = 2.d0*!pi*nu
d   = 1.d4 *!unit.pc
;--------------------------------------------------------------------------------------------
maxv_sufb = 1.e-3
minv_sufb = 1.e-5

!p.background=255
for j=start_j,end_j,step do begin
    print,j+1,'   of',end_j+1,'  ', fnames[j]
    if not keyword_set(sample) then begin
       read_amr,fnames[j],var='dens',tree=tree
       smp = max(tree.lrefine) - 7
    endif else smp = sample
    pre = dload(startn+j-start_j,var='pres',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=smp,time)
    jet = dload(startn+j-start_j,var='jet',xc=xc0,yc=yc0,zc=zc0,sample=smp)
;    d = dload(j+startn,var=var,xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=smp,time) 

    time = time - jetOn

; assume equipartition (B^2 / 8 pi = 3 pres)
    B = sqrt(3.d0*pre * 8.d0*!pi)

    E_tot = pre * jet    ;/ denxyz * !unit.mh
    CC = E_tot/!unit.me/!unit.c/!unit.c *(-2.d0+pw)/gammin^(2.d0-pw) 

; Rybicki p180
    sina = 1.d0 ; sin(a) = 1 ; maximum value
    P_tot = sqrt(3.d0)*!unit.e^3.d0*CC*B*sina/(2.d0*!pi*!unit.me*!unit.c*!unit.c*(pw+1.d0))$
          *gamma(pw/4. + 19./12.)*gamma(pw/4. - 1./12.) $
          *(!unit.me*!unit.c*wnu/(3.d0*!unit.e*B*sina))^((1.d0-pw)/2.d0)  
    dy = y[2]-y[1]
    prj_P_tot = total(P_tot,2)*dy^3.d0
    sufb = prj_P_tot / 4.d0 /!pi / dy / dy
    img0 = congrid(sufb,768,768)
;    sd = size(img0)

    loadct,3,/sil
    exestr = execute('window,0,xs=768,ys=768'+opt)
;    plot,x,y,/xst,/yst,/nodata,xr=[x[0],x[sd[1]-1]],yr=[y[0],y[sd[2]-1]],/iso,position=[0,0,sd[1],sd[2]],/dev
    tv,bytscl(alog10(img0),max=alog10(maxv_sufb),min=alog10(minv_sufb))
;    overorbit

;    legend,strtrim(time,1)+ ' s',/right,/top,text=fsc_color('yellow'),box=0
     legend,string(time/60./60.,format='(f6.2)')+' hrs',box=0,textcolor=fsc_color('yellow') $
           ,position=[600,750],/dev

    if keyword_set(png) then $
       draw,'png_radio/Jet_plt_'+string(startn+j-start_j,format='(I4.4)')
endfor
stop
end


