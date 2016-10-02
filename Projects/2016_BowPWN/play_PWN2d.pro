pro play_PWN2d,fname,var=var,sample=sample,startn=startn,endn=endn,step=step,png=png,xrange=xrange,yrange=yrange,bar=bar,nocombine=nocombine,bubble=bubble

device, decomposed=0

if not keyword_set(fname) then fname='PWN2d_hdf5_plt_cnt_'
if not keyword_set(var) then var='dens'
if not keyword_set(sample) then sample=4
if not keyword_set(step) then step=1

fnames = file_search(fname+'*')
nfiles = n_elements(fnames)

indexf = strmid(fnames,3,4,/reverse)

if keyword_set(startn) then begin
   str_start = string(startn,format='(I4.4)')
   start_j_tmp = where(strpos(fnames,str_start) ne -1)
   start_j = start_j_tmp[0]
   if start_j eq -1 then begin
      print,'out of range of startnumber'
      stop
   endif
endif else start_j = 0
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
   spawn,'mkdir png_'+var
   opt=',/pixmap' 
endif else opt=''

;yc= 6.1714e19
;zc= 6.1714e19

;maxd = 7.e-24 & mind = 9.e-25
maxd = 8.4e-24 & mind = 2.5e-30

;if not keyword_set(xrange) then xrange=[1.e18,2.4e20]

;xrange = [0.,1.2e18] 
;yrange = [-5.e17,3.5e18]
;sample = 3

if not keyword_set(xrange) then xropt = '' else xropt=',xrange=xrange'
if not keyword_set(yrange) then yropt = '' else yropt=',yrange=yrange'

loadct,0,/sil
!p.background=255

exestr1 = execute('d = loaddata(fnames[start_j],var,time=time,sample=sample,xcoords=x,ycoords=y'+xropt+yropt+')')
sd = size(d,/dimension)
d2 = fltarr(sd[0]*2,sd[1])
x2 = fltarr(sd[0]*2)

if keyword_set(nocombine) then yw = strtrim(sd[0],1) else yw = strtrim(sd[0]*2,1)
if keyword_set(nocombine) then yw2 = sd[0] else yw2 = sd[0]*2
exestr2 = execute('window,0,xs='+strtrim(sd[1],1)+',ys='+yw+opt)

; theoretical approach of bubble
vps = 6.d7
Esp = 6.66d35
d0  = 1.67d-24
dtn1 = 0.                     ; detonation time 1
dtn2 = 9.35d10 +1.d18/vps      ; detonation time 2
CC1  = (124./154./!pi)^0.2 * (Esp/d0)^0.2
CC2  = (124./154./!pi)^0.2 * (Esp/d0)^0.2

for j=start_j,end_j,step do begin
    print,j+1,'   of',end_j+1,'  ', fnames[j]
    exestr1 = execute('d = loaddata(fnames[j],var,time=time,sample=sample,xcoords=x,ycoords=y'+xropt+yropt+')')

    if not keyword_set(nocombine) then begin
      d2[sd[0]:sd[0]*2-1,*] = d
      d2[0:sd[0]-1,*] = reverse(d,1)
      x2[sd[0]:sd[0]*2-1] = x
      x2[0:sd[0]-1] = -reverse(x)
     endif else begin
      d2 = d
      x2 = x
    endelse

    if not keyword_set(bubble) then $
       tv,bytscl(alog10(transpose(d2)),min=alog10(mind),max=alog10(maxd)) $
      else begin
       tvcoord,bytscl(alog10(transpose(d2)),min=alog10(mind),max=alog10(maxd)),y,x2
       if (time ge dtn1) then begin
          bubc1 = vps*(time-dtn1)
          rad1  = CC1*(time-dtn1)^0.6
          plots,ring(bubc1,0,rad1),/data,color=0,thick=1
       endif
       if (time ge dtn2) then begin
          bubc2 = vps*(time-dtn2)
          rad2  = CC2*(time-dtn2)^0.6
          plots,ring(bubc2,0,rad2),/data,color=0,thick=1
       endif
    endelse
;    if not keyword_set(nocombine) then $
;       tv,bytscl(alog10(reverse(transpose(d),2)),min=alog10(mind),max=alog10(maxd)),1

;    legend,strtrim(time/60./60./24./365.,1)+ ' yr',position=[sd[1]-12,yw2-20],/dev,textcolor=0,box=0
;   note that time is set to be zero at contering density changed region. (ex. 1.67e18 away & 600 km/s)
    legend,string((time - dtn1)/60./60./24./365.,format='(f9.2)')+ ' yr',position=[sd[1]-130,yw2-10],/dev,textcolor=0,box=0

;    if keyword_set(bar) then color_bar,lim=[mind,maxd],/log,/right,pos=[0,0,10,200],/dev,bartitle='density' $
    if keyword_set(bar) then color_bar,lim=[mind,maxd],/log,/left,pos=[sd[1]-10,0,sd[1],200],/dev,bartitle='density' $
                            , titlegap=70
    if keyword_set(png) then $
       draw,'png_'+var+'/PWN2d_plt_'+indexf[j]
endfor
stop
end
