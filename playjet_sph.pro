pro playjet_sph,fname,var=var,startn=startn,endn=endn,step=step,png=png, zoomind=zoomind
device, decomposed=0

if not keyword_set(fname) then fname='JetSet_hdf5_plt_cnt_'
if not keyword_set(var) then var='dens'
if not keyword_set(step) then step=1
if not keyword_set(zoomind) then zoomind=-1

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
endif else begin
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
chk=''
if keyword_set(png) then begin
   dd = file_search('png_'+var)
   if (dd eq '') then $
      spawn,'mkdir png_'+var $
    else begin
Jump1:
      read,'directory is already existed. overwriting? (y/n): ',chk
      case chk of
        'y' : break 
        'n' : stop
       else : goto, Jump1
      endcase
   endelse
   opt=',/pixmap' 
endif else opt=''

xc0 = 0. & yc0 = 0. & zc0 = 0.
;maxd = 2.2d-13 & mind = 1.1d-15
;maxd = 2.76d-13 & mind = 2.9d-16
;maxd = 2.76d-13 & mind = 1.d-16
maxd = 8.e-14 & mind = 3.d-18

m1 = 4.d34
m2 = 2.d34
mtot = m1+m2
peri = 3.d12

xc1 = -peri*m1/mtot

jeton = 9.76e4

!p.background=255
nfrac=15
for j=start_j,end_j,step do begin
    print,j+1,'   of',end_j+1,'  ', fnames[j]
    if not keyword_set(sample) then begin
       read_amr,fnames[j],var='dens',tree=tree
       smp = max(tree.lrefine) - 7
    endif else smp = sample
    d = dload(startn+j-start_j,var=var,xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=smp,time) 
    img0 = congrid(reform(d[*,256,*]),768,768)
    xx   = congrid(x,768)
    zz   = congrid(z,768)
    sd = size(img0)
    loadct,0,/sil
    exestr = execute('window,0,xs=1024,ys=768'+opt)
;    plot,x,y,/xst,/yst,/nodata,xr=[x[0],x[sd[1]-1]],yr=[y[0],y[sd[2]-1]],/iso,position=[0,0,sd[1],sd[2]],/dev
    tv,bytscl(alog10(img0),max=alog10(maxd),min=alog10(mind))
;    overorbit

    if (j eq zoomind) then begin
       d1 = dload(j,var=var,xc=xc1,yc=yc0,zc=zc0,x1,y1,z1,sample=0)

       indxl_tmp = where(xx ge x1[0]) & indxl = indxl_tmp[0]       
       indxr_tmp = where(xx ge x1[511]) & indxr = indxr_tmp[0]       
       indzl_tmp = where(zz ge z1[0]) & indzl = indzl_tmp[0]       
       indzr_tmp = where(zz ge z1[511]) & indzr = indzr_tmp[0]       

       dx0 = indxr-indxl
       szfrac = findgen(nfrac)/(float(nfrac)-1.) * (512.-dx0) + dx0
       xxfrac = findgen(nfrac)/(float(nfrac)-1.) * (512.-indxl) + indxl
       zzfrac = findgen(nfrac)/(float(nfrac)-1.) * (128.-indzl) + indzl
        
       for i=0,nfrac-1 do begin
           loadct,0,/sil
           tv,bytscl(alog10(img0),max=alog10(maxd),min=alog10(mind))
           img1 = congrid(reform(d1[*,256,*]),szfrac[i],szfrac[i])
           tv,bytscl(alog10(img1),max=alog10(maxd),min=alog10(mind)),xxfrac[i],zzfrac[i]
           plots,xxfrac[i],zzfrac[i],/dev,color=fsc_color('yellow'),thick=2
           plots,xxfrac[i]+szfrac[i],zzfrac[i],/dev,color=fsc_color('yellow'),thick=2,/continue
           plots,xxfrac[i]+szfrac[i],zzfrac[i]+szfrac[i],/dev,color=fsc_color('yellow'),thick=2,/continue
           plots,xxfrac[i],zzfrac[i]+szfrac[i],/dev,color=fsc_color('yellow'),thick=2,/continue
           plots,xxfrac[i],zzfrac[i],/dev,color=fsc_color('yellow'),thick=2,/continue
           legend,string(time/60./60./24.,format='(f4.2)')+' days',box=0,textcolor=fsc_color('yellow') $
                 ,position=[600,750],/dev
           if keyword_set(png) then $
              draw,'png_'+var+'/Jet_plt_'+string(j,format='(I4.4)')+'_'+string(i,format='(I2.2)')
       endfor
    endif else if (time ge jeton) then begin
       d1 = dload(j,var=var,xc=xc1,yc=yc0,zc=zc0,sample=0) 
       img1 = reform(d1[*,256,*])
       tv,bytscl(alog10(img1),max=alog10(maxd),min=alog10(mind)),512,128
       plots,512,128,/dev,color=fsc_color('yellow'),thick=2
       plots,1023,128,/dev,color=fsc_color('yellow'),thick=2,/continue
       plots,1023,640,/dev,color=fsc_color('yellow'),thick=2,/continue
       plots,512,640,/dev,color=fsc_color('yellow'),thick=2,/continue
       plots,512,128,/dev,color=fsc_color('yellow'),thick=2,/continue
    endif

;    legend,strtrim(time,1)+ ' s',/right,/top,text=fsc_color('yellow'),box=0
   
    if (j ne zoomind) then begin
       legend,string(time/60./60./24.,format='(f4.2)')+' days',box=0,textcolor=fsc_color('yellow') $
             ,position=[600,750],/dev
       if keyword_set(png) then $
          draw,'png_'+var+'/Jet_plt_'+string(j,format='(I4.4)')
    endif
endfor
stop
end


