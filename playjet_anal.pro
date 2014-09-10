pro playjet_anal,fname,startn=startn,endn=endn,step=step,png=png
device, decomposed=0

if not keyword_set(fname) then fname='JetSet_hdf5_plt_cnt_'
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
   spawn,'mkdir png_playjet'
   opt=',/pixmap' 
endif else opt=''

xc0 = 0. & yc0 = 0. & zc0 = 0.
maxd = 2.76e-13 & mind = 2.9d-16

!p.background=255 & !p.color=0
flag=0
for j=start_j,end_j,step do begin
    print,j+1,'   of',end_j+1,'  ', fnames[j]
    d = dload(j,var='dens',xc=xc0,yc=yc0,zc=zc0,x,y,z,sample=0,time)
    sd = size(d)
    den = reform(d[*,*,sd[3]/2])
    
    v1 = dload(j,var='velx',xc=xc0,yc=yc0,zc=zc0,sample=0)   
    vx = reform(v1[*,sd[2]/2,sd[3]/2])

    if (flag eq 0) then begin
       den0=reform(den[*,sd[2]/2])
       vx0=vx
       flag=1
    endif

    loadct,0,/sil
    exestr = execute('window,0,xs='+strtrim(sd[1]+600,1)+',ys='+strtrim(sd[2],1)+opt)
    tv,bytscl(alog10(den),max=alog10(maxd),min=alog10(mind)),0,0
 
    plot,x/1.4e12,den[*,sd[2]/2],/xst,position=[sd[1]+120,300,sd[1]+580,500],/dev,/noerase,xr=[0.,7],ytitle='density',xtickformat='(a1)',/ylog
    legend,strtrim(time,1)+ ' s',/right,/top,box=0
    oplot,x/1.4e12,den0,line=2

    plot,x/1.4e12,vx/1.e5,/xst,position=[sd[1]+120,100,sd[1]+580,300],/dev,/noerase,xr=[0.,7],ytitle='velocity [km/s]',yr=[0.,3500],xtitle='r/R_s'
    oplot,x/1.4e12,vx0/1.e5,line=2
    
    if keyword_set(png) then $
       draw,'png_playjet/Jet_plt_'+string(j,format='(I4.4)')
endfor
stop
end


