pro therm1d,n=n,startn=startn,endn=endn

fnames = file_search('*'+'hdf5_plt_cnt'+'*')
nfiles = n_elements(fnames)

if not keyword_set(startn) then begin
   if (n_elements(n) ne 0) then begin
      ind = where(strmid(fnames,3,/reverse) eq string(n,format='(I4.4)'))
      startn = ind[0] & endn2=ind[0]
   endif else begin
      startn = 0 & endn2=nfiles-1
   endelse
endif else endn=nfiles-1

if keyword_set(endn) then endn2=endn

rs=1.4e12

window,0,xs=800,ys=1000
;window,0,xs=400,ys=500
dy=0.29
flag=1
yra_den = ',yst=2'
yra_pres = ',yst=2'
yra_vx = ',yst=2'
for i=startn,endn2 do begin
    den  = loaddata(fnames[i],'dens',xcoord=x,time=time)
    vx   = loaddata(fnames[i],'velx')
    pres = loaddata(fnames[i],'pres') 

    a1=execute("plot,x/rs,den ,pos=[0.16,0.1+2.*dy,0.97,0.1+3.*dy],/norm,xtickformat='(a1)',xst=2" $
              + yra_den + ",ytitle='density',/ylog")
    if (flag) then begin
       yra = !y.crange
       yra_den = ',yra=10^['+strtrim(yra[0],1)+','+strtrim(yra[1],1)+'],/yst'
       den0 = den
       vx0 = vx
       pres0 = pres
    endif 
    oplot, x/rs, den0, line=2
;    legend,fnames[i]+' ('+strtrim(time,1)+' s)',/top,/right,box=0

    a2=execute("plot,x/rs,pres,pos=[0.16,0.1+dy,0.97,0.1+2.*dy],/norm,xtickformat='(a1)',/noerase,xst=2" $
       + yra_pres + ",ytitle='pressure',/ylog")
    if (flag) then begin
       yra = !y.crange
       yra_pres = ',yra=10.^['+strtrim(yra[0],1)+','+strtrim(yra[1],1)+'],/yst'
    endif 
    oplot, x/rs, pres0, line=2

    a3=execute("plot,x/rs,vx/1.e5,pos=[0.16,0.1,0.97,0.1+dy],/norm,xtitle='r/Rs',/noerase,xst=2" $
              +",ytitle='velx [km/s]',yst=2")
    if (flag) then begin
       yra = !y.crange
       yra_vx = ',yra=['+strtrim(yra[0],1)+','+strtrim(yra[1],1)+'],/yst'
       flag = 0
    endif 
    oplot, x/rs, vx0/1.e5, line=2
endfor


stop 
end
