pro bubble2d, startn=startn, endn=endn

files = file_search('PWN2d_hdf5_plt_cnt_*')
nn    = fix(strmid(files,3,4,/reverse))

if not keyword_set(startn) then startn=nn[0]
if not keyword_set(endn) then endn=nn[n_elements(nn)-1]

st_i = where(nn eq startn)
ed_i = where(nn eq endn)

if ((st_i eq -1) or (ed_i eq -1)) then 
   print,'incorrect startn or endn'
   stop
endif

ii = 0

nfile = ed_i-st_i + 1
t = fltarr(nfile)
r1 = replicate(!values.f_nan,nfile)
r2= replicate(!values.f_nan,nfile)

vps = 6.e7
tdch = 6.03e10 + 5.e17/vps

d0 = 1.67e-25

for i=st_i,ed_i do begin
    d = loaddata(files[i],'dens',xCoord=x,yCoord=y,time=time,sample=4)
   
    t[ii] = time

    bubc1 =  vps*time
    bubc1Ind_tmp = where(y ge bubc1) & bubc1Ind = bubc1Ind_tmp[0]

    r1_tmp = x[where(d[*,bubc1Ind] eq d0)] & r1[ii] = r1_tmp[0]

    if (time ge tdch) then begin
       t2 = time - tdch
       bubc2 = vps*t2
       bubc2Ind_tmp = where(y ge bubc2) & bubc2Ind = bubc2Ind_tmp[0]-1
       r2_tmp = x[where(d[*,bubc2Ind] le 4.e-25)] & r2[ii] = r2_tmp[0]



    endif

endfor


stop
end
