pro scaling   ; Benchmark

dir='/d/d7/yoon/out_FLASH3.3_mhd/out_Jet_SphWind/sphere_1e36_rot/'
logfile='JetSet.log'

line = ''
openr,lun1,dir+logfile,/get_lun

first_read=1

cputime = 0.d0
step = 0L
phytime = 0.d0
phydt =0.d0
nblk = 0L

IO_cputime = 0.d0
IO = ''

while ~EOF(lun1) do begin
    readf,lun1,line

    line_nproc = strpos(line,'Number of processors:')
    line_step = strpos(line,'step:')
    line_blks = strpos(line,'tot blks')
    line_IO = strpos(line,'open: type=')

    if (line_nproc ne -1) then nproc = fix((strsplit(line,/extract))[3])

    if (line_blks ne -1) then begin
       strline = strsplit(line,/extract)
       nblk_ = long(strline[10])
    endif
     
    if ((line_step ne -1) or (line_IO ne -1)) then begin

       strline = strsplit(line,/extract)
       date = fix(strsplit(strline[1],'-',/extract))
       jday = julday(date[0],date[1],date[2])

       if (first_read) then begin
;          fdate = date
          jday0 = jday
          first_read = 0
       endif
       day = jday - jday0
       tt = float(strsplit(strline[2],':',/extract))
       timeInday = tt[2] + tt[1]*60. + tt[0]*60.*60.
       cputime_ = timeInDay + 60.d0*60.d0*24.d0*day
       
       if (line_step ne -1) then begin
          step_ = long((strsplit(strline[5],'=',/extract))[1])
          phytime_ = double((strsplit(strline[6],'=',/extract))[1])
          phydt_ = double((strsplit(strline[7],'=',/extract))[1])
 
          cputime = [cputime,cputime_]
          step = [step,step_]
          phytime = [phytime,phytime_]
          phydt = [phydt,phydt_]
          nblk  = [nblk,nblk_]
       endif else if (line_IO ne -1) then begin
          IO_ = (strsplit(strline[7],'=',/extract))[1]
          IO_cputime = [IO_cputime,cputime_]
          IO = [IO,IO_]
       endif
    endif
endwhile
free_lun,lun1

cputime = cputime[1:*]
step = step[1:*]
phytime = phytime[1:*]
phydt = phydt[1:*]
nblk = nblk[1:*]
IO = IO[1:*]
IO_cputime = IO_cputime[1:*]

save,file='scaling.sav',nproc,cputime,step,phytime,phydt,nblk,IO,IO_cputime
stop
end
