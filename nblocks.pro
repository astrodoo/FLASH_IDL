pro nblocks, fname

openr,lun,fname,/get_lun
line = ''

n=0L & t=0.
nlblk = 0L & tlblk = 0. & lblk = 0L
nblk  = 0L & tblk  = 0. & blk  = 0L
while ~EOF(lun) do begin
   readf, lun, line
   if (strmatch(line,'*|*')) then begin
      vars = strsplit(line,/extract,escape='|')
      if ((n_elements(vars) ne 5) or ~(strnumber(vars[0]))) then goto, BADLINE
      n = long(vars[0]) & t = float(vars[1]) 
   endif
   if (strmatch(line,'*leaf blocks*')) then begin
      vars2 = strsplit(line,'=',/extract)
      nlblk = [nlblk,n] & tlblk = [tlblk,t] 
      lblk = [lblk, long(vars2[1])]
   endif else if (strmatch(line,'*total blocks*')) then begin
      vars3 = strsplit(line,'=',/extract)
      nblk = [nblk,n] & tblk = [tblk,t] 
      blk = [blk, long(vars3[1])]
   endif
BADLINE:
endwhile

lblk  = reform(lblk[where(nlblk ne 0)])
tlblk = reform(tlblk[where(nlblk ne 0)])
nlblk = reform(nlblk[where(nlblk ne 0)])
blk  = reform(blk[where(nblk ne 0)])
tblk = reform(tblk[where(nblk ne 0)])
nblk = reform(nblk[where(nblk ne 0)])

save,filename=fname+'_blocks.dat',lblk,tlblk,nlblk,blk,tblk,nblk
end

pro draw_all

fnames = file_search('*_blocks.dat')

nfiles = n_elements(fnames)

min_n = 999999L & max_n = 0L
min_t = 999999. & max_t = 0.
min_blk = 999999L & max_blk = 0L
for i=0,nfiles-1 do begin
   restore,filename=fnames[i] 
   strexe = execute('lblk_'+strtrim(i,2)+ '= lblk')
   strexe = execute('tlblk_'+strtrim(i,2)+ '= tlblk')
   strexe = execute('nlblk_'+strtrim(i,2)+ '= nlblk')
   strexe = execute('blk_'+strtrim(i,2)+ '= blk')
   strexe = execute('tblk_'+strtrim(i,2)+ '= tblk')
   strexe = execute('nblk_'+strtrim(i,2)+ '= nblk')

   strexe = execute('min_n = min([nblk_'+strtrim(i,2)+',min_n])')
   strexe = execute('max_n = max([nblk_'+strtrim(i,2)+',max_n])')
   strexe = execute('min_t = min([tblk_'+strtrim(i,2)+',min_t])')
   strexe = execute('max_t = max([tblk_'+strtrim(i,2)+',max_t])')
   strexe = execute('min_blk = min([blk_'+strtrim(i,2)+',min_blk])')
   strexe = execute('max_blk = max([blk_'+strtrim(i,2)+',max_blk])')
endfor

window,0
plot,nblk,blk, xra=[min_n,max_n], yra=[min_blk,max_blk],/nodata,xtitle='n',ytitle='total blocks'
for i=0, nfiles-1 do begin
   strexe = execute('oplot,nblk_'+strtrim(i,2)+',blk_'+strtrim(i,2))
   strexe = execute('oplot,[nblk_'+strtrim(i,2)+'[0],nblk_'+strtrim(i,2)+'[0]],!y.crange, line=1')
endfor
window,1
plot,tblk,blk, xra=[min_t,max_t], yra=[min_blk,max_blk],/nodata,xtitle='time',ytitle='total blocks'
for i=0, nfiles-1 do begin
   strexe = execute('oplot,tblk_'+strtrim(i,2)+',blk_'+strtrim(i,2))
   strexe = execute('oplot,[tblk_'+strtrim(i,2)+'[0],tblk_'+strtrim(i,2)+'[0]],!y.crange, line=1')
endfor

stop
end
