pro showblock, tree=tree, param= param, xc=xc, yc=yc, zc=zc, color=color, cell=cell, c16=c16 

; showing a block grid

if not keyword_set(color) then color=0
if keyword_set(c16) then nc=16 else nc=8

bline = 0 & bthick = 1   & bcolor = 0     ; block line setting
cline = 0 & cthick = 0.5 & ccolor = 0     ; cell line setting

for i=0L,param.totblocks-1 do begin
;Y-Z plane at x=xc
    if (n_elements(xc) ne 0) then begin
       if ((tree[i].bndbox[0,0] le xc) and (tree[i].bndbox[1,0] ge xc)) then begin                                    
          y1 = tree[i].bndbox[0,1] & y2 = tree[i].bndbox[1,1]                                                         
          z1 = tree[i].bndbox[0,2] & z2 = tree[i].bndbox[1,2]                                                         
          oplot,[y1,y2,y2,y1,y1],[z1,z1,z2,z2,z1],color=bcolor,thick=bthick,line=bline                     
          if keyword_set(cell) then begin
             ycl = findgen(nc)/float(nc)*(y2-y1)+y1 
             zcl = findgen(nc)/float(nc)*(z2-z1)+z1 
             for j=1,nc-1 do begin
                 oplot, [ycl[j],ycl[j]], [z1,z2],color=ccolor,thick=cthick,line=cline
                 oplot, [y1,y2], [zcl[j],zcl[j]],color=ccolor,thick=cthick,line=cline
             endfor
          endif
       endif                    
;X-Z plane at y=yc
    endif else if (n_elements(yc) ne 0) then begin
       if ((tree[i].bndbox[0,1] le yc) and (tree[i].bndbox[1,1] ge yc)) then begin
          x1 = tree[i].bndbox[0,0] & x2 = tree[i].bndbox[1,0]
          z1 = tree[i].bndbox[0,2] & z2 = tree[i].bndbox[1,2]
          oplot,[x1,x2,x2,x1,x1],[z1,z1,z2,z2,z1],color=bcolor,thick=bthick,line=bline
          if keyword_set(cell) then begin
             xcl = findgen(nc)/float(nc)*(x2-x1)+x1 
             zcl = findgen(nc)/float(nc)*(z2-z1)+z1 
             for j=1,nc-1 do begin
                 oplot, [xcl[j],xcl[j]], [z1,z2],color=ccolor,thick=cthick,line=cline
                 oplot, [x1,x2], [zcl[j],zcl[j]],color=ccolor,thick=cthick,line=cline
             endfor
          endif
       endif
;X-Y plane at z=zc
    endif else if (n_elements(zc) ne 0) then begin
       if ((tree[i].bndbox[0,2] le zc) and (tree[i].bndbox[1,2] ge zc)) then begin
          x1 = tree[i].bndbox[0,0] & x2 = tree[i].bndbox[1,0]
          y1 = tree[i].bndbox[0,1] & y2 = tree[i].bndbox[1,1]
          oplot,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=bcolor,thick=bthick,line=bline
          if keyword_set(cell) then begin
             xcl = findgen(nc)/float(nc)*(x2-x1)+x1 
             ycl = findgen(nc)/float(nc)*(y2-y1)+y1 
             for j=1,nc-1 do begin
                 oplot, [xcl[j],xcl[j]], [y1,y2],color=ccolor,thick=cthick,line=cline
                 oplot, [x1,x2], [ycl[j],ycl[j]],color=ccolor,thick=cthick,line=cline
             endfor
          endif
       endif
    endif
endfor
end
