pro conicJet

x=findgen(500)-250.
y=findgen(500)-250.

coordarray,x,y,xout=xx,yout=yy

rj=200.

th0 = 30.*!dtor

r = sqrt(xx*xx+yy*yy)

sz = size(r,/dimension)
xnorm=fltarr(sz)
ynorm=fltarr(sz)
for i=0,sz[0]-1 do begin
    for j=0,sz[1]-1 do begin
        if (r[i,j] ne 0) then begin
           xnorm[i,j] = xx[i,j]/r[i,j]
           ynorm[i,j] = yy[i,j]/r[i,j]
        endif else begin
           xnorm[i,j] = 0.
           ynorm[i,j] = 0.
        endelse
    endfor
endfor
; variation of jet angle depending on the radius
th = th0 * r/rj

zz2 = cos(th)
xx2 = sin(th)*xnorm
yy2 = sin(th)*ynorm

window,0,xs=sz[0],ys=sz[1]
tvcoord,zz2,x,y,/scale
plots,ring(0,0,rj),/data

r2 = sqrt(xx2*xx2+yy2*yy2)

r1 = sqrt(xnorm*xnorm + ynorm*ynorm)

stop
end
