pro pp_jettrace,ps=ps

fname1 = '/d/d3/yoon/outputs/out_mHD_Binary_sphere/sphere_1e35/jettrace_1e35_700.sav'
fname2 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e36_high/jettrace_1e36_high_891.sav'
fname3 = '/d/d7/yoon/out_FLASH3.3_mhd/out_mHD_Binary_sphere/sphere_1e37/jettrace_1e37_490.sav'
restore,filename=fname1 & x2_35 = x2 & z2_35 = z2
restore,filename=fname2 & x2_36 = x2 & z2_36 = z2
restore,filename=fname3 & x2_37 = x2 & z2_37 = z2

loadct,39,/sil
if not keyword_set(ps) then begin
!p.background=255 & !p.color=0
window,1,xs=800,ys=700
endif else mkeps,'pp_jettrace',xs=20.,ys=17.
plot,x2_35,z2_35,xr=[-1.e13,2.e12],/xst,yr=[0.,1.e13],yst=2,psym=1,/iso,xtitle='x [cm]',ytitle='z [cm]',charsize=1.5
oplot,x2_36,z2_36,psym=1,color=50
oplot,x2_37,z2_37,psym=1,color=250

xc0=1.e12 & yc0=0.
xc1=!x.crange[0]
plotsym,3
plots,xc0,yc0,/data,psym=8,symsize=3

; fitted results
ind35 = where(x2_35 le -2.5e12)
lin35 = linfit(x2_35[ind35],z2_35[ind35])
x0=x2_35[ind35] & x0=reform(x0[0])
y0 = lin35[1]*x0+lin35[0]
oplot,[x0,!x.crange[0]],[lin35[1]*x0+lin35[0],lin35[0]+lin35[1]*!x.crange[0]]
oplot,[x0,xc1],[y0, -tan(16.69*!dtor)*(xc1-x0)+y0],line=2

ind36 = where((z2_36 ge 2.e12) and (z2_36 le 8.e12))
xx2 = x2_36[ind36] & xx2h = xx2[where(finite(xx2) eq 1)]
zz2 = z2_36[ind36] & zz2h = zz2[where(finite(xx2) eq 1)]
lin36 = linfit(xx2h,zz2h)
x0=reform(xx2h[0])
y0=lin36[0]+lin36[1]*x0
oplot,[x0,!x.crange[0]],[lin36[0]+lin36[1]*x0,lin36[0]+lin36[1]*!x.crange[0]],color=50
oplot,[x0,xc1],[y0,-tan(66.*!dtor)*(xc1-x0)+y0],color=50,line=2

ind37 = where(z2_37 ge 1.e12)
lin37 = linfit(x2_37[ind37],z2_37[ind37])
x0=x2_37[ind37] & x0=reform(x0[0])
y0=lin37[0]+lin37[1]*x0
oplot,[x0,!x.crange[0]],[lin37[0]+lin37[1]*x0,lin37[0]+lin37[1]*!x.crange[0]],color=250
oplot,[x0,xc1],[y0, -tan(78.68*!dtor)*(xc1-x0)+y0],color=250,line=2

;oplot,[xc0,xc1],[yc0, -tan(16.69*!dtor)*(xc1-xc0)+yc0],line=2
;oplot,[xc0,xc1],[yc0, -tan(66.*!dtor)*(xc1-xc0)+yc0],color=50,line=2
;oplot,[xc0,xc1],[yc0, -tan(78.68*!dtor)*(xc1-xc0)+yc0],color=250,line=2

legend,textoidl('P_{jet}=')+['1e37','1e36','1e35'],psym=1,color=[250,50,0],textcolor=[250,50,0],/right,/top,charsize=1.5

if keyword_set(ps) then epsfree
stop
end
