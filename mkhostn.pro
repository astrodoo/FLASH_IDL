pro mkhostn, startn,endn, without=without, discrete=discrete

if (n_elements(discrete) eq 0) then begin

if (n_elements(startn) eq 0) then startn=49
if (n_elements(endn) eq 0) then endn=72

nnode = endn-startn+1 
node = indgen(nnode) + startn

nwo = n_elements(without)
if (nwo ne 0) then $
   for i=0, nwo-1 do node = node[where(node ne without[i])]

endif else node = discrete    ; discrete

nnode = n_elements(node)

nodename = 'node'+string(node,format='(I2.2)')

print,'number of nodes: ', nnode
print,'number of processors: ', nnode*8
print,'output file: hostfile'
print,'nodes: ',node[*]

openw,lun1,'hostfile',/get_lun
for i=0, nnode*8-1 do printf,lun1,nodename[i mod nnode]
free_lun,lun1

stop
end
