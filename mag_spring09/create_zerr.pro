PRO create_zerr
    
    test = mrdfits('test.fits', 1)
    k = where(test.currz gt 0 and test.currz lt 1.0 and $
      test.ztype ne 3 and test.ztype ne 8 and test.zwarning eq 0)
    
    test = test[k]
    
    zmin = 0.2
    zmax = 1.0
    dz = 0.02
    
    zarray = 0.2 + findgen((zmax-zmin)/dz+1)*dz
   
    
    for i = 0, n_elements(zarray) -2 do begin
       
       keep = where(test.currz gt zarray[i] and $
         test.currz lt zarray[i+1], ct)
       
       err =  rc_hist(test[keep].zprimus-test[keep].currz, -1, 1, 0.001)
       
       if n_elements(output) eq 0 then begin
          output = {zmin : 0.0, $
            zmax : 0.0, $
            err : err.number*0.0, $
            zarray : err.midpoint}
          output = replicate(output, n_elements(zarray)-1)
       endif
       
       output[i].zmin = zarray[i]
       output[i].zmax = zarray[i+1]
       output[i].err = total(1.0*err.number, /cum)
       output[i].err = output[i].err / max(output[i].err)
       output[i].zarray = err.midpoint
       
    endfor
    
    mwrfits, output, 'zerr.fits', /create
    
    
 END

       
       
          
    
