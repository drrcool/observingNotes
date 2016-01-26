
test = mrdfits('/tmp/newdeep2.fits', 1)
    
    
    keep = where(test.newz gt 0.6 and test.newz lt 0.8 and test.magi lt 22)
    mag = test.origabsmag
      
    hogg_scatterplot, mag[1,*], mag[1,*]-mag[2,*], ps=3, xrange=[-17, -21.5], $
      yrange=[0, 1.5], /xs, /ys, /nogrey, grid=grid,  cumimage=cumimage, $
      xvec=xvec, yvec=yvec
    
    test = test[keep]
    mag = test.absmag
    
     hogg_scatterplot, mag[1,*], mag[1,*]-mag[2,*], ps=3, xrange=[-17, -21.5], $
      yrange=[0, 1.5], /xs, /ys, /nogrey, grid=grid, xnpix=30, ynpix=30
     
     hist1 = rc_hist(mag[1,*]-mag[2,*], 0, 1.5, 0.05)
     
     test = mrdfits('test.fits', 1)
     mag = -2.5*alog10(test.absmag) * (test.absmag gt 0)
     keep = where(mag[2,*] lt 22 and test.zprimus gt 0.6 and test.zprimus lt 0.8$
     and (test.zwarning and 2l^5) eq 0 and (test.zwarning and 2l^6) eq 0 and test.snr gt 12)
     test = test[keep]
     help, test
     mag = -2.5*alog10(test.absmag) * (test.absmag gt 0)
     
     hist2 = rc_hist(mag[1,*]-mag[2,*]+0.1, 0, 1.5, 0.05)
         
     hogg_scatterplot, mag[1,*]-test.distmod+0.2, $
      mag[1,*]-mag[2,*]+0.1, ps=3, xrange=[-17, -21.5], $
       yrange=[0, 1.5], /xs, /ys, /nogrey, grid=grid, xnpix=40, ynpix=40
     
    
     dfpsplot, 'test.ps', /color, /square
     loadct, 0
    display, alog10(grid+(grid eq 0))*(grid gt 0), $
      xrange=[-17, -21.5], yrange=[0, 1.5], top=250, /reverse, $
      xthick=8, ythick=8, charthick=8, xtitle=textoidl('M_{0.7 R}'), $
      ytitle=textoidl('^{0.7}(R-I)')
    
     contour, filter_image(cumimage, fwhm=1.5), xvec, yvec, $
      levels=findgen(10)*0.1, xrange=[-17, -21.5], $
      yrange=[0, 1.5], /xs, /ys, thick=8, xthick=8, ythick=8, $
      charthick=8, /over
     
     
    dfpsclose
    
    
    
    
 END
