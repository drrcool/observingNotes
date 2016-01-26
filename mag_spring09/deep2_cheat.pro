files = findfile('~/primus/masks/*deep*slits*')
delvarx, output    
    for i = 0, n_elements(files) -1 do begin
       
       obj = mrdfits(files[i], 1, columns=['ra','dec','z','magr','magi','magb'])
       
       if n_elements(output) gt 0 then output = [output, obj] else $
         output = obj
       
    endfor
    k = where(output.z gt 0 and output.z lt 0.9)
    output = output[k]
    
    mag = [[output.magb],[output.magr],[output.magi]]
    maggies = 10^(-0.4*mag) * (mag gt 0)
    err = maggies*0.1
    ivar = 1d/err^2
    k = where(mag eq 0)
    ivar[k] = 0
    z = output.z
    kcorrect, transpose(maggies), transpose(ivar), z, $
      kcorrect, band_shift=0.7, absmag=absmag, $
      filterlist=['deep_B.par','deep_R.par','deep_I.par']
    
    add = {absmag : fltarr(3)}
    add = replicate(add, n_elements(output))
    add.absmag = absmag
    
    output = struct_addtags(output, add)
    
    mwrfits, output, 'deep_kcorr.fits', /create
    
    spawn, 'gzip -vf deep_kcorr.fits'
    
    k = where(output.z gt 0.6 and output.z lt 0.8)
    output = output[k]

    djs_oplot, output.absmag[1,*], output.absmag[1,*]-output.absmag[2,*], $
      ps=3, xrange=[-17, -23], yrange=[0, 1.5]


    test = mrdfits('~/kobol/data/primus/kcorrect_test/merged.fits.gz', 1)
    k = where(test.zprimus gt 0.6 and test.zprimus lt 0.8 $
      and (test.zwarning and 2l^5) eq 0 and $
      (test.zwarning and 2l^6) eq 0 and test.snr gt 12)
    test = test[k]

    kcorrect, test.maggies[0:2,*], test.imaggies[0:2,*], test.zprimus, $
      kcorrect, absmag=absmag, filterlist=['deep_B.par','deep_R.par','deep_I.par'], band_shift=0.7

;    djs_oplot, absmag[1,*], absmag[1,*]-absmag[2,*], ps=3, color='red'
    add = {absmag : fltarr(3)}
    

    spherematch, test.ra, test.dec, output.ra, output.dec, 1.0/3600., m1, m2, dist
    output=output[m2]
    test = test[m1]




 END

    
