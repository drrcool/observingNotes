files = findfile('$PRIMUS_REDUX/1d/9998/*/*zAll*')
    for ifile = 0, n_elements(files) -1 do begin
       
       junk = strsplit(files[ifile], './-', /extract)
       rootname = junk[n_elements(junk)-4]
       primus_read_1dinputs, rootname, 9998, /nooned, /noext, slits=slits
       keep = where(slits.currz gt 0 and slits.currz lt 1.0, ct)
       if ct gt 0 then $
         primus_zvz_qaplot, mask=rootname, /coadd, rerun=9998, output=rootname+'.ps', refband = 1
       
    endfor
    
 END


       
