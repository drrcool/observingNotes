cet = '/u/rcool/kobol/data/ndwfs/cetus/cetus_j_V1.0.cat'
    spawn, "awk '{print $7,$8}' " + cet + " >/tmp/junk"
    readcol, '/tmp/junk', ra, dec
    
    files = findfile('/u/rcool/primus/masks/*slits*')
    openw, 1, '/u/rcool/match_cetus.dat'
    for ifiles = 0, n_elements(files) -1 do begin
       
       junk = mrdfits(files[ifiles], 1)
       
       spherematch, ra, dec, junk.ra, junk.dec, 1.0/3600., m1, m2, dist
       
       if n_elements(m1) gt 1 then printf, 1, files[ifiles]
       
    endfor
    
    close, 1
    
    
    
 END

