files = findfile('$PRIMUS_REDUX/1d/1500/*/????????*zAll*')

for i = 0, n_elements(files) -1 do begin

   junk = strsplit(files[i], '/.-', /extract)
   rootname = junk[n_elements(junk)-4]
   
   primus_read_1dinputs, rootname, ext=ext, slits=slits, oned=oned, $
     photoinfo=photoinfo, 1500
   
   keep = where(oned.zgrid_gal[10] gt 0)
   ext = ext[keep]
   oned = oned[keep]
   slits = slits[keep]
   photoinfo = photoinfo[keep]

   maggies = photoinfo.maggies
   maggiesivar = photoinfo.maggiesivar
   redshift = oned.zmin_gal[0]
   delvarx, vmatrix
   delvarx, lambda
   delvarx, lambda
   delvarx, coeffs
   delvarx, rmatrix

   kcorrect, maggies, maggiesivar, redshift, kcorrect, $
     filterlist=strtrim(photoinfo[0].filterlist, 2), $
     band_shift=0.7, lambda=lambda, coeffs=coeffs, rmatrix=rmatrix, $
     vmatrix=vmatrix, zvals=zvals
   
   help, vmatrix, zvals, coeffs

   k_reconstruct_maggies, coeffs, replicate(0.7, n_elements(ext)), $
     maggies, zvals=zvals, vmatrix=vmatrix, lambda=lambda, $
     filterlist=['deep_B.par','deep_R.par','deep_I.par']
   distmod = distmod(redshift)

   output = {zprimus : redshift[0], $
     absmag : maggies[*,0], $
     distmod : distmod[0], $
     zwarning : 0, $
     snr : 0.0d, $
     currz : 0.0, $
     ztype :0}
   output = replicate(output, n_elements(ext))
   output.zprimus = redshift
   output.absmag = maggies
   output.distmod = distmod
   output.zwarning = oned.zwarning
   output.snr = sqrt(ext.sn1^2+ext.sn2^2)
   output.currz = slits.currz
   output.ztype = slits.ztype
   
   
   if n_elements(global) eq 0 then global = output else $
     global= [global, output]

endfor

mwrfits, global, 'deep2_1500.fits', /create


END

