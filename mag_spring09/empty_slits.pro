PRO empty_slits, indir=indir, psfile=psfile, specfile=specfile, $
  outroot=outroot, slitsfile=slitsfile


    if not keyword_set(indir) then indir = './'
    Porig = !P
    Xorig = !X
    Yorig = !Y


    for ifile = 0, n_elements(specfile) -1 do begin
       spec = mrdfits(indir + '/' + specfile[ifile], 1)
       slits = mrdfits(slitsfile,1)
       spherematch, slits.ra, slits.dec, spec.ra, spec.dec, 1.0/3600,$
         m1, m2, dist
       spec = spec[m2]
       slits = slits[m1]

       junk = strsplit(specfile[ifile], '.', /extract)
       root = junk[n_elements(junk)-4]
       skysub = root + '_skyimage.fits.gz'
       if file_test(indir+'/'+skysub) then begin
          for ccd = 1, 8 do begin
             image = mrdfits(indir+'/'+skysub, ccd-1)          
             keep = where(spec.ccdnum eq ccd and slits.sky eq 1 , ct)
             if ct gt 0 then begin
                subspec = spec[keep]
                nx = subspec[0].xcorner[1]-subspec[0].xcorner[0]
                ny = subspec[0].ycorner[1]-subspec[0].ycorner[0]
                cutout = fltarr(nx, ny, ct)
                for jj = 1, ct -1 do begin
                   cutout[*,*,jj] = extrac(image, subspec[jj].xcorner[0],$
                     subspec[jj].ycorner[0], nx, ny)     
                endfor
                
                mediancutout = median(cutout, dimen=3)
                mwrfits, cutout, root+'_fullcutout_' + $
                  string(ccd, format='(i1)')+'.fits', /create
                spawn, 'gzip -vf ' + root+'_fullcutout_' + $
                  string(ccd, format='(i1)') +'.fits'
                
                mwrfits, mediancutout,root+'_mediancut.fits', $
                  create=(ccd eq 1)
             endif
          endfor
          spawn, 'gzip -vf ' + root +'_mediancut.fits'

       endif

    endfor












    !X = Xorig
    !Y = Yorig
    !P = Porig


 END

       

    
