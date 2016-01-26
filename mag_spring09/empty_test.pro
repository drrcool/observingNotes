Pro empty_test, specfile=specfile, psfile=psfile, rerun=rerun
    
    psfile ='empty_slits_ut061223_' + string(rerun, format='(i4.4)') + '.ps'
    specfile = findfile('$PRIMUS_REDUX/' + string(rerun, format='(i4.4)') + $
      '/ut061223/*extract*fits.gz')
    Porig = !P
    Xorig = !X
    Yorig = !Y
    specfile = reverse(specfile)
    dfpsplot, psfile, /color, /landscape                            
    !P.thick=8
    !X.thick = 8
    !Y.thick = 8
    !P.charthick=8
    !P.charsize=1.5
    
    nfile = n_elements(specfile)
    for ifile = 0, nfile -1 do begin
       
       ;;Read in the extraction file
       ext = mrdfits(specfile[ifile], 1)
       smffile = ext[0].infile
       junk = strsplit(smffile, '.', /extract)
       maskroot = junk[0]
       
       
       ;;This will issue a warning, but I really don't care
       primus_read_1dinputs, maskroot, rerun, slits=slits, /noext, /nooned
       
       if tag_exist(slits,'sky') then begin 
          ;;Quick match and trim of the files to get them to linematched
          spherematch, ext.ra, ext.dec, slits.ra, slits.dec, 1.0/3600., m1, m2, dist
          ext = ext[m1]
          slits = slits[m2]
          
          
          keep = where(slits.sky and finite(total(ext.fopt1+ext.fopt2, 1)), ct)
          Waveout = ext[100].wave1
          icronout = waveout^(-3.0)
          
          icron1 = ext.wave1^(-3.0)
          Icron2 = ext.wave2^(-3.0)
          flux1 = icron1*0.0
          flux2 = icron2*0.0
          sky1 = icron1*0.0
          sky2 = icron2*0.0
          for i = 0, ct - 1 do begin
             
             flux1[*,keep(i)] = interpol(ext[keep(i)].fopt1*ext[keep(i)].calib1, icron1[*,keep(i)], $
               icronout)
             flux2[*,keep(i)] = interpol(ext[keep(i)].fopt2*ext[keep(i)].calib2, icron2[*,keep(i)], $
               icronout)
             sky1[*,keep(i)] = interpol(ext[keep(i)].sky1*ext[keep(i)].calib1, icron1[*,keep(i)], $
               icronout)
             sky2[*,keep(i)] = interpol(ext[keep(i)].sky2*ext[keep(i)].calib2, icron2[*,keep(i)], $
               icronout)
             
          endfor
          
          k = where(finite(flux1+sky1) eq 0 or sky1 eq 0.0, ct)
          if ct gt 0 then begin
             flux1[k] = 0.0
             sky1[k] = 1.0
          endif
          
          !P.multi=[0,2,4]
          for iccd = 1, 8 do begin
             keep = where(slits.sky and finite(total(ext.fopt1+ext.fopt2, 1)) and $
               ext.ccdnum eq iccd, ct)
             
             case iccd of
                1 : !P.position = [0.1, 0.5, 0.3, 0.9]
                2 : !P.position = [0.3, 0.5, 0.5, 0.9]
                3 : !P.position = [0.5, 0.5, 0.7, 0.9]
                4 : !P.position = [0.7, 0.5, 0.9, 0.9]
                6 : !P.position = [0.1, 0.1, 0.3, 0.5]
                5 : !P.position = [0.3, 0.1, 0.5, 0.5]
                8 : !P.position = [0.5, 0.1, 0.7, 0.5]
                7 : !P.position = [0.7, 0.1, 0.9, 0.5]
             endcase
             
             ytickformat='(A1)'
             if iccd eq 6 or iccd eq 1 then ytickformat=''
             xtickformat='(A1)'
             if iccd gt 4 then xtickformat=''
             
             xtitle=''
             ytitle='' 
             if iccd gt 4 then xtitle='Wavelength (\AA)'
             if iccd eq 1 or iccd eq 6 then ytitle='Sky Subtracted Empty Slits'
             title=''
             junk = strsplit(specfile[ifile], './', /extract)
             root = junk[n_elements(junk)-4]
             if iccd eq 2 then title=root
             if iccd eq 3 then title="CCD " +  string(iccd, format='(i1)')

             
             if ct gt 10 then begin          
                
                yval1 = waveout*0.0
                err1 = waveout*0.0
                yval2 = yval1
                err2 = err1
                
                for ii = 0, n_elements(waveout) -1 do begin
                   djs_iterstat, flux1[ii,keep], median=median, sigma=sigma
                   yval1[ii] = median
                   err1[ii] = sigma
                   djs_iterstat, flux2[ii,keep], median=median, sigma=sigma
                   yval2[ii] = median
                   err2[ii] = sigma
                endfor
                
                
                
                Rjc_ploterr, waveout, yval1, $
                  xrange=[4500, 9500], ps=2, yrange=[-0.02, 0.02],  /xs, /ys, $
                  xtickformat=xtickformat, ytickformat=ytickformat, xtitle=xtitle, ytitle=ytitle, $
                  title=title
                
                range = .04
                djs_xyouts, 5000, -0.02+0.9*range, 'CCD ' + string(iccd, format='(i1)')
                
                xxval = findgen(500)*10+4500
                yyval = interpol(yval1, waveout, xxval)
                eeval = interpoL(err1, waveout, xxval)
                for ii = 0, n_elements(xxval) -1 do $
                  djs_oplot, xxval[ii]*[1,1], yyval[ii]+eeval[ii]*[-1,1], thick=1, color='grey'
                djs_oplot, waveout, yval, ps=2
                
                
                k = where(finite(flux2+sky2) eq 0 or sky2 eq 0.0, ct)
                if ct gt 0 then begin
                   flux2[k] = 0.0
                   sky2[k] = 1.0
                endif
                
                
                
                xxval = findgen(500)*10+4500
                yyval = interpol(yval2, waveout, xxval)
                eeval = interpoL(err2, waveout, xxval)
                for ii = 0, n_elements(xxval) -1 do $
                  djs_oplot, xxval[ii]*[1,1], yyval[ii]+eeval[ii]*[-1,1], color='pink', thick=1
                rjc_oploterr, waveout, yval2, color='red', ps=4
                djs_oplot, waveout, yval1, ps=2
                djs_oplot, [-100, 100000], [0,0], line=3
                
                djs_oplot, waveout, median(sky1[*,keep], dimen=2)*0.01, ps=10, color='dark blue'
                
                
             endif
             
             
          endfor
       endif
       
    endfor
    dfpsclose
    !P = Porig
    !X = Xorig
    !Y = Yorig
    
 END
