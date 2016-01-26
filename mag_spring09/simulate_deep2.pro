

dat = mrdfits('~/kobol/data/deep2/zcat*fits.gz',1, $
  columns=['magr','magi','magb','z'])

zerr = mrdfits('zerr.fits', 1)

sort = sort(dat.z)
dat = dat[sort]

k = where(dat.z gt 0.3 and dat.z lt 0.9)
dat = dat[k]

newz = dat.z*0.0

for i = 0, n_elements(dat) -1 do begin
   k = where(zerr.zmin lt dat[i].z and $
     zerr.zmax gt dat[i].z)
   
   newz[i] = dat[i].z + interpol(zerr[k].zarray, $
     zerr[k].err, randomu(seed))
   
endfor

add = replicate( {newz:0.0}, n_elements(dat))
add.newz = newz

dat = struct_addtags(dat, add)

mags = [[dat.magb],[dat.magr],[dat.magi]]
maggies = 10.0^(-0.4*mags)
maggiesivar = 1d/(maggies*0.1)^2 *(mags gt 0)

filterlist = ['deep_B.par', 'deep_R.par', 'deep_I.par']
kcorrect, transpose(maggies), transpose(maggiesivar), newz, $
  band_shift=0.7, kcorrect, filterlist=filterlist, absmag=absmag

kcorrect, transpose(maggies), transpose(maggiesivar), $
  dat.z, band_shift=0.7, kcorrect, filterlist=filterlist, absmag=origabsmag

add = replicate( {absmag : fltarr(3), $
  origabsmag : fltarr(3)}, n_elements(dat))

add.absmag = absmag
add.origabsmag = origabsmag

dat = struct_addtags(dat, add)




mwrfits, dat, '/tmp/newdeep2.fits', /create


END
