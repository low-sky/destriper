#!/usr/bin/env python
import astropy.io.fits as fits
import numpy as np
#import matplotlib.pyplot as plt
from spectral_cube.spectral_cube import SpectralCube, LazyMask, BooleanArrayMask, Projection
#import spectral_cube.io
#from signal_id.noise import Noise
#from radio_beam import beam
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import scipy.stats as ss
#import scipy.signal as sg
from scipy import fftpack
import scipy.ndimage as ndimage
import os
import sys
import warnings
from astropy.wcs import WCS

def fft_plot(image):
    # Take the fourier transform of the image.
    F1 = fftpack.fft2(image)
    # Now shift the quadrants around so that low spatial frequencies are in
    # the center of the 2D fourier transformed image.
    F2 = fftpack.fftshift( F1 )
    return F2

def frequency_calc(image):
    # calculate cubes of frequency
    freqsx0 = np.fft.fftfreq(image.shape[0])
    freqsx00 = np.fft.fftshift(freqsx0)
    freqsy0 = np.fft.fftfreq(image.shape[1])
    freqsy00 = np.fft.fftshift(freqsy0)
    fxx0,fyy0 = np.meshgrid(freqsx00,freqsy00,indexing='ij')
    
    return fxx0,fyy0

def find_cutoff(image):
    cutoff = np.percentile(np.abs(image[np.isfinite(image)]),95.4499736)
    print 'peak cutoff: ',np.abs(cutoff)
    return cutoff

def find_peaks(fftimage,fxx,fyy):
    
    #mask inner area of emission. These need to be flipped depending on the
    #orientation of the noise in fft image. 

    # determine central emission region around frequency valye of 0 in x and y
    cutoff = np.abs(np.sqrt(fxx**2+fyy**2)) < 0.06

    # trim along axis that shows noise;  hardwired to be the y-axis here.
    cutoff2 = (np.abs(fyy) < 0.004) & (np.abs(fxx) > 0.03)
    
    #create 2D mask

    fftmask = np.ones(fftimage.shape)
    fftmask[cutoff] = 0.
    fftmask[cutoff2] = 1.
    fftmask[(fftmask == 0.)] = np.nan

    mask_image = np.copy(fftimage)
    maskpeak = mask_image*fftmask
    #determine scaling factor
    threshold = find_cutoff(maskpeak)
    #find peaks
    maxmap = ndimage.filters.maximum_filter(np.abs(maskpeak).real,size=[3,3])
    #maxima = (np.abs(maskpeak) == maxmap)
    data_min = ndimage.filters.median_filter(np.abs(maskpeak).real, size=[9,9])
    diffmap = (maxmap - data_min)
    diff = (diffmap > threshold)
   
    invc = np.invert(diff)
    invc1 = np.ones(invc.shape)
    invc1[~invc] = np.nan

#    print 'repairing ',invalids[0].shape[0],' values'
    return invc1

def inverse_transform(repaired_image):
    inverseshift0 = fftpack.ifftshift(repaired_image)
    inverse0 = fftpack.ifft2(inverseshift0)
    return inverse0.real

def fillbadcube(cube):
    for i in np.arange(cube.shape[0]):
        cube[i,:,:] = fillbad(cube[i,:,:])
    return(cube)

def fillbad(blanked):
    from scipy import interpolate
    valid_mask = ~np.isnan(blanked)
    coords = np.array(np.nonzero(valid_mask)).T
    values = blanked[valid_mask]
    it = interpolate.LinearNDInterpolator(coords, values,fill_value=0)
    output = it(list(np.ndindex(blanked.shape))).reshape(blanked.shape)
    return(output)

def mask_make(mask,image):
    # create map excluding central emission region
    noise_removed = image*mask
    coord_change = np.where(np.isnan(mask))
    blanked = np.copy(image)
    blanked *= mask

    # create median filter
    from scipy import interpolate
    valid_mask = ~np.isnan(blanked)
    coords = np.array(np.nonzero(valid_mask)).T
    values = np.log(np.abs(blanked[valid_mask]))

    it = interpolate.LinearNDInterpolator(coords, values, fill_value=0)
    output = np.exp(it(list(np.ndindex(image.shape))).reshape(image.shape))

    amplitude = output[coord_change[0],coord_change[1]]
    phase = np.random.rand(coord_change[0].shape[0])
    real = amplitude.real*np.cos(2*np.pi*phase)
    imag = amplitude.real*np.sin(2*np.pi*phase)
    #replace noise values, return new map
    newvals = np.vectorize(np.complex)(real,imag)
    noise_removed[coord_change[0],coord_change[1]] = newvals[:]

    return noise_removed,output

def fitswrite(initial_image, final_image, name):
    #header = initial_image.wcs.to_header()
    header = fits.getheader(initial_image)
    hdu = fits.PrimaryHDU(data = final_image, header=header)
    hdu.writeto(name,clobber=True)
    
def fftclean(InputFile,OutputFile=None,SaveDiagnostics=False,MomentOnly=False):
    print 'starting clean..'
    cubes = SpectralCube.read(InputFile,format='fits')
    root = InputFile.replace('.fits', '_fftclean')
    momentmap = cubes.moment0(axis=0)

    if SaveDiagnostics:
        hdu = fits.PrimaryHDU(data=momentmap.value)
        hdu.writeto(root+'_moment.fits',clobber=True)

    moment0 = np.ndarray.astype(momentmap.value,dtype=float)
    moment0 = fillbad(moment0)
    fft_image0 = fft_plot(moment0)

    if SaveDiagnostics:
        hdu = fits.PrimaryHDU(data=moment0)
        hdu.writeto(root+'_moment_fillbad.fits',clobber=True)
        hdu = fits.PrimaryHDU(data=np.abs(fft_image0))
        hdu.writeto(root+'_fftmoment.fits',clobber=True)

    # find frequency cubes
    fxx,fyy = frequency_calc(fft_image0)

    # determine map of noise spikes
    mask = find_peaks(fft_image0,fxx,fyy)
    print 'found peaks'

    # repair spikes in moment map
    repaired, output = mask_make(mask,fft_image0)

    if SaveDiagnostics:
        hdu = fits.PrimaryHDU(data=mask)
        hdu.writeto(root+'_fftmask.fits',clobber=True)
        hdu = fits.PrimaryHDU(data=np.abs(output))
        hdu.writeto(root+'_fftmedian.fits',clobber=True)
        hdu = fits.PrimaryHDU(data=np.abs(repaired))
        hdu.writeto(root+'_fftrepaired.fits',clobber=True)
        hdu = fits.PrimaryHDU(data=inverse_transform(repaired))
        hdu.writeto(root+'_repaired.fits',clobber=True)


    if MomentOnly:
        return

    datacube = fits.getdata(InputFile)
    datacube = np.nan_to_num(datacube)
    datacube = np.ndarray.astype(datacube,dtype=float)
    datacube = fillbadcube(datacube)



    fft_cube = np.ndarray(shape=datacube.shape,dtype=complex)
    for i in range(datacube.shape[0]):
        fft_cube[i,:,:] = fft_plot(datacube[i,:,:])
        
    print 'fft complete'


    # use mask to sample and repair noise regions in each channel of fft cube
    cube_repair = np.ndarray(shape=fft_cube.shape,dtype=complex)
    inverse_cube = np.ndarray(shape=fft_cube.shape)
    mediancube = np.ndarray(shape=fft_cube.shape)
    for i in range(fft_cube.shape[0]):
        new_vals, output = mask_make(mask,fft_cube[i,:,:])
        cube_repair[i,:,:] = new_vals
        mediancube[i,:,:] = output
        inverse_cube[i,:,:] = inverse_transform(cube_repair[i,:,:])
        
    head = fits.getheader(InputFile)
    cubewcs = WCS(head)
    OutputCube = SpectralCube(data=inverse_cube.astype(np.float32),wcs=cubewcs)
#    import pdb; pdb.set_trace()

    if isinstance(OutputFile,basestring):
        OutputCube.write(OutputFile)
        return(OutputCube)

    else:
        return(OutputCube)

#    import pdb ; pdb.set_trace()

if __name__=="__main__":
    InputFile = sys.argv[1]
    output_file = sys.argv[2]
    fftclean(InputFile,OutputFile = output_file)
