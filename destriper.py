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
    fxx0,fyy0 = np.meshgrid(freqsx00,freqsy00)
    
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
    data_min = ndimage.filters.median_filter(np.abs(maskpeak).real, size=[3,3])
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

def mask_make(mask,image):

    # create map excluding central emission region
    noise_removed = image*mask
    coord_change = np.where( np.isnan(mask))

    # create median filter
    output = ndimage.filters.median_filter(np.absolute(image).real,size=[1,15])

    # determine amplitude and random phase values to replace noise spikes
    amplitude = output[coord_change[0],coord_change[1]]
    phase = np.random.rand(coord_change[0].shape[0])
    real = amplitude.real*np.cos(2*np.pi*phase)
    imag = amplitude.real*np.sin(2*np.pi*phase)
    #replace noise values, return new map
    newvals = np.vectorize(np.complex)(real,imag)
#    newvals = np.complex(real,imag)
    noise_removed[coord_change[0],coord_change[1]] = newvals[:]
        
    return noise_removed, output

def fitswrite(initial_image, final_image, name):
    #header = initial_image.wcs.to_header()
    header = fits.getheader(initial_image)
    hdu = fits.PrimaryHDU(data = final_image, header=header)
    hdu.writeto(name)
    
def fftclean(InputFile,OutputFile=None):
    print 'starting clean..'
    cubes = SpectralCube.read(InputFile,format='fits')
    momentmap = cubes.moment0(axis=0)
    image00 = np.nan_to_num(momentmap.value)

    image000 = np.ndarray.astype(image00,dtype=float)

    image3 = fits.getdata(InputFile)
    image30 = np.nan_to_num(image3)
    image300 = np.ndarray.astype(image30,dtype=float)
            
    fft_image0 = fft_plot(image000)
   # hdu = fits.PrimaryHDU(data=np.abs(fft_image0))
   # hdu.writeto('a24mom0fft.fits')
    fft_cube30 = np.ndarray(shape=image300.shape,dtype=complex)
    for i in range(image300.shape[0]):
        fft_cube30[i,:,:] = fft_plot(image300[i,:,:])
        
    print 'fft complete'

    # find frequency cubes
    fxx,fyy = frequency_calc(fft_image0)

    # determine map of noise spikes
    mask = find_peaks(fft_image0,fxx,fyy)
    print 'found peaks'

    # repair spikes in moment map
    repaired, output = mask_make(mask,fft_image0)
    hdu = fits.PrimaryHDU(data=mask)
#    hdu.writeto('a24mom0mask.fits')
   # hdu = fits.PrimaryHDU(data=np.abs(output))
   # hdu.writeto('a24mom0fftmedian.fits')
    # use mask to sample and repair noise regions in each channel of fft cube
    cube30_repair = np.ndarray(shape=fft_cube30.shape,dtype=complex)
    inverse_cube = np.ndarray(shape=fft_cube30.shape)
    mediancube = np.ndarray(shape=fft_cube30.shape)
    for i in range(fft_cube30.shape[0]):
        new_vals, output = mask_make(mask,fft_cube30[i,:,:])
        cube30_repair[i,:,:] = new_vals
        mediancube[i,:,:] = output
        inverse_cube[i,:,:] = inverse_transform(cube30_repair[i,:,:])
        
    head = fits.getheader(InputFile)
    cubewcs = WCS(head)
    OutputCube = SpectralCube(data=inverse_cube.astype(np.float32),wcs=cubewcs)

    if isinstance(OutFile,basestring):
        OutputCube.write(OutFile)
        return(OutputCube)

    else:
        return(OutputCube)

    import pdb ; pdb.set_trace()

if __name__=="__main__":
    InputFile = sys.argv[1]
    output_file = sys.argv[2]
    fftclean(InputFile,OutputFile = output_file)
