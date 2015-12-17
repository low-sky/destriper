    # plotting final images
    
    # plot fft
    #p.figure(1,figsize=(8,8))
    #p.clf()
    #p.jet()
    #axes = p.gca()
    #im = axes.imshow(np.log10(np.abs(fft_image0)),extent=[fxx.min(),fxx.max(),fyy.min(),fyy.max()],vmin=0,vmax=6.5)
    #divider = make_axes_locatable(axes)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #p.colorbar(im, cax=cax)
    # plot repaired fft
    #p.figure(2,figsize=(8,8))
    #p.clf()
    #p.jet()
    #axes = p.gca()
    #im = axes.imshow(np.log10(np.abs(cube30_repair.sum(axis=0))),extent=[fxx.min(),fxx.max(),fyy.min(),fyy.max()])
    #divider = make_axes_locatable(axes)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #p.colorbar(im, cax=cax)
    #p.show()
    #i_repaired = inverse_transform(repaired)
    #inv_repaired = i_repaired.astype(np.float32)
    #cutoff = find_cutoff(fft_image0)
    #peaks_x,peaks_y = find_peaks(cutoff,fft_image0,fxx,fyy)
   
    #newvals = np.zeros(peaks_x.size,dtype=complex)

    #cube30_repair = np.ndarray(shape=fft_cube30.shape,dtype=complex)
    #inverse_cube = np.ndarray(shape=fft_cube30.shape,dtype=complex)
    #for i in range(fft_cube30.shape[0]):
   #     new_vals = mask_make(peaks_x,peaks_y,fft_cube30[i,:,:])
   #     cube30_repair[i,:,:] = new_vals

   #     shift = fftpack.ifftshift(cube30_repair[i,:,:])
   #     inverse_cube[i,:,:] = fftpack.ifft2(shift)

    # cleaned final
    #p.figure(3,figsize=(8,8))
    #p.clf()
    #p.copper()
    #axes = p.gca()
    #im = axes.imshow(inv_repaired.real,vmin=-75000,vmax=125000)
    #divider = make_axes_locatable(axes)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #p.colorbar(im, cax=cax)

    # mask
   # p.figure(4,figsize=(8,8))
   # p.clf()
   # p.jet()
   # axes = p.gca()
   # im = axes.imshow(mask,extent=[fxx.min(),fxx.max(),fyy.min(),fyy.max()])

    # moment map initial
   # p.figure(5,figsize=(8,8))
   # p.clf()
   # p.copper()
   # axes = p.gca()
   # im = axes.imshow(momentmap.value,vmin=-75000,vmax=125000)
   # divider = make_axes_locatable(axes)
   # cax = divider.append_axes("right", size="5%", pad=0.05)
  #  p.colorbar(im, cax=cax)
  #  p.show()

    # median cube
   # p.figure(6,figsize=(8,8))
   # p.clf()
   # p.jet()
   # axes = p.gca()
   # im = axes.imshow(np.log10(np.abs(mediancube[0,:,:])),extent=[fxx.min(),fxx.max(),fyy.min(),fyy.max()])
    #divider = make_axes_locatable(axes)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #p.colorbar(im, cax=cax)
    #p.show()
    print 'image repaired'
#axes = p.gca()
    #im = axes.imshow(np.log10(np.abs(repaired)),extent=[fxx.min(),fxx.max(),fyy.min(),fyy.max()])
    #divider = make_axes_locatable(axes)
   # cax = divider.append_axes("right", size="5%", pad=0.05)
   # p.colorbar(im, cax=cax)
    #write(input_file,inverse_cube2,output_file)

    print 'output file created'
    
    # Now plot up both
    """
    p.figure(3,figsize=(8,8))
    p.clf()
    p.imshow( image00 )
    p.colorbar()
    
    p.figure(4,figsize=(8,8))
    p.clf()
    p.imshow( inverse_cube.real.sum(axis=0) )
    p.colorbar()
    p.show()"""
