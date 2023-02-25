# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 13:07:48 2022

@author: berej
"""

import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from matplotlib import image
from astropy.io import fits
from lmfit import Model, fit_report
import glob2
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import median_filter

single_image = 0
threshold = 5
path = os.getcwd()
os.chdir(path+'\\mrk304')
file_list = os.listdir()

print('======================================================')
print('The current version of the script doesn\'t dedark the flat or debias the science frames because of a problem with the QHY600 bias frames')
print('======================================================')

# =============================================================================
# if 'master_bias.fits' in file_list:
#     print('Master bias already exists, skipping...')
# else:
#     list_bias = glob2.glob('bias/*fits')
#     nb_bias = len(list_bias)
#     
#     hdu = fits.open(list_bias[0])
#     hdu0 = hdu[0]
#     data = hdu[0].data
#     hdu.close()
#     
#     telescope = hdu0.header.get('TELESCOP').strip()
#     acq_type = hdu0.header.get('ACQTYPE').strip()
#     epoch = hdu0.header.get('DATE-OBS')
#     exptime = hdu0.header.get('EXPTIME')
#     size_x  = hdu0.header.get('NAXIS1')
#     size_y  = hdu0.header.get('NAXIS2')
#     
#     print(f'{acq_type} ({size_x}x{size_y}) taken with {telescope} on {epoch} ({exptime}s exposure).')
#     
#     print('Creating master bias via np.median...')
#     
#     biases = np.empty((size_x, size_y, nb_bias))
#     for i in range(nb_bias):
#         print(list_bias[i])
#         
#         # Read FITS and store first hdu into hdu variable
#         hdu = fits.open(list_bias[i])[0]
#         
#         # Each slice of the darks matrix is a dark frame
#         biases[:,:,i] = hdu.data.T
#         
#         
#     master_bias = np.median(biases, axis = 2)
#     
#     print('Saving to disk...')
#     
#     # Save to disk
#     hdu.data = master_bias.T
#     hdu.writeto('master_bias.fits', overwrite=True)
# =============================================================================

if 'master_dark.fits' in file_list:
    print('Master dark already exists, skipping...')
    
else:
    
    list_dark = glob2.glob('darks/*fits')
    nb_dark = len(list_dark)
    
    hdu = fits.open(list_dark[0])
    hdu0 = hdu[0]
    data = hdu[0].data
    hdu.close()
    
    telescope = hdu0.header.get('TELESCOP').strip()
    acq_type = hdu0.header.get('ACQTYPE').strip()
    epoch = hdu0.header.get('DATE-OBS')
    exptime = hdu0.header.get('EXPTIME')
    size_x  = hdu0.header.get('NAXIS1')
    size_y  = hdu0.header.get('NAXIS2')
    
    print(f'{acq_type} ({size_x}x{size_y}) taken with {telescope} on {epoch} ({exptime}s exposure).')
    
    print('Creating master dark via np.median...')
    
    darks = np.empty((size_x, size_y, nb_dark))
    for i in range(nb_dark):
        print(list_dark[i])
        
        # Read FITS and store first hdu into hdu variable
        hdu = fits.open(list_dark[i])[0]
        
        # Each slice of the darks matrix is a dark frame
        darks[:,:,i] = hdu.data.T
        
        
    master_dark = np.median(darks, axis = 2)
    
    print('Calculating hot pixels...')
    
    threshold = 5
    
    bkg_mean, bkg_median, bkg_sigma = sigma_clipped_stats(master_dark, sigma=3.0)
    
    hot_pixel = np.where(master_dark > bkg_median + threshold * bkg_sigma)
        
    print( f'Number of pixel in the dark: {len(master_dark.flatten()):8d}')
    print( f'Number of hot pixels       : {len(master_dark[hot_pixel]):8d}')
    print( f'Fraction of hot pixels (%) : {100*len(master_dark[hot_pixel])/len(master_dark.flatten()):.2f}' )
    
    # =============================================================================
    # First method, we smooth the hot pixels with surrounding pixels
    # =============================================================================
    
    print('Smoothing dark via 5x5 median filter...')
    
    # Create a smoothed version of the frame
    smoothed_dark = median_filter(master_dark, (5,5))
    
    # Replace the value of hot pixels in the frame by the value in the smoothed-version of the frame
    master_dark[hot_pixel] = smoothed_dark[hot_pixel]
    
    # Create a binary map: 0:Ok, 1:Hot pixel
    map_hot = master_dark*0
    map_hot[hot_pixel] = 1
    
    print('Saving to disk...')
    
    # Save to disk
    hdu.data = master_dark.T
    hdu.writeto('master_dark.fits', overwrite=True)
    
    hdu.data = map_hot.T
    hdu.writeto('bad_pixels_hot.fits', overwrite=True)
    
if 'master_flat.fits' in file_list:
    print('Master flat already exists, skipping...')
    
else:
    
    list_flat = glob2.glob('flats/*fits')
    
    hdu = fits.open(list_flat[0])
    hdu0 = hdu[0]
    data = hdu[0].data
    hdu.close()
    
    telescope = hdu0.header.get('TELESCOP').strip()
    acq_type = hdu0.header.get('ACQTYPE').strip()
    epoch = hdu0.header.get('DATE-OBS')
    exptime = hdu0.header.get('EXPTIME')
    size_x  = hdu0.header.get('NAXIS1')
    size_y  = hdu0.header.get('NAXIS2')
    
    print(f'{acq_type} ({size_x}x{size_y}) taken with {telescope} on {epoch} ({exptime}s exposure).')
    
    print('Creating master flat via np.median...')
    
    nb_flat = len(list_flat)
    
    hdu_dark = fits.open('master_dark.fits')
    dark = hdu_dark[0]
    t_dark = dark.header['EXPTIME']
    
    flats = np.empty((size_x, size_y, nb_flat))
    for i in range(nb_flat):
        print(list_flat[i])
        
        # Read FITS and store first hdu into hdu variable
        hdu = fits.open(list_flat[i])[0]
        
        # Each slice of the flats matrix is a flat frame
        #flats[:,:,i] = hdu.data.T - dark.data.T * (hdu.header['EXPTIME']/dark.header['EXPTIME'])
        
# =============================================================================
#         dark_ratio = hdu.header['EXPTIME']/t_dark
#         if dark_ratio < 2 and mode == 'photo':
#             new_dark = (dark.data.T - bias.data.T) * hdu.header['EXPTIME']/t_dark + bias.data.T
#             flats[:,:,i] = hdu.data.T - new_dark
#             flats[:,:,i] /= np.median(flats[:,:,i])
#         else:
#             print('Ratio between flat exp time and dark exp time too large, skipping dark correction...')
# =============================================================================
        flats[:,:,i] = hdu.data.T / np.median(hdu.data.T)
        
    #hdu_dark.close()
    master_flat = np.median(flats, axis=2)
    
    print('Calculating dead pixels...')
    
    bkg_flat_mean, bkg_flat_median, bkg_flat_sigma = sigma_clipped_stats(master_flat, sigma=3.0)
    
    dead_pixel = np.where(master_flat < bkg_flat_median-threshold*bkg_flat_sigma)
    
    print( f'Number of flat pixels       : {len(master_flat.flatten()):8d}')
    print( f'Number of dead pixels       : {len(master_flat[dead_pixel]):8d}')
    print( f'Fraction of dead pixels (%) : {100*len(master_flat[dead_pixel])/len(master_flat.flatten()):.2f}')
    
    print('Smoothing flat via 5x5 median filter...')
    
    # Create a smoothed version of the frame
    smoothed_flat = median_filter(master_flat, (5,5))
    
    # Replace the value of hot pixels in the frame by the value in the smoothed-version of the frame
    master_flat[dead_pixel] = smoothed_flat[dead_pixel]
    
    # Create a binary map: 0:Ok, 1:Dead pixel
    map_dead = master_flat*0
    map_dead[dead_pixel] = 1
    
    print('Saving to disk...')
    
    # Save to disk
    hdu.data = master_flat.T
    hdu.writeto('master_flat.fits', overwrite=True)
    
    hdu.data = map_dead.T
    hdu.writeto('bad_pixels_dead.fits', overwrite=True)

print('Reducing science image(s)...')

if 'output' not in os.listdir():
    os.mkdir('output')

hdu = fits.open(glob2.glob('science/*fits')[0])
hdu0 = hdu[0]
data = hdu[0].data
hdu.close()

telescope = hdu0.header.get('TELESCOP').strip()
acq_type = hdu0.header.get('ACQTYPE').strip()
epoch = hdu0.header.get('DATE-OBS')
exptime = hdu0.header.get('EXPTIME')
size_x  = hdu0.header.get('NAXIS1')
size_y  = hdu0.header.get('NAXIS2')

print(f'{acq_type} ({size_x}x{size_y}) taken with {telescope} on {epoch} ({exptime}s exposure).')

darkfits = fits.open('master_dark.fits')[0]
master_dark = darkfits.data
t_dark = darkfits.header['EXPTIME']
flatfits = fits.open('master_flat.fits')[0]
master_flat = flatfits.data
t_flat = flatfits.header['EXPTIME']

hp_map = fits.open('bad_pixels_hot.fits')[0]
dp_map = fits.open('bad_pixels_dead.fits')[0]
hot_pixel = np.where(hp_map.data == 1)
dead_pixel = np.where(dp_map.data == 1)

list_science = glob2.glob('science/*fits')
nb_science = len(list_science)  

if single_image:
    print('Single image flag active, looping through...')
    print('This may take some time...')
    
    for i in range(nb_science):

        print(list_science[i])
        
        hdu = fits.open(list_science[i])[0]
        
        #reduced = ((hdu.data - master_bias) - (hdu.header['EXPTIME']/t_dark)*(master_dark - master_bias)) /     master_flat
        
        reduced = (hdu.data - master_dark) / master_flat
        
        smoothed = median_filter(reduced, (5,5))
        
        reduced[hot_pixel] = smoothed[hot_pixel]
        reduced[dead_pixel] = smoothed[dead_pixel]
        
        hdu_reduced = hdu.copy()
        hdu_reduced.data = reduced
        hdu_reduced.writeto(r'.\output\reduced_'+list_science[i].lstrip('science\\'), overwrite = True)
        
        print(i+1, '/', nb_science, 'reduced')
            
else:
    print('Creating a Z projected stack...')
    
    finals = np.empty((size_x, size_y, nb_science))
    
    for i in range(nb_science):
        
        print(list_science[i])
        
        hdu = fits.open(list_science[i])[0]
        
        finals[:,:,i] = hdu.data.T
        
    final = np.median(finals, axis = 2)
    
    print('Reducing...')
    
    #reduced = ((final.T - master_bias) - (hdu.header['EXPTIME']/t_dark)*(master_dark - master_bias)) / master_flat
    
    reduced = (final.T - master_dark) / master_flat
    
    smoothed = median_filter(reduced, (5,5))
    
    reduced[hot_pixel] = smoothed[hot_pixel]
    reduced[dead_pixel] = smoothed[dead_pixel]
    
    hdu_reduced = hdu.copy()
    hdu_reduced.data = reduced
    hdu_reduced.writeto(r'.\output\reduced_'+list_science[0].lstrip('science\\').split('_')[0]+'.fits', overwrite = True)
    
    
print('Done!')