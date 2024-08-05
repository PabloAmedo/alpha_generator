# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:53:49 2024

@author: diego
"""
print('Running...')
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

im1path = '../tiffs/July/highgain_compare/100ms/ss_single_1.tiff'
im2path = '../tiffs/July/highgain_compare/100ms/ss_single_2.tiff'
im3path = '../tiffs/July/highgain_compare/100ms/ss_single_3.tiff'
im4path = '../tiffs/July/highgain_compare/100ms/ss_single_4.tiff'
im5path = '../tiffs/July/highgain_compare/100ms/ss_single_5.tiff'
im6path = '../tiffs/July/highgain_compare/100ms/ss_single_6.tiff'
im7path = '../tiffs/July/highgain_compare/100ms/ss_single_7.tiff'
im8path = '../tiffs/July/highgain_compare/100ms/ss_single_8.tiff'
im9path = '../tiffs/July/highgain_compare/100ms/ss_single_9.tiff'
im10path = '../tiffs/July/highgain_compare/100ms/ss_single_10.tiff'


im1 = np.array(Image.open(im1path))
im2 = np.array(Image.open(im2path))
im3 = np.array(Image.open(im3path))
im4 = np.array(Image.open(im4path))
im5 = np.array(Image.open(im5path))
im6 = np.array(Image.open(im6path))
im7 = np.array(Image.open(im7path))
im8 = np.array(Image.open(im8path))
im9 = np.array(Image.open(im9path))
im10 = np.array(Image.open(im10path))

sum_im = im1 + im2 + im3 + im4 + im5 + im6 +im7 + im8 + im9 +im10

image_final = Image.fromarray(sum_im)
image_final.save('../tiffs/July/combinned_image.tiff')


plt.figure()
plt.imshow(image_final)

