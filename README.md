# Moire

## Introduction
Final Year Project in investigating the Moire Effect. The idea is adopted from Professor Thierry Blu, EE, CUHK.

## main.m
Main program, investigating the removal of moire patterns.

'lena' is  encoded into moire patterns with phase shifting.

The moire patterns are retrieved and multiplied with the encoded image to decode.

Yet, the decoded image is full of moire patterns.

Different filters are tried to remove the patterns and obtain the original image.

## rgb_and_pca.m
Splitting the image into rgb channels to see if it can reduce moire patterns. Principal component analysis is also used.

## imfrest2.m and immax.m
Finding moire patterns of an image.

## ft2D
Trying to remove moire patterns of an image with multiple methods.