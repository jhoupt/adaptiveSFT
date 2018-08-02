from __future__ import division
def Lab2RGB(L, a, b):
#LAB2RGB Convert an image from CIELAB to RGB
#Lab2RGB takes L, a, and b double matrices and returns an image in the RGB color space.  
#Values for L are in the range [0,100] while a* and b* are roughly in the range [-110,110].

#R, G, and B are the values will be returned as doubles in the range [0,1]
#RGB is be uint8s in the range [0,255].

#This transform is based on ITU-R Recommendation BT.709 using the D65
#white point reference. The error in transforming RGB -> Lab -> RGB is
#approximately 10^-5.  
#Elizabeth Fox - Translated from MATLAB function lab2rgb.m (22 October 2017)


  # Import important packages/functions
  import numpy as np 
  from numpy import arange, empty 
  from numpy.random import randint
  
  
  ## Thresholds
  T1 = 0.008856
  T2 = 0.206893
  
  if type(L) is tuple or type(L) is np.ndarray:
    M = L.shape[0]
    N = L.shape[1]
    RGB = empty([3,M,N], float)
  else :
    M = 1
    N = 1
    RGB = empty(3, float)
  s = M * N
  
  L =np.array(L).flatten() * 1.
  a =np.array(a).flatten() * 1.
  b =np.array(b).flatten() * 1.
  
  # Compute Y
  fY = pow((L + 16) /116, 3)
  YT = fY > T1
  fY = (~YT) * (L / 903.3) + YT * fY
  Y = fY
  
  #Alter fY slightly for further calculations
  fY = YT * pow(fY,1/3) + (~YT) * (7.787 * fY + 16/116)
  
  # Compute X
  fX = a / 500 + fY
  XT = fX > T2
  X = (XT * pow(fX, 3) + (~XT) * ((fX - 16/116) / 7.787))
  
  # Compute Z
  fZ = fY - b / 200
  ZT = fZ > T2
  Z = (ZT * pow(fZ, 3) + (~ZT) * ((fZ - 16/116) / 7.787))
  
  # Normalize for D65 white point
  X = X * 0.950456
  Z = Z * 1.088754
  
  #XYZ to RGB
  MAT = np.matrix([[3.240479,-1.537150,-0.498535],[-0.969256,1.875992,0.041556],[0.055648,-0.204043,1.057311]])
  
  XYZ = np.vstack((X,Y,Z))
  
  tempRGB = np.maximum(np.minimum(np.dot(MAT, XYZ),1), 0)
  
  R = tempRGB[0].reshape(M,N)
  G = tempRGB[1].reshape(M,N)
  B = tempRGB[2].reshape(M,N)
  
  
  RGB[0] = R
  RGB[1] = G
  RGB[2] = B
  RGB = np.array(np.round(RGB*255), dtype=np.uint8)
  return(RGB)
