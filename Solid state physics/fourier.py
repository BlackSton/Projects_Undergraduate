import numpy as np
import cv2
import matplotlib.pyplot as plt

im = cv2.imread('6.png',cv2.IMREAD_GRAYSCALE)
#img = cv2.bitwise_not(im)
f=np.fft.fft2(im)
fshift = np.fft.fftshift(f)
plt.imshow(np.abs(fshift),cmap=plt.cm.gray)
plt.axis('off')
fig = plt.gcf()
fig.savefig('6F.png',dpi=500,transparent= True)
plt.show()
