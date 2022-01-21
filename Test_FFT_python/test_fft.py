import numpy as np

data = np.loadtxt('data.txt')
#data = np.ones(100).reshape(10,10)
dim = np.shape(np.fft.fft(data))
data = np.fft.fft(data)
#data = np.fft.fft(data.T)

print(data)
