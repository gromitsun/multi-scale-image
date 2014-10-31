import numpy as np
import matplotlib.pyplot as plt
def tukey(alpha,N):
    w=np.empty(N)
    for n in range(N):
        if 0<=n<alpha*(N-1)/2.:
            w[n]=1/2.*(1+np.cos(np.pi*(2*n/alpha/(N-1.)-1)))
        elif alpha*(N-1)/2.<=n<=(N-1)*(1-alpha/2.):
            w[n]=1
        else:
            w[n]=1/2.*(1+np.cos(np.pi*(2*n/alpha/(N-1.)-2./alpha+1)))
    return w


def tukey2d(alpha,M,N):
    x=tukey(alpha,N)
    y=tukey(alpha,M)
    X,Y=np.meshgrid(x,y)
    return X*Y

plt.imshow(tukey2d(0.5,100,100),cmap='gray')
plt.xticks([0,99],['0','N-1'])
plt.yticks([0,99],['0','N-1'])
plt.colorbar()
plt.show()
