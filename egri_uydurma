import numpy as np
import matplotlib.pyplot as plt

# verilen değerlere uyacak şekilde doğru veya eğri uydurma
#
# matris yapısı:
    #
    # [sk]*[an]=[ti]
    #
    #
        #      m   k              m   k
        # sk = Σ  x          tk = Σ [x  * f(x )]
        #     i=0  i             i=0  i      i
        #
            # |s0 s1   s2   ... sn   ||a0|   |t0|  
            # |s1 s2   s3   ... sn+1 ||a1|   |t1|
            # |. . . . . . . . . . . ||..| = |..|
            # |. . . . . . . . . . . ||..|   |..|
            # |sn sn+1 sn+2 ... s2n  ||an|   |tn|
                #
                #   istenen j. dereceden eğri/doğru için 
                #   elde edilecek denklem
                #
                #                        2           j
                #   g(x) = a  + a x + a x + ... + a x
                #           0    1     2           j
                #
                
def arange(x1,x2,dx):
    l = list()
    while x1 < x2:
        l.append(float('{:.5}'.format(float(x1))))
        x1 += dx
    return l
                
def sk(x, k):
    _sk = 0.0
    for _ in x:
        _sk += _**k
    return _sk

def tk(x, fx, k):
    _tk = 0.0
    for _,__ in zip(x, fx):
        _tk += __ * _**k
    return _tk

def gx(a,n):
    _gx = ''
    for _ in range(n + 1):
        _gx += f'{a[_][0]}*x**{_}'
        if _ < n:
            _gx += '+'
    return _gx.replace('[','').replace(']','')

x  = [1, 4,  6, 9, 10]                          # x    değerlerinin girilişi
fx = [4, 9, 15, 7,  3]                          # f(x) değerlerinin girilişi

n = int(input('>>> Uydurulacak eğri derecesi :  '))

__sk = list()
__tk = list()

for _ in range(n + 1):
    ___sk = list()
    for __ in range(n + 1):
        ___sk.append(sk(x,_ + __))
    __sk.append(___sk)
    
for _ in range(n + 1):
    ___tk = list()
    for __ in range(1):
        ___tk.append(tk(x,fx,_ + __))
    __tk.append(___tk)
    
_sk = np.matrix(__sk)
_tk = np.matrix(__tk)

_a  = np.dot(_sk.I,_tk)

_gx = lambda x: eval(gx(_a,n))

val = list()

for _ in arange(min(x),max(x),.005):
    val.append(_gx(_))


plt.plot(x, fx, 'bo',arange(min(x),max(x),.005), val,'--')
plt.show()
