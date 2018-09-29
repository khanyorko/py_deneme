# İkinci dereceden (ya da istenen) bir fonksiyonun istenen aralıkta 
# dikdörtgenlerden yararlanarak belli bir doğrulukta integralinin
# hesaplanmasını sağlayacaktır
#
# x1 : Aralık başlangıcı
# x2 : Aralık sonu (dahil)
# d  : İstenen bölme sayısı
# linspace(0,10,5) ---> [0.0, 2.5, 5.0, 7.5, 10.0]
#
def linspace(x1, x2, d):
    try:
        return [x1 + ((x2 - x1) / (d - 1)) * (_) for _ in range(d)]
    except:
        return [x1, x2]
#        
# f  : İntegrali istenen denklem
# x1 : İntegral başlangıç değeri
# x2 : İntegral bitiş değeri
# d  : İstenen hassasiyete göre değişen dikdörtgen sayısı
#
#
#
#
def intg(f, x1, x2, d):
    toplam = 0
    aralık = linspace(x1, x2, d)
    for _ in range(len(aralık) - 1):
        toplam += (aralık[_ + 1] - aralık[_])*f(aralık[_])
    return toplam
#
#
# f(x) = 3x^2 + x + 2 
# 
# intg(f) = 1070
#
# intg(fonk,0,10,50  ) :  1038.5755935027073
# intg(fonk,0,10,150 ) :  1038.5755935027073
# intg(fonk,0,10,200 ) :  1062.223681220171
# intg(fonk,0,10,250 ) :  1063.783164787665
# intg(fonk,0,10,500 ) :  1066.895795599214
# intg(fonk,0,10,1000) :  1068.4489494499503
# intg(fonk,0,10,2500) :  1069.3798319647997
# intg(fonk,0,10,5000) :  1069.6899579955984
#
print('>>> ax^2 + bx + c')
a, b, c = list(map(float,input('>>> Katsayılar: ').split(' ')))
fonk = lambda x: a*x**2 + b*x + c
