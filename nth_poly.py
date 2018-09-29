import matplotlib.pyplot as plt

def strcat(args):
    ans = ''
    for _ in args:  
        if _ is args[-1]:                                       # köklerin tekrarlanması durumuna karşı >is< kullanıldı (!)
            ans += _                                            # gönderilen string'lerin aralarına * konur,
            break                                               # eval komutunda kullanılmak üzere birleştirilir
        ans += _ + '*'                                          # son elemana gelindiğinde * konmaz
    return ans

def arange(x1,x2,dx):
    l = list()                                                  # aralığın elemanlarının tutulacağı boş bir liste oluşturulur
    while x1 < x2:                                              # x1 sayısı x2 sayısından küçük olduğu sürece
        l.append(float('{:.5}'.format(float(x1))))              # istenen dx miktarına göre liste elemanlarla doldurulur
        x1 += dx
    return l

def drawPol(roots):
    intv    = arange(min(roots)-.3,max(roots)+.3,.1)            # çizilecek aralık arange(min-.3, max+.3, .1) şeklinde ayarlanır
    val_s   = [[f'({x}-{_})' for _ in roots] for x in intv]     # arange elemanları teker teker '(x0-x1)','(x0-x2)',..,'(x0-xn)' 
    val_ss  = list(map(strcat,val_s))                           # şeklinde val_s listesinde stringler olarak saklanır.
    val     = list(map(eval,val_ss))                            # val_s'in elemanları '(x0-x1)*(x0-x2)*..*(x0-xn)' şeklinde
    plt.plot(intv, val, roots,[0]*len(roots), 'bo')             # birleşir ve her birine eval uygulanarak val listesinde
    plt.show()                                                  # saklanır, son olarak pyplot ile polinom ve ayrıca 
                                                                # kökleri nokta olarak ekrana basılır

while True:
    try:
        roots = list(map(float,list(input('[Kökler]>>> ').split())))
        drawPol(roots)
    except KeyboardInterrupt:
        break
    except:
        print('>>> Kökler sayı olmalı ve boşluk ile ayrılmalı...')
pass
