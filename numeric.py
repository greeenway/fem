import math
from math import sqrt as sqrt

#gauss quadrature weights
alpha = [5./18, 8./18, 5./18]
xi = [(5-math.sqrt(15))/10., 1/2., (5+math.sqrt(15))/10.]

def gauss_2D_3(func):
    """Applies 2D gaussian quadrature to approximate the integrand <func> on [0]1)x[0]1) """
    res = 0.
    for i in range(0,3):
        for j in range(0,3):
            res += alpha[i] * alpha[j] * func(xi[i], xi[j])
    
    return res


#triangle
#page 313
alpha_t = [(155. - sqrt(15)) / 2400,
           (155. - sqrt(15)) / 2400,
           (155. - sqrt(15)) / 2400,
           (155. + sqrt(15)) / 2400,
           (155. + sqrt(15)) / 2400,
           (155. + sqrt(15)) / 2400,
           9./80,]
xi_t = [ ( (6-sqrt(15))/21., (6-sqrt(15))/21. ),
         ( (9+2*sqrt(15.))/21., (6-sqrt(15.))/21. ),
         ( (6-sqrt(15))/21., (9+2*sqrt(15))/21. ),
         ( (6+sqrt(15))/21., (9-2*sqrt(15))/21. ),
         ( (6+sqrt(15))/21., (6+sqrt(15))/21. ),
         ( (9-2*sqrt(15))/21., (6+sqrt(15))/21. ),
         ( (1./3., 1./3.) ),
        ] 

def gauss_triangle_6(func):
    """Applies 2D gaussian quadrature to approximate the integrand <func> on (0,0)-(1,0)-(0,1) """
    res = 0.

    for i in range(0,7):
        #for j in range(0,7):
        res += alpha_t[i] * func( xi_t[i][0], xi_t[i][1] )

    return res


#gauss 1d quadrature weights (equal to 2d)
alpha1 = [5./18, 8./18, 5./18]
xi1 = [(5-math.sqrt(15))/10., 1/2., (5+math.sqrt(15))/10.]

def gauss_1D_6(func):
    """Applies 1D gaussian quadrature to approximate the integrand <func> on (0,0) """
    res = 0.

    for i in range(0,3):
        res += alpha[i] * func(xi[i])
    
    return res




if __name__ == '__main__':
    #
    f = lambda x,y: x**2 +y**2 -2*x*y + x

    res = gauss_triangle_6(f)
    print(res)
    print(1/4.)

    f1D = lambda x: x**2-x+3
    print(gauss_1D_6(f1D))
    print(17./6)






