def is_gamma_equiv(cusp1,cusp2,level):
    """
    return flag, matrix g such that
    flag = if cusp 1 is equivalent to cusp 2
    modulo gamma(level), and g(cusp1) = cusp2
    """
    N = level
    a1,c1 = cusp1.numerator(),cusp1.denominator()
    a2,c2 = cusp2.numerator(),cusp1.denominator()
    if (a1 + a2 ) % N == 0 and (c1+c2) % N == 0:
        a2 = - a2
        c2 = -c2

    if (a1 - a2 ) % N == 0 and (c1 -c2) % N == 0:
        _,d,b = xgcd(a1,c1)
        assert a1*d+c1*b == 1
        A = matrix([[a1,-b],[c1,d]])
        verbose('A = %s'%A)
        cusp3 = cusp2.apply(A.inverse().list()) # cusp3 =  A^{-1}(cusp2).
        a3,c3 = cusp3.numerator(),cusp3.denominator()
        verbose('a3,c3 = %s,%s'%(a3,c3))
        _,x,y = xgcd(a3,c3*N)
        B = matrix([[a3,-N*y],[c3,x]])
        verbose('B = %s'%B) # B(oo) = C3, and B \in Gamma(N)

        D = A*B*A.inverse()
        assert cusp1.apply(D.list()) == cusp2
        if D[0][0] < 0:
            D = -D
        return (True,D)

    else:
        return False