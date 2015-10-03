def char_sum(p,order, alg = 'complex-analytic'):
    """
    p -- a prime
    psi = psi_b is an additive character of F_p.
    psi_b(x) = exp(2 pi i Tr(bx)).
    """
    print 'order = %s'%order
    if order == 1:
        raise NotImplementedError
    if Mod(p+1, order) != 0:
        raise ValueError('p+1 = %s, order=  %s'%(p+1, order))
    m = order
    K.<zeta> = CyclotomicField(p*m)
    zetap, zetam = zeta**m, zeta**p

    prec = max(10*p,200)

    F.<a> = GF(p**2)
    gen = F.multiplicative_generator()
    Fcross = [gen**i for i in range(p^2-1)] # all elements in F^\times.

    if alg == 'complex-analytic':
        v = list(range(1,p))
        result = 1
        phi = K.complex_embedding(prec)
        zetap, zetam = phi(zeta**m), phi(zeta**p)
        for sigma in v:
            # print('sigma = %s'%sigma)
            result *= sum([ zetap**(sigma*(ZZ(Fcross[i] + Fcross[i]**p + Fcross[i]**(p+1))))*zetam**(i) for i in range(p^2-1)]) # conjugates.
        return round(result.real_part())
    else:
        result = sum([ zetap**(Fcross[i] + Fcross[i]**p + Fcross[i]**(p+1))*zetam**(i) for i in range(p^2-1)]) # conjugates.
        f = result.minpoly()
        return ZZ(f[0])


# 1. How large is T? (p^2-1)p





# to-do: fix the character sum so that it works for any character \psi_b.