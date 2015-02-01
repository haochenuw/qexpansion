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


def prime_onemod(N):
    """
    return the smallest prime that is 1 modulo N
    """
    p = next_prime(N)
    while p % N != 1:
        p = next_prime(p)
    return p


def periodgC(f,gg,phi,prec,test = False,conjugate = False):
    """
    given a newform f attached to elliptic
    curve E, an element gg in Gamma1(N),
    computes the period
    2*pi*i*\int_{a}^{gg(a)} f(z) dz
    by truncating the q-expansion of f to a polynomial
    with degree = terms depending on the lower left entry of g, a is the optimal point in [Cremona97].
    """
    a=gg[0,0]
    b=gg[0,1]
    c=gg[1,0]
    d=gg[1,1]

    if c == 0: # parabolic matrices have zero period, stated in Cremona.
        return 0
    terms = ZZ((prec*abs(c)/10).floor()) # temporary
    if test:
        terms = terms // 5

    f = f.q_expansion(terms)
    C = ComplexField(prec)

    ff = f.shift(-1).integral()
    v = list(ff.polynomial())
    K = v[0].parent()
    if K is not QQ:
        #phi = K.complex_embeddings()[0]
        q = var('q')
        if not conjugate:
            w = [phi(t) for t in v]
        else:
            w = [phi(t).conjugate() for t in v]
        ffpC = C[q](w)


    else:
        ffpc = ff.polynomial()

    tau=-(d/c)+(C(I)/abs(c))

    gammatau=(a*tau+b)/(c*tau+d)
    return ffpC.substitute(exp(2*C(pi)*C(I)*gammatau)) - ffpC.substitute(exp(2*C(pi)*C(I)*tau))


def numerical_global_constant(f,phi,prec = 53):
    t = cputime()
    N = f.level()
    C = ComplexField(prec)

    v = [matrix(2,2,g.matrix().list()) for g in list(Gamma1(N).gens()) if abs(g[1][0]) < C(N*10)]
    cand1 = []
    for g in v:
        pg = periodgC(f,g,phi,prec,test = True,conjugate = False)
        if abs(pg) > C(1.0/10):
            cand1.append(g)
            verbose('g= %s'%g)
            verbose('pg =  %s'%pg)
    verbose('length1 = %s'%len(cand1))

    cand2 = []

    for g in cand1:
        pgnew = periodgC(f,g,phi,prec,test = False,conjugate = False)
        if abs(pgnew) > C(1.0)/10:
            cand2.append((g,pgnew))
    verbose('length2 = %s'%len(cand2))


    cand2.sort(key = lambda a: abs(a[0][0][1]))

    ratios = []
    for g,pgnew in cand2[:3]:
        gmodified = modify(g)
        print g,gmodified,pgnew
        wNg = wN(gmodified,N)
        print 'wNg = %s'%wNg
        p1 = periodgC(f,wNg,prec,test = False,conjugate = True)
        print 'ratio = %s'%C(pgnew/p1)
        ratios.append((pgnew/p1))
    mu = mean(ratios)
    deviation = sum([abs(ratio-mu) for ratio in ratios])
    if deviation < C(1.0/1000):
        verbose('time = %s'%cputime(t))
        return ratio
    else:
        raise ValueError('did not stabilize')

def fricke(N):
    return matrix([[0,-1],[N,0]])


def wN(g,N):
    try:
        g  = g.matrix()
    except:
        pass
    return fricke(N)*(g*~fricke(N))