load('PeriodMapping.sage')

def generate_points(N,number_of_points = 5,prec = 53):
    """
    return some random points on the arc |z| = 1/sqrt(N), im(z) > 0.
    """
    C = ComplexField(prec)
    try:
        N = ZZ(N)
    except Exception:
        print N.parent()
        raise ValueError('N(=%s) must be an integer.'%N)

    thetas = [C(0.4*random() + 0.1) for _ in range(number_of_points)]
    verbose('A list of random angles theta computed')
    zs = [C(1/sqrt(N))*C(exp(C(pi)*C(I)*theta)) for theta in thetas]
    verbose('upper half plane points obtained.')
    return zs


def period(f,phi,a,b,terms = None, prec =53, exponent = 0):
    """
    helper funcion. Evaluate the following integral

    \int_a^b 2 \pi i f(z) z^{exponent} dz

    using integration by parts.

    """
    try:
        e = ZZ(exponent)
    except:
        raise ValueError
    if e < 0:
        raise ValueError('exponent must be non-negative')
    if e == 0:
        return _period(f,phi,a,b,terms =terms,prec= prec)
    if terms is None:
        terms = prec*50 # Linear growth is good, see William Stein: chow heegner paper.
    try:
        # if f is not a q-expansion.
        f = f.q_expansion(terms)
    except:
        pass

    # now use integration by parts. Let f1 be the integral of f
    f1 = f.shift(-1).integral()
    v = list(f1.polynomial())
    K = v[0].parent()
    if K is not QQ:
        w = [phi(t) for t in v]
    else:
        w = v
    q = var('q')
    C = ComplexField(prec)

    ffC = C[q](w)

    endpt = C(exp(2*C(pi)*C(I)*C(b)))
    startpt = C(exp(2*C(pi)*C(I)*C(a)))

    partone = ffC.substitute(endpt)*C(b)**e - ffC.substitute(startpt)*C(a)**e
    parttwo =  e*period(f1,phi,a,b,terms =terms,prec= prec, exponent = e - 1)/C(2*pi*I)

    #verbose('partone = %s'%partone)
    # verbose('parttwo = %s'%parttwo)

    return  partone -parttwo



def _period(f,phi,a,b,terms = None, prec =53,weight = 2):
    """
    return an approximation of the  integral \int_a^b 2\pi i f(z)dz

    when weight = 4, return the integral
    \int_a^b 2 \pi i f(z) z dz.
    """
    if terms is None:
        terms = prec*50 # Linear growth is good, see William Stein: chow heegner paper.
    try:
        # if f is not a q-expansion.
        f = f.q_expansion(terms)
    except:
        pass

    verbose('got here on compuing periods, terms = %s'%f.prec())
    C = ComplexField(prec)


    ff = f.shift(-1).integral()
    v = list(ff.polynomial())
    K = v[0].parent()
    if K is not QQ:
        q = var('q')
        w = [phi(t) for t in v]
        ffpC = C[q](w)
    else:
        q = var('q')
        ffpC = C[q](ff.polynomial())
    endpt = C(exp(2*C(pi)*C(I)*C(b)))
    startpt = C(exp(2*C(pi)*C(I)*C(a)))
    verbose('two end points computed')

    if weight == 2:
        result = ffpC.substitute(endpt) - ffpC.substitute(startpt)
        verbose('substitution made')
        return result

    elif weight == 4:
        result1 = ffpC.substitute(endpt)*C(b) - ffpC.substitute(startpt)*C(a)

        secondPoly = ffpC.shift(-1).integral()
        result2 = secondPoly.substitute(endpt) - secondPoly.substitute(startpt)

        return result1 - result2*(1/C(2*pi*I))

    else:
        raise NotImplementedError

def global_constant(f,phi = None,points = 5,prec = 53,terms = 1000,level = None, weight = 2):
    """
    f -- some newform on Gamma1(N) with coeff in number field K
    phi -- a complex embedding of K
    OUTPUT:
        w(f) -- the atkin-lehner pseudo-eigenvalue of f.
        Or fail.
    """
    if level is None:
        level = f.level()
    if phi is None:
        phi = QQ.embeddings(QQ)[0]
    # print 'f = %s'%f
    pts = generate_points(level,number_of_points = points, prec= prec)
    results = []

    terms = max(terms, prec*50)

    C = ComplexField(prec)

    base_point = C(I/sqrt(level))
    for i in range(len(pts)):
        if weight %2  == 0:
            e = weight//2 - 1
            c = period(f,phi,base_point,pts[i],prec = prec,terms = terms, exponent = e)
        else:
            # if weight is odd. we need to use the other choice of P, which is
            # P(x,y) = -sqrt{N}x^{r+1}y^r + x^ry^{r+1}
            # where r = (k-3)/2
            N = level
            sqrtN = CC(N**(0.5))
            r = (weight - 3)//2
            c1 = period(f,phi,base_point,pts[i],prec = prec,terms = terms, exponent = r+1)
            c2 = period(f,phi,base_point,pts[i],prec = prec,terms = terms, exponent = r)
            c = -sqrtN*c1 + c2

        if abs(c) > 1e-3:
            results.append(c/c.conjugate())
    if len(results) == 0:
        raise ValueError('no nonzero integrals, please give more points.')
    else:
        import numpy as np
        std = np.std(results)
        if std < 1e-3:
            return C(mean(results))
        else:
            raise ValueError('results ( = %s) have a large standard deviation = %s, please give a higher precision.'%(results,std))
    return None