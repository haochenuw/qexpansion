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


def period(f,phi,a,b,terms = None, prec =53):
    """
    return an approximation of the  integral \int_a^b f
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

    result = ffpC.substitute(endpt) - ffpC.substitute(startpt)

    verbose('substitution made')

    return result

def global_constant(f,phi,points = 5,prec = 53,terms = 1000,level = None):
    """
    f -- some newform on Gamma1(N) with coeff in number field K
    phi -- a complex embedding of K
    OUTPUT:
        w(f) -- the atkin-lehner pseudo-eigenvalue of f.
        Or fail.
    """
    if level is None:
        level = f.level()
    # print 'f = %s'%f
    pts = generate_points(level,number_of_points = points, prec= prec)
    results = []

    terms = max(terms, prec*50)

    C = ComplexField(prec)

    base_point = C(I/sqrt(level))
    for i in range(len(pts)):
        c = period(f,phi,base_point,pts[i],prec = prec,terms = terms)
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