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
    return zs


def period(f,phi,a,b,terms = None, prec =53):
    """
    return an approximation of the  integral \int_a^b f
    """
    if terms is None:
        terms = prec*20 # Linear model is good, see William Stein, chow heegner paper.

    f = f.q_expansion(terms)
    C = ComplexField(prec)

    ff = f.shift(-1).integral()
    v = list(ff.polynomial())
    K = v[0].parent()
    if K is not QQ:
        #phi = K.complex_embeddings()[0]
        q = var('q')
        w = [phi(t) for t in v]
        ffpC = C[q](w)
    else:
        verbose('field is rational')
        ffpC = ff.polynomial()

    return ffpC.substitute(exp(2*C(pi)*C(I)*b)) - ffpC.substitute(exp(2*C(pi)*C(I)*a))


def global_constant(f,phi, points = 5,prec = 53):
    """
    f -- some newform on Gamma1(N) with coeff in number field K
    phi -- a complex embedding of K
    OUTPUT:
        w(f) -- the atkin-lehner pseudo-eigenvalue of f.
        Or output fail.
    """
    N = f.level()
    pts = generate_points(N,number_of_points = points)
    results = []

    C = ComplexField(prec)

    base_point = C(I/sqrt(N))
    for i in range(len(pts)):
        c = period(f,phi,base_point,pts[i],prec = prec)
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
            raise ValueError('results have a large standard deviation = %s, please give a higher precision.'%std)
    return None