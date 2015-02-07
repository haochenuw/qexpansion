load('PeriodMapping.sage')

def generate_points(N,number_of_points = 10):
    """
    return some random points on the arc |z| = 1/sqrt(N), im(z) > 0.
    """
    try:
        N = ZZ(N)
    except Exception:
        print N.parent()

    thetas = [CC(0.4*random() + 0.1) for _ in range(number_of_points)]
    verbose('thetas computed')
    print N.parent()
    zs = [CC(1/sqrt(N))*CC(exp(CC(pi)*CC(I)*theta)) for theta in thetas]
    return zs


def period(f,phi,a,b,terms = 1000,prec =53):
    """
    return an approximation of the  integral \int_a^b f
    """
    f = f.q_expansion(terms)
    C = ComplexField(prec)

    ff = f.shift(-1).integral()
    v = list(ff.polynomial())
    K = v[0].parent()
    verbose('K = %s'%K)
    verbose('heads of v = %s'%v[:5])
    if K is not QQ:
        #phi = K.complex_embeddings()[0]
        q = var('q')
        w = [phi(t) for t in v]
        ffpC = C[q](w)
    else:
        verbose('field is rational')
        ffpC = ff.polynomial()

    return ffpC.substitute(exp(2*C(pi)*C(I)*b)) - ffpC.substitute(exp(2*C(pi)*C(I)*a))


def global_constant(f,phi, points =10):
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
    base_point = CC(I/sqrt(N))
    for i in range(len(pts)):
        c = period(f,phi,base_point,pts[i])
        if abs(c) > 1e-3:
            results.append(c/c.conjugate())
    if len(results) == 0:
        raise ValueError('no nonzero integrals, please give more points.')
    else:
        import numpy as np
        std = np.std(results)
        if std < 1e-3:
            return CC(mean(results))
        else:
            raise ValueError('results have a large standard deviation = %s, please give a higher precision.'%std)
    return None