class QExpComputer(object):
    def __init__(self,E,p):
        #if not isinstance(E,EllipticCurve):
        #    raise ValueError('The first argument must be an elliptic curve')
        self.E = E
        self.p = p
        self.N = E.conductor()
        if not is_prime_power(p):
            raise ValueError('the second argument must be a prime power.')
        elif self.N % p:
            raise ValueError('the second argument must divide the conductor of E.')

        self.f = E.modular_form()

    def __repr__(self):
        return 'Computer for the q-expansion of the newform attached to the elliptic cuvrve %s at the cusp 1/%s'%(self.curve().label(), self.denom())

    def curve(self):
        return self.E

    def denom(self):
        return self.N // self.p

    def characters(self):
        """
        return the set of Dirichlet characters of modulus N with conductor dividing p,
        which is relevanet of the computation of expansions.
        """
        return list(DirichletGroup(self.p))

    def newform_twist(self,chi):
        """
        some work should be put here
        """

    def pseudo_eigenvalue_numeric(self,chi):
        pass


    def fchiwN(self,chi):
        """
        a simple formula
        """
        gchi, Nchi = self.newform_twist(chi)
        return

    def expansion_data(self):
        result =[]
        for chi in self.characters():
            if chi.conductor() == p: # non-trivial conductor
                Tchi = TwistedNewform(self.f,chi)
                result.append((chi, Tchi.constant(),Tchi.newform()))

    def expansion_numerical(self,n,prec = 100):
        # return the n-th coefficient of the expansion of self.f at 1/self.denom()
        # only implemented when self.p = p is a prime and v_p(N) == 2.
        f = self.f
        wf = f.atkin_lehner_eigenvalue()
        p = self.p
        N = self.N
        an = self.E.an(n)

        C = ComplexField(prec)
        if (not p.is_prime() or N.valuation(p) != 2):
            raise NotImplementedError

        result = 0
        for chi, wchi, gchi, phi in expansion_data(self):
            gausschi = chi.bar().gauss_sum().complex_embeddings()[0]
            Nchi = gchi.level()
            verbose('Nchi = %s'%Nchi)
            m = Nchi.valuation(p)
            if m == 0:
                result += an*gausschi*wchi*C(chi.bar()(n))
            elif m == 1:
                bp = gchi.qexp(p+1)[-1]
                result +=  (gausschi*wchi*C(chi.bar()(n)))*(-bp/p)
                if n % p == 0:
                    result += p*self.E.an(n//p)*wchi*C(chi.bar()(n//p))

            elif m == 2:
                result += gausschi*wf*wchi*chi(n)*f.an(n)
            else:
                raise ValueError('the valuation at p is wrong. Pleaes debug.')
        # finally we add the last term -f
        return result - an


    def compute_expansion(self,prec):
        pass








def all_coeffs_newforms(level,weight,prec):
    """
    return a matrix of coefficients of cusp forms in the canonical basis of S_k(N)
    up to prec, and(!) the corresponding matrix of the same forms under under w_N.
    See the documentation for more details
    """
    indexes = []
    result = []
    result_wNtwisted = []
    N = level
    k = weight
    v = N.divisors()
    for M in v:
        for d in (ZZ(N//M)).divisors():
            newformlist = CuspForms(M,k).newforms('a')
            for newform in newformlist:
                result.append(expansion(newform,prec,d,N))
                twist_const = QQ(newform.atkin_lehner_eigenvalue()*((QQ(N/(d**2*M)))**2)) # !
                print newform, M,d,twist_const
                result_wNtwisted.append([a*twist_const for a in expansion(newform,prec,N//(d*M),N)])
    return (result, result_wNtwisted)


def expansion(f,prec,d,N):
    """
    return the first prec terms of f|B_d = f(q^d)
    where f is newform of level N' and N'd divides N.
    """
    Nprime = f.level()
    assert N % (Nprime*d) == 0

    fq = f.q_expansion(prec+1)
    v = list(fq)
    verbose('len(v) = %s'%len(v))
    result = [0 for _ in range(prec)]
    for i in range(prec):
        if i % d == 0:
            result[i] = v[i//d]
    return result

def solve_over_nf(A,v):
    K = v[0].parent()
    print 'K = %s'%K
    Anp = np.array(A)
    vnp = np.array([v])
    B1 = np.concatenate((vnp,Anp),axis = 0)
    return matrix(B1).left_kernel()


def period_ratio(f,phi,gg,terms = 2000,prec = 100):
    C = ComplexField(prec)
    N = f.level()
    wNgg = wN(gg.matrix(),N)
    f = f.q_expansion(terms)
    ff = f.shift(-1).integral()
    v = list(ff.polynomial())
    K = v[0].parent()
    verbose('field = %s'%K)
    if K is not QQ:
        # phi = K.complex_embeddings()[0]
        q = var('q')
        ffpC = C[q]([phi(a) for a in v])
        ffpCbar = C[q]([phi(a).conjugate() for a in v])
    else:
        q = var('q')
        ffpC = C[q](v)
        ffpCbar = ffpC

    r1 = periodgC(ffpC,wNgg,terms,prec)

    if abs(r1) < 1e-6:
        raise ValueError('denominator of ratio too small = %s'%abs(r1))
    r2 = periodgC(ffpCbar,gg,terms,prec)
    print 'r1 = %s, r2 = %s'%(r1,r2)
    if abs(r2/r1-1) < 1e-3:
        return r2/r1
    else:
        raise ValueError('Too far away from 1 = %s'%abs(r2/r1-1))




def fricke(N):
    return matrix([[0,-1],[N,0]])


def wN(g,N):
    return fricke(N)*(g*~fricke(N))