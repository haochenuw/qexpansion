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

class TwistedNewform(object):

    def __init__(self,f,chi,check =False):
        self.f = f
        self.chi = chi
        N = f.level()
        Q = chi.conductor()
        self.Q = Q
        self.N = N
        self.newformdata = None

    def __repr__(self):
        return 'Twist of Newform %s by chi = %s'%(self.f.elliptic_curve().label(),self.chi)

    def newform(self):
        """
        returns a tuple(g,phi)
        where g is a newform of some level corr to f_chi
        and phi is a embedding of K_g into K_chi
        """
        if self.newformdata is not None:
            return self.newformdata
        f = self.f
        N = self.N
        chi = self.chi
        Q = self.Q
        if Q == 1:
            print 'twist is trivial'
            return (f,QQ.embeddings(QQ)[0])

        B = self.B_bound()
        verbose('B-bound = %s'%B)

        fqB = list(f.qexp(B+10))
        chisquare = (chi^2).primitive_character()
        Q1 = chisquare.conductor()

        Kchi = chi.base_ring()
        for M in N.divisors():
            if M % Q1 == 0:
                verbose('Working on level M = %s'%M)
                glist = CuspForms(chisquare.extend(M),2).newforms('a')
                for g in glist:
                    verbose('Working with one g = %s'%g)
                    Kg = g.base_ring()
                    v = Kg.embeddings(Kchi)
                    if len(v) > 0:
                        gqB  = list(g.qexp(B+10))
                        for phi in v:
                            print 'phi = %s'%phi
                            print 'fqB = %s'
                            if all([phi(gqB[n]) == fqB[n]*chi(n) for n in range(1,B+1) if gcd(n,Q) == 1]):
                                print 'all checked out '
                                self.newformdata = (g,phi) # save the result
                                return (g,phi)

        raise ValueError('newform not found! Please debug.')

    def numerical_global_constant(self):
        chisquare = (self.chi)**2
        if chisquare.conductor() == 1: # the twist is again a form on Gamma0(N)
            return self.f.atkin_lehner_eigenvalue()
        else:
            gchi,phi = self.newform()
            G = Gamma0(gchi.level())
            for gg in G.gens():
                try:
                    w = period_ratio(gchi,phi,gg)
                    print 'gg = %s'%gg
                    return w
                except ValueError, msg:
                    verbose('error = %s'%msg)
                    verbose('gg = %s giving virtually 0 period. '%gg)
        raise ValueError('did not find eigenvalue')


    def B_bound(self):
        N = self.f.level()
        k = self.f.weight()
        return 2*CuspForms(N,k).dimension()


    def expansion_data(self):
        """
        return a FormalSum
        representing the expansion of f|W_NR_\chi(p)W_N.
        (c,d) such that
        f|W_N R_\chi(p) W_N. = w(f)w(g_chi)* \sum c_i\bar(g(q^d_i))
        """
        p = self.Q
        if p == 1:
            return (1,1)
        else:
            g, phi = self.newform()
            m = g.level().valuation(p)
            print 'got here!'
            print 'm = %s'%m
            bp = g.qexp(p+1)[-1]
            if m == 2:
                print 'got here'
                return [(1,1)]
            elif m == 1:
                return [(CC(-bp/p), 1),(p,p)]
            elif m == 0:
                return [(CC(1/p),1),(-bp,p),(p^2,p^2)]

    def expansion_string(self):
        lst = self.expansion_data()
        return FormalSum([(CC(a),'B%s'%b) for a,b in lst],CC)






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


def periodgC(ffpC,gg,terms,prec):
    """
    given a newform f attached to elliptic
    curve E, an element gg in Gamma0(N),
    computes the period
    2*pi*i*\int_{a}^{gg(a)} f(z) dz
    by truncating the q-expansion of f to a polynomial
    with degree = terms, a is the optimal point in [Cremona97].
    """
    C = ComplexField(prec)
    a=gg[0,0]
    b=gg[0,1]
    c=gg[1,0]
    d=gg[1,1]

    if c == 0: # parabolic matrices have zero period, stated in Cremona.
        return 0
    else:
        tau=-(d/c)+(C(I)/abs(c))
        gammatau=(a*tau+b)/(c*tau+d)
        return ffpC.substitute(exp(2*C(pi)*C(I)*gammatau)) - ffpC.substitute(exp(2*C(pi)*C(I)*tau))


def fricke(N):
    return matrix([[0,-1],[N,0]])


def wN(g,N):
    return fricke(N)*(g*~fricke(N))