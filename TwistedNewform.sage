load('GlobalConstant.sage')

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

    def character(self):
        """
        return the character of the twisted form. Note that we assumed f has trivial
        character. Then the result should be chi^2
        """
        return self.chi**2

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

        chi = chi.primitive_character()
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
                            verbose('phi = %s'%phi)
                            if all([phi(gqB[n]) == fqB[n]*chi(n) for n in range(1,B+1) if gcd(n,Q) == 1]):
                                verbose('all checked out')
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
            bp = g.qexp(p+1)[-1]
            if m == 2:
                return [(1,1)]
            elif m == 1:
                return [(CC(-bp/p), 1),(p,p)]
            elif m == 0:
                return [(CC(1/p),1),(-bp,p),(p^2,p^2)]

    def expansion_string(self):
        lst = self.expansion_data()
        return FormalSum([(CC(a),'B%s'%b) for a,b in lst],CC)

    def al_eigenvalue(self):
        """
        return the atkin-lehner eigenvalue of the original(untwisted) form.
        """
        return self.f.atkin_lehner_eigenvalue()

    def pseudo_eigenvalue(self):
        g,phi = self.newform()
        return global_constant(g,phi)

    def expansion(self,terms):
        """
        returns the numerical q-expansion of f|W_NR_chiW_N
        """
        chi = self.chi.primitive_character()
        chibar = chi.bar()
        gauss = chibar.gauss_sum_numerical()
        wf = self.al_eigenvalue()
        wg = self.pseudo_eigenvalue()
        const = wf*wg*gauss
        print 'const = %s'%const
        result = []
        v = self.expansion_data()
        max_degree = max([d for c,d in v])
        print 'max degree = %s'%max_degree
        g,phi = self.newformdata
        gg = g.qexp(terms*max_degree)
        vggp = list(gg.polynomial())
        vggpC = [phi(a).conjugate() for a in vggp] # Don't forget to take complex conjugate!
        result = [0 for _ in range(terms)]
        for c, d in v:
            newcoeffs = vggpC[::ZZ(d)][:terms]
            result = [a+c*b for a,b in zip(result,newcoeffs)]
        q = var('q')
        return CC[[q]]([const*t for t in vggpC]).add_bigoh(terms)




