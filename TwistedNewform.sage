load('GlobalConstant.sage')
load('Decompose-Characters.sage')
load('KnOperator.sage')

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
        return (self.chi**2).primitive_character()

    def is_minimal(self):
        """
        self is associated to an elliptic curve.
        """
        Q = self.Q
        if not is_prime(Q):
            raise ValueError('conductor of the twisting character must be prime.')
        from sage.modular.local_comp.type_space import TypeSpace
        from sage.modular.modform.element import Newform
        E = self.f.elliptic_curve()
        fnew = Newform(CuspForms(E.conductor()), E.modular_symbol_space(),'a')
        T = TypeSpace(fnew,Q)
        return T.is_minimal()


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
            verbose('twist is trivial')
            self.newformdata = (f,QQ.embeddings(QQ)[0])
            return self.newformdata

        B = self.B_bound()
        verbose('B-bound = %s'%B)

        fqB = list(f.qexp(B+10))
        chisquare = (chi^2).primitive_character()
        Q1 = chisquare.conductor()

        chi = chi.primitive_character()
        Kchi = chi.base_ring()
        name = 'chi'
        for a in chi.element():
            name += str(a)

        # the coprime part
        Nprime = N//gcd(N,Q)
        print 'Nprime = %s'%Nprime

        for M in N.divisors():
            if M % (lcm(Q1,Nprime)) == 0:
                verbose('Working on level M = %s'%M)
                glist = CuspForms(chisquare.extend(M),2).newforms('a'+name+'M%s'%M)
                verbose('Number of conjugacy classes of newforms on level %s is %s'%(M,len(list(glist))))
                for g in glist:
                    verbose('Working with the newform g = %s'%g)
                    Kg = g.base_ring()
                    if Kg is not QQ:
                        L = Kg.absolute_field('c')
                        from_L, to_L = L.structure()
                    else:
                        L = Kg
                        to_L = Kg.embeddings(QQ)[0]
                    v = L.embeddings(Kchi)

                    if len(v) > 0:
                        gqB  = list(g.qexp(B+10))
                        for phi in v:
                            verbose('phi = %s'%phi)
                            if all([phi(to_L(gqB[n])) == fqB[n]*chi(n) for n in range(1,B+1) if gcd(n,Q) == 1]):
                                verbose('found the twisted newform.')
                                self.newformdata = (g,phi) # save the result
                                return (g,phi)

        raise ValueError('newform not found for chi = %s! Please debug.'%self.chi)

    def numerical_global_constant(self):
        # deprecated method. should be deleted.
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

    def _expansion_data(self):
        """
        Warning: deprecated. Do not need this anymore.
        The helper function used in expansion_data
        returns a list of tuples (c,d)
        such that f_chi = \sum c_ig(q^d_i), where g is the newform
        associated to f_chi.
        """
        Q = self.chi.conductor()
        if Q == 1:
            return [(1,1)]

        g,_ = self.newform()

        v = Q.prime_divisors()

        # check that the conductor is a prime power.
        if len(v) > 1:
            raise NotImplementedError
        p = v[0]

        m = g.level().valuation(p)
        n = self.f.level().valuation(p)

        # simplest case: twist is a newform.
        if m == n:
            return [(1,1)]


        else:
            # second case where twist has non-trivial p-level
            # f_chi = g - ap(g)*g|B_p
            bp = g.coefficients([p])[0]
            if m > 0:
                return [(1,1),(-bp,p)]

            # Last case, where twist has level not-divisible by p.
            # f_chi = g + a_p(g)g|B_p + chi^2(p)pg|B_p^2
            else: # m = 0
                chisquare = self.character()
                if chisquare.conductor() > 1:
                    raise ValueError('This is impossible. Please debug.')
                else:
                    return [(1,1),(-bp,p),(p,p^2)]

    def _delaunay_factors(self,M,prec = 53):
        """
        return chi', c, such that f|R_chi(M) = c f_chi'.
        We require that cond'(chi) = M.
        Use Delaunay's formula on page 74.
        """

        chi = self.chi
        Mchi = chi.modulus()
        if M % Mchi:
            raise ValueError('Mchi must divide M.')
        if conductor_prime(chi) != M:
            raise ValueError('cond \' (chi) must equal M')

        C = ComplexField(prec)
        chi = chi.extend(M)
        chint = factorization(chi) # loaded function from DC.sage

        # use chint to obtain info on tr, and ...
        nt = chint.modulus()
        tr = ZZ(M // nt)
        i = len(tr.prime_divisors())
        if nt > 1: # if chint is not the trivial character
            return (chint, C((-1)**i*chint(tr)*chint.bar().gauss_sum_numerical()))
        else: # if chint is the trivial character. We need to separate out this case
        # because the gauss sum for trivial character modulo 1 in sage is 0 (?) but I want it to be 1.
            return (chint, (-1)**i)


    def _expansion_data_comp(self):
        """
        Replaces _expansion_data. Work for composite modulus.
        Let g be the twisted newform associated to self.
        Returns a list of tuples (c_j,d_j) such that
        f_chi = f_chint = \sum_j c_j g(q^d_j) = g|K_n.
        """
        Kn = self._kn()
        return Kn.exp_data()

    def _kn(self):
        chint,_ = self._delaunay_factors(self.chi.modulus())
        n = chint.modulus().radical()
        g,_ = self.newform()
        return KnOperator(g,n)

    def expansion_data_comp(self):
        """
        Returns a list of tuples (c,d) such that
        representing the expansion of f|W_NR_\chi(M)W_N.

        f|W_N R_\chi(M) W_N.
                             = A * f_chi |W_N
                             = A * \sum (c_i * \bar(g(q^d_i)))

        """
        N = self.f.level()
        Ng = self.newform()[0].level()
        return [(coeff*QQ(N/(Ng*degree**2)),ZZ(N/(degree*Ng))) for coeff,degree in self._expansion_data_comp()]

    def constant(self,M,prec = 53):
        """
        returns the constant A mentioned in expansion_data_comp
        """
        # return w(f)w(g)*(the second output of delaunay_factors)
        wf = self.al_eigenvalue()
        wg = self.pseudo_eigenvalue(prec = prec)
        alconst = wf*wg
        verbose('alconst = %s'%alconst)
        delaunayconst = self._delaunay_factors(M,prec = prec)[1]
        return alconst*delaunayconst

    def _const(self,M):
        """
        deprecated.
        Being replaced by _delaunay_factors.
        return the constant c_{chi,M} such that
        f|R_chi(M) = c_{chi,M} f_chi.
        where M = some prime power divisible by cond(chi). We assume M > 1.
        """
        chi = self.chi
        Mchi = chi.modulus()

        if M % Mchi:
            raise ValueError('M must be divisible by the modulus of chi')

        chi = chi.extend(M)
        Q = chi.conductor()

        vchi = chi.decomposition()
        result = 1

        vchi_nondegenerate = []
        # first we deal with some degeneracies
        for chip in vchi:
            mod = chip.modulus()
            cond = chip.conductor()
            if 1 < cond and cond < mod:
                return 0
            elif cond == 1:
                if not mod.is_prime():
                    return 0
                else:
                    result *= -1
            else: # conductor = modulus
                result *= chip.bar().gauss_sum()
        # Delaunay's other constants should also be taken into account!!!
        # Why don't I use Delaunay's formula??
        v = M.prime_divisors()



        # check that M is a prime power.
        if len(v) > 1:
            raise NotImplementedError

        p = v[0]
        d = M.valuation(p)

        if Q == 1: # chi is trivial.
            if d == 1:
                return -1
            else: return 0
        elif Q < M:
                return 0 # use the fact that ord_p(N_f) >= 2 => a(pn) = 0 for all n.
        else:
            # the standard cse where M = cond(chi).
            # Then f|R_chi(M) = g(\chibar)f_chi
            return self.chi.bar().gauss_sum()



    def expansion_data(self,M):
        """
        WARNING: DEPRECATED. USE EXPANSION_DATA_COMP INSTEAD.
        representing the expansion of f|W_NR_\chi(M)W_N.
        (c,d) such that
        f|W_N R_\chi(M) W_N. = w(f)w(g_chi)* \sum (c_i * \bar(g(q^d_i)))

        New content: now I modify it to work for prime power.
        i.e. the expansion data for
            f|W_NR_\chi(M)W_N
        where M is a prime power.

        """
        print 'M = %s'%M
        v = M.prime_divisors()
        if len(v) > 1:
            raise NotImplementedError

        p = v[0]
        Q = self.chi.conductor() # conductor of chi.
        if M % Q:
            raise ValueError('Conductor of chi must divide M')

        N = self.N
        c = self._const(M)

        Ng = self.newform()[0].level()
        return [(c*coeff*QQ(N/(Ng*degree**2)),ZZ(N/(degree*Ng))) for coeff,degree in self._expansion_data()]


        """
        # First we use delaunay's formula to deal with the degenerate case
        # where condchi is trivial.
        if condchi == 1:
            if M.valuation(p) >= 2:
                return [(0,1)] # 0
            else:
                return [(-1,1)] # -f

        # Now we deal with the second degenerate case, where
        # condchi = p^a, M = p^b, and a < b.
        elif condchi < M:
            return [(0,1)] # zero. Note I can only show this works for Gamma0(N) newforms.
            # However, by commutativity, I can move a non-primitive sum to front and
            # always get zero.

        # The standard case, where M = cond(chi).
        else:
            g, phi = self.newform()
            m = g.level().valuation(p)
            bp = g.coefficients([p])[0]
            if m == 2:
                return [(1,1)]
            elif m == 1:
                return [(CC(-bp/p), 1),(p,p)]
            elif m == 0:
                return [(CC(1/p),1),(-bp,p),(p^2,p^2)]
        """


    def al_eigenvalue(self):
        """
        return the atkin-lehner eigenvalue of the original(untwisted) form.
        """
        return self.f.atkin_lehner_eigenvalue()

    def pseudo_eigenvalue(self,prec = 53):
        g,phi = self.newform()
        #First we deal with some easy cases:
        if g.character().conductor() == 1: # twist is on Gamma0(N).
            return g.atkin_lehner_eigenvalue()

        return global_constant(g,phi,prec = prec)

    def expansion(self,M,terms = 15,prec = 53):
        """
        returns the numerical q-expansion of f|W_NR_chi(M)W_N.
        """
        chi = self.chi


        # Degenerate case first.
        if M != conductor_prime(chi):
            verbose('in the zero case')
            return 0

        C = ComplexField(prec)

        # computing the constants in front,
        const = self.constant(M)

        # computing the q-expansion of gbar
        g,phi = self.newform()
        gg = g.qexp(terms+10)
        q = gg.parent().gen()
        #vggp = list(gg.polynomial())
        #vggpC = [phi(a).conjugate() for a in vggp] # Don't forget to take complex conjugate!

        # add up the result. V contains the constants plus the B_p information.
        v = self.expansion_data_comp()

        verbose('v = %s'%v)



        result = [0 for _ in range(terms)]
        for coeff, d in v:
            ggd = gg(q = q^d)
            vggp = list(ggd.polynomial())
            vggpC = [phi(a).conjugate() for a in vggp] # Don't forget to take complex conjugate!
            newcoeffs = vggpC[:terms]
            result = [a+C(coeff)*b for a,b in zip(result,newcoeffs)]
        q = var('q')
        return C(const)*C[[q]](result).add_bigoh(terms)




