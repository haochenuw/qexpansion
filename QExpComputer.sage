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
        return [(cchi,) for chi in self.characters()]


    def compute_expansion(self,prec):
        pass

class TwistedNewform(object):

    def __init__(self,f,chi,check =False):
        self.f = f
        self.chi = chi
        N = f.level()
        Q = chi.conductor()
        if N % Q:
            raise ValueError('Q( = %s) has to divide N( = %s)'%(Q,N))
        self.gchi = None
        self.Q = Q
        self.N = N

    def __repr__(self):
        return 'Twist of Newform %s by chi = %s'%(self.f.elliptic_curve().label(),self.chi)

    def newform(self):
        """
        returns a tuple(g,phi)
        where g is a newform of some level corr to f_chi
        and phi is a embedding of K_g into K_chi
        """
        f = self.f
        N = self.N
        Q = self.Q
        if Q == 1:
            print 'twist is trivial'
            return (f,QQ.embeddings(QQ)[0])

        B = self.B_bound()
        verbose('B-bound = %s'%B)

        fqB = list(f.qexp(B+1))
        chisquare = (chi^2).primitive_character()
        Q1 = chisquare.conductor()

        Kchi = chi.base_ring()

        for M in N.divisors():
            if M % Q1 == 0:
                glist = CuspForms(chisquare.extend(M),2).newforms('a')
                for g in glist:
                    Kg = g.base_ring()
                    v = Kg.embeddings(Kchi)
                    if len(v) > 0:
                        gqB  = list(g.qexp(B+1))
                        for phi in v:
                            if all([phi(gqB[n]) == fqB[n]*chi(n) for n in range(1,B+1) if gcd(n,Q) == 1]):
                                self.gchi = g # save the result
                                return (g,phi)

        raise ValueError('newform not found! Please debug.')

    def B_bound(self):
        N = self.f.level()
        k = self.f.weight()
        return 2*CuspForms(N,k).dimension()

    def twist_level(self):
        if self.gchi == None:
            return self.newform()[0].level()
        else:
            return self.gchi.level()

    def Qth_power_of_global_constant(self):
        prec = 100
        Q = self.Q # we are going to raise to Q-th power.
        gchi = self.gchi
        vgchi  = list(gchi.q_expansion(prec+1)**Q)[:prec]
        A,B = all_coeffs_newforms(self.level(),self.weight()*Q,prec)
        C = solve_over_nf(A,vgchi)
        sol = -C.basis()[0][1:]

        print '%s-th power of w is: '%Q
        return solv*matrix(B)[Q+1]

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