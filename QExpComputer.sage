load('TwistedNewform.sage')



class QExpComputer(object):

    def __init__(self,E,d):
        self.E = E
        self.d = d
        self.N = E.conductor()
        if self.N % (d^2):
            raise ValueError('the second argument squared must divide the conductor of E.')
        self.f = E.modular_form()

    def __repr__(self):
        return 'Computer for the q-expansion of the newform attached to the elliptic cuvrve %s at the cusp 1/%s'%(self.curve().label(), self.denom())

    def curve(self):
        return self.E

    def denom(self):
        return self.N // self.d

    def characters(self):
        """
        return the set of Dirichlet characters of modulus d with
        cond'(chi) == d in the sense of Winnie Li's twist paper.
        """
        return [chi for chi in list(DirichletGroup(self.d)) if conductor_prime(chi) == self.d]


    def expansion_numerical(self,numerator = 1,terms = 15,prec = 53):
        """
        return the first (terms) terms of the q-expansion of self.f at the cusp
        c = numerator/self.denom(). The default numerator is 1.
        """
        f = self.f
        d = self.d
        if gcd(numerator,self.denom()) > 1:
            raise ValueError('numerator must be coprime to the cusp denom.')
        try:
            numerator = ZZ(numerator)
        except:
            raise ValueError('numerator must be an integer.')
        a = numerator.inverse_mod(self.denom())

        C = ComplexField(prec)
        q = var('q')

        result = CC[[q]](0)
        for chi in self.characters():
            verbose('chi = %s'%chi)
            Tchi = TwistedNewform(f,chi)
            exp_chi = Tchi.expansion(d,terms)
            verbose('exp_chi = %s'%exp_chi)
            result += exp_chi*CC(chi(-a))
        return result/euler_phi(d)





"""

def all_coeffs_newforms(level,weight,prec):
    #return a matrix of coefficients of cusp forms in the canonical basis of S_k(N)
    #up to prec, and(!) the corresponding matrix of the same forms under under w_N.
    #See the documentation for more details
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
"""