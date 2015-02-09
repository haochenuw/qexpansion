class KnOperator(object):
    """
    The Kn operator acting on a newform
    Kn = \prod p\mid n K_p, where

    K_p = Id - U_pB_p.
    """
    def __init__(self,f,n):
        try:
            n = ZZ(n)
        except:
            raise ValueError('n must be an integer')
        if not n.is_squarefree():
            raise ValueError('n must be square free')
        self.f = f
        self.n = n

    def __repr__(self):
        return 'The K%s operator acting on the newform f = %s'%(self.n,self.f)

    def character(self):
        return self.f.character()


    def kn_polynomial(self):
        try:
            K = self.f.hecke_eigenvalue_field()
        except:
            K = self.f.base_ring()
        chi = self.character()
        n = self.n
        vn = n.prime_divisors()
        R = PolynomialRing(K,['B%s'%p for p in vn])
        variables = R.gens()
        verbose('R = %s'%R)
        vf = self.f.padded_list(n+1)
        result = R(1)
        for i in range(len(vn)):
            factor =  (1 - vf[vn[i]]*R(variables[i]) + chi(p)*p*R(variables[i])**2)
            result *= factor
        return result

    def expansion_data(self):
        P = self.kn_polynomial()
        vn = (self.n).prime_factors()
        result = []
        for a in P.monomials():
            coeff = P.monomial_coefficient(a)
            verbose('a = %s'%a)
            result.append((coeff,a(vn)))
        return result
