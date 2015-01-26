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
    def fchiwN(self,chi):
        """
        a simple formula
        """
        gchi, Nchi = self.newform_twist(chi)
        return

