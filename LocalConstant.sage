def epsilon_factor(p, order, conjugate = False, multiplicative_order = None, twist =  0):
    """
    p -- a prime.
    order -- the order of a multiplicative character of F_p^2.
    twist_order -- Twist by a Dirichlet character.

    Computes the epsilon factor at p of a supercuspidal representation of conductor p^2.
    twisted by a Dirichlet character.
    """
    d = order
    if d not in [3,4,6]:
        raise ValueError
    if multiplicative_order == None:
        multiplicative_order = p-1
    d1 = multiplicative_order
    K.<zeta> = CyclotomicField(p*d*d1)
    zetap, zetad, zetad1 = zeta**(d*d1), zeta**(p*d1), zeta**(p*d)
    F.<a> = GF(p**2)
    gen = F.multiplicative_generator()

    m = p^2-1
    Fcross = [gen**i for i in range(m)]

    Fp = GF(p)
    genp = Fp.multiplicative_generator()

    # Defining a multiplicative character
    chi = dict([(genp**i, zetad1**(twist*i)) for i in range(p)])
    #verbose('additive twist = %s'%b)

    if conjugate:
        result = sum([ zetap**((Fcross[i] + Fcross[i]**p))*(zetad**(i))*chi[Fp(gen**(i*(p+1)))] for i in range(m)])
    else:
        result = sum([ zetap**((Fcross[i] + Fcross[i]**p))*(zetad**(-i))*chi[Fp(gen**(i*(p+1)))] for i in range(m)])


    try:
        return QQ(-result/p)
    except:
        return -result/p


# Formula for first term.
def first_term(p, order, conjugate = False, multiplicative_order = None, twist =  0):
    """
    p -- a prime.
    order -- the order of a multiplicative character of F_p^2.
    twist_order -- Twist by a Dirichlet character.

    Computes the first term of the expansion of the modular form f attached to E/Q with additive reduction at prime p.
    """
    d = order
    if d not in [3,4,6]:
        raise ValueError
    if multiplicative_order == None:
        multiplicative_order = p-1
    d1 = multiplicative_order
    K.<zeta> = CyclotomicField(p*d*d1)
    zetap, zetad, zetad1 = zeta**(d*d1), zeta**(p*d1), zeta**(p*d)
    F.<a> = GF(p**2)
    gen = F.multiplicative_generator()

    m = p^2-1
    Fcross = [gen**i for i in range(m)]

    Fp = GF(p)
    genp = Fp.multiplicative_generator()

    # Defining a multiplicative character
    chi = dict([(genp**i, zetad1**(twist*i)) for i in range(p)])
    #verbose('additive twist = %s'%b)

    if conjugate:
        denom = sum([ zetap**((Fcross[i] + Fcross[i]**p))*(zetad**(i))*chi[Fp(gen**(i*(p+1)))] for i in range(m)])
        num = sum([ zetap**((Fcross[i] + Fcross[i]**p + Fcross[i]**(p+1)))*zetad**(i) for i in range(m)])
    else:
        denom = sum([ zetap**((Fcross[i] + Fcross[i]**p))*(zetad**(-i))*chi[Fp(gen**(i*(p+1)))] for i in range(m)])
        num = sum([ zetap**((Fcross[i] + Fcross[i]**p + Fcross[i]**(p+1)))*zetad**(-i) for i in range(m)])


    return num/denom

# to-do: fix the character sum so that it works for any character \psi_b.