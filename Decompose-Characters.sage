def conductor_prime(chi):
    """
    return the 'cond(chi)' defined by Li and
    also seen in Delaunay's thesis of a character chi.
    """
    vchi = chi.decomposition()
    result = 1
    for chip in vchi:
        cond = chip.conductor()
        if cond > 1:
            result *= cond
        else:
            result*= chip.modulus().radical() # must be a prime
    return result


def factorization(chi):
    """
    factors a Dirichlet character chi into trivial and nontrivial part and return the
    non-trivial part.
    For the notations, see Delaunay's thesis, page 73-74.
    """
    verbose('chi being factored = %s'%chi)
    vchi = chi.decomposition()
    verbose('vhi = %s'%vchi)
    tr = 1
    vtr = []
    vnt = []
    for chip in vchi:
        cond = chip.conductor()
        if cond == 1:
            tr*= chip.modulus()
            vtr.append(chip)
        else:
            vnt.append(chip)

    nt = ZZ(chi.modulus()//tr)

    verbose('got here, tr = %s, nt = %s'%(tr,nt))
    # Now we compute the characters.
    if tr == 1:
        chitr = list(DirichletGroup(1))[0]
    else:
        chitr = prod([sigma.extend(tr) for sigma in vtr])
    if nt == 1:
        chint = list(DirichletGroup(1))[0]
    else:
        chint = prod([tau.extend(nt) for tau in vnt])

    return chint
