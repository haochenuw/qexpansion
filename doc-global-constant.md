2/1:

The method.

    global_constant(f,phi)

What it does is computing the pesudo-eigenvalue of the newform f with an embedding $\phi$ from its field of coefficients into complex numbers.

If f is a $\Gamma_0(N)$-newform, then should simply compute numerically the Atkin-Lehner eigenvalue of $W_N$, which is
$\pm 1$. We test it here:

    sage: load('GlobalConstant.sage')
    sage: f = EllipticCurve('37a').modular_form()
    sage: phi = QQ.embeddings(CC)[0]
    sage: global_constant(f,phi)
    1.00000000000000 - 6.79564560804387e-18*I

This confirms that the elliptic curve 37a has odd analytic rank.

Next we test an example where $f$ is a twist of some newform of $\Gamma_0(N)$. Let f be the newform associated with the elliptic curve E = 150a, and $\chi$ be the character of conductor 5 that maps [2] to $\zeta_4$. Then $g = f \otimes \chi \in \mathcal{N}_2(30,\chi^2)$, we have its q-expansion
\[
g = q + a_{0}q^{2} - a_{0}q^{3} - q^{4} + \left(a_{0} - 2\right)q^{5} + O(q^{6})
\]
where $a_0 = \zeta_4$. Since the 5-th coefficient of g is nonzero, we can directly compute the global constant $w(g)$
using [Winnie Li: Theorem 2.1], that
\[
    w(g) = \frac{\mathfrak{g}(\chi^2)}{a_q(g)}
\]
,where $\mathfrak{g}$ denotes the Gauss sum function. Hence we computed

    sage: chi = DirichletGroup(5).0**2
    sage: g = Newforms(chi.extend(30),names = 'a')[0]
    sage: phi = g.base_ring().complex_embeddings()[1]
    sage: chi0.gauss_sum_numerical()/phi(g.qexp(10)[5])
    -0.894427190999916 - 0.447213595499958*I
    sage: global_constant(g,phi) # takes about 10 seconds.
    -0.894427190999916 - 0.447213595499958*I


So our global constant method agrees with the theoretical result.


