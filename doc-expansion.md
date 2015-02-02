2/1: Doc for the function

    TwistedNewform.expansion(terms)

    Input: terms -- a positive integer
    Output: the q-expansion of f|W_NR_chiW_N, upto q^{terms},
    with coefficients in the double precision complex field

Note that we assumed the `TwistedNewform` class is constructed
with $f, \chi$, where $f$ is a newform on $\Gamma_0(N)$ attached to some elliptic curve defined over $\mathbb{Q}$, and
$\chi$ is a Dirichlet character of prime conductor $p \geq 5$, with $p^2 || N(f)$.

Examples:

1. (not implemented yet)Let $E$ be the elliptic curve '48a' and $\chi_4$ be the unique Dirichlet character of conductor 4. Then we know $f \otimes \chi_4 = g$, where $g$ is the newform associated with the curve `24a', seen as follows:

        sage: E = EllipticCurve('48a');F= EllipticCurve('24a')
        sage: f, g = E.modular_form();F.modular_form()
        sage: load('TwistedNewform.sage')
        sage: chi = DirichletGroup(4).0;
        sage: Tf = TwistedNewform(f,chi)
        sage:...? wrong answer

2. Let $\chi_5$ be the unique character of conductor 5 and order 2. Then the form '50b' is the twist of '50a' by $\chi_5$. Hence we have


        sage: E = EllipticCurve('50a');F= EllipticCurve('50b')
        sage: f, g = E.modular_form(),F.modular_form()
        sage: load('TwistedNewform.sage')
        sage: chi = DirichletGroup(5).0\**2;
        sage: Tf = TwistedNewform(f,chi)
        sage: Tf.newform()[0]
        q + q^2 - q^3 + q^4 + O(q^6)
        sage: g # the same form.
        q + q^2 - q^3 + q^4 + O(q^6)