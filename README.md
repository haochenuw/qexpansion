### 1. Background.

This project aims at solving the following problem: Let $E$ be an
elliptic curve defined over $\mathbb{Q}$ and let $f$ be the normalized newform of weight 2 and level $N$ attached to $E$. Suppose we have access to the $q$-expansion of $f$ at the cusp $\infty$.

> Problem: compute the expansion of $f$ at __all__ cusps of $\Gamma_0(N)$.

First, note that the expansion is well-defined only when the cusp satsifies certain properties. (See [Diamond,Shurman pp 17-18]).
We make a definition

**Definition** Let $G = \Gamma_0(N)$. A _big cusp_ of $G$ is a cusp
$z$ such that $z$ is equivalent to a cusp of form $[\frac{c}{d}]$,
with $d \mid N$ and $N \mid d^2$.

**Definition** Let $z$ be a big cusp of $\Gamma_0(N)$, take
any matrix $\alpha \in SL_2(\mathbb{Z})$ with $\alpha(\infty) = z$. The _Fourier expansion_ of f at $z$ is
$$
    f(q;z) = f|[\alpha](q) = \sum_{n =1}^{\infty} a_n(f;z)q^n \in \mathbb{C}[[q]].
$$
The Fourier expansion is well-defined up to Gal(Qbar/Q) conjugates.

In this project, we will compute the coefficients of $f(q;z)$ both numerically and exactly.

### 2. Cloning.

To use these code, it's advised that you have a Sagemathcloud account.

1. Clone this repository.
(In your project, do the following):

        git clone https://github.com/haochenuw/qexpansion

2. Create a worksheet.

3. Load the main file by doing

        load('QExpComputer.sage')

### 3. Usages and Examples.

We use the class "QExpComputer" to compute expansions. It takes an elliptic curve and an integer $d'$ such that $d'^2 \mid N$. Then the expansion_numerical() method computes the expansion at the cusp $1/(N/d')$.

1. We start with a sanity-check example: computing the $q$-expansion at $\infty$. In your worksheet, do

        sage: f = EllipticCurve('37a').modular_form();
        sage: Comp = QExpComputer(f,1)
        sage: Comp.expansion_numerical()
        0.000000000000000 + 1.00000000000000*q - 2.00000000000000*q^2 - 3.00000000000000*q^3 + 2.00000000000000*q^4 - 2.00000000000000*q^5 + 6.00000000000000*q^6 - 1.00000000000000*q^7 + 0.000000000000000*q^8 + 6.00000000000000*q^9 + 4.00000000000000*q^10 - 5.00000000000000*q^11 - 6.00000000000000*q^12 - 2.00000000000000*q^13 + 2.00000000000000*q^14 + O(q^15)
        sage: f.qexp(15)
        q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + 4*q^10 - 5*q^11 - 6*q^12 - 2*q^13 + 2*q^14 + O(q^15)

    We see these two expansions match up. Note that the default option computes up to $q^{15}$, but one can change the keyword 'terms' to get arbitrary number of terms.

2. As a second exmple, we compute with the curve '48a' and cusp $z = 1/12$, the smallest power of $q$ in this  expansion is known to be $q^2$.

        sage: f = EllipticCurve('48a').modular_form()
        sage: Q = QExpComputer(f,4)
        sage: Q.expansion_numerical()
        -0.000000000000000 - 0.000000000000000*q + (-2.44929359829471e-16 - 2.00000000000000*I)*q^2 - 0.000000000000000*q^3 - 0.000000000000000*q^4 - 0.000000000000000*q^5 + (2.44929359829471e-16 + 2.00000000000000*I)*q^6 - 0.000000000000000*q^7 - 0.000000000000000*q^8 - 0.000000000000000*q^9 + (4.89858719658941e-16 + 4.00000000000000*I)*q^10 + O(q^15)

    so we have $f(q;[\frac{1}{12}]) = -2iq^2 + O(q^3)$.

3. As a non-trivial example, we compute the expansion at a cusp of denominator 7 on $X_0(49)$.

        sage: f = EllipticCurve('49a').modular_form()
        sage: Q = QExpComputer(f,7)
        sage: Q.expansion_numerical()
        0.000000000000000 + (0.623489801858731 - 1.29468991410431*I)*q + (-0.222520933956313 + 0.177454523299420*I)*q^2 + 0.000000000000000*q^3 + (0.900968867902419 + 0.205640264727405*I)*q^4 + 0.000000000000000*q^5 + 0.000000000000000*q^6 + 0.000000000000000*q^7 + (-1.87046940557619 + 3.88406974231292*I)*q^8 + (0.667562801868939 - 0.532363569898259*I)*q^9 + 0.000000000000000*q^10 + (-3.60387547160967 - 0.822561058909619*I)*q^11 + O(q^15)


