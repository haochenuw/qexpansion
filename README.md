###1. Background.

This project aims at solving the following problem: Let $E$ be an
elliptic curve over $\mathbb{Q}$ and let $f$ be the normalized newform
of weight 2 and level $N$ attached to $E$. Suppose we have access to the $q$-expansion of $f$ at the cusp $\infty$. Then how do we compute the expansion of $f$ at other cusps of $\Gamma_0(N)$?

First, note that the expansion is well-defined only when the cusp satsifies certain properties. (See [Diamond,Shurman pp 17-18]).
We make a definition

**Definition** Let $G = \Gamma_0(N)$. A big cusp of $G$ is a cusp
$z$ such that z is equivalent to a cusp of form $[\frac{c}{d}]$,
with $d \mid N$ and $N \mid d^2$.

**Definition** Let $z$ be a big cusp of $\Gamma_0(N)$, take
any matrix $\alpha \in SL_2(\mathbb{Z})$ with $\alpha(\infty) = z$. The
Fouier expansion of f at $z$ is defined as
$$
    f(q;z) = f|[\alpha](q) = \sum_{n =1}^{\infty} a_n(f;z)q^n \in \mathbb{C}[[q]].
$$
The code in this project will compute the coefficients of $f(q;z)$ numerically.

###2. Preparations.

To use these codes, it's advised that you have a sagemathcloud project.

1. Clone this repository.
(In your project, do the following):

        git clone https://github.com/haochenuw/qexpansion

2. Create a worksheet.

3. Load the main file.

        load('QExpComputer.sage')

###3. Example.

We use the class "QExpComputer" to compute expansions. It takes an elliptic curve and an integer $d'$ such that $d'^2 \mid N$. Then the expansion_numerical() method computes the expansion at the cusp $1/(N/d')$.

1. We start with a sanity-check example: computing the $q$-expansion at $\infty$. In your worksheet, do

        sage: E = EllipticCurve('37a'); Comp = QExpComputer(E,1)
        sage: Comp.expansion_numerical()
        0.000000000000000 + 1.00000000000000*q - 2.00000000000000*q^2 - 3.00000000000000*q^3 + 2.00000000000000*q^4 - 2.00000000000000*q^5 + 6.00000000000000*q^6 - 1.00000000000000*q^7 + 0.000000000000000*q^8 + 6.00000000000000*q^9 + 4.00000000000000*q^10 - 5.00000000000000*q^11 - 6.00000000000000*q^12 - 2.00000000000000*q^13 + 2.00000000000000*q^14 + O(q^15)
        sage: E.modular_form().qexp(15)
        q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + 4*q^10 - 5*q^11 - 6*q^12 - 2*q^13 + 2*q^14 + O(q^15)

    We see these two expansions match up. Note that the default option computes up to $q^{15}$, but one can change the keyword 'terms' to get arbitrary number of terms.

2. As a second exmple, we compute with the curve '48a' and cusp $z = 1/12$, the smallest power of $q$ in this  expansion is known to be $q^2$.

        sage: E = EllipticCurve('48a);
        sage: Comp = QExpComputer(E,4)
        sage: Comp.expansion_numerical(terms=20)
        -0.000000000000000 - 0.000000000000000*q + (-2.44929359829471e-16 - 2.00000000000000*I)*q^2 - 0.000000000000000*q^3 - 0.000000000000000*q^4 - 0.000000000000000*q^5 + (2.44929359829471e-16 + 2.00000000000000*I)*q^6 - 0.000000000000000*q^7 - 0.000000000000000*q^8 - 0.000000000000000*q^9 + (4.89858719658941e-16 + 4.00000000000000*I)*q^10 - 0.000000000000000*q^11 - 0.000000000000000*q^12 - 0.000000000000000*q^13 - 0.000000000000000*q^14 - 0.000000000000000*q^15 - 0.000000000000000*q^16 - 0.000000000000000*q^17 + (-2.44929359829471e-16 - 2.00000000000000*I)*q^18 + O(q^20)

    so we have $f(q;[\frac{1}{12}]) = -2iq^2 + O(q^3)$.





