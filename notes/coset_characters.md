# Coset Character Calculation

We have an expression for evaluating the characters of affine Lie Algebras. What would be nice would be to obtain the characters for the quotient of Lie algebras. In some sense if we have two simple Lie algebras $\mathfrak{g}, \mathfrak{h}$ we want to calculate the characters of $\hat{\mathfrak{g}}/\hat{\mathfrak{h}}$. The main tool for doing this would be to calculate the *branching* functions. 

Let's start by studying the characters of $\mathfrak{su}(2)$ in more detail. 

**<u>Proposition:</u>** The normalized reduced character of an $\hat{\mathfrak{su}}(2)$ representation with weight $\hat l = l + k\omega$ is given by
$$
\chi_{\hat l}(q) = \sum_{m = -k + 1}^{k} \sigma_{l + m\omega}(q) q^{-\frac{m^2}{4k}} \Theta_{m + k\omega}(q) = \sum_{m= - k+1}^k c_{m}^l(q) \Theta_{m}^k(q).
$$
where $c^{l}_m$ is the normalized string function. 

Notice that these generalized Theta functions have very cute expansions
$$
\Theta_{m}^{k}(z,\tau) = \sum_{n \in \mathbb{Z}}x^{k\left(n+\frac{m}{2k}\right)}q^{k\left( n + \frac{m}{2k} \right)^2},
$$
where $x=e^{-2\pi i z}, q =e^{2\pi i \tau}$ which is quote nice. Now we want to find an identity of string function products. Let's start with an identity of the $\Theta$'s. We do this by 
$$
\Theta_{n}^{k}(z,\tau) \Theta_{m}^l(z,\tau) = \sum_{r,s \in \mathbb{Z}} x^{k\left(r+\frac{n}{2k}\right) + l\left(s+\frac{m}{2l}\right)}q^{k\left( r + \frac{n}{2k} \right)^2 + l\left( s + \frac{m}{2l} \right)^2}
$$
Let's work this out. In particular consider the power of $q$ since this is more restrictive

