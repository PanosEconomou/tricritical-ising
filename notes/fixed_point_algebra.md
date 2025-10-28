# Fixed Point algebra

Here we calculate the Fixed point algebra $\mathcal{W}$ of $\mathcal{A} = \text{Vir}\oplus \text{Vir}$ by the $\mathbb{Z}_2^2$ generated through exchange $S$ and the diagonal $\eta$ action $H = \eta\boxtimes \eta$. We do this by calculating the decomposition of the characters and then eventually putting everything together to find the S-matrix.

[toc]

# Character Decompositions

We will find the irreducible modules of the fixed point algebra $\mathcal{W}$ by decomposing the irreps of $\mathcal{A}$ into eigenspaces by $S,H$ and then calculate the characters of these modules by identities of Virasoro characters.

**<u>Lemma:</u>** Each irreducible highest weight unitary representation $A$ of $\mathcal{A}$ can be written as a product of two unique irreducible highest weight representations $i,j$ of $\text{Vir}$
$$
A = A_{i,j} \coloneqq v=i\otimes j.
$$
We will use this a lot in our decompositions.

**<u>Corollary:</u>** The character of the $\mathcal{A}$-representation $A_{ij}$ is given by
$$
\chi_{ij}^{\mathcal{A}}(q) = \chi_i(q) \chi_{j}(q),
$$
where $\chi(q)$ denote the $\text{Vir}$ characters. 

## Fixed point under $S$

We will calculate the fixed point algebra in two steps. The first one is to calculate the intermediate algebra $\mathcal{W}\subset \mathcal{A}^S \subset \mathcal{A}$ where $\mathcal{A}^S$ is the fixed point algebra under $S$.

**<u>Definition:</u>** Given a representation $A_{ij}$ of $\mathcal{A}$, the **exchange map** $S: A_{ij} \to A_{ji}$ is defined such that  if $v_i \in i$ is the vacuum of the representation $i$ 
$$
Sv_i\otimes v_j = v_{j} \otimes v_{i}.
$$
and it commutes with the representation of the diagonal Virasoro algebra in $\mathcal{A}$. 

This lifts to a Lie algebra automorphism on $\mathcal{A}$ such that if  $l_n^k$ are the generators of the $k^{\text{th}}$ copy of $\text{Vir}$ in $\mathcal{A}$ then $l_n^1 S = Sl_n^1$ and $S^2 = 0$.

Namely we can work in a different basis in $\mathcal{A}$ given by $L^{\pm}_n = l_n \otimes 1 \pm 1\otimes l_n$, therefore we see that
$$
\begin{align*}
SL_n^{\pm} = \pm L_n^{\pm}S. 
\end{align*}
$$
**<u>Definition:</u>** The **fixed point subalgebra** $\mathcal{A}^S \subset \mathcal{A}$ is defined by
$$
\mathcal{A}^S = \{L \in \mathcal{A} \mid SLS = L\}.
$$
Now we can actually start working out how the representations of $\mathcal{A}$ decompose into representations of $\mathcal{A}^S$.

**<u>Proposition:</u>** Let $i$ be an irreducible unitary highest weight representation of $\text{Vir}$ with highest weight vector $v_i$. Then the representation $A_{ii} = i\otimes i$ of $\mathcal{A}$ splits into two representations $A_{ii}^{\pm}$ of $\mathcal{A}^S$, i.e. $A_{ii} = A_{ii}^+ \oplus A_{ii}^-$ such that 
$$
\begin{align*}
A_{ii}^+ = \mathcal{A}^S v_i\otimes v_i && A_{ii}^- = \mathcal{A}^S L^-_k v_i\otimes v_i,
\end{align*}
$$
where $k \in \mathbb{Z}$ is such that $L_k^- \in \mathcal{A}$ is the first generator where $L_{k}^-v_i\otimes v_i \neq 0$. In addition, $v_i\otimes v_j$ and its descendant are the highest weight vectors of the corresponding representations.

***Note:*** Usually $k=1$ unless $i$ is the vacuum module in which case $k=2$. 

***Proof:*** We know that since $v_i$ is cyclic in $i$, for each $u \in A_{ii}$ in the level basis of $A_{ii}$ there a string of generators $U \in \mathcal{A}$ such that $u = U v_i\otimes v_i$. Then if $U \in \mathcal{A}^S$ we have that $u \in A_{ii}^+$. If $U\notin \mathcal{A}^S$ we can proceed as follows. Take $L^-_k$ be the generator with the property in the proposition, and consider $\alpha_k \in \mathbb{C}$ defined by $[(L^-_k)^\dagger,L_{k^-}] v_i\otimes v_i = \alpha_k v_i\otimes v_i$. Then we define $\hat U = \frac{1}{\alpha_k} U (L_{k}^{-})^\dagger$. Then, we can see that $\hat U \in \mathcal{A}^S$ and that the $u = \hat U L_{k}^-v_i\otimes v_i$. 

To show that they are irreducible we see that they are cyclic subspaces of a vector space with an $\mathcal{A}^S$ invariant nondegenerate Hermitian inner product, therefore they are highest weight cyclic representations which implies that they are irreducible. Otherwise the inner product would be degenerate.
$$
\begin{equation}\tag*{$\Box$}\end{equation}
$$
But these are not the only irreps that appear here! Let's see what happens in the other case.

**<u>Proposition:</u>** Let $i,j$ be irreducible unitary highest weight representations of $\text{Vir}$ then there exists a unique, up to isomorphism, irreducible representation $A_{ij}^S$ of $\mathcal{A}^S$ such that it is isomorphic to $A_{ij}$ as a Virasoro module.

***Proof:*** We will consider the representation $A_{ij}\oplus A_{ji}$, and show that it decomposes into $\mathcal{A}^S$ irreducible representations that are isomorphic to $A_{ij}$ as Virasoro modules. Let $v_i,v_j$ be the highest weight vectors of $i,j$ and consider $v_{ij}^{\pm} = v_i \otimes v_j \pm v_j \otimes v_i$ and consider the modules
$$
A_{ij}^{\pm} = \mathcal{A}^S v_{ij}^{\pm}.
$$
Notably, $A_{ij}^{\pm}$ are no longer $\mathcal{A}$ modules because $L_0^{-}$ permutes them. However, by construction, they are $\mathcal{A}^S$ modules. We want to show that they are highest weight. To do so we notice that any generator $L_k^+$ will, either annihilate or increase the $L_0^+$ weight of $v_{ij}^{\pm}$, and so will any generator in $\mathcal{A}^S$ formed by a string of $L^-_k$. However, we need to show that any positive mode of $\mathcal{A}^S$ as a VOA annihilates $v_{ij}^{\pm}$. This is true because we have just showed that $v_{ij}^{\pm}$ has the smallest conformal weight, so if $T$ annihilates everything, then so do the rest of the positive modes of any current since they must satisfy $[L_0,J_n] = -nJ_n \implies L_0 J_nv_{ij}^{\pm} = (h_i + h_j -n)J_0v_{ij}^{\pm}$ which would be less for positive $n$. 

Finally, we show that the two modules are isomorphic. Since they are cyclic, highest weight, unitary modules, they are irreducible. They also have the same highest weight, as well as the action of any generator of $\mathcal{A}^S$ does not depend on $\pm$. Such a module is therefore isomorphic as a Virasoro module to $A_{ij}$.
$$
\begin{equation}\tag*{$\Box$}\end{equation}
$$
This leads us to the following, description of the simple objects in $\text{Rep}(\mathcal{A}^S)$.

**<u>Theorem:</u>** Given a set of simple objects $I \subset \text{Rep}(\text{Vir})$ a set of simple objects $I^S \subset \text{Rep}(\mathcal{A}^S)$ is given by
$$
I^S = \{A_{ii}^{\pm} \in \text{Rep}(\mathcal{A}) \mid i \in I\} \cup \{A_{ij}^S \mid i\neq j \in I \text{ and } (i,j) \sim (j,i)\}.
$$
 ***Proof:*** We know that a complete set of simple objects should appear in a rational CFT with chiral algebra $\mathcal{A}$ and the $A_{ij}^{\pm,S}$ are the ones that appear when we express the same theory as a rational CFT with algebra $\mathcal{A}^S$. 
$$
\begin{equation}\tag*{$\Box$}\end{equation}
$$
Awesome, so we were able to make a nontrivial statement about the representation theory of $\mathcal{A}^S$. Now let's calculate the characters of these representations. Before we proceed, let's prove this useful property. 

**<u>Lemma:</u>** Let $i,j$ be irreducible representations of $\text{Vir}$ then for $L_0 \in \mathcal{A}^S$ and $q \in \mathbb{C}$
$$
\text{tr}_{A_{ij}} Sq^{L_0} = \delta_{ij} \text{tr}_{i} q^{2l_0}.
$$
***Proof:*** Let $B_{ij}$ be the a basis of $A_{ij}$. Then we have that since $v\in A_{ij} \implies Sv \in A_{ji}$ then $\langle Sv,q^{L_0}v \rangle = 0$ since the two modules are orthogonal. Therefore
$$
\begin{align*}
\text{tr}_{A_{ij}}S q^{L_0} = \sum_{v \in B_{ij}} \langle Sv,q^{L_0} v \rangle =  \delta_{ij} \text{tr}_{A_{ii}}S q^{L_0}
\end{align*}
$$
Now we just need to calculate the trace for $A_{ii} = i\otimes i$. We can decompose $i = \bigoplus_{n=0}^\infty i^n$, where $i^n$ is the $n^{\text{th}}$ eigenspace of $l_0 \in \text{Vir}$. As a result, we have that
$$
A_{ii} = \bigoplus_{n=0}^{\infty} \bigoplus_{m=0}^{n} i^{m}\otimes i^{n-m},
$$
when we try to decompose $A_{ii}$ into eigenspaces of $L_0 \in \mathcal{A}$. So we can write the trace as
$$
\text{tr}_{A_{ii}} Sq^{L_0} = \sum_{n=0}^{\infty} q^{2h_i +n} \sum_{m=0}^{n} \text{tr}_{i^m\otimes i^{n-m}} S.
$$
 To calculate the innermost trace we can do a similar trick! Notice that vectors in different eigenspaces are also orthogonal. So if $v \in i^{m}\otimes i^{n-m}$ then $Sv \in i^{n-m}\otimes i^{m}$ which implies that
$$
\tr_{i^m \otimes i^{n-m}} S = \delta_{m,n-m} \tr_{i^m\otimes i^m} S.
$$
We can simplify this further, since if we pick an orthogonal basis $B$ for $i^m$ then
$$
\begin{align*}
\text{tr}_{i^m \otimes i^m} S = \sum_{v,u \in B} \langle v\otimes u,u\otimes v \rangle = \sum_{v,u \in B} \delta_{u,v} = \text{dim\,} i^m.
\end{align*} 
$$
Putting everything together we have that
$$
\text{tr}_{A_{ij}} Sq^{L_0} = \delta_{ij} \sum_{n=0}^{\infty} q^{2h_i +n} \sum_{m=0}^{n} \delta_{n-2m,0} \text{dim\,}i^m = \delta_{ij} \sum_{m=0}^{\infty} q^{2h_i +2m} \text{dim\,}i^m = \delta_{ij}\text{tr}_{i}q^{2l_0}.
$$

$$
\begin{equation}\tag*{$\Box$}\end{equation} 
$$



We will use this property a lot in the upcoming calculations of characters.

**<u>Proposition:</u>** The $\mathcal{A}^S$ characters of $A_{ij}^\pm$ are given by
$$
\begin{align*}
\chi_{A_{ii}^{\pm}}(q) &= \frac{1}{2} \left[ \chi_i(q)^2  \pm \chi_{i}(q^2)\right]\\
\chi_{A_{ij}^S}(q) &= \chi_i(q) \chi_j(q),
\end{align*} 
$$
where $\chi_i$ is a Virasoro character. 

***Proof:*** We use the previous lemma and the first lemma to write the following system of equations
$$
\begin{align*}
\chi_{i}(q)^2 &=\chi_{A_{ii}}(q) = \text{tr}_{A_{ii}^+\oplus A_{ii}^{-}} q^H = \text{tr}_{A_{ii}^+} q^H + \text{tr}_{A_{ii}^-} q^H = \chi_{A_{ii}^+}(q) + \chi_{A_{ii}^-}(q)\\
\chi_{i}(q^2) &= \text{tr}_{A_{ii}^+\oplus A_{ii}^{-}} Sq^H = \text{tr}_{A_{ii}^+} q^H - \text{tr}_{A_{ii}^-} q^H = \chi_{A_{ii}^+}(q) - \chi_{A_{ii}^-}(q).
\end{align*}
$$
Then we solve. For the characters of $A_{ij}^S$ we have that they representations are isomorphic as Virasoro representations to $A_{ij}$ so they must have the same Virasoro characters.
$$
\begin{equation}\tag*{$\Box$}\end{equation} 
$$


## Fixed point under $H$

We want to now find the fixed point subalgebra $\mathcal{W} \subset \mathcal{A}^S$ which is also invariant under $H$. Conceptually we repeat the same idea as before. We only need to find out how the modules of $\mathcal{A}^S$ split under the action of $H$. 

First, let's understand the action of $H$ on $\mathcal{A}^S$. 

**<u>Definition:</u>** Let $\mathcal{A}^S$ be the fixed point algebra of $\mathcal{A}$ under the automorphism $S$. Then the fixed point algebra under $S,H$ of $\mathcal{A}$ is given by
$$
\mathcal{W}=\{L \in \mathcal{A}^S \mid HLH = L\}.
$$
**<u>Lemma:</u>** $\mathcal{W}= \mathcal{A}^S$.

***Proof:*** As a Verlinde line $\eta_i$ commutes with all the generators $l_n^i$ for $i=1,2$. Therefore $H = \eta_1\eta_2$ commutes with $L_n^{\pm}$ which means that it must commute with $L \in \mathcal{A}^S$.
$$
\begin{equation}\tag*{$\Box$}\end{equation}
$$
So we are done!

**<u>Proposition:</u>** Let $i,j$ be irreducible highest weight unitary representations of $\text{Vir}$, then for any representation $A \in I^S$ generated by $i,j$ the map $H:A\to A$ is $H = h_{ij} \text{id}_A$ where
$$
h_{ij} = \frac{s_{i0}s_{j0}}{s_{00}^2},
$$
where $s$ is the $\text{Vir}$ $S$-matrix.



# S-Matrix of Fixed Point Algebra

All we need to do now is to calculate the $S$-matrix for the fixed point algebra irreducible modules. This is straightforward since we know the $S$-matrix for the Virasoro characters, and these are given in terms of them.

**<u>Proposition:</u>** Let $s$ be the Virasoro $S$-Matrix. Then the $S$-Matrix for $\mathcal{W}$ is given by 
$$
\begin{align*}
S_{A_{ij}^S,A_{kl}^S} &= s_{ij} s_{kl} + s_{ik}s_{jl}\\
S_{A_{ij}^S,A_{kk}^\pm} &= s_{ik}s_{jk}
\end{align*} 
$$












