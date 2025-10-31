# Ising Invariance

One of the main properties of the fixed point we are looking for is that it commutes with Ising generators. Here we will interpret that condition on the boundary state that corresponds to the defect itself.

[toc]

# Diagonal Action

Let $D$ be a conformal defect along the real line on the sphere, and consider a topological defect $L$ parallel to it such that $LD = DL$. When folding along $D$ we obtain a boundary $\psi_D$ in the folded theory, therefore that boundary state must satisfy
$$
\left(L\otimes 1 - 1\otimes L\right) \psi_D = 0.
$$
Let's call the operator $L^- \coloneqq L\otimes 1 - 1\otimes L$. Then we have that if $\psi_D$ originates via pullback from the $\mathbb{Z}_2$ cyclic orbifold, it must be invariant under the exchange map $S$, i.e. $S\psi_D = \psi_D$. However, we find that $SL^- = - L^-S$. Now consider the projection operator $P^{\pm} = \frac{1 + S}{2}$, which leads to the decomposition of $L^-$ into 4 blocks
$$
\begin{align*}
L^- = P^{+}L^- P^+ + P^{-}L^- P^+ + P^{+}L^- P^- + P^{-}L^- P^- = P^{-}L^- P^+  + P^{-}L^- P^-,
\end{align*}
$$
since $P^+L^- = 0$. Finally, what we need to calculate is
$$
L^- \psi_D = P^- L^- P^+ \psi_D = 0,
$$
so we only really need to calculate $P^-L^-P^+$.

## Action of Ising Verlinde Lines

Let's now figure out what the action of the Verlinde lines is on irreducible representations of the fixed point algebra $\mathcal{V}^S$. Here is a useful result.

**<u>Proposition:</u>** The operators associated to Verlinde lines in a CFT with VOA $\mathcal{V}$ act like multiples of the identity in the irreducible representations of the fixed point algebra $\mathcal{V}^S$. 

***Proof:*** Since the fixed point algebra $\mathcal{V}^S$ is a subalgebra of $\mathcal{V}$ the Verlinde lines commute with it, therefore the result is guaranteed by Schur's lemma.
$$
\begin{equation}\tag*{$\Box$}\end{equation}
$$
Now all we need to do is to find these numbers. The good thing is though that we only need to find how they act on the highest weight vector, and since (for most of them) their highest weight vectors are in Tricritical Ising, we can figure it out relatively simply.

**<u>Proposition:</u>** The representation of the Ising fusion category in Tricritical Ising is given by $\eta,N$ such that for the highest vector $v_i$ of the Virasoro irreducible unitary highest weight representation $i$
$$
\begin{align*}
\eta v_i = \frac{S_{\eta i}}{S_{\eta 1}}v_i && Nv_i = \frac{S_{Ni}}{s_{N1}}v_i,
\end{align*}
$$
where the corresponding irreducible highest weight Virasoro representations have weights $h_{\eta} = \frac{3}{2}$ and $h_N = \frac{7}{16}$.

***Proof:*** Oh nooo, just work out the fusion rules of the representations of $\text{Vir}_{c=\frac{7}{10}}$.
$$
\begin{equation}\tag*{$\Box$}\end{equation} 
$$






