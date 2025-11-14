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
\eta v_i\otimes \bar v_i = \eta_i v_i \otimes \bar v_i= \frac{S_{\eta i}}{S_{\eta 1}}v_i\otimes \bar v_i && Nv_i\otimes \bar v_i = N_i v_i\otimes \bar v_i= \frac{S_{Ni}}{S_{N1}}v_i \otimes \bar v_i,
\end{align*}
$$
where the corresponding irreducible highest weight Virasoro representations have weights $h_{\eta} = \frac{3}{2}$ and $h_N = \frac{7}{16}$.

***Proof:*** Oh nooo, just work out the fusion rules of the representations of $\text{Vir}_{c=\frac{7}{10}}$.
$$
\begin{equation}\tag*{$\Box$}\end{equation}
$$

Before we proceed, let's state a result that we have derived elsewhere so that we can fix notation.

**<u>Lemma:</u>** Consider a rational CFT obtained from folding a unitary minimal model with central charge $c$. It has chiral algebra $\mathcal{V}=\text{Vir}_c\otimes \text{Vir}_c$ and a discrete symmetry $S$ inherited by the outer automorphism of $\mathcal{V}$ that exchanges the two copies. If $I$ is the set of irreducible highest weight unitary representations of $\text{Vir}_{c}$ the the bulk Hilbert space decomposes as
$$
\mathcal{H}= \bigoplus_{i,j \in I} (i\otimes \bar i) \otimes (j \otimes \bar j).
$$
If the fixed point algebra under $S$ is given by $\mathcal{V}^S$ then the same Hilbert space can be written as
$$
\mathcal{H}= \left[\bigoplus_{\substack{i\in I\\a,b = \pm}} V_{ii}^a\otimes \overline{V_{ii}^b}\right] \oplus \left[ \bigoplus_{ij \in I} V_{ij}^S\otimes \overline{V_{ij}^S}  \right]
$$
 where $V_{ij}^s \cong V_{ji}^S$ as irreducible $\mathcal{V}^S$ modules.



**<u>Proposition:</u>** Let $L$ be a simple TDL and $L^- = L\otimes 1 -  1\otimes L$. Then the action of $L^-$  on the irreducible representations of $\mathcal{V}^S$ in the bulk Hilbert space on the theory with maximal chiral algebra $\mathcal{V}$ is given by
$$
\begin{align*}
\left.L^-\right|_{V_{ii}^{\pm}} = 0 && \left.L^-\right|_{V_{ij}^S} = (L_i - L_j) \text{Id}_{V_{ji}^{S}},
\end{align*}
$$
where $V_{ij}^S \otimes V_{ji}^S = V_{ij} \otimes V_{ji}$ with $V_{ij}^S \cong V_{ji}^S$.

***Proof:*** The highest weight vectors of $V_{ii}^\pm$ are descendants of $v_i\otimes v_i$, therefore $L^-v_{i}\otimes v_{i} = (h_i - h_i)v_{i}\otimes v_i$. For $V_{ij}^S$ the highest weight vector in given by $P^\pm v_i \otimes v_j$ where $P^{\pm} = \frac{1}{\sqrt{2}}(1 \pm S)$ where the $\pm$ corresponds to weather we are considering the left or right copy in the decomposition of $V_{ij} \otimes V_{ji}$. From here we have that $L^- P^{\pm} = P^{\mp} L^-$ and that $L^- v_i\otimes v_j = (L_i - L_j) v_i \otimes v_j$.
$$
\begin{equation}\tag*{$\Box$}\end{equation}
$$

So we see that one of them acts trivially, but the other exchanges the irreducible representations, which is slightly nontrivial. What we need to figure out, however, is how to induce the action of the lines in the rest of the irreducible representations of $\mathcal{V}^S$. 









