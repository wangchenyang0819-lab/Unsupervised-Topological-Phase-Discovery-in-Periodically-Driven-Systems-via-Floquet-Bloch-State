# Unsupervised Classification of Floquet Topological Phases

This repository contains the complete MATLAB implementation for the paper:  
**“Unsupervised Classification of Topological Phases in Periodically Driven Systems via Floquet-Bloch States”**  
*Machine Learning: Science and Technology* (2026)

The code demonstrates an unsupervised learning framework that clusters Floquet topological phases using a kernel defined directly in momentum‑time \((k,t)\) space, without requiring prior knowledge of topological invariants. The method is validated on three symmetry classes: **A (2D)**, **AIII (1D)**, and **D (1D)**.

## Repository Structure
.
├── README.md # This file
├── class_A/ # 2D Floquet topological insulator (class A)
│ ├── README.md # Detailed instructions for class A
│ ├── TwoAdiffusion.m # Main clustering script
│ ├── blochA.m # Bloch vector visualisation
│ ├── twoDAinvariants.m # Theoretical invariant calculation
│ └── fig2.m # Figure reproduction script (not provided)
├── class_AIII/ # 1D Chiral symmetric Floquet system (class AIII)
│ ├── README.md # Detailed instructions for class AIII
│ ├── AIIIinvariantsdiffusion.m # Clustering + invariants
│ ├── blochAIII.m # Bloch vector visualisation
│ └── fig1.m # Figure reproduction script (not provided)
└── class_D/ # 1D Floquet topological superconductor (class D)
├── README.md # Detailed instructions for class D
├── dinvariantsdiffusion.m # Clustering + theoretical phase diagram
├── blochd.m # Bloch vector visualisation + Chern number
└── fig3.m # Figure reproduction script (not provided)

text

## Two Ways to Compute Floquet‑Bloch States

All scripts support two equivalent representations of the Floquet‑Bloch states:

1. **Frequency‑momentum (Floquet) representation**  
   - The effective Hamiltonian is constructed in an extended Floquet space including multiple frequency copies (photon sectors).  
   - Diagonalisation yields the quasi‑energy spectrum and the Floquet eigenstates \(|\phi_n(\boldsymbol{k},t)\rangle\) expanded in Fourier harmonics.  
   - This approach is used in all `*diffusion.m` and `bloch*.m` scripts.

2. **Direct time evolution**  
   - Alternatively, the time‑ordered evolution operator over one period \(U(\boldsymbol{k},T)\) can be computed by numerical integration.  
   - The Floquet eigenstates are then obtained from the eigenvectors of \(U(\boldsymbol{k},T)\).  
   - This method is used in `twoDAinvariants.m` for invariant calculation.

Both representations are physically equivalent when the frequency‑domain truncation is sufficiently large.

## Key Numerical Considerations

### Frequency‑domain truncation (number of photon sectors)

- The Floquet Hamiltonian is truncated to include **2×omg + 1** frequency copies (ranging from `-omg` to `+omg`).  
- **Sufficiently large `omg`** (e.g., 10–30) is required to ensure that the Floquet‑Bloch states are properly normalised and that the quasi‑energies converge.  
- Insufficient `omg` leads to artificial couplings and incorrect topological signatures.  
- In the provided scripts, `omg` is chosen high enough that further increases do not change the clustering results (convergence tested).

### Discretisation of momentum and time

- The Brillouin zone and one driving period are discretised into `N` (momentum) and `Nt` (time) grid points.  
- **The grid must be fine enough to resolve all band‑inversion points** (spin‑flip points) that separate topologically distinct samples.  
- Typically, `N = 40–60` and `Nt = 40–60` are sufficient.  
- For 2D systems (`N × N` momentum grid) the total number of points is larger; careful balancing between resolution and computational cost is required.

### Kernel hyperparameter ε

- The similarity kernel is defined as \(K_{ij} = \prod_{\boldsymbol{k},t} \left(1 - e^{-|\det(\cdots)|^2/\epsilon^2}\right)\).  
- **ε must be small enough** that topologically identical samples are not split into separate clusters.  
- At the same time, ε should not be so small that numerical noise causes artificial distinctions.  
- The optimal ε is class‑dependent:  
  - Class A: ε = 0.1  
  - Class AIII: ε = 0.5  
  - Class D: ε = 0.03  

In all cases, ε is chosen such that the kernel yields eigenvalues of the diffusion map that clearly separate the spectral gaps (eigenvalues close to 1 indicate the number of clusters).

## Clustering Pipeline (Common to All Diffusion Scripts)

Each `*diffusion.m` script implements the following three‑step workflow:

1. **Kernel construction**  
   - For each sample (parameter point), compute the flattened Floquet operator (FFO) from the Floquet‑Bloch states.  
   - Evaluate the kernel matrix \(K_{ij}\) using the determinant of summed projection operators.

2. **Diffusion map**  
   - Normalise the kernel to a Markov matrix \(P\).  
   - Perform eigenvalue decomposition.  
   - The number of eigenvalues close to 1 (spectral gap) determines the embedding dimension.

3. **Clustering**  
   - Take the first \(d\) non‑trivial eigenvectors as low‑dimensional coordinates.  
   - Apply **hierarchical clustering** (single linkage) with the number of clusters set to the number of eigenvalues near 1.  
   - Map the cluster labels back to parameter space to obtain the unsupervised phase diagram.

## Topological Invariants for Validation

For each symmetry class, separate scripts compute the theoretical invariants:

- **Class A** – Chern number \(C\) and Pontryagin index \(\nu\) (3D winding number).  
- **Class AIII** – Winding numbers \(W_0\) (0‑gap) and \(W_\pi\) (π‑gap).  
- **Class D** – Analytic phase diagram based on trigonometric functions + Chern number from Bloch vectors.

These invariants are used **only for validation**, never as input to the unsupervised clustering. They confirm that the clusters correspond exactly to the distinct topological phases.

## Bloch Vector Visualisation

The `bloch*.m` scripts generate the Bloch vector components \((x,y,z)\) of the Floquet‑Bloch state on the \((k,t)\) torus.  
- These data (saved as `*.mat` files) are used to create **linking‑point visualisations** (e.g., red/blue points in the paper).  
- They also provide an intuitive geometric picture of why the kernel groups states with the same topological texture.

## Reproducibility

All parameters (grid sizes, ε, number of samples, harmonic truncation order, etc.) are explicitly set in the respective scripts.  
The code is self‑contained and produces exactly the same figures as in the paper when the reproduction scripts (`fig1.m`, `fig2.m`, `fig3.m`) are executed.

For details specific to each symmetry class (parameter ranges, driving protocols, invariant formulas, figure generation), please refer to the README inside each subfolder.

## Citation

If you use this code in your own research, please cite:
[Authors], “Unsupervised Classification of Topological Phases in Periodically Driven Systems via Floquet-Bloch States”,
