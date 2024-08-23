# Perturb-STNet
\texttt{Perturb-STNet} is a generic framework for spatial and spatiotemporally resolved single-cell perturbation data, and it can also accommodate non-spatially resolved data by using the high-dimensional gene or protein expression profile to derive cell embeddings used as pseudo-spatial cell coordinates to identify differential variable regulators as demonstrated in the EMT validation analysis. In addition, \texttt{Perturb-STNet} is not limited to two-dimensional spatial coordinates; it supports three-dimensional spatial data by creating polyhedron instead of polygon meshes for network nodes \citep{wang2018three,guo2023three}. \texttt{Perturb-STNet} quantifies nested regulator variabilities and can validate different perturbation hypotheses or trajectories by specifying control and perturbed timelines. Unlike prediction-focused methods
