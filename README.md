# ***Perturb-STNet***: Fusion of Spatiotemporal and Network Models to Prioritize Nested Effects in Single-Cell Perturbations
## Description
_Perturb-STNet_ is a generic framework for analyzing spatial and spatiotemporally resolved single-cell perturbation data. It identifies spatial and temporally differentially expressed regulators arising from the impact of perturbations. _Perturb-STNet_ offers significant advantages over existing techniques. One key benefit is its ability to relax the traditional assumption of matching pairs between control and perturbed conditions in differential analysis. This flexibility allows for testing multiple hypotheses more effectively. The framework achieves this through a method called Differential Nested Effects (DNE), which compares regulator variability in cell-neighborhood network graphs between control conditions and combined control and perturbed conditions.

_Perturb-STNet_ analysis pipeline can be broken down into two major phases. The first phase is to estimate a cell-neighborhood biological network and adjust for spatial misalignment across multiple time points and samples. The second part uses the estimated network as leverage in the spatiotemporal statistical modeling of regulators to detect differentially expressed regulators (e.g. genes, proteins) across space and time. 

The analysis pipeline is summarized in the figure below. 

![image](https://github.com/user-attachments/assets/624a5741-b688-4034-9f90-ac688d1a18b7)

<p align="justify">From the Figure,  (a) shows the introduction of external intervention (drugs) on grouped subjects across different days, followed by tissue extraction and protein imaging or spatial transcriptomics. (b) This process results in spatially resolved image and expression data with discrepancies in spatial coordinates over time, and possibly samples (not shown). (c) The spatially indexed image is partitioned using a mesh, where each polygon represents a collection of cells in the same neighborhood. (d) Each cell type within a polygon forms a node, and edges connecting nodes in the network are established using a minimum spanning tree algorithm on the gene or protein expression profiles at each time point, which is then combined to obtain the background network consistent across time points. (e) In this illustration, the network at a given time point (e.g. day 0) is considered the control, while other days are considered perturbed because, in this example, the interest includes identifying the regulators (e.g. genes) whose expression changes over space and time compared to the time at day 0. The background network under control and perturbed conditions are represented as nested precision matrices. (f) With the precision matrices and the expression data, spatio-temporal statistical models are fitted and the signal-to-noise ratio (SNR) is estimated for the control and perturbed conditions. (e) The framework generates interpretable outputs such as regulator ranking, regulatory networks, and annotations of cellular microenvironments. </p>

## Usage
Detailed instructions or examples on how to use the tool. This can include code snippets, command-line examples, or screenshots.
## Parameters
If your tool has configurable options or parameters, describe them here.
## Examples
Provide examples or case studies showing how your tool can be used in practice.
## Acknowledgments
Any acknowledgments for libraries, contributors, or inspiration for your project?
## Contact Information
How users can reach you for questions or feedback.
