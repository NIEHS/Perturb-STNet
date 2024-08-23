# ***Perturb-STNet***: Fusion of Spatiotemporal and Network Models to Prioritize Nested Effects in Single-Cell Perturbations
## Description
_Perturb-STNet_ is a generic framework for analyzing spatial and spatiotemporally resolved single-cell perturbation data. It identifies spatial and temporally differentially expressed regulators arising from the impact of perturbations. _Perturb-STNet_ offers significant advantages over existing techniques. One key benefit is its ability to relax the traditional assumption of matching pairs between control and perturbed conditions in differential analysis. This flexibility allows for testing multiple hypotheses more effectively. The framework achieves this through a method called Differential Nested Effects (DNE), which compares regulator variability in cell-neighborhood network graphs between control conditions and combined control and perturbed conditions.

_Perturb-STNet_ analysis pipeline can be broken down into two major phases. The first phase is to estimate a cell-neighborhood biological network and adjust for spatial misalignment across multiple time points and samples. The second part uses the estimated network as leverage in the spatiotemporal statistical modeling of regulators to detect differentially expressed regulators (e.g. genes, proteins) across space and time. 

The analysis pipeline is summarized in the figure below. 

![image](https://github.com/user-attachments/assets/12c9bbd2-aa28-4e8d-ad85-4688dd120c5e)

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
