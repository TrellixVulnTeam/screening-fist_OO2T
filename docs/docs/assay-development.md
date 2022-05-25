# Assay Development

## Summary

This page describes development of the high throughput P450-compound binding assay used in this project.
The technique on which the assay is based is also described.

## Aim

The initial aims of this development work were:

- Develop a high throughput P450-ligand binding assay based on established biophysical characterisation techniques.
- Develop necessary software for design and analysis of each assay.
- Compare the precision and accuracy of the assay to existing techniques.

## Basis: UV-Visible Spectroscopy for Monitoring Cytochrome P450-Ligand Binding

The assay is based on a technique for quantifying P450-ligand interactions based on UV-visible photospectroscopy.
The technique consists of the purified Cytochrome P450 heme domain in question in a neutral buffer at around 5-10 uM in a optically clear cuvette.
Since only the heme-containing domain of the P450 is used, no chemical reactions are expected to take place which removes time-sensitivity from the assay.



The UV-visible light absorbance of the sample is typically measured for all wavelengths between 200 and 800 nm, which for a P450 without a ligand bound in the active site should show a large and defined absorbance peak at around 420 nm.

After an initial absorbance measurement of the ligind-free P450, the compound of interest can be titrated into the sample.
On binding to the ligand, the absorbance profile of the P450 changes such that the absorbance at 420 nm ($A_{420}$) decreases and absorbance at 390 nm ($A_{390}$) increases.

The change in $A_{420}$ and $A_{390}$ in response to change in ligand concentration can be quantified and used to derive metrics that indicate affinity between the ligand and P450 using Michaelis-Menten kinetics models.

The original Michaelis-Menten model of enzyme kinetics states:

$$ v = V_{max} \frac{[S]}{[S] + K_M} $$

where $v$ is the reaction velocity - the rate of an enzymatic reaction. 
$V_{max}$ is the maximum possible $v$ for this particular enzyme-substrate pair, $[S]$ is the concentration of the substrate and $K_M$ is the $[S]$ at which $v = V_{max}$.

$V_{max}$ and $K_M$ are useful metrics for quantifying the binding interaction between enzyme and substrate, where low $K_M$ indicates a tight binding interaction and a high $V_{max}$ indicates a large magnitude of response.

Important assumptions in the Michaelis-Menten model of kinetics are:

1. The concentration of enzyme is $< K_d$ 
2. The rate of reaction is directly proportional to the concentration of the substituents
3. The reaction is at chemical equilibrium at the time of measurement
4. The interaction is reversible

A variant of this model is applied to Cytochrome P450 photospectroscopy assays.
$v$ is substituted for $\Delta A_{390} - \Delta A_{420}$ - the magnitude of the P450 response and $K_M$ is substituted for $K_d$ - the dissociation constant between the enzyme and ligand.
This yields the formula:

$$\Delta A_{390} - \Delta A_{420} = V_{max} \frac{[S]}{[S] + K_d} $$


!!! todo
	- reference this
	- figure
	- apply to example?

## Early Iterations

- compound dispensing
- buffer conditions
- echo
- `echo`

## Refinement

