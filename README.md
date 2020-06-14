# Logarithmic_Strain_Space-dealii
C++ code using the deal.ii library for the transformation into the logarithmic strain space

## References and Literature
* Miehe papers
* ...

## The goal
@todo formulate
expand small strain  material models to finite strain
easy to use pre- and post-processing
const. time for pre- and postprocessing, more efficient especially for complicated material models

## Background
@todo add some equations
The transformation consists of three steps. First, we transform the deformation gradient `F` into the logarithmic strain space (ln-space) and obtain the Hencky strain (preprocessing). Secondly, the standard small strain material model is called using the Hencky strain and the usual history. The outcome is the stress measure `T` and the fourth order tangent `C`. Finally, we transform the stress and tangent modulus back from the logarithmic strain space to obtain e.g. the Second Piola-Kirchhoff stress tensor `S`and the Lagrangian tangent modulus `L` (postprocessing).

## How to/Interface
* Add the ln_space header files (.h) and auxiliary functions from this repository to your working directory.
* Include the header "ln_space.h" in your code
* Create an instance of the ln_space class
```
    	ln_space<dim> ln_space;
```
* Execute the preprocessing step using the deformation gradient as input
```
    	ln_space.pre_ln(F);
```
* Get the Hencky strain `hencky_strain` from the ln_space and call your material model (e.g. elastoplasticity(...)) with the Hencky strain as the input strain
```
        SymmetricTensor<2,dim> hencky_strain = ln_space.hencky_strain;
        SymmetricTensor<4,3> tangent = elastoplasticity(/*input->*/ hencky_strain, history, /*output->*/ T_stress);
```
* Execute the postprocessing using the computed `T_stress` and tangent as input
```
        ln_space.post_ln(T_stress,tangent);
```
* Extract the Second Piola-Kirchhoff stress tensor `stress_S` and the Lagrangian tangent modulus `L`
```
    	SymmetricTensor<2,dim> stress_S = ln_space.second_piola_stress_S;
    	SymmetricTensor<4,dim> L = ln_space.C;
```

## ToDo
* rename variables, e.g. 'C' should be lagrangian tangent ...
* avoid directly accessing the member variables, use access functions instead, e.g. with const return values
* merge the for-loops in the 2D pre_ln
* check which header files are needed
* pack all needed external functions into an own auxiliary script (maybe requires namespace to avoid overloading)
* add a remark on how to handle the 2D case (extract_dim, ...)
* add an intro text
* complete the documentation (comments)
* maybe add an example of a basic code calling the matmod with the ln-space (e.g. deal.ii example expanded to finite strains)

