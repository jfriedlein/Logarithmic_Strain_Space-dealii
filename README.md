# Logarithmic_Strain_Space-dealii
C++ code using the deal.ii library for the transformation into the logarithmic strain space

## Background
The transformation consists of three steps. First, we transform the deformation gradient F into the logarithmic strain space and obtain the Hencky strain (preprocessing). Secondly, the standard small strain material model is called using the Hencky strain and the usual history. The outcome is the stress measure T and the tangent C. Finally, we transform the stress and tangent modulus back from the logarithmic strain space to obtain e.g. the Second Piola-Kirchhoff stress tensor and the Lagrangian tangent modulus (postprocessing).

## How to/Interface
* Add the ln_space header files and auxiliary functions to your working directory.
* Include the header "ln_space.h"
* Create an instance of the ln_space class
    	ln_space<dim> ln_space;
* Execute the preprocessing step using the deformation gradient as input
    	ln_space.pre_ln(F);
* Call your material model (e.g. elastoplasticity(*)) with the Hencky strain as the input strain
        SymmetricTensor<2,dim> hencky_strain = ln_space.hencky_strain;
        SymmetricTensor<4,3> tangent = elastoplasticity(/*input->*/ hencky_strain, history, /*output->*/ T_stress);
* Execute the postprocessing using the computed T_stress and tangent as input
        	ln_space.post_ln(T_stress,tangent);
* Extract the Second Piola-Kirchhoff stress tensor stress_S and the Lagrangian tangent modulus C
    	SymmetricTensor<2,dim> stress_S = ln_space.second_piola_stress_S;
    	SymmetricTensor<4,dim> C = ln_space.C;

## ToDo
* ADD the special outer_product_sym functions !!!
* avoid directly accessing the member variables, use access functions instead, e.g. with const return values
* merge the for-loops in the 2D pre_ln
* check which header files are needed
* pack all needed external functions into an own auxiliary script (maybe requires namespace to avoid overloading)
* think about merging the pre_ln 2D and 3D cases
* add a remark on how to handle the 2D case (extract_dim, ...)
* add an intro text
* complete the documentation (comments)
* maybe add an example of a basic code calling the matmod with the ln-space
* maybe combine all three steps into a single call, e.g. with function pointer to material model

