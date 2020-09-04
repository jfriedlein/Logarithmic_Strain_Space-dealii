# Logarithmic_Strain_Space-dealii
C++ code using the deal.ii library for the transformation into the logarithmic strain space

(also available for Fortran, but not yet documented and uploaded)

## The goal/When to use this code
The logarithmic strain space (herein often abbreviated as ln-space) is a very simple (in terms of the application) way to apply small strain material models to finite strains. So, in case you have a small strain model that you would like to use for applications exposed to large deformations/finite strains, the ln-space it probably the easiest way to achieve this.

All you need are three steps as schematically (simplified) shown in the figure:
1. Pre-processing from the world of finite strains into the logarithmic strain space
2. Calling your small strain model with the Hencky strain computed in the pre-processing
3. Post-processing the computed stress and tangent(s) by transforming them from the ln-space back into the real world.
And the best is, steps 1 and 3 are always the same (copy-paste) and the second step is the same as done for the small strain implementation.

<img src="https://github.com/jfriedlein/Logarithmic_Strain_Space-dealii/blob/master/images/ln-space%20-%20overview.png" width="700">

Drawbacks?

Especially step 3, the post-processing is expensive independent of your material model. So, as far as efficiency is concerned, a simple material model utilising the ln-space is most certainly slower, than its derived finite strain equivalent (a model that was developed for finite strains). However, for complicated material models it can be faster to use the small strain model in the ln-space, instead of a finite strain equivalent (and it requires no further derivations/development to get from small to finite strains). Another disadvantage, is the limitation to small elastic strains. The latter is usually satisfied for metal forming and similar applications, where elastic strain are small and large plastic strain occur.

## Background
@todo add some equations

The transformation consists of three steps. First, we transform the deformation gradient `F` into the logarithmic strain space (ln-space) and obtain the Hencky strain (preprocessing). Secondly, the standard small strain material model is called using the Hencky strain and the usual history. The outcome is the stress measure `T` and the fourth order tangent `C`. Finally, we transform the stress and tangent modulus back from the logarithmic strain space to obtain e.g. the Second Piola-Kirchhoff stress tensor `S`and the Lagrangian tangent modulus `L` (postprocessing).

## Interface/How to use this code
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

## The code/ A look under the hood
The documentation for the ln-space code can be found here https://jfriedlein.github.io/Logarithmic_Strain_Space-dealii/html/classln__space.html.

@todo Create a proper mainpage instead of this namespace-thing

## References and Literature
* papers on the algorithms and application:
    * Miehe, C. and Lambrecht, M. (2001), Algorithms for computation of stresses and elasticity moduli in terms of Seth–Hill's family of generalized strain tensors. Commun. Numer. Meth. Engng., 17: 337-353. doi:10.1002/cnm.404
    * Miehe, C. & Apel, N. & Lambrecht, M.. (2002). Anisotropic additive plasticity in the logarithmic strain space: Modular kinematic formulation and implementation based on incremental minimization principles for standard materials. Computer Methods in Applied Mechanics and Engineering - COMPUT METHOD APPL MECH ENG. 191. 5383-5425. 10.1016/S0045-7825(02)00438-3.
    * Apel, N.: Approaches to the description of anisotropic material behaviour at finite
elastic and plastic deformations. Stuttgart, Universität Stuttgart, 2004. (more comprehensive than Miehe et al papers)
    * Schmaltz, S. (2015). Inverse Materialparameteridentifikation von Blechwerkstoffen für ein anisotropes elasto-plastisches Materialmodell bei finiten Deformationen. [PDF](urn:nbn:de:bvb:29-opus4-58089) (great outline of algorithms and implementation)


## Credits
A huge thanks to Christian Burkhardt for implementing the algorithms outlined in the papers by Miehe et al. into C++ (and making sense of the immense (and partly insane) index notation). So I just had to organise the code into a class and add the 2D case.

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

