# Logarithmic_Strain_Space-dealii
C++ code using the deal.ii library for the transformation into the logarithmic strain space

## ToDo
* combine the three steps (pre, call to small strain material model, post) into as single call and e.g. put the second step into a function pointer (don't know how in combination with template, ...)
e.g. by a function pointer that contains the call to the material model
	ln_space::transform (  /*input->*/ DeformationGradient, alpha_n, eps_p_n, func_pointer, /*output->*/ stress_S, Tangent);
This would eradicate all the unnecessary declarations and arguments and make the ln-space a one-liner.
Current Issues in this regard: No clue of function pointers for template functions or workarounds for this
