/**

@{

\brief Implementation of transformation into the logarithmic strain space in deal.ii

\author jfriedlein

@tableofcontents

@mainpage Transformation into the logarithmic strain space

@section intro Introduction
@todo refer to papers by Miehe


@subsection subsec_overview Overview
This overview shall give you a first impression what to expect ...


@subsection subsec_interface Interface - How to incorporate the transformation into your code
@todo add a simple example code that calls the pre and post with e.g. an elastic matmod as coding example in section code


@subsection subsec_more_basics Some more basics



@subsection subsec_resources Some resources/links



@section code The commented program

\code
#ifndef ln_space_H
#define ln_space_H
 
\endcode
@section includes Include Files
The data type SymmetricTensor and some related operations, such as trace, symmetrize, deviator, ... for tensor calculus
\code
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/exceptions.h>
 
\endcode
@todo Check whether the following three headers are needed at all
\code
#include <iostream>
#include <fstream>
#include <cmath>
 
#include "functions.h"
 
using namespace dealii;
 
namespace ln_space
{
//	template<int dim>
//	void transform (  /*input->*/ Tensor<2, dim> &F, double &alpha_n, SymmetricTensor<2,dim> &eps_p_n, void (*func_pointer)(auto),
//			 	 	  /*output->*/ SymmetricTensor<2,dim> &stress_S, SymmetricTensor<4,dim> &C, SymmetricTensor<4,dim> &Lambda, double &stress_vM )
//	{
//		 Vector<double> eigenvalues(dim);
//		 std::vector< Tensor<1, dim> > eigenvector(dim);
//		 std::vector< SymmetricTensor<2, dim> > eigenbasis(dim);
//		 Vector<double> ea(dim);
//		 Vector<double> da(dim);
//		 Vector<double> fa(dim);
//		 SymmetricTensor<2, dim> hencky_strain;
\endcode

\code
//		 pre_ln ( /*input->*/ F,
//				  /*output->*/ hencky_strain, ea, da, fa, eigenvector, eigenvalues, eigenbasis );
\endcode

\code
//		 // We only create this variable to emphasise the difference between the stress \a T from the
//		 // small strain material modell and the 2nd Piola Kirchhoff stress \a S transformed back by the post processing
//		  SymmetricTensor<2,dim> T_stress;
\endcode

\code
//		 func_pointer( /*input->*/ hencky_strain, alpha_n, eps_p_n, /*output->*/ T_stress, stress_vM, C, Lambda );
\endcode

\code
//		 // ToDo: We also have to transform Lambda, stress_vM
\endcode

\code
//		 // also transforms the T_stress to the 2nd Piola-Kirchhoff stress
//		  post_ln ( /*input->*/ ea, da, fa, eigenvalues, eigenbasis,
//				 				  /*output->*/ T_stress, C );
\endcode

\code
//		 stress_S = T_stress;
//	}
 
 
	template<int dim>
	void pre_ln ( 	/*input->*/ Tensor<2, dim> &F, 
					/*output->*/ SymmetricTensor<2,dim> &hencky_strain, Vector<double> &ea, Vector<double> &da, Vector<double> &fa,
								 std::vector<Tensor<1,dim> > &eigenvector, Vector<double> &eigenvalues, std::vector< SymmetricTensor<2,dim> > &eigenbasis ) 
	{
\endcode
Following "Algorithms for computation of stresses and elasticity moduli in terms of Seth–Hill’s family of generalized strain tensors" by Miehe&Lambrecht \n
Table I. Algorithm A
\code
		/*
		 * 1. Eigenvalues, eigenvalue bases and diagonal functions:
		 */
 
\endcode
Get the symmetric right cauchy green tensor
\code
		 SymmetricTensor<2, dim> right_cauchy_green_sym = symmetrize( transpose(F) * F);//Physics::Elasticity::Kinematics::C(F);
 
\endcode
Compute Eigenvalues, Eigenvectors and Eigenbasis
\code
		 {
\endcode
Get Eigenvalues and Eigenvectors from the deal.ii function \a eigenvectors(*)
\code
			for (unsigned int i = 0; i < dim; ++i) {
				eigenvalues[i] = eigenvectors(right_cauchy_green_sym)[i].first;
				eigenvector[i] = eigenvectors(right_cauchy_green_sym)[i].second;
			}
 
\endcode
Check if the found eigenvectors are perpendicular to each other
\code
			if ((std::fabs(eigenvalues(0) - 1) > 1e-10)
				&& (std::fabs(eigenvalues(1) - 1) > 1e-10)
				&& (std::fabs(eigenvalues(2) - 1) > 1e-10)) {
				for (unsigned int i = 0; i < dim; ++i) {
					for (unsigned int j = i + 1; j < dim; ++j) {
						AssertThrow( (std::fabs(eigenvector[i] * eigenvector[j]) < 1e-12),
									 ExcMessage("ln-space<< Eigenvectors are not perpendicular to each other.") );
					}
				}
			}
 
\endcode
Compute eigenbasis Ma: eigenbasis = eigenvector \otimes eigenvector
\code
			for (unsigned int i = 0; i < dim; ++i) {
				eigenbasis[i] = symmetrize( outer_product(eigenvector[i], eigenvector[i]) );
				AssertThrow( eigenvalues(i) >= 0.0,
							 ExcMessage("ln-space<< Eigenvalue is negativ. Check update_qph.") );
			}
		 }
 
\endcode
Compute diagonal function \a ea and its first and second derivate \a da and \a fa \n
\code
		 for (unsigned int i = 0; i < dim; ++i)
		 {
			ea(i) = 0.5 * std::log( std::abs(eigenvalues(i)) );	// diagonal function
			da(i) = std::pow(eigenvalues(i), -1.0);				// first derivative of diagonal function ea
			fa(i) = -2.0 * std::pow(eigenvalues(i), -2.0);			// second derivative of diagonal function ea
			AssertThrow( ea(i) == ea(i),
						 ExcMessage( "ln-space<< Ea is nan due to logarithm of negativ eigenvalue. Check update_qph.") );
			AssertThrow( da(i) > 0.0,
						 ExcMessage( "ln-space<< First derivative da of diagonal function is "+std::to_string(da(i))+" < 0.0 . Check update_qph.") );
		 }
 
\endcode
Compute the Hencky strain
\code
		 for (unsigned int a = 0; a < dim; ++a)
			hencky_strain += ea(a) * eigenbasis[a];
 
\endcode
Output-> SymmetricTensor<2, dim> hencky_strain, Vector<double> ea, da, fa,
\code
		//			std::vector<Tensor<1, dim>> eigenvector, Vector<double> eigenvalues, std::vector< SymmetricTensor<2, dim> > eigenbasis
	}
 
 
	template<int dim>
	void post_ln ( 	/*input->*/ Vector<double> &ea, Vector<double> &da, Vector<double> &fa,
			 	 	 	 	 	Vector<double> &eigenvalues, std::vector< SymmetricTensor<2,dim> > &eigenbasis,
					/*output->*/ SymmetricTensor<2,dim> &second_piola_stress_S, SymmetricTensor<4,dim> &elasto_plastic_tangent	)
	{
		const double comp_tolerance = 1e-8;
		SymmetricTensor<2,dim> stress_measure_T_sym = second_piola_stress_S; // the output argument is abused as an input argument
 
		/*
		 * 3. Set up coefficients \a theta, \a xi and \a eta
		 */
 
\endcode
Compute the coefficients based on the eigenvalues, eigenvectors and ea,da,fa
\code
		 Tensor<2, dim> theta;
		 Tensor<2, dim> xi;
		 double eta = 999999999.0;
 
\endcode
For three different eigenvalues \f$ \lambda_a \neq \lambda_b \neq \lambda_c \f$
\code
		 if (
				 ( !(std::fabs(eigenvalues(0) - eigenvalues(1)) < comp_tolerance) )
				 &&
				 ( !(std::fabs(eigenvalues(0) - eigenvalues(2)) < comp_tolerance) )
				 &&
				 ( !(std::fabs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance) )
			 )
		 {
			eta = 0.0;
			for (unsigned int a = 0; a < dim; ++a)
				for (unsigned int b = 0; b < dim; ++b)
					if (a != b)
					{
						theta[a][b] = (ea(a) - ea(b))
									  / (eigenvalues(a) - eigenvalues(b));
						xi[a][b] = (theta[a][b] - 0.5 * da(b))
								   / (eigenvalues(a) - eigenvalues(b));
 
						for (unsigned int c = 0; c < dim; ++c)
							if ((c != a) && (c != b))
							{
								eta +=
										ea(a)
										/ (2.0
										   * (eigenvalues(a)
											  - eigenvalues(b))
										   * (eigenvalues(a)
											  - eigenvalues(c)));
							}
					}
		 }
		//	For three equal eigenvalues \f$ \lambda_a = \lambda_b = \lambda_c \f$
		 else if ( (std::fabs(eigenvalues(0) - eigenvalues(1)) < comp_tolerance)
					&&
				   (std::fabs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance) )
		 {
			eta = 0.0;
			for (unsigned int a = 0; a < dim; ++a)
			{
				for (unsigned int b = 0; b < dim; ++b)
					if (a != b)
					{
						theta[a][b] = 0.5 * da(0);
						xi[a][b] = (1.0 / 8.0) * fa(0);
					}
			}
			eta = (1.0 / 8.0) * fa(0);
		 }
 
\endcode
For two equal eigenvalues a and b: \f$ \lambda_a = \lambda_b \neq \lambda_c \f$
\code
		 else if ( (std::fabs(eigenvalues(0) - eigenvalues(1)) < comp_tolerance)
				   &&
				   ( !(std::fabs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance) ) )
		 {
			eta = 0.0;
			for (unsigned int a = 0; a < dim; ++a)
				for (unsigned int b = 0; b < dim; ++b)
					if ((a != b) && ((a == 2) || (b == 2)))
					{
						theta[a][b] = (ea(a) - ea(b))
									  / (eigenvalues(a) - eigenvalues(b));
						xi[a][b] = (theta[a][b] - 0.5 * da(b))
								   / (eigenvalues(a) - eigenvalues(b));
					}
 
			theta[0][1] = 0.5 * da(0);
			theta[1][0] = theta[0][1];
			xi[0][1] = (1.0 / 8.0) * fa(0);
			xi[1][0] = xi[0][1];
			eta = xi[2][0];
		 }
\endcode
For two equal eigenvalues a and c: \f$ \lambda_a = \lambda_c \neq \lambda_b \f$
\code
		 else if ( (std::fabs(eigenvalues(0) - eigenvalues(2)) < comp_tolerance)
				   &&
				   (!(std::fabs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance)) )
		 {
			eta = 0.0;
			for (unsigned int a = 0; a < dim; ++a)
				for (unsigned int b = 0; b < dim; ++b)
					if ( (a != b) && ((a == 1) || (b == 1)) )
					{
						theta[a][b] = (ea(a) - ea(b))
									  / (eigenvalues(a) - eigenvalues(b));
						xi[a][b] = (theta[a][b] - 0.5 * da(b))
								   / (eigenvalues(a) - eigenvalues(b));
					}
 
			theta[0][2] = 0.5 * da(0);
			theta[2][0] = theta[0][2];
			xi[0][2] = (1.0 / 8.0) * fa(0);
			xi[2][0] = xi[0][2];
			eta = xi[1][0];
		 }
\endcode
For two equal eigenvalues b and c: \f$ \lambda_b = \lambda_c \neq \lambda_a \f$
\code
		 else if ( (std::fabs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance)
				   &&
				   (!(std::fabs(eigenvalues(0) - eigenvalues(1)) < comp_tolerance)) )
		 {
			eta = 0.0;
			for (unsigned int a = 0; a < dim; ++a)
				for (unsigned int b = 0; b < dim; ++b)
					if ( (a != b) && ((a == 0) || (b == 0)) )
					{
						theta[a][b] = (ea(a) - ea(b))
									  / (eigenvalues(a) - eigenvalues(b));
						xi[a][b] = (theta[a][b] - 0.5 * da(b))
								   / (eigenvalues(a) - eigenvalues(b));
					}
 
			theta[1][2] = 0.5 * da(1);
			theta[2][1] = theta[1][2];
			xi[1][2] = (1.0 / 8.0) * fa(1);
			xi[2][1] = xi[1][2];
			eta = xi[0][1];
		 }
		 else
		 {
			deallog << "ln-space<< eigenvalues:0: " << eigenvalues[0] << std::endl;
			deallog << "ln-space<< eigenvalues:1: " << eigenvalues[1] << std::endl;
			deallog << "ln-space<< eigenvalues:2: " << eigenvalues[2] << std::endl;
			AssertThrow( false,
						 ExcMessage("ln-space<< Eigenvalue case not possible, check update_qph!") );
		 }
 
\endcode
Ensure that \a eta was initialised in one of the cases
\code
		 AssertThrow( (eta < 9999999), ExcMessage("Eta in update_qph not initialised") );
 
 
		/*
		 * 4. Lagrangian stresses and elasticity moduli
		 */
 
\endcode
Compute projection tensor P
\code
		 Tensor<4, dim> projection_tensor_P;
		 for (unsigned int a = 0; a < dim; ++a)
		 {
			projection_tensor_P += da(a) * (Tensor<4,dim> ) outer_product(eigenbasis[a],eigenbasis[a]);
			for (unsigned int b = 0; b < dim; ++b)
				if (b != a)
					projection_tensor_P += theta[a][b] * get_tensor_operator_G(eigenbasis[a],eigenbasis[b]);
		 }
 
\endcode
Check whether the projecton tensor is symmetric and store it into a \a SymmetricTensor
\code
		 AssertThrow( symmetry_check(projection_tensor_P), ExcMessage( "ln-space<< Projection tensor P is not symmetric") );
		 SymmetricTensor<4,dim> projection_tensor_P_sym = symmetrize(projection_tensor_P);
 
\endcode
Compute the double contraction of T and L
\code
		 Tensor<4, dim> projection_tensor_T_doublecon_L;
		 for (unsigned int a = 0; a < dim; ++a)
		 {
			projection_tensor_T_doublecon_L += fa(a)
											   * (stress_measure_T_sym * eigenbasis[a])
											   * (Tensor<4, dim> ) (outer_product(eigenbasis[a], eigenbasis[a]));
			for (unsigned int b = 0; b < dim; ++b)
				if (b != a)
				{
					projection_tensor_T_doublecon_L += 2.0 * xi[a][b]
													   * (
															get_tensor_operator_F_right( eigenbasis[a], eigenbasis[b], eigenbasis[b], stress_measure_T_sym )
														  + get_tensor_operator_F_left(  eigenbasis[a], eigenbasis[b], eigenbasis[b], stress_measure_T_sym )
														  + get_tensor_operator_F_right( eigenbasis[b], eigenbasis[b], eigenbasis[a], stress_measure_T_sym )
														  + get_tensor_operator_F_left(  eigenbasis[b], eigenbasis[b], eigenbasis[a], stress_measure_T_sym )
														  + get_tensor_operator_F_right( eigenbasis[b], eigenbasis[a], eigenbasis[b], stress_measure_T_sym )
														  + get_tensor_operator_F_left(  eigenbasis[b], eigenbasis[a], eigenbasis[b], stress_measure_T_sym )
														 );
 
					for (unsigned int c = 0; c < dim; ++c)
						if ( (c != a) && (c != b) )
						{
							projection_tensor_T_doublecon_L += 2.0 * eta
															   * (
																	  get_tensor_operator_F_right( eigenbasis[a], eigenbasis[b], eigenbasis[c], stress_measure_T_sym )
																	+ get_tensor_operator_F_left(  eigenbasis[b], eigenbasis[c], eigenbasis[a], stress_measure_T_sym )
																  );
						}
				}
		 }
 
\endcode
Check whether the tensor is symmetric and store it into a \a SymmetricTensor
\code
		 AssertThrow( symmetry_check(projection_tensor_T_doublecon_L),
					  ExcMessage("ln-space<< Projection tensor T:L is not symmetric") );
		 SymmetricTensor<4,dim> projection_tensor_T_doublecon_L_sym = symmetrize(projection_tensor_T_doublecon_L);
 
\endcode
Compute the retransformed values
\code
		 second_piola_stress_S = stress_measure_T_sym * projection_tensor_P_sym;
 
		 elasto_plastic_tangent = projection_tensor_P_sym * elasto_plastic_tangent * projection_tensor_P_sym // Note: we work on the input argument \a elasto_plastic_tangent
								  + projection_tensor_T_doublecon_L_sym;
	}
}
 
 
#endif // ln_space_H
\endcode

@section END The End

Hosted via GitHub according to https://goseeky.wordpress.com/2017/07/22/documentation-101-doxygen-with-github-pages/ \n
Design of the documentation inspired by the deal.ii tutorial programs.

@}
*/
