#ifndef ln_space_H
#define ln_space_H

/*
 *
 * Authors: Christian Burkhardt and Johannes Friedlein,
 * 	    FAU Erlangen-Nuremberg, 2019/2020
 *
 */

// @section includes Include Files
// The data type SymmetricTensor and some related operations, such as trace, symmetrize, deviator, ... for tensor calculus
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/exceptions.h>

// @todo Check whether the following three headers are needed at all
#include <iostream>
#include <fstream>
#include <cmath>

#include "functions.h"

#include "../MA-Code/auxiliary_functions/StandardTensors.h"
#include "../MA-Code/auxiliary_functions/tensor_operators.h"
#include "../MA-Code/auxiliary_functions/auxiliary_functions.h"
#include "../2D_axial-symmetry_plane-strain_dealii/handling_2D.h"

using namespace dealii;

/**
 * @note CONVENTION: We explicitly don't use the Lagrangian tangent defined as \f$ C = 2 \frac{\partial S}{\partial C} \f$,
 * but only the derivative \f$ \frac{\partial S}{\partial C} \f$. Similarly, additional tangent contributions as
 * \f$ \frac{\partial (\bullet)}{\partial C} \f$ are computed and provided without (!) the factor of two. Hence, make
 * sure that your assembly routine incorporates the factor of 2 for instance in the linearisation of RCG.
 */
template <int dim>
class ln_space
{
  public:
	ln_space();

	SymmetricTensor<2,dim> hencky_strain; // only dim components
	SymmetricTensor<2,3> hencky_strain_3D; // full 3D components

	// Contain dim components from the 3D quantities
	 SymmetricTensor<2,3> second_piola_stress_S;
	 SymmetricTensor<4,3> dS_dC_3D;

	bool processing_failed = false;

	void pre_ln ( /*input->*/ const Tensor<2,3> &F
				  /*output->*/ );

	void post_ln ( /*input->*/ SymmetricTensor<2,3> &stress_measure_T_sym, SymmetricTensor<4,3> &elasto_plastic_tangent
				   /*output->second_piola_stress_S, C*/ );

	SymmetricTensor<2,3> post_transform_strain ( /*input->*/ SymmetricTensor<2,3> &ln_tensor);
	SymmetricTensor<2,3> post_transform_stress ( /*input->*/ SymmetricTensor<2,3> &ln_tensor);

	SymmetricTensor<2,3> plastic_right_cauchy_green_AS (SymmetricTensor<2,3> plastic_hencky_strain);

  private:
  	 Vector<double> eigenvalues;
  	 std::vector< Tensor<1,3> > eigenvector;
  	 std::vector< SymmetricTensor<2,3> > eigenbasis;
  	 Vector<double> ea;
  	 Vector<double> da;
  	 Vector<double> fa;

  	SymmetricTensor<4,3> projection_tensor_P_sym;

 	const double comp_tolerance = 1e-8;
};


template <int dim>
ln_space<dim>::ln_space( )
:
eigenvalues(3),
eigenvector(3),
eigenbasis(3),
ea(3),
da(3),
fa(3)
{
}


// @section 3D Tranformation for 3D
/*
 * #################################################################### 3D ##############################################################
 */

// 2D and 3D
template<int dim>
void ln_space<dim>::pre_ln ( /*input->*/ const Tensor<2,3> &F /*output->hencky_strain, eigenvalues, eigenvector, eigenbasis, ea, da, fa*/ )
{
	// Following
	// "Algorithms for computation of stresses and elasticity moduli in terms of Seth–Hill’s family of generalized strain tensors"
	// by Miehe&Lambrecht \n
	// Table I. Algorithm A
	/*
	 * 1. Eigenvalues, eigenvalue bases and diagonal functions:
	 */

	// Get the symmetric right cauchy green tensor
	 SymmetricTensor<2,3> right_cauchy_green_sym = symmetrize( transpose(F) * F);//Physics::Elasticity::Kinematics::C(F);

	// Compute Eigenvalues, Eigenvectors and Eigenbasis
	 {
		 // array of pairs with EW double, EV tensor, 3 array entries
		  std::array< std::pair< double, Tensor<1, dim> >, 3 > eigenlist;

		 // Get Eigenvalues and Eigenvectors from the deal.ii function \a eigenvectors(*)
		 // determine the eigenvalues and eigenvectors all at once and save both entries of the pair
		  try
		  {
			  eigenlist = eigenvectors(right_cauchy_green_sym);
		  }
		  catch (std::exception &exc)
		  {
			  std::cerr << std::endl
						<< "----------------------------------------------------"
						<< std::endl;
			  std::cerr << "GG<< pre_ln(*) failed to find eigenvalues or eigenvectors for the given tensor. "
						   "Aborting the current load step and trying a restart with a reduced load increment." << std::endl
						<< "----------------------------------------------------"
						<< std::endl;

			  processing_failed = true;

			  return;
		  }

		// Store the values in a nicer format
		 for (unsigned int i = 0; i < 3; ++i)
		 {
			eigenvalues[i] = eigenlist[i].first;
			eigenvector[i] = eigenlist[i].second;
		 }

		// The deal.ii function \a eigenvectors return the EWe in descending order, but in Miehe et al. the C33 EW
		// shall belong to lambda_3, hence we search for the EW that equals C33
		// and move this to the end of the list of eigenvalues and eigenvectors
		 // ToDo-optimize: if we find the EW at i=2, we can leave it there, hence the loop could be limited to (i<2)
		// @todo-optimize: Check the use of the std::swap function
		if ( dim==2 )
			for (unsigned int i = 0; i < 3; ++i)
				 if ( std::abs(eigenvalues[i]-right_cauchy_green_sym[2][2]) < comp_tolerance )
				 {
					 double tmp_EW = eigenvalues[2];
					 eigenvalues[2] = eigenvalues[i]; // truely lambda_3
					 eigenvalues[i] = tmp_EW;

					 Tensor<1,3> tmp_EV = eigenvector[2];
					 eigenvector[2] = eigenvector[i];
					 eigenvector[i] = tmp_EV;
				 }

		// Check if the found eigenvectors are perpendicular to each other
//		if ((std::fabs(eigenvalues(0) - 1) > 1e-10)
//			&& (std::fabs(eigenvalues(1) - 1) > 1e-10)
//			&& (std::fabs(eigenvalues(2) - 1) > 1e-10)) {
//			for (unsigned int i = 0; i < 3; ++i) {
//				for (unsigned int j = i + 1; j < 3; ++j) {
//					Assert( (std::fabs(eigenvector[i] * eigenvector[j]) < 1e-12),
//								 ExcMessage("ln-space<< Eigenvectors are not perpendicular to each other.") );
//				}
//			}
//		}

		// Compute eigenbasis Ma: eigenbasis = eigenvector \otimes eigenvector
		 for (unsigned int i = 0; i < 3; ++i)
		 {
			eigenbasis[i] = outer_product_sym(eigenvector[i]);
			Assert( eigenvalues(i) >= 0.0,
						 ExcMessage("ln-space<< Eigenvalue is negativ. Check update_qph.") );
		 }
	 }

	// Compute diagonal function \a ea and its first and second derivate \a da and \a fa \n
	 for (unsigned int i = 0; i < 3; ++i)
	 {
		ea(i) = 0.5 * std::log( std::abs(eigenvalues(i)) );	// diagonal function
		da(i) = std::pow(eigenvalues(i), -1.0);				// first derivative of diagonal function ea
		fa(i) = -2.0 * std::pow(eigenvalues(i), -2.0);			// second derivative of diagonal function ea
		Assert( ea(i) == ea(i),
					 ExcMessage( "ln-space<< Ea is nan due to logarithm of negativ eigenvalue. Check update_qph.") );
		Assert( da(i) > 0.0,
					 ExcMessage( "ln-space<< First derivative da of diagonal function is "+std::to_string(da(i))+" < 0.0 ."
							 	 "Check update_qph.") );
	 }

	// Compute the Hencky strain
	 for (unsigned int a = 0; a < 3; ++a)
		hencky_strain_3D += ea(a) * eigenbasis[a];

	 hencky_strain = extract_dim<dim> (hencky_strain_3D);

	// Output-> SymmetricTensor<2, dim> hencky_strain, Vector<double> ea, da, fa,
	//			std::vector<Tensor<1, dim>> eigenvector, Vector<double> eigenvalues,
	//			std::vector< SymmetricTensor<2, dim> > eigenbasis
}



// 3D
template<>
void ln_space<3>::post_ln ( /*output->*/ SymmetricTensor<2,3> &stress_measure_T_sym, SymmetricTensor<4,3> &elasto_plastic_tangent )
{
	/*
	 * 3. Set up coefficients \a theta, \a xi and \a eta
	 */

	// Compute the coefficients based on the eigenvalues, eigenvectors and ea,da,fa
	 Tensor<2,3> theta;
	 Tensor<2,3> xi;
	 double eta = 999999999.0;

	// For three different eigenvalues \f$ \lambda_a \neq \lambda_b \neq \lambda_c \f$
	 if (
			 ( !(std::fabs(eigenvalues(0) - eigenvalues(1)) < comp_tolerance) )
			 &&
			 ( !(std::fabs(eigenvalues(0) - eigenvalues(2)) < comp_tolerance) )
			 &&
			 ( !(std::fabs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance) )
		 )
	 {
		eta = 0.0;
		for (unsigned int a = 0; a < 3; ++a)
			for (unsigned int b = 0; b < 3; ++b)
				if (a != b)
				{
					theta[a][b] = (ea(a) - ea(b))
								  / (eigenvalues(a) - eigenvalues(b));
					xi[a][b] = (theta[a][b] - 0.5 * da(b))
							   / (eigenvalues(a) - eigenvalues(b));

					for (unsigned int c = 0; c < 3; ++c)
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
		for (unsigned int a = 0; a < 3; ++a)
		{
			for (unsigned int b = 0; b < 3; ++b)
				if (a != b)
				{
					theta[a][b] = 0.5 * da(0);
					xi[a][b] = (1.0 / 8.0) * fa(0);
				}
		}
		eta = (1.0 / 8.0) * fa(0);
	 }

	// For two equal eigenvalues a and b: \f$ \lambda_a = \lambda_b \neq \lambda_c \f$
	 else if ( (std::fabs(eigenvalues(0) - eigenvalues(1)) < comp_tolerance)
			   &&
			   ( !(std::fabs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance) ) )
	 {
		eta = 0.0;
		for (unsigned int a = 0; a < 3; ++a)
			for (unsigned int b = 0; b < 3; ++b)
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
	// For two equal eigenvalues a and c: \f$ \lambda_a = \lambda_c \neq \lambda_b \f$
	 else if ( (std::fabs(eigenvalues(0) - eigenvalues(2)) < comp_tolerance)
			   &&
			   (!(std::fabs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance)) )
	 {
		eta = 0.0;
		for (unsigned int a = 0; a < 3; ++a)
			for (unsigned int b = 0; b < 3; ++b)
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
	// For two equal eigenvalues b and c: \f$ \lambda_b = \lambda_c \neq \lambda_a \f$
	 else if ( (std::fabs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance)
			   &&
			   (!(std::fabs(eigenvalues(0) - eigenvalues(1)) < comp_tolerance)) )
	 {
		eta = 0.0;
		for (unsigned int a = 0; a < 3; ++a)
			for (unsigned int b = 0; b < 3; ++b)
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

	// Ensure that \a eta was initialised in one of the cases
	 //AssertThrow( (eta < 9999999), ExcMessage("Eta in update_qph not initialised") );


	/*
	 * 4. Lagrangian stresses and elasticity moduli
	 */

	 std::vector< SymmetricTensor<4,3> > Ma_x_Ma (3);
	// Compute projection tensor P
	 Tensor<4,3> projection_tensor_P;
	 for (unsigned int a = 0; a < 3; ++a)
	 {
		Ma_x_Ma[a] = outer_product_sym(eigenbasis[a]);
		projection_tensor_P += da(a) * (Tensor<4,3> ) Ma_x_Ma[a];
		for (unsigned int b = 0; b < 3; ++b)
			if (b != a)
				projection_tensor_P += theta[a][b] * get_tensor_operator_G(eigenbasis[a],eigenbasis[b]);
	 }

	// Check whether the projecton tensor is symmetric and store it into a \a SymmetricTensor
	 Assert( symmetry_check(projection_tensor_P), ExcMessage( "ln-space<< Projection tensor P is not symmetric") );
	 projection_tensor_P_sym = symmetrize(projection_tensor_P);

	// Compute the double contraction of T and L
	 Tensor<4,3> projection_tensor_T_doublecon_L;
	 for (unsigned int a = 0; a < 3; ++a)
	 {
		projection_tensor_T_doublecon_L += fa(a)
										   * (stress_measure_T_sym * eigenbasis[a])
										   * (Tensor<4,3> ) Ma_x_Ma[a];
		for (unsigned int b = 0; b < 3; ++b)
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

				for (unsigned int c = 0; c < 3; ++c)
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

	// Check whether the tensor is symmetric and store it into a \a SymmetricTensor
	 Assert( symmetry_check(projection_tensor_T_doublecon_L),
				  ExcMessage("ln-space<< Projection tensor T:L is not symmetric") );
	 SymmetricTensor<4,3> projection_tensor_T_doublecon_L_sym = symmetrize(projection_tensor_T_doublecon_L);

	// Compute the retransformed values
	 second_piola_stress_S = stress_measure_T_sym * projection_tensor_P_sym;

	// Factor 0.5 to get the sole dS_dC derivative
	 dS_dC_3D = 0.5 * ( projection_tensor_P_sym * elasto_plastic_tangent * projection_tensor_P_sym
							  + projection_tensor_T_doublecon_L_sym );
}


// @section 2D Tranformation for 2D
/*
 * #################################################################### 2D ##############################################################
 */
// One might think, that calling the 3D transformation with an expanded deformation gradient removes the need for a seperate 2D
// transformation. But the 3D code seems to become singular e.g. for a plane strain deformation gradient.
// @todo Does this also happen for the more general axial symmetry

// 2D
template<>
void ln_space<2>::post_ln ( /*input->*/ SymmetricTensor<2,3> &stress_measure_T_sym, SymmetricTensor<4,3> &elasto_plastic_tangent )
{
	/*
	 * 3. Set up coefficients \a theta, \a xi and \a eta
	 */

	// Compute the coefficients based on the eigenvalues, eigenvectors and ea,da,fa
	 Vector<double> gamma (3);
	 Vector<double> kappa (2);
	 double beta = 999999999.0;

	// For two different eigenvalues \f$ \lambda_a \neq \lambda_b \f$
	 if (std::fabs(eigenvalues(0) - eigenvalues(1)) > comp_tolerance)
	 {
		beta = 2. * (ea[0]-ea[1]) / (eigenvalues(0) - eigenvalues(1));
		for (unsigned int a = 0; a < 3; ++a)
			gamma[a] = da[a] - beta;
		for (unsigned int alpha = 0; alpha < 2; ++alpha)
			kappa[alpha] = 2. * gamma[alpha] / (eigenvalues(0) - eigenvalues(1));
	 }
	// For two equal eigenvalues a and b: \f$ \lambda_a = \lambda_b \f$
	 else if (std::fabs(eigenvalues(0) - eigenvalues(1)) <= comp_tolerance)
	 {
		beta = da[0];
		for (unsigned int a = 0; a < 3; ++a)
			gamma[a] = da[a] - da[0];
		for (unsigned int alpha = 0; alpha < 2; ++alpha)
			kappa[alpha] = (1.5 - (alpha+1)) * fa[0]; // index alpha starts at 0, hence (alpha+1)
	 }
	 else
	 {
		deallog << "ln-space<< eigenvalues:0: " << eigenvalues[0] << std::endl;
		deallog << "ln-space<< eigenvalues:1: " << eigenvalues[1] << std::endl;
		deallog << "ln-space<< eigenvalues:2: " << eigenvalues[2] << std::endl;
		AssertThrow( false,
					 ExcMessage("ln-space<< Eigenvalue case not possible, check update_qph!") );
	 }

	// Ensure that \a eta was initialised in one of the cases
	 Assert( (beta < 9999999), ExcMessage("Beta in update_qph not initialised") );


	/*
	 * 4. Lagrangian stresses and elasticity moduli
	 */

	// zeta_a = T : M_a
	 Vector<double> zeta_a (3);
	 for ( unsigned int a=0; a<3; ++a )
		 zeta_a(a) = stress_measure_T_sym * eigenbasis[a];

	// Compute projection tensor P
	 Tensor<4,3> projection_tensor_P = beta * Tensor<4,3> (StandardTensors::II<3>());

	 std::vector< SymmetricTensor<4,3> > Ma_x_Ma (3);
	 // ToDo-optimize: check whether storing (M_a x M_a) saves some computation time, because we need it three times
	 for (unsigned int a = 0; a < 3; ++a)
	 {
		 Ma_x_Ma[a] = outer_product_sym(eigenbasis[a]);
		 projection_tensor_P += gamma[a] * (Tensor<4,3>) Ma_x_Ma[a];
	 }

	 projection_tensor_P_sym = symmetrize(projection_tensor_P);

	// Compute the double contraction of T and L
	 // ToDo-optimize: merge the for-loops
	 Tensor<4,3> projection_tensor_T_doublecon_L;
	 for (unsigned int a = 0; a < 3; ++a)
	 {
		projection_tensor_T_doublecon_L += fa(a)
										   * zeta_a[a]
										   * (Tensor<4,3>) Ma_x_Ma[a];
	 }
	 for (unsigned int alpha = 0; alpha < 2; ++alpha)
	 {
		projection_tensor_T_doublecon_L += std::pow(-1, (alpha+1)+1)
										   * kappa[alpha]
										   * (
											   zeta_a[alpha] * Tensor<4,3> (StandardTensors::II<3>())
											   + outer_product_sym(stress_measure_T_sym,eigenbasis[alpha])
											 );
	 }
	 for (unsigned int a = 0; a < 3; ++a)
		 for (unsigned int alpha = 0; alpha < 2; ++alpha)
		 {
			 projection_tensor_T_doublecon_L
				  += std::pow(-1,alpha+1) * kappa[alpha]
					 * (
						   zeta_a[alpha] * Tensor<4,3> (Ma_x_Ma[a])
						   + zeta_a[a] * outer_product_sym(eigenbasis[alpha],eigenbasis[a])
						);
		 }

	// Check whether the tensor is symmetric and store it into a \a SymmetricTensor
	 Assert( symmetry_check(projection_tensor_T_doublecon_L),
				  ExcMessage("ln-space<< Projection tensor T:L is not symmetric") );
	 SymmetricTensor<4,3> projection_tensor_T_doublecon_L_sym = symmetrize(projection_tensor_T_doublecon_L);

	// Compute the retransformed values
	 second_piola_stress_S = stress_measure_T_sym * projection_tensor_P_sym;

	 dS_dC_3D = 0.5 * ( projection_tensor_P_sym * elasto_plastic_tangent * projection_tensor_P_sym
			 	 	 	+ projection_tensor_T_doublecon_L_sym );
}


/*
 * Transform the given second order tensor from the ln-space into the real world
 * works for e.g. @todo fin
 * \f$ \frac{\partial (\bullet)}{\partial C} \f$
 */
template<int dim>
SymmetricTensor<2,3> ln_space<dim>::post_transform_strain ( /*input->*/ SymmetricTensor<2,3> &ln_tensor)
{
	// The factor 0.5 ensures that we get \f$ \frac{\partial (\bullet)}{\partial C} \f$ and not "2*"
	return 0.5 * ln_tensor * projection_tensor_P_sym;
}


/*
 * Transform the given second order tensor from the ln-space into the real world
 * works for e.g. @todo fin
 * \f$ \frac{\partial T}{\partial (\bullet)} \f$
 */
template<int dim>
SymmetricTensor<2,3> ln_space<dim>::post_transform_stress ( /*input->*/ SymmetricTensor<2,3> &ln_tensor)
{
	return 1. * ln_tensor * projection_tensor_P_sym;
}

// @todo-assure Still not verified, How?
template <int dim>
SymmetricTensor<2,3> ln_space<dim>::plastic_right_cauchy_green_AS (SymmetricTensor<2,3> plastic_hencky_strain)
{
//	AssertThrow( false, ExcMessage("plastic_right_cauchy_green_AS<< Sorry. Even though the algorithm "
//								   "to compute the plastic RCG tensor was tested, the implementation into this "
//								   "ln-space class has not been tested at all. So either you have enough faith to "
//								   "simply remove this AssertThrow in the code or you do some testing to validate the function yourself."));

	// Compute the eigenvalues and eigenvectors of the plastic hencky strain
	 Vector<double> eigenvalues_pl(3);
	 std::vector< Tensor<1,3> > eigenvector_pl(3);
	 for (unsigned int i = 0; i < 3; ++i)
	 {
		eigenvalues_pl[i] = eigenvectors(plastic_hencky_strain)[i].first;
		eigenvector_pl[i] = eigenvectors(plastic_hencky_strain)[i].second;
	 }

	// Check if the found eigenvectors are perpendicular to each other
	// @todo Change the \a false to a criterion detecting whether the code is run in debug or release mode
	 if ( false /*no debugging*/)
	 {
		if ((fabs(eigenvalues_pl(0) - 1) > 1e-10)
				&& (fabs(eigenvalues_pl(1) - 1) > 1e-10)
				&& (fabs(eigenvalues_pl(2) - 1) > 1e-10))
			for (unsigned int i = 0; i < 3; ++i)
				for (unsigned int j = i + 1; j < 3; ++j)
					Assert( (fabs(eigenvector[i] * eigenvector[j]) < 1e-12),
							ExcMessage("Eigenvectors are not perpendicular to each other") );
	 }

	// Compute eigenbasis
	 std::vector< SymmetricTensor<2,3> > eigenbasis_pl(3);
	 for (unsigned int i = 0; i < 3; ++i)
		eigenbasis_pl[i] = outer_product_sym( eigenvector_pl[i] );

	// Compute the values of \a ea
	 Vector<double> ea(3);
	 for (unsigned int i = 0; i < 3; ++i)
		ea(i) = exp(2.0* eigenvalues_pl(i));			// diagonal

	// Finally compute the plastic right Cauchy-Green tensor
	 SymmetricTensor<2,3> plastic_right_cauchy_green;
	 for (unsigned int a = 0; a < 3; ++a)
		plastic_right_cauchy_green += ea(a) * eigenbasis_pl[a];

	return plastic_right_cauchy_green;
}


#endif // ln_space_H
