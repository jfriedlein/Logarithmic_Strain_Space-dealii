#ifndef ln_space_H
#define ln_space_H

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
#include "../MA-Code/handling_2D.h"

using namespace dealii;

template <int dim>
class ln_space
{
  public:
	ln_space();

	SymmetricTensor<2,dim> hencky_strain; // only dim components
	SymmetricTensor<2,3> hencky_strain_3D; // full 3D components

	// Contain dim components from the 3D quantities
	 SymmetricTensor<2,dim> second_piola_stress_S;
	 SymmetricTensor<4,dim> C;
	 SymmetricTensor<4,3> C_3D;

	void pre_ln ( /*input->*/ Tensor<2,3> &F
				  /*output->*/ );

//	void pre_ln ( /*input->*/ Tensor<2,2> &F_2D
//				  /*output->*/ );

	void post_ln ( /*input->*/ SymmetricTensor<2,3> &stress_measure_T_sym, SymmetricTensor<4,3> &elasto_plastic_tangent
				   /*output->second_piola_stress_S, C*/ );

	// TESTING:
	double eps_p_22 = 0.;

  private:
  	 Vector<double> eigenvalues;
  	 std::vector< Tensor<1,3> > eigenvector;
  	 std::vector< SymmetricTensor<2,3> > eigenbasis;
  	 Vector<double> ea;
  	 Vector<double> da;
  	 Vector<double> fa;

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

// 3D
template<int dim>
void ln_space<dim>::pre_ln ( /*input->*/ Tensor<2,3> &F /*output->hencky_strain, eigenvalues, eigenvector, eigenbasis, ea, da, fa*/ )
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
		// Get Eigenvalues and Eigenvectors from the deal.ii function \a eigenvectors(*)
		for (unsigned int i = 0; i < 3; ++i) {
			eigenvalues[i] = eigenvectors(right_cauchy_green_sym)[i].first;
			eigenvector[i] = eigenvectors(right_cauchy_green_sym)[i].second;
		}

		// The deal.ii function \a eigenvectors return the EWe in descending order, but in Miehe et al. the C33 EW
		// shall belong to lambda_3, hence we search for the EW that equals C33
		// and move this to the end of the list of eigenvalues and eigenvectors
		 // ToDo-optimize: if we find the EW at i=2, we can leave it there, hence the loop could be limited to (i<2)
		 // ToDo-optimize: Besides the following loop everything else is the same as for 3D maybe merge this
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
	// ToDo-clean: Only added the deformation gradient to differentiate between the 3D and truely 2D functions

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
	 SymmetricTensor<4,3> projection_tensor_P_sym = symmetrize(projection_tensor_P);

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

	 C = projection_tensor_P_sym * elasto_plastic_tangent * projection_tensor_P_sym
							  + projection_tensor_T_doublecon_L_sym;
}


//// 2D
//template<>
//void ln_space<2>::pre_ln ( /*input->*/ Tensor<2,2> &F_2D
//			  	  	  	   /*output->*/ )
//{
//	// Following "Algorithms for computation of stresses and elasticity moduli in terms of Seth–Hill’s family of generalized strain tensors" by Miehe&Lambrecht \n
//	// Table I. Algorithm A
//	/*
//	 * 1. Eigenvalues, eigenvalue bases and diagonal functions:
//	 */
//
//	// Expand the 2D DefoGradient to pseudo 3D
//	 Tensor<2,3> F;
//	 for ( unsigned int i=0; i<2; ++i )
//	  for ( unsigned int j=0; j<2; ++j )
//		  F[i][j] = F_2D[i][j];
//
//	 // ToDo-assure: @q: We set the third component to 1, because the gradient of u_3/x_3
//	 // shall be zero when there is no deformation in the out-of plane direction (plane strain state)
//	  F[2][2] = 1.;
//
//	  if ( true /*axisymmetric*/){
//			F[2][2]/*theta*/ += eps_p_22; // abusing the zz-component of eps_p_n for u/r
//			eps_p_22 = 0.;
//	  }
//
//	// Get the symmetric right cauchy green tensor and expand it to pseudo 3D
//	 SymmetricTensor<2,3> right_cauchy_green_sym = symmetrize( contract<1,0>(transpose(F),F) );
//
//	// Compute Eigenvalues, Eigenvectors and Eigenbasis
//	 {
//		// Get Eigenvalues and Eigenvectors from the deal.ii function \a eigenvectors(*)
//		 for (unsigned int i = 0; i < 3; ++i) {
//			eigenvalues[i] = eigenvectors(right_cauchy_green_sym)[i].first;
//			eigenvector[i] = eigenvectors(right_cauchy_green_sym)[i].second;
//		 }
//
//		// The deal.ii function \a eigenvectors return the EWe in descending order, but in Miehe et al. the C33 EW
//		// shall belong to lambda_3, hence we search for the EW that equals C33
//		// and move this to the end of the list of eigenvalues and eigenvectors
//		 // ToDo-optimize: if we find the EW at i=2, we can leave it there, hence the loop could be limited to (i<2)
//		 // ToDo-optimize: Besides the following loop everything else is the same as for 3D maybe merge this
//		 for (unsigned int i = 0; i < 3; ++i)
//			 if ( std::abs(eigenvalues[i]-right_cauchy_green_sym[2][2]) < comp_tolerance )
//			 {
//				 double tmp_EW = eigenvalues[2];
//				 eigenvalues[2] = eigenvalues[i]; // truely lambda_3
//				 eigenvalues[i] = tmp_EW;
//
//				 Tensor<1,3> tmp_EV = eigenvector[2];
//				 eigenvector[2] = eigenvector[i];
//				 eigenvector[i] = tmp_EV;
//			 }
//
//		// Check if the found eigenvectors are perpendicular to each other
////		 if ((std::fabs(eigenvalues(0) - 1) > 1e-10)
////			&& (std::fabs(eigenvalues(1) - 1) > 1e-10)
////			&& (std::fabs(eigenvalues(2) - 1) > 1e-10)) {
////			for (unsigned int i = 0; i < 3; ++i) {
////				for (unsigned int j = i + 1; j < 3; ++j) {
////					AssertThrow( (std::fabs(eigenvector[i] * eigenvector[j]) < 1e-12),
////								 ExcMessage("ln-space<< Eigenvectors are not perpendicular to each other.") );
////				}
////			}
////		 }
//
//		// Compute eigenbasis Ma: eigenbasis = eigenvector \otimes eigenvector
//		 for (unsigned int i = 0; i < 3; ++i)
//		 {
//			eigenbasis[i] = outer_product_sym(eigenvector[i]);
//			Assert( eigenvalues(i) >= 0.0,
//						 ExcMessage("ln-space<< Eigenvalue is negativ. Check update_qph.") );
//		 }
//	 }
//
//	// Compute diagonal function \a ea and its first and second derivate \a da and \a fa
//	 for (unsigned int i = 0; i < 3; ++i)
//	 {
//		ea(i) = 0.5 * std::log( std::abs(eigenvalues(i)) );	// diagonal function
//		da(i) = std::pow(eigenvalues(i), -1.0);				// first derivative of diagonal function ea
//		fa(i) = -2.0 * std::pow(eigenvalues(i), -2.0);			// second derivative of diagonal function ea
//		Assert( ea(i) == ea(i),
//					 ExcMessage( "ln-space<< Ea is nan due to logarithm of negativ eigenvalue. Check update_qph.") );
//		Assert( da(i) > 0.0,
//					 ExcMessage( "ln-space<< First derivative da of diagonal function is "+std::to_string(da(i))+" < 0.0 ."
//							 	 "Check update_qph.") );
//	 }
//
//	// Compute the Hencky strain
//	 for (unsigned int a = 0; a < 3; ++a)
//		 hencky_strain_3D += ea(a) * eigenbasis[a];
//
//	 hencky_strain = extract_dim<2> (hencky_strain_3D);
//
//	// Output-> SymmetricTensor<2, dim> hencky_strain, Vector<double> ea, da, fa,
//	//			std::vector<Tensor<1, dim>> eigenvector, Vector<double> eigenvalues,
//	//			std::vector< SymmetricTensor<2, dim> > eigenbasis
//}



// @section 2D Tranformation for 2D
/*
 * #################################################################### 2D ##############################################################
 */
// One might think, that calling the 3D transformation with an expanded deformation gradient removes the need for a seperate 2D
// transformation. But the 3D code seems to become singular for a plane strain deformation gradient.

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

	 SymmetricTensor<4,3> projection_tensor_P_sym = symmetrize(projection_tensor_P);

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
	 {
		 SymmetricTensor<2,3> stress_S_3D = stress_measure_T_sym * projection_tensor_P_sym;
		 second_piola_stress_S = extract_dim<2> ( stress_S_3D );
	 }
	 {
		 C_3D = projection_tensor_P_sym * elasto_plastic_tangent * projection_tensor_P_sym
				+ projection_tensor_T_doublecon_L_sym;

		 C = extract_dim<2> ( C_3D );
	 }
}


#endif // ln_space_H
