#ifndef CONSTITUTIVE_LAW_THERMO_ELASTO_PLASTIC
#define CONSTITUTIVE_LAW__THERMO_ELASTO_PLASTIC


/*
---------------------------------------------------
---------------------------------------------------
---------------------------------------------------
Constitutive law for geometrically linear and nonlinear
thermo elasto plastic material behaviour
using the logarithmic strain space
in an incremental formulation
iwtm96
--------------------------------------------------
---------------------------------------------------
---------------------------------------------------
*/


#include <map>
#include <vector>
#include <math.h>


#include "parameters/AllParameters.h"
#include <deal.II/base/exceptions.h>
#include "constitutive_laws/base/Constitutive_Law.h"
#include "assembly/MISCFunctions.h"
#include "materials/mechanical/base/MechanicalMaterial.h"
#include "materials/mechanical/derivatives/MechanicalMaterial_TiAl6V4.h"
#include "assembly/qph/PointHistoryMechanicalLargeStrainIncremental.h"


namespace Constitutive_Laws
{


    template<int dim>
    class Thermo_Elasto_Plastic : public Constitutive_Law<dim>
    {
    public:
        Thermo_Elasto_Plastic(const Parameters::AllParameters<dim> &prm);
        ~Thermo_Elasto_Plastic() = default;

        inline void return_stress_tangent_int_variables_for_small_strain                                        (unsigned int material_id,
                                                                                                                 SymmetricTensor<2,dim> e_n1,
                                                                                                                 double alpha_n,
                                                                                                                 SymmetricTensor<2,dim> e_n,
                                                                                                                 SymmetricTensor<2,dim> ep_n,
                                                                                                                 SymmetricTensor<2,dim> eth_n,
                                                                                                                 SymmetricTensor<2,dim> cauchy_stress_n,
                                                                                                                 double db_temperature_previous_time_step,
                                                                                                                 unsigned int int_material_id_previous_time_step,
                                                                                                                 double db_current_temperature,
                                                                                                                 double db_current_time_step_length,
                                                                                                                 SymmetricTensor<2,dim> &cauchy_stress,
                                                                                                                 SymmetricTensor<4,dim> &elasto_plastic_tangent_ss,
                                                                                                                 QPH::internal_variables_incremental<dim> &updated_int_variables);


        inline void return_stress_tangent_int_variables_for_large_strain                                         (unsigned int material_id,
                                                                                                                  Tensor<2, dim> deformation_gradient,
                                                                                                                  double alpha_n,
                                                                                                                  SymmetricTensor<2, dim> e_n, // total strain n
                                                                                                                  SymmetricTensor<2, dim> ep_n, // plastic strain n
                                                                                                                  SymmetricTensor<2, dim> eth_n, // thermal strain n
                                                                                                                  SymmetricTensor<2, dim> cauchy_stress_n,
                                                                                                                  double db_temperature_previous_time_step,
                                                                                                                  unsigned int int_material_id_previous_time_step,
                                                                                                                  double db_current_temperature,
                                                                                                                  double db_current_time_step_length,
                                                                                                                  SymmetricTensor<2,dim> &second_piola_stress,
                                                                                                                  SymmetricTensor<4,dim> &elasto_plastic_tangent,
                                                                                                                  QPH::internal_variables_incremental<dim> &updated_int_variables);

    protected:

    private:

        MechanicalMaterials::MechanicalMaterial<dim>*                                                           ptr_mech_material;

        inline SymmetricTensor<2,dim> return_delta_e_th						                            		(unsigned int material_id,
                                                                                                                    unsigned int int_material_id_previous_time_step,
                                                                                                                    double db_temperature_previous_time_step,
                                                                                                                    double db_current_temperature);

        inline double return_delta_lambda								                                        (unsigned int material_id,
                                                                                                                  double yield_function,
                                                                                                                  double db_delta_time_step,
                                                                                                                  double db_current_temperature);

        // yield function for now Von Mises
        inline double return_yield_function		                                        						(unsigned int material_id,
                                                                                                                    SymmetricTensor<2,dim> sigma_dev_trial,
                                                                                                                    double hardening_stress_R_trial,
                                                                                                                    double db_current_temperature);

        inline double return_hardening_stress_R_trial			                                        		(unsigned int material_id_n,
                                                                                                                  unsigned int material_id_n1,
                                                                                                                  double alpha_n,
                                                                                                                  double alpha_n1,
                                                                                                                  double db_T_n,
                                                                                                                  double db_current_temperature);

        inline SymmetricTensor<4,dim> return_elastic_stress_strain_tensor	                                    (unsigned int material_id,
                                                                                                                  double db_current_temperature);

        inline SymmetricTensor<4,dim> return_elasto_plastic_stress_strain_tensor	                            (unsigned int material_id,
                                                                                                                 double sigma_dev_trial_norm,
                                                                                                                 SymmetricTensor<2,dim> normal,
                                                                                                                 double delta_lambda,
                                                                                                                 double db_delta_timestep,
                                                                                                                 double db_current_temperature);

        inline SymmetricTensor<2,dim> return_delta_cauchy_stress_deviator		                            	(unsigned int material_id,
                                                                                                                  SymmetricTensor<2,dim> e_delta,
                                                                                                                  SymmetricTensor<2,dim> ep_delta,
                                                                                                                  SymmetricTensor<2,dim> e_n,
                                                                                                                  SymmetricTensor<2,dim> ep_n,
                                                                                                                  double db_T_delta,
                                                                                                                  double db_current_temperature);

        inline SymmetricTensor<2,dim> return_delta_cauchy_stress_vol			                            	(unsigned int material_id,
                                                                                                                 SymmetricTensor<2,dim> emech_delta,
                                                                                                                 SymmetricTensor<2,dim> emech_n,
                                                                                                                 double db_T_delta,
                                                                                                                 double db_current_temperature);

        const SymmetricTensor<2,dim> I = unit_symmetric_tensor<dim>();
        const SymmetricTensor<4,dim> IxI = outer_product(I, I);
        const SymmetricTensor<4,dim> deviatoric_identity = deviator_tensor<dim>();

        const unsigned int material_id_powder;
    };



    template<int dim>
    inline void Thermo_Elasto_Plastic<dim>::return_stress_tangent_int_variables_for_small_strain(
            unsigned int material_id,
            SymmetricTensor<2,dim> e_n1,
            double alpha_n,
            SymmetricTensor<2,dim> e_n,
            SymmetricTensor<2,dim> ep_n,
            SymmetricTensor<2,dim> eth_n,
            SymmetricTensor<2,dim> cauchy_stress_n,
            double db_temperature_previous_time_step,
            unsigned int int_material_id_previous_time_step,
            double db_current_temperature,
            double db_current_time_step_length,
            SymmetricTensor<2,dim> &cauchy_stress,
            SymmetricTensor<4,dim> &elasto_plastic_tangent_ss,
            QPH::internal_variables_incremental<dim> &updated_int_variables){

        double T_delta = db_current_temperature - db_temperature_previous_time_step;
        SymmetricTensor<2, dim> ep_delta;
        SymmetricTensor<2, dim> eth_delta = return_delta_e_th(material_id,
                                                              int_material_id_previous_time_step,
                                                              db_temperature_previous_time_step,
                                                              db_current_temperature);
        SymmetricTensor<2, dim> eth_n1 = eth_n + eth_delta;

        SymmetricTensor<2, dim> e_delta = e_n1 - e_n;
        SymmetricTensor<2, dim> sigma_dev_trial;


        if(cauchy_stress_n.norm() >= 1e-10){
            sigma_dev_trial = deviator(cauchy_stress_n)
                              + return_delta_cauchy_stress_deviator(material_id,
                                                                    e_delta,
                                                                    ep_delta,
                                                                    e_n,
                                                                    ep_n,
                                                                    T_delta,
                                                                    db_current_temperature);

        }
        else{
            sigma_dev_trial = return_delta_cauchy_stress_deviator(material_id,
                                                                  e_delta,
                                                                  ep_delta,
                                                                  e_n,
                                                                  ep_n,
                                                                  0.0,
                                                                  db_current_temperature);
        }


        SymmetricTensor<2, dim> emech_n = e_n - eth_n /* - ep_n */ ; // the volumetric part of the plastic strains is suppost to be 0 so the " - ep_n" part is not necessary
        SymmetricTensor<2, dim> emech_delta = e_delta - eth_delta; // the volumetric part of the plastic strains is suppost to be 0 so the " - ep_delta" part is not necessary

        SymmetricTensor<2, dim> sigma_vol;
        if(cauchy_stress_n.norm() >= 1e-10){
            sigma_vol = return_volumetric_tensor(cauchy_stress_n)
                        + return_delta_cauchy_stress_vol(material_id,
                                                         emech_delta,
                                                         emech_n,
                                                         T_delta,
                                                         db_current_temperature);
        }
        else{
            sigma_vol = return_delta_cauchy_stress_vol(material_id,
                                                       emech_delta,
                                                       emech_n,
                                                       0.0,
                                                       db_current_temperature);
        }

        // compute trial hardening stress
        double R_trial = return_hardening_stress_R_trial(int_material_id_previous_time_step,
                                                         material_id,
                                                         alpha_n,
                                                         alpha_n,
                                                         db_temperature_previous_time_step,
                                                         db_current_temperature);

        //compute the yield function (von Mises)
        double phi_trial = return_yield_function(material_id,
                                                 sigma_dev_trial,
                                                 R_trial,
                                                 db_current_temperature);


        if ( phi_trial <= 0  /* elastic case */   ) {
            cauchy_stress = (sigma_vol + sigma_dev_trial);
            elasto_plastic_tangent_ss = return_elastic_stress_strain_tensor(material_id, db_current_temperature);

            updated_int_variables.alpha = alpha_n;
            updated_int_variables.total_strain = e_n1;
            updated_int_variables.plastic_strain = ep_n;
            updated_int_variables.thermal_strain = eth_n1;
            updated_int_variables.cauchy_stress = cauchy_stress;
            updated_int_variables.temperature = db_current_temperature;
            updated_int_variables.material_id = material_id;
            updated_int_variables.hardening_stress = R_trial;

            return;
        }
        else { /* plastic case */
            //deallog << "PLASTIC CASE " << std::endl;
            double sigma_dev_trial_norm = sigma_dev_trial.norm();
            SymmetricTensor<2, dim> normal = sigma_dev_trial / sigma_dev_trial_norm;

            double delta_lambda = return_delta_lambda(material_id, phi_trial, db_current_time_step_length, db_current_temperature);
            double alpha_n1 = alpha_n + delta_lambda * sqrt(2.0 / 3.0);
            double delta_R = - delta_lambda * sqrt(2.0 / 3.0) * ptr_mech_material->return_hardening_modulus(material_id, db_current_temperature);


            SymmetricTensor<2, dim> ep_n1 = ep_n + delta_lambda * normal;
            ep_delta = ep_n1 - ep_n;

            SymmetricTensor<2, dim> sigma_dev;
            if(cauchy_stress_n.norm()>= 1e-10){
                sigma_dev = deviator(cauchy_stress_n)
                            + return_delta_cauchy_stress_deviator(material_id,
                                                                  e_delta,
                                                                  ep_delta,
                                                                  e_n,
                                                                  ep_n,
                                                                  T_delta,
                                                                  db_current_temperature);
            }
            else{
                sigma_dev = return_delta_cauchy_stress_deviator(material_id,
                                                                e_delta,
                                                                ep_delta,
                                                                e_n,
                                                                ep_n,
                                                                0.0,
                                                                db_current_temperature);
            }

            cauchy_stress = sigma_vol + sigma_dev;
            elasto_plastic_tangent_ss = return_elasto_plastic_stress_strain_tensor(material_id,
                                                                                   sigma_dev_trial_norm,
                                                                                   normal,
                                                                                   delta_lambda,
                                                                                   db_current_time_step_length,
                                                                                   db_current_temperature);

            updated_int_variables.alpha = alpha_n1;
            updated_int_variables.total_strain = e_n1;
            updated_int_variables.plastic_strain = ep_n1;
            updated_int_variables.thermal_strain = eth_n1;
            updated_int_variables.cauchy_stress = cauchy_stress;
            updated_int_variables.temperature = db_current_temperature;
            updated_int_variables.material_id = material_id;
            updated_int_variables.hardening_stress = R_trial + delta_R;

            return;
        }
    }


    template<int dim>
    inline void Thermo_Elasto_Plastic<dim>::return_stress_tangent_int_variables_for_large_strain
		(
    		/*input->*/
            unsigned int material_id,
            Tensor<2, dim> deformation_gradient,
            double alpha_n,
            SymmetricTensor<2, dim> e_n,
            SymmetricTensor<2, dim> ep_n,
            SymmetricTensor<2, dim> eth_n,
            SymmetricTensor<2, dim> cauchy_stress_n,
            double db_temperature_previous_time_step,
            unsigned int int_material_id_previous_time_step,
            double db_current_temperature,
            double db_current_time_step_length,
			/*output->*/
            SymmetricTensor<2,dim> &second_piola_stress,
            SymmetricTensor<4,dim> &elasto_plastic_tangent,
            QPH::internal_variables_incremental<dim> &updated_int_variables
		)
	{
        const double comp_tolerance = 1e-8;

    	// Following "Algorithms for computation of stresses and elasticity moduli in terms of Seth–Hill’s family of generalized strain tensors" by Miehe&Lambrecht \n
    	// Table I. Algorithm A
    	/*
    	 * 1. Eigenvalues, eigenvalue bases and diagonal functions:
    	 */

    	// Get the symmetric right cauchy green tensor
         SymmetricTensor<2, dim> right_cauchy_green_sym = Physics::Elasticity::Kinematics::C(deformation_gradient);

        // Compute Eigenvalues, Eigenvectors and Eigenbasis
         Vector<double> eigenvalues(dim);
         std::vector<Tensor<1, dim>> eigenvector(dim);
         std::vector<SymmetricTensor<2, dim>> eigenbasis(dim);
         {
			// Get Eigenvalues and Eigenvectors from the deal.ii function \a eigenvectors(*)
			for (unsigned int i = 0; i < dim; ++i) {
				eigenvalues[i] = eigenvectors(right_cauchy_green_sym)[i].first;
				eigenvector[i] = eigenvectors(right_cauchy_green_sym)[i].second;
			}

			// Check if the found eigenvectors are perpendicular to each other
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

			// Compute eigenbasis Ma: eigenbasis = eigenvector \otimes eigenvector
			for (unsigned int i = 0; i < dim; ++i) {
				eigenbasis[i] = symmetrize( outer_product(eigenvector[i], eigenvector[i]) );
				AssertThrow( eigenvalues(i) >= 0.0,
							 ExcMessage("ln-space<< Eigenvalue is negativ. Check update_qph.") );
			}
         }

        // Compute diagonal function \a ea and its first and second derivate \a da and \a fa \n
         Vector<double> ea(dim);
         Vector<double> da(dim);
         Vector<double> fa(dim);
         for (unsigned int i = 0; i < dim; ++i)
         {
            ea(i) = 0.5 * std::log( std::abs(eigenvalues(i)) );	// diagonal function
            da(i) = std::pow(eigenvalues(i), -1.0);				// first derivative of diagonal function ea
            fa(i) = -2.0 * std::pow(eigenvalues(i), -2.0);			// second derivative of diagonal function ea
            AssertThrow( ea(i) == ea(i),
                         ExcMessage( "ln-space<< Ea is nan due to logarithm of negativ eigenvalue. Check update_qph.") ;
            AssertThrow( da(i) > 0.0,
                         ExcMessage( "ln-space<< First derivative da of diagonal function is "+std::to_string(da(i))+" < 0.0 . Check update_qph.") );
         }

        // Compute the Hencky strain
         SymmetricTensor<2, dim> hencky_strain;
         for (unsigned int a = 0; a < dim; ++a)
            hencky_strain += ea(a) * eigenbasis[a];

        // Call the small strain material model with the hencky_strain and history
         // Declare the return arguments
          SymmetricTensor<2, dim> stress_measure_T_sym;
          SymmetricTensor<4,dim> elasto_plastic_tangent_ss;

         return_stress_tangent_int_variables_for_small_strain
		 	 (
				/*input->*/
				material_id,
				hencky_strain,
				alpha_n,
				e_n,
				ep_n,
				eth_n,
				cauchy_stress_n,
				db_temperature_previous_time_step,
				int_material_id_previous_time_step,
				db_current_temperature,
				db_current_time_step_length,
				/*output->*/
				stress_measure_T_sym,
				elasto_plastic_tangent_ss,
				updated_int_variables
			);


        /*
         * 3. Set up coefficients \a theta, \a xi and \a eta (step 2 was bypassed)
         */

        // Compute the coefficients based on the eigenvalues, eigenvectors and ea,da,fa
         Tensor<2, dim> theta;
         Tensor<2, dim> xi;
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

        // For two equal eigenvalues a and b: \f$ \lambda_a = \lambda_b \neq \lambda_c \f$
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
        // For two equal eigenvalues a and c: \f$ \lambda_a = \lambda_c \neq \lambda_b \f$
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
        // For two equal eigenvalues b and c: \f$ \lambda_b = \lambda_c \neq \lambda_a \f$
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

        // Ensure that \a eta was initialised in one of the cases
         AssertThrow( (eta < 9999999), ExcMessage("Eta in update_qph not initialised") );


        /*
         * 4. Lagrangian stresses and elasticity moduli
         */

        // Compute projection tensor P
         Tensor<4, dim> projection_tensor_P;
         for (unsigned int a = 0; a < dim; ++a)
         {
            projection_tensor_P += da(a) * (Tensor<4, dim> ) outer_product(eigenbasis[a],eigenbasis[a]);
            for (unsigned int b = 0; b < dim; ++b)
                if (b != a)
                    projection_tensor_P += theta[a][b] * get_tensor_operator_G(eigenbasis[a],eigenbasis[b]);
         }

        // Check whether the projecton tensor is symmetric and store it into a \a SymmetricTensor
         AssertThrow( symmetry_check(projection_tensor_P), ExcInternalError( "ln-space<< Projection tensor P is not symmetric") );
         SymmetricTensor<4, dim> projection_tensor_P_sym = symmetrize(projection_tensor_P);

        // Compute the double contraction of T and L
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
                                                    	    get_tensor_operator_F_right(
                                                    	    	eigenbasis[a],
																eigenbasis[b], eigenbasis[b],
																stress_measure_T_sym )
                                                          + get_tensor_operator_F_left(
																eigenbasis[b], eigenbasis[b],
																eigenbasis[a],
																stress_measure_T_sym )
														  + get_tensor_operator_F_right(
																eigenbasis[b], eigenbasis[a],
																eigenbasis[b],
																stress_measure_T_sym )
														  + get_tensor_operator_F_left(
																eigenbasis[a], eigenbasis[b],
																eigenbasis[b],
																stress_measure_T_sym )
														  + get_tensor_operator_F_right(
																eigenbasis[b], eigenbasis[b],
																eigenbasis[a],
																stress_measure_T_sym )
														  + get_tensor_operator_F_left(
																eigenbasis[b], eigenbasis[a],
																eigenbasis[b],
																stress_measure_T_sym )
														   );

                    for (unsigned int c = 0; c < dim; ++c)
                        if ( (c != a) && (c != b) )
                        {
                            projection_tensor_T_doublecon_L += 2.0 * eta
                                                               * (
                                                            		  get_tensor_operator_F_right(
																		eigenbasis[a], eigenbasis[b],
																		eigenbasis[c],
																		stress_measure_T_sym )
																	+ get_tensor_operator_F_left(
																		eigenbasis[b], eigenbasis[c],
																		eigenbasis[a],
																		stress_measure_T_sym )
																  );
                        }
                }
         }

        // Check whether the tensor is symmetric and store it into a \a SymmetricTensor
         AssertThrow( symmetry_check(projection_tensor_T_doublecon_L),
                      ExcInternalError("ln-space<< Projection tensor T:L is not symmetric") );
         SymmetricTensor<4, dim> projection_tensor_T_doublecon_L_sym = symmetrize(projection_tensor_T_doublecon_L);

        // Compute the retransformed values
         SymmetricTensor<2, dim> second_piola_stress_S;
         second_piola_stress_S = stress_measure_T_sym
                                * projection_tensor_P_sym;

         elasto_plastic_tangent = projection_tensor_P_sym * elasto_plastic_tangent_ss * projection_tensor_P_sym
                                  + projection_tensor_T_doublecon_L_sym;
         second_piola_stress = second_piola_stress_S;

         updated_int_variables.cauchy_stress_ls = (second_piola_stress * (1 / determinant(deformation_gradient)) );
    }

//////////////////////////////
//// private functions   ////
/////////////////////////////

    template<int dim>
    SymmetricTensor<2, dim>
    Thermo_Elasto_Plastic<dim>::return_delta_e_th(unsigned int material_id,
                                                  unsigned int int_material_id_previous_time_step,
                                                  double db_temperature_previous_time_step,
                                                  double db_current_temperature) {

        SymmetricTensor<2,dim> delta_eth_integral;
        if(material_id == material_id_powder){
            return delta_eth_integral;
        }
        double db_T_upper_bound = db_current_temperature;
        double db_T_lower_bound = db_temperature_previous_time_step;
        if(int_material_id_previous_time_step == material_id_powder){
            db_T_lower_bound = ptr_mech_material->return_reference_temperature();
        }

        double integral_of_CTE = ptr_mech_material->return_expansion_coefficient_integral_form(material_id,db_T_upper_bound,db_T_lower_bound);
        delta_eth_integral  = integral_of_CTE * I;
        return delta_eth_integral;
    }

    template<int dim>
    double
    Thermo_Elasto_Plastic<dim>::return_yield_function(unsigned int material_id,
                                                      SymmetricTensor<2, dim> sigma_dev_trial,
                                                      double hardening_stress_R_trial,
                                                      double db_current_temperature) {
        double tmp;
        double db_yield_stress_at_current_temp = ptr_mech_material->return_yield_stress(material_id, db_current_temperature);
        tmp = sigma_dev_trial.norm() -  sqrt(2.0/3.0) * ( db_yield_stress_at_current_temp - hardening_stress_R_trial);

        return tmp;
    }

    template<int dim>
    double Thermo_Elasto_Plastic<dim>::return_delta_lambda(unsigned int material_id,
                                                           double yield_function,
                                                           double db_delta_time_step,
                                                           double db_current_temperature) {
        double tmp;
        double db_shear_modulus_at_current_temp = ptr_mech_material->return_shear_modulus(material_id, db_current_temperature);
        double db_hardening_modulus_at_current_temp = ptr_mech_material->return_hardening_modulus(material_id, db_current_temperature);
        double db_viscosity_at_current_temp = ptr_mech_material->return_viscosity(material_id, db_current_temperature);

        tmp = yield_function
              /
              (   2.0*db_shear_modulus_at_current_temp
                  + (2.0/3.0) * db_hardening_modulus_at_current_temp
                  + (db_viscosity_at_current_temp/db_delta_time_step)
              );

        return tmp;
    }

    template<int dim>
    double Thermo_Elasto_Plastic<dim>::return_hardening_stress_R_trial(unsigned int material_id_n,
                                                                       unsigned int material_id_n1,
                                                                       double alpha_n,
                                                                       double alpha_n1,
                                                                       double db_T_n,
                                                                       double db_current_temperature) {
        double tmp;

        double db_hardening_modulus_at_current_temp = ptr_mech_material->return_hardening_modulus(material_id_n1, db_current_temperature);
        double db_hardening_modulus_at_T_n = ptr_mech_material->return_hardening_modulus(material_id_n, db_T_n);
        double db_dev_hardening_modulus_wrt_temp = ptr_mech_material->return_dev_hardening_modulus_wrt_temp(material_id_n1, db_current_temperature);
        double delta_alpha = alpha_n1 - alpha_n;
        double delta_temp = db_current_temperature - db_T_n;

        tmp = - alpha_n * db_hardening_modulus_at_T_n
              - (  db_hardening_modulus_at_current_temp * delta_alpha
                   + db_dev_hardening_modulus_wrt_temp * alpha_n  * delta_temp );

        return tmp;
    }

    template<int dim>
    SymmetricTensor<4, dim>
    Thermo_Elasto_Plastic<dim>::return_elastic_stress_strain_tensor(unsigned int material_id,
                                                                    double db_current_temperature) {
        SymmetricTensor<4,dim> tmp;

        double db_bulk_modulus_at_current_temp = ptr_mech_material->return_bulk_modulus(material_id, db_current_temperature);
        double db_shear_modulus_at_current_temp = ptr_mech_material->return_shear_modulus(material_id, db_current_temperature);

        tmp = 	  db_bulk_modulus_at_current_temp * IxI
                   + 2.0 * db_shear_modulus_at_current_temp *deviatoric_identity;

        return tmp;
    }

    template<int dim>
    SymmetricTensor<2, dim> Thermo_Elasto_Plastic<dim>::return_delta_cauchy_stress_deviator(unsigned int material_id,
                                                                                            SymmetricTensor<2, dim> e_delta,
                                                                                            SymmetricTensor<2, dim> ep_delta,
                                                                                            SymmetricTensor<2, dim> e_n,
                                                                                            SymmetricTensor<2, dim> ep_n,
                                                                                            double db_T_delta,
                                                                                            double db_current_temperature) {
        SymmetricTensor<2,dim> tmp;
        double db_shear_modulus_at_current_temp = ptr_mech_material->return_shear_modulus(material_id, db_current_temperature);
        double db_dev_shear_modulus_wrt_temp = ptr_mech_material->return_dev_shear_modulus_wrt_temp(material_id, db_current_temperature);

        tmp = 2.0 * db_shear_modulus_at_current_temp * (deviator(e_delta) - ep_delta)
              + 2.0 * (deviator (e_n) - ep_n) * db_dev_shear_modulus_wrt_temp * db_T_delta;

        return tmp;
    }

    template<int dim>
    SymmetricTensor<2, dim> Thermo_Elasto_Plastic<dim>::return_delta_cauchy_stress_vol(unsigned int material_id,
                                                                                       SymmetricTensor<2, dim> emech_delta,
                                                                                       SymmetricTensor<2, dim> emech_n,
                                                                                       double db_T_delta,
                                                                                       double db_current_temperature) {
        SymmetricTensor<2,dim> tmp;
        double db_bulk_modulus_at_current_temp = ptr_mech_material->return_bulk_modulus(material_id, db_current_temperature);
        double db_dev_bulk_modulus_wrt_temp = ptr_mech_material->return_dev_bulk_modolus_wrt_temp(material_id, db_current_temperature);

        tmp= 3.0 * db_bulk_modulus_at_current_temp * (emech_delta - deviator(emech_delta))
             + 3.0 * (emech_n - deviator(emech_n) ) * db_dev_bulk_modulus_wrt_temp * db_T_delta;

        return tmp;
    }

    template<int dim>
    SymmetricTensor<4, dim>
    Thermo_Elasto_Plastic<dim>::return_elasto_plastic_stress_strain_tensor(unsigned int material_id,
                                                                           double sigma_dev_trial_norm,
                                                                           SymmetricTensor<2,dim> normal,
                                                                           double delta_lambda,
                                                                           double db_delta_timestep,
                                                                           double db_current_temperature)
    {
        SymmetricTensor<4,dim> tmp;
        double db_bulk_modulus_at_current_temp = ptr_mech_material->return_bulk_modulus(material_id, db_current_temperature);
        double db_shear_modulus_at_current_temp = ptr_mech_material->return_shear_modulus(material_id, db_current_temperature);
        double db_hardening_modulus_at_current_temp = ptr_mech_material->return_hardening_modulus(material_id, db_current_temperature);
        double db_viscosity_at_current_temp = ptr_mech_material->return_viscosity(material_id, db_current_temperature);

        tmp = db_bulk_modulus_at_current_temp * IxI
              + 2.0  * db_shear_modulus_at_current_temp * deviatoric_identity
              - 2.0 * db_shear_modulus_at_current_temp *
                (
                        (   (2.0 * db_shear_modulus_at_current_temp)
                            / (2.0 * db_shear_modulus_at_current_temp
                               + (2.0/3.0) * db_hardening_modulus_at_current_temp
                               + (db_viscosity_at_current_temp / db_delta_timestep)
                            )
                        ) * outer_product(normal,normal)
                        +
                        (   (2.0 * db_shear_modulus_at_current_temp * delta_lambda )
                            / sigma_dev_trial_norm
                        ) * ( deviatoric_identity - outer_product(normal,normal) )
                );

        return tmp;
    }


} //END Namespace

#endif // CONSTITUTIVE_LAW_THERMO_ELASTO_PLASTIC
