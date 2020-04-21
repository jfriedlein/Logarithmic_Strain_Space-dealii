/* ---------------------------------------------------------------------
 *
 * The first coding assignment to get familiar with tensor calculus related
 * deal.II class templates. This includes:
 * 
 * - Tensor<1,dim>
 * - Tensor<2,dim>
 * - Tensor<4,dim>
 * - Vector<double>
 * - FullMatrix<double>
 * 
 * dim is a template variable which allows to vary e.g. between the two- and
 * three dimensional case. As described in the brief repetition of the 
 * essentials in C++, the respective tensor class templates allow different
 * dimensions as input, i.e. dim=1,dim=2,dim=3 for each respective rank
 *
 * ---------------------------------------------------------------------
 */


#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/timer.h>

#include <iostream>
#include <vector>

using namespace dealii;


template<int dim>
Tensor<4,dim> get_tensor_operator_G(const SymmetricTensor<2,dim> &Ma, const SymmetricTensor<2,dim> &Mb)
{
	Tensor<4,dim> tmp; // has minor symmetry of indices k,l

	for(unsigned int i=0; i<dim; ++i)
		for(unsigned int j=0; j<dim; ++j)
			for(unsigned int k=0; k<dim; ++k)
				for(unsigned int l=0; l<dim; ++l)
					tmp[i][j][k][l] = Ma[i][k] * Mb[j][l] + Ma[i][l] * Mb[j][k];

	return tmp;
}


template<int dim>
Tensor<4,dim> get_tensor_operator_G_1(const SymmetricTensor<2,dim> &Ma, const SymmetricTensor<2,dim> &Mb)
{
	Tensor<4,dim> tmp; // has minor symmetry of indices k,l

	for(unsigned int i=0; i<dim; ++i)
		for(unsigned int j=0; j<dim; ++j)
			for(unsigned int k=0; k<dim; ++k)
				for(unsigned int l=k; l<dim; ++l)
				{
					double tmp_scalar = Ma[i][k] * Mb[j][l] + Ma[i][l] * Mb[j][k];
					tmp[i][j][k][l] = tmp_scalar;
					tmp[i][j][l][k] = tmp_scalar;
				}

	return tmp;
}


// F right F_{a(bc)}
template<int dim>
Tensor<4,dim> get_tensor_operator_F_right(const SymmetricTensor<2,dim> &Ma,
										  const SymmetricTensor<2,dim> &Mb,
										  const SymmetricTensor<2,dim> &Mc,
										  const SymmetricTensor<2,dim> &T )
{
	Tensor<4,dim> tmp;

	Tensor<2,dim> temp_tensor = contract<1,0>((Tensor<2,dim>)T, (Tensor<2,dim>)Mc);
	Tensor<2,dim> MbTMc = contract<1,0>((Tensor<2,dim>)Mb,temp_tensor);

	for(unsigned int i=0; i<dim; ++i)
		for(unsigned int j=0; j<dim; ++j)
			for(unsigned int k=0; k<dim; ++k)
				for(unsigned int l=0; l<dim; ++l)
					tmp[i][j][k][l] = Ma[i][k] * MbTMc[j][l] + Ma[i][l] * MbTMc[j][k];

	return tmp;
}


// F right F_{a(bc)}
template<int dim>
Tensor<4,dim> get_tensor_operator_F_right_1(const SymmetricTensor<2,dim> &Ma,
										  const SymmetricTensor<2,dim> &Mb,
										  const SymmetricTensor<2,dim> &Mc,
										  const SymmetricTensor<2,dim> &T )
{
	Tensor<4,dim> tmp; // has minor symmetry of indices k,l

	Tensor<2,dim> temp_tensor = contract<1,0>((Tensor<2,dim>)T, (Tensor<2,dim>)Mc);
	Tensor<2,dim> MbTMc = contract<1,0>((Tensor<2,dim>)Mb,temp_tensor);

	for(unsigned int i=0; i<dim; ++i)
		for(unsigned int j=0; j<dim; ++j)
			for(unsigned int k=0; k<dim; ++k)
				for(unsigned int l=k; l<dim; ++l)
				{
					double tmp_scalar = Ma[i][k] * MbTMc[j][l] + Ma[i][l] * MbTMc[j][k];
					tmp[i][j][k][l] = tmp_scalar;
					tmp[i][j][l][k] = tmp_scalar;
				}

	return tmp;
}


// F right F_{(ab)c}
template<int dim>
Tensor<4,dim> get_tensor_operator_F_left(const SymmetricTensor<2,dim> &Ma,
										 const SymmetricTensor<2,dim> &Mb,
										 const SymmetricTensor<2,dim> &Mc,
										 const SymmetricTensor<2,dim> &T){
	Tensor<4,dim> tmp;

	Tensor<2,dim> temp_tensor = contract<1,0>((Tensor<2,dim>)T, (Tensor<2,dim>)Mb);
	Tensor<2,dim> MaTMb = contract<1,0>((Tensor<2,dim>)Ma,temp_tensor);

	for(unsigned int i=0; i<dim; ++i)
		for(unsigned int j=0; j<dim; ++j)
			for(unsigned int k=0; k<dim; ++k)
				for(unsigned int l=0; l<dim; ++l)
					tmp[i][j][k][l] = MaTMb[i][k] * Mc[j][l] + MaTMb[i][l] * Mc[j][k];

	return tmp;
}

// F right F_{(ab)c}
template<int dim>
Tensor<4,dim> get_tensor_operator_F_left_1(const SymmetricTensor<2,dim> &Ma,
										 const SymmetricTensor<2,dim> &Mb,
										 const SymmetricTensor<2,dim> &Mc,
										 const SymmetricTensor<2,dim> &T){
	Tensor<4,dim> tmp;

	Tensor<2,dim> temp_tensor = contract<1,0>((Tensor<2,dim>)T, (Tensor<2,dim>)Mb);
	Tensor<2,dim> MaTMb = contract<1,0>((Tensor<2,dim>)Ma,temp_tensor);

	for(unsigned int i=0; i<dim; ++i)
		for(unsigned int j=0; j<dim; ++j)
			for(unsigned int k=0; k<dim; ++k)
				for(unsigned int l=k; l<dim; ++l)
				{
					double tmp_scalar = MaTMb[i][k] * Mc[j][l] + MaTMb[i][l] * Mc[j][k];
					tmp[i][j][k][l] = tmp_scalar;
					tmp[i][j][l][k] = tmp_scalar;
				}

	return tmp;
}

//----------------------------------------------------

int main ()
{
    /* For now the three dimensional case is considered.
     * This information is essential for deal.II since
     * is suited for arbitrary dimensions due to its
     * template character. Use this variable for all
     * deal.II templates that will be created in the
     * sequent
     * 
     * It is always helpful to read the manual and see
     * how functions are implemented, i.e. what is the
     * return value, how is the function called or what
     * is the input.
     * 
     * Further some functions may be declared 
     * "DEPRECATED" which means they are still usable
     * but will be removed in future releases of the 
     * library -> not recommended to use those
     */
    
    
    const int dim=3;
    //------------------------------------------------
    /* START IMPLEMENTATION HERE                    */
    //------------------------------------------------
    
    
    //------------------------------------------------
    //          EX - 1
    /* Create two tensors of rank one and name them
     * u and v respectively and print them to the screen.
     * Therefore consider the available documentation
     * and manual on the deal.II website
     */
    TimerOutput timer (std::cout, TimerOutput::summary,
                       TimerOutput::cpu_times);
    
    Tensor<1,dim> u;
    Tensor<1,dim> v;
    Tensor<1,dim> w;
    
    //BEGIN - INSERT YOUR CODE HERE
    u[0]=1.5; u[1]=0.61; u[2]=3.63;        
    v[0]=4.5; v[1]=5; v[2]=6.8;
    w[0]=7; w[1]=8.98; w[2]=9.35;
    //END - INSERT YOUR CODE HERE
    
    std::vector< SymmetricTensor<2,dim> > Ma (3);
    Ma[0] = symmetrize(outer_product(u,u));
    Ma[1] = symmetrize(outer_product(v,v));
    Ma[2] = symmetrize(outer_product(w,w));

    SymmetricTensor<2,dim> T;
    T[0][0] = 5.654;
    T[1][1] = 1.97;
    T[2][2] = 0.61;
    T[0][1] = 3.651;
    T[0][2] = 7.125;
    T[1][2] = 4.99;

    std::cout << Ma[0] << std::endl;
    std::cout << Ma[1] << std::endl;
    std::cout << Ma[2] << std::endl;

    Tensor<4,dim> F;
    timer.enter_subsection("F");
    for ( unsigned int i=0; i<10000; i++)
    	F = get_tensor_operator_F_left(Ma[0],Ma[1],Ma[2],T);
    timer.leave_subsection("F");
    
    Tensor<4,dim> F1;
    timer.enter_subsection("F1");
    for ( unsigned int i=0; i<10000; i++)
    	F1 = get_tensor_operator_F_left_1(Ma[0],Ma[1],Ma[2],T);
    timer.leave_subsection("F1");
    
    std::cout << "F=" << std::endl;
	for ( unsigned int i=0; i<dim; ++i )
    	for ( unsigned int j=0; j<dim; ++j )
        	for ( unsigned int k=0; k<dim; ++k )
            	for ( unsigned int l=0; l<dim; ++l )
            		std::cout << i<<j<<k<<l<< ": " << F[i][j][k][l] << std::endl;
	
    std::cout << "G_for1=" << std::endl;
	for ( unsigned int i=0; i<dim; ++i )
    	for ( unsigned int j=0; j<dim; ++j )
        	for ( unsigned int k=0; k<dim; ++k )
            	for ( unsigned int l=0; l<dim; ++l )
            		std::cout << i<<j<<k<<l<< ": " << F1[i][j][k][l] << std::endl;
	
    double error_dfor = 0.;
    double error_dfor2 = 0.;
//    double error_dfor3 = 0.;
	for ( unsigned int i=0; i<dim; ++i )
    	for ( unsigned int j=0; j<dim; ++j )
        	for ( unsigned int k=0; k<dim; ++k )
            	for ( unsigned int l=0; l<dim; ++l )
            	{
            		error_dfor += std::abs(F[i][j][k][l]-F1[i][j][k][l]);
//            		error_dfor2 += std::abs(D[i][j][k][l]-D_for2[i][j][k][l]);
//            		error_for3 += std::abs(Ma_x_Ma[i][j][k][l]-Ma_x_Ma_for3[i][j][k][l]);
            	}
	std::cout << "error 1=" << error_dfor << std::endl;
}

//----------------------------------------------------