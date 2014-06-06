// ==========================================================================
// Copyright (C) 2007 Ben Sunshine-Hill
// ==========================================================================

#ifndef CONJGRAD_H
#define CONJGRAD_H

#pragma warning(disable: 4244 4267 4996)

#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include <boost/numeric/ublas/io.hpp>

#include <utility>


/** Use the conjugate gradient method to solve a symmetric, positive semidefinite system of linear equations.
  *
  * This is an iterative method which attempts to converge on a solution.
  * The method will exit when either the number of iterations reaches max_iter,
  * or the 2-norm of the residual goes below tol.
  *
  * @param A         The system of linear equations to solve. This must be symmetric and positive definite.
  * @param b         The constant terms of the system.
  * @param x         Where the output is stored. This must be a vector of the appropriate size. It is overwritten by this method.
  * @param max_iter  The maximum number of iterations the solver will make.
  * @param tol       The maximum tolerable error, below which the solver will consider the solution to have converged.
  */
template<typename Matrix, typename Vector>
bool cg_solve(const Matrix &A, const Vector &b, Vector &x, int max_iter, double tol)
{
	std::fill(x.begin(), x.end(), 0);

	Vector r = b;
	Vector p = b;
	Vector r_old;
	Vector temp;

	// CG loop
	for(int niter = 0; niter < max_iter; niter++)
	{
		temp = prod(A,p);
		double alpha = inner_prod(r,r)/inner_prod(p,temp);
		x += (p*alpha);
		r_old = r;
		r -= (temp*alpha);
		double residn = norm_2(r);
		if(residn < tol){
         //std::cout << "numiters: "<< niter <<std::endl;
			return true;
		}
		double beta = inner_prod(r,r)/inner_prod(r_old,r_old);      
		p = r + p*beta;
	}

   std::cout << "WARNING: cg_solve did not converge" << std::endl;
	return false;
}

template<typename Matrix, typename Vector>
bool cg_psolve(const Matrix &A, const Matrix &Cinv, const Vector &b, Vector &x, int max_iter, double tol)
{
	std::fill(x.begin(), x.end(), 0);
   Vector r = b;
   Vector z = prod(Cinv, r); 
   Vector d = z;
   Vector temp;
   double resign;

	for(int niter = 0; niter < max_iter; niter++)
	{
      temp = prod(A, d);
      double dotzr = inner_prod(z, r);
      double alpha = dotzr/inner_prod(d, temp);
      x += alpha*d;

      r -= alpha*temp;
      resign = norm_2(r);
      if(resign < tol) 
      { 
        /*std::cout << "numiters: "<< niter <<std::endl;*/ 
        return true; 
      }
      
      z = prod(Cinv, r);
      double beta = inner_prod(z, r)/dotzr;
      d = z + beta*d;
   }

   std::cout << "WARNING: cg_psolve did not converge: " << resign << std::endl;
	return false;
}

#endif