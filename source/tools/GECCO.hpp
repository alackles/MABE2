/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  GECCO.hpp
 *  @brief This file provides code to build GECCO Niching Landscapes.
 *  @note This file is a modified composite of cfunction.h and cfunction.cpp in CEC2013: https://github.com/mikeagn/CEC2013/ 
 * 
 *  
 */

#ifndef MABE_TOOL_GECCO_H
#define MABE_TOOL_GECCO_H

#include <map>
#include <string>
#include <cassert>

#include "emp/base/assert.hpp"
#include "emp/base/Ptr.hpp"
#include "emp/base/vector.hpp"

#include "emp/io/File.hpp"
#include "emp/math/math.hpp"
#include "emp/math/Random.hpp"
#include "emp/base/vector.hpp"
#include "emp/bits/BitVector.hpp"
#include "emp/math/constants.hpp"

namespace mabe {

  /// The GECCO Niching Competition provides popular simple functions on which to test evolution of bitstrings.

  /* Composition Functions framework */
  class CFunction {
  	//non-copyable
  	CFunction(const CFunction&);
  	CFunction& operator=(const CFunction &);
  public:
  	CFunction();
  	CFunction(const int &dim, const int &nofunc);
  	virtual ~CFunction();

  	virtual double evaluate(const emp::vector<double> x) = 0;
  	double get_lbound(const int &ivar) const { return lbound_[ivar]; } 
  	double get_ubound(const int &ivar) const { return ubound_[ivar]; } 

  protected:
  	int dimension_;
  	int nofunc_;
  	double C_;
  	double f_bias_;
  	emp::vector<double> lambda_;
  	emp::vector<double> sigma_;
  	emp::vector<double> bias_;
  	emp::vector<double> weight_;
  	emp::vector<double> lbound_;
  	emp::vector<double> ubound_;
  	emp::vector<double> fi_;
  	emp::vector<double> z_;
  	emp::vector<double> fmaxi_;
  	emp::vector<double> tmpx_;
  	emp::vector<emp::vector<double>> O_;
    emp::vector<emp::vector<emp::vector<double>>> M_;
    emp::Random rng_;


  	/* Inner help functions */
  	void init_rotmat_identity();
  	void init_optima_rand();
  	void load_optima(const std::string &filename);
  	void load_rotmat(const std::string &filename);
  	void calculate_weights(const emp::vector<double> x);
  	void transform_to_z(const emp::vector<double> x, const int &index);
  	void transform_to_z_noshift(const emp::vector<double> x, const int &index);
  	void calculate_fmaxi();
  	double evaluate_inner_(const emp::vector<double> x);
  	std::vector< std::vector<double> > get_copy_of_goptima() const;
  };

  /* Basic Benchmark functions */

  /******************************************************************************
    * F1: Five-Uneven-Peak Trap 
    * Variable ranges: x in [0, 30
    * No. of global peaks: 2
    * No. of local peaks:  3. 
    *****************************************************************************/
  double five_uneven_peak_trap(const emp::vector<double> x, const int &dim) {
    double result=-1.0;
    if (x[0]>=0 && x[0]< 2.5) {
    result = 80*(2.5-x[0]);
    } else if (x[0] >= 2.5 && x[0] < 5.0) {
    result = 64*(x[0]-2.5);
    } else if (x[0] >= 5.0 && x[0] < 7.5) {
    result = 64*(7.5-x[0]);
    } else if (x[0] >= 7.5 && x[0] < 12.5) {
    result = 28*(x[0]-7.5);
    } else if (x[0] >= 12.5 && x[0] < 17.5) {
    result = 28*(17.5-x[0]);
    } else if (x[0] >= 17.5 && x[0] < 22.5) {
    result = 32*(x[0]-17.5);
    } else if (x[0] >= 22.5 && x[0] < 27.5) {
    result = 32*(27.5-x[0]);
    } else if (x[0] >= 27.5 && x[0] <= 30) {
    result = 80*(x[0]-27.5);
    }
    return result;
  }

  /******************************************************************************
  * F2: Equal Maxima
  * Variable ranges: x in [0, 1]
  * No. of global peaks: 5
  * No. of local peaks:  0. 
  *****************************************************************************/
  double equal_maxima(const emp::vector<double> x, const int &dim) {
    double s = sin(5.0 * emp::PI * x[0]);
    return pow(s, 6);
  }

  /******************************************************************************
  * F3: Uneven Decreasing Maxima
  * Variable ranges: x in [0, 1]
  * No. of global peaks: 1
  * No. of local peaks:  4. 
  *****************************************************************************/
  double uneven_decreasing_maxima(const emp::vector<double> x, const int &dim) {
    double tmp1 = -2*log(2)*((x[0]-0.08)/0.854)*((x[0]-0.08)/0.854);
    double tmp2 = sin( 5*emp::PI*(pow(x[0],3.0/4.0)-0.05) );
    return exp(tmp1) * pow(tmp2, 6);
  }

  /******************************************************************************
  * F4: Himmelblau
  * Variable ranges: x, y in [−6, 6
  * No. of global peaks: 4
  * No. of local peaks:  0.
  *****************************************************************************/
  double himmelblau(const emp::vector<double> x, const int &dim) {
    return 200 - (x[0]*x[0] + x[1] - 11)*(x[0]*x[0] + x[1] - 11) - 
    (x[0] + x[1]*x[1] - 7)*(x[0] + x[1]*x[1] - 7);
  }

  /******************************************************************************
  * F5: Six-Hump Camel Back
  * Variable ranges: x in [−1.9, 1.9]; y in [−1.1, 1.1]
  * No. of global peaks: 2
  * No. of local peaks:  2.
  *****************************************************************************/
  double six_hump_camel_back(const emp::vector<double> x, const int &dim) {
    return -( (4 - 2.1*x[0]*x[0] + pow(x[0],4.0)/3.0)*x[0]*x[0] + 
    x[0]*x[1] + (4*x[1]*x[1] -4)*x[1]*x[1] );
  }

  /******************************************************************************
  * F6: Shubert
  * Variable ranges: x_i in  [−10, 10]^n, i=1,2,...,n
  * No. of global peaks: n*3^n
  * No. of local peaks: many
  *****************************************************************************/
  double shubert(const emp::vector<double> x, const int &dim) {
    double result(1), sum(0); 
    for (int i=0; i<dim; i++) {
    sum=0;
    for (int j=1; j<6; j++) {
      sum = sum + j * cos((j+1) * x[i] + j);
    }
    result = result * sum;
    }
    return -result;
  }

  /******************************************************************************
  * F7: Vincent
  * Variable range: x_i in [0.25, 10]^n, i=1,2,...,n
  * No. of global optima: 6^n
  * No. of local optima:  0.
  *****************************************************************************/
  double vincent(const emp::vector<double> x, const int &dim) {
    double result(0);
    for (int i=0; i<dim; i++){
    if (x[i]<=0){
      std::cerr << "Illegal value: " << x[i] << std::endl;
      exit(-1);
    }
    result = result + sin(10 * log(x[i]));
    }
    return result/dim;
  }

  /******************************************************************************
  * F8: Modified Rastrigin - All Global Optima
  * Variable ranges: x_i in [0, 1]^n, i=1,2,...,n
  * No. of global peaks: \prod_{i=1}^n k_i
  * No. of local peaks:  0.
  *****************************************************************************/
  /* Modified Rastrigin -- All Global Optima */
  constexpr static double MPPF92[2] = {3, 4};
  constexpr static double MPPF98[8] = {1, 2, 1, 2, 1, 3, 1, 4};
  constexpr static double MPPF916[16] = {1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 4};

  double modified_rastrigin_all(const emp::vector<double> x, const int &dim)
  {
    double result(0);
    for (int i=0; i<dim; i++){
    if (dim == 2)  { result = result + 10+ 9*cos(2*emp::PI*MPPF92[i]*x[i]); }
    if (dim == 8)  { result = result + 10+ 9*cos(2*emp::PI*MPPF98[i]*x[i]); }
    if (dim == 16) { result = result + 10+ 9*cos(2*emp::PI*MPPF916[i]*x[i]); }
    }
    return -result;
  }

  /******************************************************************************
  * Basic functions for composition 
  *****************************************************************************/
  /* Ackley's function */
  double FAckley(const emp::vector<double> x, const int &dim) {
    double sum1(0.0), sum2(0.0), result;
    for (int i=0; i<dim; ++i) {
    sum1 += x[i]*x[i];
    sum2 += cos(2.0*emp::PI*x[i]);
    }
    sum1 = -0.2*sqrt(sum1/dim);
    sum2 /= dim;
    result = 20.0 + emp::E - 20.0*exp(sum1) - exp(sum2);
    return result;
  }

  /* Rastrigin's function */
  double FRastrigin(const emp::vector<double> x, const int &dim) {
    double result(0.0);
    for (int i=0; i<dim; ++i) {
    result += (x[i]*x[i] - 10.0*cos(2.0*emp::PI*x[i]) + 10.0);
    }
    return result;
  }

  /* Weierstrass's function */
  double FWeierstrass(const emp::vector<double> x, const int &dim) {
    double result(0.0), sum(0.0), sum2(0.0), a(0.5), b(3.0);
    int k_max(20);

    for (int j=0; j<=k_max; ++j) {
    sum2 += pow(a,j)*cos(2.0*emp::PI*pow(b,j)*(0.5));
    }
    for (int i=0; i<dim; ++i) {
    sum = 0.0;
    for (int j=0; j<=k_max; ++j) {
      sum += pow(a,j)*cos(2.0*emp::PI*pow(b,j)*(x[i]+0.5));
    }
    result += sum;
    }
    return result - sum2*dim;
  }

  /* Griewank's function */
  double FGriewank(const emp::vector<double> x, const int &dim) {
    double sum(0.0), prod(1.0), result(0.0);

    for (int i=0; i<dim; ++i) {
    sum  += x[i]*x[i]/4000.0;
    prod *= cos( x[i]/sqrt(double(1.0+i)) );
    }
    result = 1.0 + sum - prod;
    return result;
  }

  /* Sphere function */
  double FSphere(const emp::vector<double> x, const int &dim) {
    double result(0.0);
    for (int i=0; i<dim; ++i) {
    result += x[i]*x[i];
    }
    return result;
  }

  /* Schwefel's function */
  double FSchwefel(const emp::vector<double> x, const int &dim) {
    double sum1(0.0), sum2(0.0);

    for (int i=0; i<dim; ++i) {
    sum2 = 0.0;
    for (int j=0; j<=i; ++j) {
      sum2 += x[j];
    }
    sum1 += sum2*sum2;
    }
    return sum1;
  }

  /* Rosenbrock's function */
  double FRosenbrock(const emp::vector<double> x, const int &dim) {
    double result(0.0);

    for (int i=0; i<dim-1; ++i) {
    result += 100.0*pow((x[i]*x[i]-x[i+1]),2.0) + 1.0*pow((x[i]-1.0),2.0);
    }
    return result;
  }

  /* FEF8F2 function */
  double FEF8F2(const emp::vector<double> xx, const int &dim) {
    double result(0.0);
    double x(0), y(0), f(0), f2(0);

    for (int i=0; i<dim-1; ++i) {
    x = xx[i]   +1;
    y = xx[i+1] +1;

    f2 = 100.0*(x*x - y)*(x*x - y) + (1.0 - x)*(1.0 - x);
    f  = 1.0 + f2*f2/4000.0 - cos(f2);

    result += f;
    }
    /* do not forget the (dim-1,0) case! */
    x = xx[dim-1] +1;
    y = xx[0]     +1;

    f2 = 100.0*(x*x - y)*(x*x - y) + (1.0 - x)*(1.0 - x);
    f  = 1.0 + f2*f2/4000.0 - cos(f2);

    result += f;

    return result;
  }

  /* Interfaces for Composition functions */
  class CF1 : public CFunction {
    //non-copyable
    CF1(const CF1 &);
    CF1& operator=(const CF1&);
  public:
    CF1(const int dim);
    double evaluate(const emp::vector<double> x);
  };

  class CF2 : public CFunction {
    //non-copyable
    CF2(const CF2 &);
    CF2& operator=(const CF2&);
  public:
    CF2(const int dim);
  double evaluate(const emp::vector<double> x);
  };

  class CF3 : public CFunction {
    //non-copyable
    CF3(const CF3 &);
    CF3& operator=(const CF3&);
  public:
    CF3(const int dim);
    double evaluate(const emp::vector<double> x);
  };

  class CF4 : public CFunction {
    //non-copyable
    CF4(const CF4 &);
    CF4& operator=(const CF4&);
  public:
    CF4(const int dim);
    double evaluate(const emp::vector<double> x);
  };

  /******************************************************************************
  * Composition Functions
  *****************************************************************************/
  /* Constructors */
  CFunction::CFunction() : 
    dimension_(-1), nofunc_(-1), C_(-1), lambda_(NULL), sigma_(NULL),
    bias_(NULL), O_(NULL), M_(NULL), weight_(NULL), lbound_(NULL),
    ubound_(NULL), fi_(NULL), z_(NULL), f_bias_(0), fmaxi_(NULL), 
    tmpx_(NULL), function_(NULL), rng_(-1) 
  {
  }

  CFunction::CFunction(const int & dim, const int & nofunc, emp::Random & random) : 
    dimension_(dim), nofunc_(nofunc), C_(2000.0), lambda_(NULL), 
    sigma_(NULL), bias_(NULL), O_(NULL), M_(NULL), weight_(NULL),
    lbound_(NULL), ubound_(NULL), fi_(NULL), z_(NULL), f_bias_(0),
    fmaxi_(NULL), tmpx_(NULL), function_(NULL), rng_(random)
  {
  }

  /* Destructor */
  CFunction::~CFunction()
  {
  }

  // Set the random number generator for this landscape
  // Useful when you want to use the composite functions
  void SetSeed(emp::Random & random) {
    
  }

  void CFunction::calculate_weights(const emp::vector<double> x)
  {
    double sum(0), maxi(emp::MIN_INT), maxindex(0);
    for (int i=0; i<nofunc_; ++i) {
      sum = 0.0;
      for (int j=0; j<dimension_; ++j) {
        sum += ( x[j] - O_[i][j] ) * ( x[j] - O_[i][j] );
      }
      weight_[i] = exp( -sum/(2.0 * dimension_ * sigma_[i] * sigma_[i]) );
      if (i==0) { maxi = weight_[i]; }
      if (weight_[i] > maxi) {
        maxi = weight_[i];
        maxindex = i;
      }
    }
    sum = 0.0;
    for (int i=0; i<nofunc_; ++i) {
      //if (weight_[i] != maxi) {
      if (i != maxindex) {
        weight_[i] *= (1.0 - pow(maxi, 10.0));
      }
      sum += weight_[i];
    }
    for (int i=0; i<nofunc_; ++i) {
      if (sum == 0.0) {
        weight_[i] = 1.0/(double)nofunc_;
      } else {
        weight_[i] /= sum;
      }
    }
  }

  // Load specified optima from an outside .dat file
  void CFunction::load_optima(const std::string &filename)
  {
    std::fstream file;
    file.open(filename.c_str(), std::fstream::in);
    if (!file.is_open()) {
      std::cerr<< "Error: Can not open file: " << filename << std::endl;
      exit(0);
    }
    double tmp;
    for (int i=0; i< nofunc_; ++i) {
      for (int j=0; j< dimension_; ++j) {
        file >> tmp; 
        O_[i][j] = tmp;
      }
    }
    file.close();
  }

  // Load specified rotation matrix from an outside .dat file
  void CFunction::load_rotmat(const std::string &filename)
  {
    std::fstream file;
    file.open(filename.c_str(), std::fstream::in);
    if (!file.is_open()) {
      std::cerr<< "Error: Can not open file: " << filename << std::endl;
      exit(0);
    }
    double tmp(-1);
    for (int i=0; i<nofunc_; ++i) {
      for (int j=0; j<dimension_; ++j) {
        for (int k=0; k<dimension_; ++k) {
          file >> tmp; 
          M_[i][j][k] = tmp;
        }
      }
    }
    file.close();
  }
 
 // 3-Dimensional identity matrix 
 // One "layer" per function to be composited
 // Each "layer" of the matrix is an ID matrix
  void CFunction::init_rotmat_identity() {
    for (int i=0; i<nofunc_; ++i) {
      for (int j=0; j<dimension_; ++j) {
        for (int k=0; k<dimension_; ++k) {
          M_[i][j][k] = (j==k ? 1 : 0 );
        }
      }
    }	
  }

  // 2-dimensional random optima matrix
  // Each row represents a function
  // Each column represents a random optima for each dimension of the function
  void CFunction::init_optima_rand()
  {	
    for (int i=0; i< nofunc_; ++i) {
      for (int j=0; j< dimension_; ++j) {
        O_[i][j] = lbound_[j] + (ubound_[j] - lbound_[j]) * rng_.GetDouble();
      }
    }
  }

  void CFunction::transform_to_z(const emp::vector<double> x, const int &index)
  {
    /* Calculate z_i = (x - o_i)/\lambda_i */
    for (int i=0; i<dimension_; ++i) {
      tmpx_[i] = (x[i] - O_[index][i])/lambda_[index];
    }
    /* Multiply z_i * M_i */
    for (int i=0; i<dimension_; ++i) {
      z_[i] = 0;
      for (int j=0; j<dimension_; ++j) {
        z_[i] += M_[index][j][i] * tmpx_[j];
      }
    }
  }

  void CFunction::transform_to_z_noshift(const emp::vector<double> x, const int &index)
  {
    /* Calculate z_i = (x - o_i)/\lambda_i */
    for (int i=0; i<dimension_; ++i) {
      tmpx_[i] = (x[i])/lambda_[index];
    }
    /* Multiply z_i * M_i */
    for (int i=0; i<dimension_; ++i) {
      z_[i] = 0;
      for (int j=0; j<dimension_; ++j) {
        z_[i] += M_[index][j][i] * tmpx_[j];
      }
    }
  }

  void CFunction::calculate_fmaxi()
  {
    /* functions */
    for (int i=0; i<nofunc_; ++i) assert(function_[i] != NULL);
    double *x5 = new double[dimension_];
    for (int i=0; i<dimension_; ++i) { x5[i] = 5 ; }

    for (int i=0; i<nofunc_; ++i) {
      transform_to_z_noshift(x5, i);
      fmaxi_[i] = (*function_[i])(z_, dimension_);
    }
    delete [] x5;
  }

  double CFunction::evaluate_inner_(const emp::vector<double> x)
  {
    double result(0);
    calculate_weights(x);
    for (int i=0; i<nofunc_; ++i) {
      transform_to_z(x, i);
      fi_[i] = (*function_[i])(z_, dimension_);
    }
    for (int i=0; i<nofunc_; ++i) {
      result += weight_[i]*( C_ * fi_[i] / fmaxi_[i] + bias_[i] );
    }

    return -result + f_bias_;
  }

  std::vector< std::vector<double> > CFunction::get_copy_of_goptima() const
  {
    assert(O_ != NULL && "O_ == NULL");
    std::vector< std::vector<double> > OO;

    for (int i=0; i< nofunc_; ++i) {
      std::vector<double> kk;
      for (int j=0; j< dimension_; ++j) {
        kk.push_back(O_[i][j]);
      }
      OO.push_back(kk);
    }
    return OO;
  }

  CF1::CF1(const int dim) : CFunction(dim, 6)
  {
    for (int i=0; i<nofunc_; ++i) {
      sigma_[i] = 1;
      bias_[i]  = 0;
      weight_[i]= 0;
    }
    lambda_[0] = 1.0; 
    lambda_[1] = 1.0; 
    lambda_[2] = 8.0; 
    lambda_[3] = 8.0; 
    lambda_[4] = 1.0/5.0; 
    lambda_[5] = 1.0/5.0;
    /* Lower/Upper Bounds */
    for (int i=0; i<dimension_; ++i) {
      lbound_[i] = -5.0;
      ubound_[i] = 5.0;
    }
    /* load optima */
    if (dimension_ == 2 || dimension_ == 3 || dimension_ == 5 
        || dimension_ == 10 || dimension_ == 20 ) {
      std::string fname;
      fname = "data/CF1_M_D" + std::to_string(dim) + "_opt.dat";
      load_optima(fname);
    } else { 
      init_optima_rand();
    }
    /* M_ Identity matrices */
    init_rotmat_identity();
    /* Initialize functions of the composition */
    function_[0] = function_[1] = &FGriewank;
    function_[2] = function_[3] = &FWeierstrass;
    function_[4] = function_[5] = &FSphere;
    calculate_fmaxi();
  }

  double CF1::evaluate(const emp::vector<double> x)
  {
    return evaluate_inner_(x);
  }

  CF2::CF2(const int dim) : CFunction(dim, 8)
  {
    for (int i=0; i<nofunc_; ++i) {
      sigma_[i] = 1.0;
      bias_[i]  = 0.0;
      weight_[i]= 0.0;
    }
    lambda_[0] = 1.0; 
    lambda_[1] = 1.0; 
    lambda_[2] = 10.0; 
    lambda_[3] = 10.0; 
    lambda_[4] = 1.0/10.0; 
    lambda_[5] = 1.0/10.0;
    lambda_[6] = 1.0/7.0;
    lambda_[7] = 1.0/7.0;
    /* Lower/Upper Bounds */
    for (int i=0; i<dimension_; ++i) {
      lbound_[i] = -5.0;
      ubound_[i] = 5.0;
    }
    /* load optima */
    if (dimension_ == 2 || dimension_ == 3 || dimension_ == 5 
        || dimension_ == 10 || dimension_ == 20 ) {
      std::string fname;
      fname = "data/CF2_M_D" + std::to_string(dim) + "_opt.dat";
      load_optima(fname);
    } else { 
      init_optima_rand();
    }
    /* M_ Identity matrices */
    init_rotmat_identity();

    /* Initialize functions of the composition */
    function_[0] = function_[1] = &FRastrigin;
    function_[2] = function_[3] = &FWeierstrass;
    function_[4] = function_[5] = &FGriewank;
    function_[6] = function_[7] = &FSphere;
    calculate_fmaxi();
  }

  double CF2::evaluate(const emp::vector<double> x)
  {
    return evaluate_inner_(x);
  }

  CF3::CF3(const int dim) : CFunction(dim, 6)
  {
    for (int i=0; i<nofunc_; ++i) {
      bias_[i]  = 0.0;
      weight_[i]= 0.0;
    }
    sigma_[0] = 1.0;
    sigma_[1] = 1.0;
    sigma_[2] = 2.0;
    sigma_[3] = 2.0;
    sigma_[4] = 2.0;
    sigma_[5] = 2.0;
    lambda_[0] = 1.0/4.0; 
    lambda_[1] = 1.0/10.0; 
    lambda_[2] = 2.0; 
    lambda_[3] = 1.0; 
    lambda_[4] = 2.0; 
    lambda_[5] = 5.0;
    /* Lower/Upper Bounds */
    for (int i=0; i<dimension_; ++i) {
      lbound_[i] = -5.0;
      ubound_[i] = 5.0;
    }
    /* load optima */
    if (dimension_ == 2 || dimension_ == 3 || dimension_ == 5 
        || dimension_ == 10 || dimension_ == 20 ) {
      std::string fname;
      fname = "data/CF3_M_D" + std::to_string(dim) + "_opt.dat";
      load_optima(fname);
      fname = "data/CF3_M_D" + std::to_string(dim) + ".dat";
      load_rotmat(fname);
    } else { 
      init_optima_rand();
      /* M_ Identity matrices */
      init_rotmat_identity();
    }
    /* Initialize functions of the composition */
    function_[0] = function_[1] = &FEF8F2;
    function_[2] = function_[3] = &FWeierstrass;
    function_[4] = function_[5] = &FGriewank;
    calculate_fmaxi();
  }

  double CF3::evaluate(const emp::vector<double> x)
  {
    return evaluate_inner_(x);
  }

  CF4::CF4(const int dim) : CFunction(dim, 8)
  {
    for (int i=0; i<nofunc_; ++i) {
      sigma_[i] = 1.0;
      bias_[i]  = 0.0;
      weight_[i]= 0.0;
    }
    sigma_[0] = 1.0;
    sigma_[1] = 1.0;
    sigma_[2] = 1.0;
    sigma_[3] = 1.0;
    sigma_[4] = 1.0;
    sigma_[5] = 2.0;
    sigma_[6] = 2.0;
    sigma_[7] = 2.0;
    lambda_[0] = 4.0; 
    lambda_[1] = 1.0; 
    lambda_[2] = 4.0; 
    lambda_[3] = 1.0; 
    lambda_[4] = 1.0/10.0; 
    lambda_[5] = 1.0/5.0;
    lambda_[6] = 1.0/10.0;
    lambda_[7] = 1.0/40.0;
    /* Lower/Upper Bounds */
    for (int i=0; i<dimension_; ++i) {
      lbound_[i] = -5.0;
      ubound_[i] = 5.0;
    }
    /* load optima */
    if (dimension_ == 2 || dimension_ == 3 || dimension_ == 5 
        || dimension_ == 10 || dimension_ == 20) {
      std::string fname;
      fname = "data/CF4_M_D" + std::to_string(dim) + "_opt.dat";
      load_optima(fname);
      fname = "data/CF4_M_D" + std::to_string(dim) + ".dat";
      load_rotmat(fname);
    } else {
      init_optima_rand();
      /* M_ Identity matrices */
      init_rotmat_identity();
    }
    /* Initialize functions of the composition */
    function_[0] = function_[1] = &FRastrigin;
    function_[2] = function_[3] = &FEF8F2;
    function_[4] = function_[5] = &FWeierstrass;
    function_[6] = function_[7] = &FGriewank;
    calculate_fmaxi();
  }

  double CF4::evaluate(const emp::vector<double> x)
  {
    return evaluate_inner_(x);
  }

};

#endif