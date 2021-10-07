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
#include <functional>

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

  using comp_func_t = std::function<double(const emp::vector<double> &, const size_t dims)>;

  /* Composition Functions framework */
  class CFunction {
  protected:
    size_t dim;
  	size_t numfunc;
    emp::Random rng;
  	
    double C;
  	
  	emp::vector<double> lambda;
  	emp::vector<double> sigma;
  	emp::vector<double> bias;
   
  	emp::vector<emp::vector<double>> O;
    emp::vector<emp::vector<emp::vector<double>>> M;
   
  	emp::vector<double> weight;
  	emp::vector<double> lbound;
  	emp::vector<double> ubound;
  	emp::vector<double> fi;
  	emp::vector<double> z;
    double fbias;
  	emp::vector<double> fmaxi;
  	emp::vector<double> tmpx;
   
    emp::vector<comp_func_t> function;

    // internal helper functions
  	void init_rotmat_identity();
  	void init_optima_rand(emp::Random & random);
  	void load_optima(const std::string &filename);
  	void load_rotmat(const std::string &filename);
  	void calculate_weights(const emp::vector<double> x);
  	void transform_to_z(const emp::vector<double> x, const int &index);
  	void transform_to_z_noshift(const emp::vector<double> x, const int &index);
  	void calculate_fmaxi();

    // initialize composite functions
    comp_func_t FAckley(const emp::vector<double> x, const size_t & dim);
    comp_func_t FRastrigin(const emp::vector<double> x, const size_t & dim);
    comp_func_t FWeierstrass(const emp::vector<double> x, const size_t & dim);
    comp_func_t FGriewank(const emp::vector<double> x, const size_t & dim);
    comp_func_t FSphere(const emp::vector<double> x, const size_t & dim);
    comp_func_t FSchwefel(const emp::vector<double> x, const size_t & dim);
    comp_func_t FRosenbrock(const emp::vector<double> x, const size_t & dim);
    comp_func_t FEF8F2(const emp::vector<double> x, const size_t & dim);
  
  public:
  	CFunction() 
    : dim(-1), numfunc(-1), rng(-1)
    , C(-1)
    , lambda(0), sigma(0), bias(0)
    , O(0), M(0)
    , weight(0), lbound(0), ubound(0)
    , fi(0), z(0), fbias(0)
    , fmaxi(0), tmpx(0), function(0)
    { ; }
    CFunction(const CFunction &) = default;
    CFunction(CFunction &&) = default;

    // dimensions are the number of dimensions each function should have
    // numfunc is the number of functions to composite
    // rng is the random number generator used to generate this landscape
    // I actually do not know what C is?
  	CFunction(size_t _dim, size_t _numfunc, emp::Random & random) 
    : dim(_dim), numfunc(_numfunc), rng(random)
    , C(2000.0)
    , lambda(8), sigma(8), bias(8)
    , O(_dim * 8), M(_dim * _dim * 8)
    , weight(8), lbound(8), ubound(8)
    , fi(8), z(8), fbias(8)
    , fmaxi(8), tmpx(8), function(8)
    { ; }
    CFunction & operator=(const CFunction &) = delete;
    CFunction & operator=(CFunction &&) = default;

  	double GetLower(const int &ivar) const { return lbound[ivar]; } 
  	double GetUpper(const int &ivar) const { return ubound[ivar]; }
    emp::vector<emp::vector<double>> GetMaxima() const { return O; } 

    void Config(const size_t _dim, const size_t _numfunc, emp::Random & _rng) {
      dim = _dim;
      rng = _rng;
      numfunc = _numfunc;
      lbound.assign(dim, -5.0);
      ubound.assign(dim, 5.0);
    }

    double GetFitness(const emp::vector<double> x) {
      double result = 0;
      calculate_weights(x);
      for (size_t i=0; i<numfunc; ++i) {
        transform_to_z(x, i);
        fi[i] = (function[i])(z, dim);
      }
      for (size_t i=0; i<numfunc; ++i) {
        result += weight[i]*( C * fi[i] / fmaxi[i] + bias[i] );
      }
      return -result + fbias;
    }

  };

  /* Interfaces for Composition functions */
  class CF1 : public CFunction {
  public:
    CF1();
    void Config(const size_t _dim, const size_t _numfunc, emp::Random & _rng) {
      CFunction::Config(_dim, _numfunc, _rng);
      sigma.assign(numfunc, 1.0);
      lambda = {1.0, 1.0, 8.0, 8.0, 1.0/5.0, 1.0/5.0};
      /* load optima */
      if (dim == 2 || dim == 3 || dim == 5 
          || dim == 10 || dim == 20 ) {
        std::string fname;
        fname = "DataGECCO/CF1_M_D" + std::to_string(dim) + "_opt.dat";
        load_optima(fname);
      } else { 
        init_optima_rand(rng);
      }
      /* M_ Identity matrices */
      init_rotmat_identity();
      /* Initialize functions of the composition */
      //function[0] = function[1] = FGriewank;
      //function[2] = function[3] = FWeierstrass;
      //function[4] = function[5] = FSphere;
      //calculate_fmaxi();
    };
  };

  class CF2 : public CFunction {
  public:
    CF2();
    void Config(const size_t _dim, const size_t _numfunc, emp::Random & _rng) {
      CFunction::Config(_dim, _numfunc, _rng);
      sigma.assign(numfunc, 1.0);
      lambda = {1.0, 1.0, 10.0, 10.0, 1.0/10.0, 1.0/10.0, 1.0/7.0, 1.0/7.0};
      /* load optima */
      if (dim == 2 || dim == 3 || dim == 5 
          || dim == 10 || dim == 20 ) {
        std::string fname;
        fname = "DataGECCO/CF2_M_D" + std::to_string(dim) + "_opt.dat";
        load_optima(fname);
      } else { 
        init_optima_rand(rng);
      }
      /* M_ Identity matrices */
      init_rotmat_identity();

      /* Initialize functions of the composition */
      //function[0] = function[1] = FRastrigin;
      //function[2] = function[3] = FWeierstrass;
      //function[4] = function[5] = FGriewank;
      //function[6] = function[7] = FSphere;
      //calculate_fmaxi();
    }
  };

  class CF3 : public CFunction {
  public:
    CF3();
    // Set up the composition
    void Config(const size_t _dim, const size_t _numfunc, emp::Random & _rng) {
      CFunction::Config(_dim, _numfunc, _rng);
      sigma = {1.0, 1.0, 2.0, 2.0, 2.0, 2.0};
      lambda = {1.0/4.0, 1.0/10.0, 2.0, 1.0, 2.0, 5.0};
      /* load optima */
      if (dim == 2 || dim == 3 || dim == 5 
          || dim == 10 || dim == 20 ) {
        std::string fname;
        fname = "DataGECCO/CF3_M_D" + std::to_string(dim) + "_opt.dat";
        load_optima(fname);
        fname = "DataGECCO/CF3_M_D" + std::to_string(dim) + ".dat";
        load_rotmat(fname);
      } else { 
        init_optima_rand(rng);
        /* M_ Identity matrices */
        init_rotmat_identity();
      }
      /* Initialize functions of the composition */
      //function[0] = function[1] = FEF8F2;
      //function[2] = function[3] = FWeierstrass;
      //function[4] = function[5] = FGriewank;
      //calculate_fmaxi();
    }
  };

  class CF4 : public CFunction {
  public:
    CF4();
    void Config(const size_t _dim, const size_t _numfunc, emp::Random & _rng) {
      CFunction::Config(_dim, _numfunc, _rng);
      sigma = {1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
      lambda = {4.0, 1.0, 4.0, 1.0, 1.0/10.0, 1.0/5.0, 1.0/10.0, 1.0/40.0};
      /* load optima */
      if (dim == 2 || dim == 3 || dim == 5 
          || dim == 10 || dim == 20) {
        std::string fname;
        fname = "DataGECCO/CF4_M_D" + std::to_string(dim) + "_opt.dat";
        load_optima(fname);
        fname = "DataGECCO/CF4_M_D" + std::to_string(dim) + ".dat";
        load_rotmat(fname);
      } else {
        init_optima_rand(rng);
        /* M_ Identity matrices */
        init_rotmat_identity();
      }
      /* Initialize functions of the composition */
      function[0] = FRastrigin; 
      //function[1] = FRastrigin;
      //function[2] = function[3] = FEF8F2;
      //function[4] = function[5] = FWeierstrass;
      //function[6] = function[7] = FGriewank;
      //calculate_fmaxi();
    };
  };

  /* Basic Benchmark functions */
  // NOTE: Each can take a dimension as a second argument, but this dimension does nothing.
  // It is just to conform the benchmark function arguments to the composite function arguments.

  /******************************************************************************
    * F1: Five-Uneven-Peak Trap 
    * Variable ranges: x in [0, 30
    * No. of global peaks: 2
    * No. of local peaks:  3. 
    *****************************************************************************/
  double five_uneven_peak_trap(const emp::vector<double> x, const size_t) {
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
  double equal_maxima(const emp::vector<double> x, const size_t) {
    double s = sin(5.0 * emp::PI * x[0]);
    return pow(s, 6);
  }

  /******************************************************************************
  * F3: Uneven Decreasing Maxima
  * Variable ranges: x in [0, 1]
  * No. of global peaks: 1
  * No. of local peaks:  4. 
  *****************************************************************************/
  double uneven_decreasing_maxima(const emp::vector<double> x, const size_t) {
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
  double himmelblau(const emp::vector<double> x, const size_t) {
    return 200 - (x[0]*x[0] + x[1] - 11)*(x[0]*x[0] + x[1] - 11) - 
    (x[0] + x[1]*x[1] - 7)*(x[0] + x[1]*x[1] - 7);
  }  	

  /******************************************************************************
  * F5: Six-Hump Camel Back
  * Variable ranges: x in [−1.9, 1.9]; y in [−1.1, 1.1]
  * No. of global peaks: 2
  * No. of local peaks:  2.
  *****************************************************************************/
  double six_hump_camel_back(const emp::vector<double> x, const size_t) {
    return -( (4 - 2.1*x[0]*x[0] + pow(x[0],4.0)/3.0)*x[0]*x[0] + 
    x[0]*x[1] + (4*x[1]*x[1] -4)*x[1]*x[1] );
  }

  /******************************************************************************
  * F6: Shubert
  * Variable ranges: x_i in  [−10, 10]^n, i=1,2,...,n
  * No. of global peaks: n*3^n
  * No. of local peaks: many
  *****************************************************************************/
  double shubert(const emp::vector<double> x, const size_t &dim) {
    double result(1), sum(0); 
    for (size_t i=0; i<dim; i++) {
      sum=0;
      for (size_t j=1; j<6; j++) {
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
  double vincent(const emp::vector<double> x, const size_t &dim) {
    double result(0);
    for (size_t i=0; i<dim; i++){
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

  double modified_rastrigin_all(const emp::vector<double> x, const size_t &dim)
  {
    double result(0);
    for (size_t i=0; i<dim; i++) {
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
  comp_func_t FAckley = [](const emp::vector<double> x, const size_t &dim) {
    double sum1(0.0), sum2(0.0), result;
    for (size_t i=0; i<dim; ++i) {
      sum1 += x[i]*x[i];
      sum2 += cos(2.0*emp::PI*x[i]);
    }
    sum1 = -0.2*sqrt(sum1/dim);
    sum2 /= dim;
    result = 20.0 + emp::E - 20.0*exp(sum1) - exp(sum2);
    return result;
  }; 

  /* Rastrigin's function */
  comp_func_t FRastrigin = [](const emp::vector<double> x, const size_t &dim) {
    double result(0.0);
    for (size_t i=0; i<dim; ++i) {
      result += (x[i]*x[i] - 10.0*cos(2.0*emp::PI*x[i]) + 10.0);
    }
    return result;
  };

  /* Weierstrass's function */
  comp_func_t FWeierstrass = [](const emp::vector<double> x, const size_t &dim) {
    double result(0.0), sum(0.0), sum2(0.0), a(0.5), b(3.0);
    int k_max(20);
    for (int j=0; j<=k_max; ++j) {
      sum2 += pow(a,j)*cos(2.0*emp::PI*pow(b,j)*(0.5));
    }
    for (size_t i=0; i<dim; ++i) {
      sum = 0.0;
    for (int j=0; j<=k_max; ++j) {
      sum += pow(a,j)*cos(2.0*emp::PI*pow(b,j)*(x[i]+0.5));
    }
    result += sum;
    }
    return result - sum2*dim;
  };

  /* Griewank's function */
  comp_func_t FGriewank = [](const emp::vector<double> x, const size_t &dim) {
    double sum(0.0), prod(1.0), result(0.0);
    for (size_t i=0; i<dim; ++i) {
      sum  += x[i]*x[i]/4000.0;
      prod *= cos( x[i]/sqrt(double(1.0+i)) );
    }
    result = 1.0 + sum - prod;
    return result;
  };

  /* Sphere function */
  comp_func_t FSphere = [](const emp::vector<double> x, const size_t &dim) {
    double result(0.0);
    for (size_t i=0; i<dim; ++i) {
      result += x[i]*x[i];
    }
    return result;
  };

  /* Schwefel's function */
  comp_func_t FSchwefel = [](const emp::vector<double> x, const size_t &dim) {
    double sum1(0.0), sum2(0.0);
    for (size_t i=0; i<dim; ++i) {
      sum2 = 0.0;
      for (size_t j=0; j<=i; ++j) {
        sum2 += x[j];
      }
      sum1 += sum2*sum2;
    }
    return sum1;
  };

  /* Rosenbrock's function */
  comp_func_t FRosenbrock = [](const emp::vector<double> x, const size_t &dim) {
    double result(0.0);
    for (size_t i=0; i<dim-1; ++i) {
      result += 100.0*pow((x[i]*x[i]-x[i+1]),2.0) + 1.0*pow((x[i]-1.0),2.0);
    }
    return result;
  };

  /* FEF8F2 function */
  comp_func_t FEF8F2 = [](const emp::vector<double> xx, const size_t &dim) {
    double result(0.0);
    double x(0), y(0), f(0), f2(0);
    for (size_t i=0; i<dim-1; ++i) {
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
  };


  /******************************************************************************
  * Composition Function Helpers
  *****************************************************************************/
  void CFunction::calculate_weights(const emp::vector<double> x) {
    double sum(0), maxi(emp::MIN_INT), maxindex(0);
    for (size_t i=0; i<numfunc; ++i) {
      sum = 0.0;
      for (size_t j=0; j<dim; ++j) {
        sum += ( x[j] - O[i][j] ) * ( x[j] - O[i][j] );
      }
      weight[i] = exp( -sum/(2.0 * dim * sigma[i] * sigma[i]) );
      if (i==0) { maxi = weight[i]; }
      if (weight[i] > maxi) {
        maxi = weight[i];
        maxindex = i;
      }
    }
    sum = 0.0;
    for (size_t i=0; i<numfunc; ++i) {
      //if (weight[i] != maxi) {
      if (i != maxindex) {
        weight[i] *= (1.0 - pow(maxi, 10.0));
      }
      sum += weight[i];
    }
    for (size_t i=0; i<numfunc; ++i) {
      if (sum == 0.0) {
        weight[i] = 1.0/(double)numfunc;
      } else {
        weight[i] /= sum;
      }
    }
  }

  // Load specified optima from an outside .dat file
  void CFunction::load_optima(const std::string &filename) {
    std::fstream file;
    file.open(filename.c_str(), std::fstream::in);
    if (!file.is_open()) {
      std::cerr<< "Error: Can not open file: " << filename << std::endl;
      exit(0);
    }
    double tmp;
    for (size_t i=0; i< numfunc; ++i) {
      for (size_t j=0; j< dim; ++j) {
        file >> tmp; 
        O[i][j] = tmp;
      }
    }
    file.close();
  }

  // Load specified rotation matrix from an outside .dat file
  void CFunction::load_rotmat(const std::string &filename) {
    std::fstream file;
    file.open(filename.c_str(), std::fstream::in);
    if (!file.is_open()) {
      std::cerr<< "Error: Can not open file: " << filename << std::endl;
      exit(0);
    }
    double tmp(-1);
    for (size_t i=0; i<numfunc; ++i) {
      for (size_t j=0; j<dim; ++j) {
        for (size_t k=0; k<dim; ++k) {
          file >> tmp; 
          M[i][j][k] = tmp;
        }
      }
    }
    file.close();
  }
 
 // 3-Dimensional identity matrix 
 // One "layer" per function to be composited
 // Each "layer" of the matrix is an ID matrix
  void CFunction::init_rotmat_identity() {
    for (size_t i=0; i<numfunc; ++i) {
      for (size_t j=0; j<dim; ++j) {
        for (size_t k=0; k<dim; ++k) {
          M[i][j][k] = (j==k ? 1 : 0 );
        }
      }
    }	
  }

  // 2-dimensional random optima matrix
  // Each row represents a function
  // Each column represents a random optima for each dimension of the function
  void CFunction::init_optima_rand(emp::Random & random) {	
    for (size_t i=0; i< numfunc; ++i) {
      for (size_t j=0; j< dim; ++j) {
        O[i][j] = lbound[j] + (ubound[j] - lbound[j]) * random.GetDouble();
      }
    }
  }

  void CFunction::transform_to_z(const emp::vector<double> x, const int &index) {
    /* Calculate z_i = (x - o_i)/\lambdai */
    for (size_t i=0; i<dim; ++i) {
      tmpx[i] = (x[i] - O[index][i])/lambda[index];
    }
    /* Multiply z_i * M_i */
    for (size_t i=0; i<dim; ++i) {
      z[i] = 0;
      for (size_t j=0; j<dim; ++j) {
        z[i] += M[index][j][i] * tmpx[j];
      }
    }
  }

  void CFunction::transform_to_z_noshift(const emp::vector<double> x, const int &index) {
    /* Calculate z_i = (x - o_i)/\lambdai */
    for (size_t i=0; i<dim; ++i) {
      tmpx[i] = (x[i])/lambda[index];
    }
    /* Multiply z_i * M_i */
    for (size_t i=0; i<dim; ++i) {
      z[i] = 0;
      for (size_t j=0; j<dim; ++j) {
        z[i] += M[index][j][i] * tmpx[j];
      }
    }
  }

  void CFunction::calculate_fmaxi() {
    /* functions */
    for (size_t i=0; i<numfunc; ++i) emp_assert(function[i] != NULL);
    emp::vector<double> x5(dim);
    for (size_t i=0; i<dim; ++i) { 
      x5[i] = 5 ;
    }
    for (size_t i=0; i<numfunc; ++i) {
      transform_to_z_noshift(x5, i);
      fmaxi[i] = (function[i])(z, dim);
    }
  }

};

#endif