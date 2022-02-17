/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  EvalGECCO.hpp
 *  @brief MABE Evaluation module for (a subset of) the GECCO Niching Competition Landscapes. Find the landscapes here: http://epitropakis.co.uk/gecco2021/
 */

#ifndef MABE_EVAL_GECCO_H
#define MABE_EVAL_GECCO_H

#include "../../core/MABE.hpp"
#include "../../core/Module.hpp"
#include "../../tools/GECCO.hpp"

#include "emp/datastructs/reference_vector.hpp"

namespace mabe {

  class EvalGECCO : public Module {
  private:
    size_t dims;
    CF3 comp3;
    CF4 comp4;
    mabe::Collection target_collect;
    std::string fcn_name;
    std::string vals_trait;
    std::string fitness_trait;
    std::string genome_file;
    std::string dat_path;

  public:
    EvalGECCO(mabe::MABE & control,
           const std::string & name="EvalGECCO",
           const std::string & desc="Module to evaluate bitstrings on one of the GECCO Niching Competition 3D landscapes.",
           const size_t & _dims=3,
           const std::string & _fname="Shubert",
           const std::string & _vtrait="vals", const std::string & _ftrait="fitness",
           const std::string & _gfile="genome.csv",
           const std::string & _dpath="./../source/tools/DataGECCO/")
      : Module(control, name, desc)
      , target_collect(control.GetPopulation(0))
      , fcn_name(_fname)
      , dims(_dims)
      , vals_trait(_vtrait)
      , fitness_trait(_ftrait)
      , genome_file(_gfile)
      , dat_path(_dpath)
    {
      SetEvaluateMod(true);
    }
    ~EvalGECCO() { }

    void SetupConfig() override {
      LinkCollection(target_collect, "target", "Which population(s) should we evaluate?");
      LinkVar(fcn_name, "fcn_name", "Which function should we use? [Shubert, Vincent, CF3, CF4]");
      LinkVar(dims, "dims", "How many dimensions should the fitness function have? [2, 3]");
      LinkVar(vals_trait, "vals_trait", "Which trait stores the 3-tuple to evaluate?");
      LinkVar(fitness_trait, "fitness_trait", "Which trait should we store fitness in?");
      LinkVar(genome_file, "genome_file", "Where should we output the genome?");
      LinkVar(dat_path, "dat_path", "Where do we pull the .dat files for the composite functions?");
    }

    void SetupModule() override {
      // Setup the traits.
      AddRequiredTrait<emp::vector<double>>(vals_trait);
      AddOwnedTrait<double>(fitness_trait, "Landscape fitness value", 0.0);
      
      // set up composite functions
      if (fcn_name == "CF3") {
        comp3.SetDataPath(dat_path);
        comp3.Config(dims, control.GetRandom());
      }
      if (fcn_name == "CF4") {
        comp3.SetDataPath(dat_path);
        comp4.Config(dims, control.GetRandom());
      }
      
    }
  
    void OnUpdate(size_t /* update */) override {
      emp_assert(control.GetNumPopulations() >= 1);

      // Loop through the population and evaluate each organism.
      double max_fitness = 0.0;
      emp::vector<double> max_val(dims);
      emp::Ptr<Organism> max_org = nullptr;
      mabe::Collection alive_collect( target_collect.GetAlive() );
      for (Organism & org : alive_collect) {
        org.GenerateOutput();
        const auto & val = org.GetTrait<emp::vector<double>>(vals_trait);
        double fitness = 0.0;
        if (fcn_name == "Shubert") {
          fitness = gecco::shubert(val, dims);
        } else if (fcn_name == "Vincent") {
          fitness = gecco::vincent(val, dims);
        } else if (fcn_name == "CF3") {
          fitness = comp3.GetFitness(val);
        } else if (fcn_name == "CF4") {
          fitness = comp4.GetFitness(val);
        } else {
          std::cout << "Invalid function name:" << fcn_name << std::endl;
        }

        // Store the count on the organism in the fitness trait.
        org.SetTrait<double>(fitness_trait, fitness);
        
        if (fitness > max_fitness || !max_org) {
          max_fitness = fitness;
          max_org = &org;
          max_val = val;
        }
      }
      std::cout << "Max " << fitness_trait << " = " << max_fitness << std::endl;
    }
  };

  MABE_REGISTER_MODULE(EvalGECCO, "Evaluate bitstrings on a 3D GECCO niching lanscape.");
}

#endif
