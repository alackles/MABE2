/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  EvalAddMult.hpp
 *  @brief MABE Evaluation module to compare additive and multiplicative models of epistasis. 
 */

#ifndef MABE_EVAL_ADDMULT_H
#define MABE_EVAL_ADDMULT_H

#include "../../core/MABE.hpp"
#include "../../core/Module.hpp"

#include "emp/datastructs/reference_vector.hpp"

namespace mabe {

  class EvalAddMult : public Module {
  private:
    size_t N;
    size_t n_add;

    std::string vals_trait;
    std::string fitness_trait;
    std::string mutant_file;

  public:
    EvalAddMult(mabe::MABE & control,
           const std::string & name="EvalAddMult",
           const std::string & desc="Module to evaluate bitstrings on an NK Fitness Lanscape",
           size_t _N=100, size_t _nadd=50, const std::string & _vtrait="vals", const std::string & _ftrait="fitness",
           const std::string & _mfile="mutants.csv")
      : Module(control, name, desc)
      , N(_N)
      , n_add(_nadd)
      , vals_trait(_vtrait)
      , fitness_trait(_ftrait)
      , mutant_file(_mfile)
    {
      SetEvaluateMod(true);
    }
    ~EvalAddMult() { }

    // Setup member functions associated with this class.
    static void InitType(emplode::TypeInfo & info) {
      info.AddMemberFunction("EVAL",
                             [](EvalAddMult & mod, Collection list) { return mod.Evaluate(list); },
                             "Evaluate all orgs in an OrgList.");
    }

    void SetupConfig() override {
      LinkVar(N, "N", "Number of values required in output");
      LinkVar(n_add, "n_add", "Number of values in output to add (vs. multiply)");
      LinkVar(vals_trait, "vals_trait", "Which trait stores the integer sequence to evaluate?");
      LinkVar(fitness_trait, "fitness_trait", "Which trait should we store NK fitness in?");
      LinkVar(mutant_file, "mutant_file", "Where should we store the mutants file?");
    }

    void SetupModule() override {
      // Setup the traits.
      AddRequiredTrait<emp::vector<int>>(vals_trait);
      AddOwnedTrait<int>(fitness_trait, "Landscape fitness value", 0);
    }

    int CalcFitness(const emp::vector<int> & vals, const int & num_add) {
      int sum = std::accumulate(vals.begin(), vals.begin() + num_add, 0);
      int prod = std::accumulate(vals.begin() + num_add, vals.end(), sum, std::multiplies<int>());
      int result = sum + prod;
      return result;
    }

    emp::vector<int> max_vals;

    double Evaluate(const Collection & orgs) {
      // Loop through the population and evaluate each organism.
      double max_fitness = 0.0;
      emp::Ptr<Organism> max_org = nullptr;
      mabe::Collection alive_orgs( orgs.GetAlive() );
      for (Organism & org : alive_orgs) {
        org.GenerateOutput();
        const auto & vals = org.GetTrait<emp::vector<int>>(vals_trait);
        if (vals.size() != N) {
          emp::notify::Error("Org returns ", vals.size(), " bits, but ",
                             N, " bits needed for NK landscape.",
                             "\nOrg: ", org.ToString());
        }
        int mult = vals.size() - n_add; // How many values to multiply by each other 
        double fitness = CalcFitness(vals, n_add); // add up the first N/2 values and multiply in the next N/2 values 
        org.SetTrait<double>(fitness_trait, fitness);

        if (fitness > max_fitness || !max_org) {
          max_fitness = fitness;
          max_org = &org;
          max_vals = vals;
        }
      }

      return max_fitness;
    }

    // If a population is provided to Evaluate, first convert it to a Collection.
    double Evaluate(Population & pop) { return Evaluate( Collection(pop) ); }

    // If a string is provided to Evaluate, convert it to a Collection.
    double Evaluate(const std::string & in) { return Evaluate( control.ToCollection(in) ); }
    
    void BeforeExit() override {
      // use the existing landscape to evaluate and output the mutants
      std::ofstream mutFile(mutant_file);
      mutFile << "org_ID,pos_REF,pos_MUT,score_REF,score_MUT,\n";
      int org_id = 0;
      const auto & vals = max_vals;
      for (size_t i = 0; i < N ; ++i) {
        int pos_ref = i;
        auto genome = vals;
        for (int i = 1; i < 10; ++i) {
          // mutate genome to single mutant
          // TOGGLE HERE
          // get fitness of org with single mutation (i)
          int fitness_ref = CalcFitness(genome, n_add);
          for (size_t j = 0 ; j < N ; ++j) {
            if (j != i) {
              int pos_mut = j;
              //mutate genome to double mutant
              // - toggle here
              // get fitness of org with dual mutations (i and j)
              int fitness_mut = CalcFitness(genome, n_add);
              //mutate genome back to original single mutant
              // - toggle here
              mutFile << org_id << "," << pos_ref << "," << pos_mut << "," << fitness_ref << "," << fitness_mut << "," << "\n";
            }
          }
          // mutate genome back to original 
          // toggle here 
      }
      mutFile.close();
    }
    }
  };

  MABE_REGISTER_MODULE(EvalAddMult, "Evaluate vectors of integers on a simple additive & multiplicative landscape.");
}

#endif
