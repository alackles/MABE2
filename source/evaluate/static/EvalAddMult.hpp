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

    std::string vals_trait;
    std::string fitness_trait;

  public:
    EvalAddMult(mabe::MABE & control,
           const std::string & name="EvalAddMult",
           const std::string & desc="Module to evaluate bitstrings on an NK Fitness Lanscape",
           size_t _N=100, const std::string & _vtrait="vals", const std::string & _ftrait="fitness")
      : Module(control, name, desc)
      , N(_N)
      , vals_trait(_vtrait)
      , fitness_trait(_ftrait)
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
      LinkVar(vals_trait, "vals_trait", "Which trait stores the integer sequence to evaluate?");
      LinkVar(fitness_trait, "fitness_trait", "Which trait should we store NK fitness in?");
    }

    void SetupModule() override {
      // Setup the traits.
      AddRequiredTrait<emp::vector<int>>(vals_trait);
      AddOwnedTrait<int>(fitness_trait, "Landscape fitness value", 0);
    }

    int CalcFitness(const emp::vector<int> & vals, const int & n_add, const int & n_mult) {
      int result = 1;
      return result;
    }

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
        int add = vals.size()/2; // How many values to add together
        int mult = vals.size() - add; // How many values to multiply by each other 
        double fitness = CalcFitness(vals, add, mult); // add up the first N/2 values and multiply in the next N/2 values 
        org.SetTrait<double>(fitness_trait, fitness);

        if (fitness > max_fitness || !max_org) {
          max_fitness = fitness;
          max_org = &org;
        }
      }

      return max_fitness;
    }

    // If a population is provided to Evaluate, first convert it to a Collection.
    double Evaluate(Population & pop) { return Evaluate( Collection(pop) ); }

    // If a string is provided to Evaluate, convert it to a Collection.
    double Evaluate(const std::string & in) { return Evaluate( control.ToCollection(in) ); }
  };

  MABE_REGISTER_MODULE(EvalAddMult, "Evaluate vectors of integers on a simple additive & multiplicative landscape.");
}

#endif
