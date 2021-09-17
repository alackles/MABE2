/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  EvalNKRank.hpp
 *  @brief MABE Evaluation module for NK Landscapes which hardcodes metrics from rank epistasis. Hopefully soon depreciable. 
 */

#ifndef MABE_EVAL_NK_RANK_H
#define MABE_EVAL_NK_RANK_H

#include "../../core/MABE.hpp"
#include "../../core/Module.hpp"
#include "../../tools/NK.hpp"

#include "emp/datastructs/reference_vector.hpp"

namespace mabe {

  class EvalNKRank : public Module {
  private:
    size_t N;
    size_t K;    
    NKLandscape landscape;
    mabe::Collection target_collect;

    std::string bits_trait;
    std::string fitness_trait;
    std::string knockout_file;

  public:
    EvalNKRank(mabe::MABE & control,
           const std::string & name="EvalNKRank",
           const std::string & desc="Module to evaluate bitstrings on an NK Fitness Lanscape WITH rank epistasis baked in.",
           size_t _N=100, size_t _K=3, const std::string & _btrait="bits", const std::string & _ftrait="fitness", const std::string & _kfile="knockouts.csv")
      : Module(control, name, desc)
      , N(_N), K(_K)
      , target_collect(control.GetPopulation(0))
      , bits_trait(_btrait)
      , fitness_trait(_ftrait)
      , knockout_file(_kfile)
    {
      SetEvaluateMod(true);
    }
    ~EvalNKRank() { }

    void SetupConfig() override {
      LinkCollection(target_collect, "target", "Which population(s) should we evaluate?");
      LinkVar(N, "N", "Number of bits required in output");
      LinkVar(K, "K", "Number of bits used in each gene");
      LinkVar(bits_trait, "bits_trait", "Which trait stores the bit sequence to evaluate?");
      LinkVar(fitness_trait, "fitness_trait", "Which trait should we store NK fitness in?");
      LinkVar(knockout_file, "knockout_file", "Where should we save the information about knockouts?");
    }

    void SetupModule() override {
      // Setup the traits.
      AddRequiredTrait<emp::BitVector>(bits_trait);
      AddOwnedTrait<double>(fitness_trait, "NK fitness value", 0.0);

      // Setup the fitness landscape.
      landscape.Config(N, K, control.GetRandom());  // Setup the fitness landscape.
    }

    emp::BitVector max_bits;
    
    void OnUpdate(size_t /* update */) override {
      emp_assert(control.GetNumPopulations() >= 1);

      // Loop through the population and evaluate each organism.
      mabe::Collection alive_collect( target_collect.GetAlive() );
      double max_fitness = 0.0;
      emp::Ptr<Organism> max_org = nullptr;
      for (Organism & org : alive_collect) {
        org.GenerateOutput();
        const auto & bits = org.GetTrait<emp::BitVector>(bits_trait);
        if (bits.size() != N) {
          AddError("Org returns ", bits.size(), " bits, but ",
                   N, " bits needed for NK landscape.",
                   "\nOrg: ", org.ToString());
        }
        double fitness = landscape.GetFitness(bits);
        org.SetTrait<double>(fitness_trait, fitness);

        if (fitness > max_fitness || !max_org) {
          max_fitness = fitness;
          max_org = &org;
          max_bits = bits;
        }
      }
      std::cout << "Max " << fitness_trait << " = " << max_fitness << std::endl;
    }

    // Rank epistasis analysis on the final population
    void BeforeExit() override {
      std::ofstream kfileout(knockout_file);
      kfileout << "org_ID,mt_pos,ko_pos,score_MT,score_KO,\n";
      int org_id = 0;
      std::cout << "exit test " << std::endl;
      const auto & bits = max_bits;
      for (int i = 0; i < N ; ++i) {
        int mt_pos = i;
        auto knockout = bits;
        knockout.Toggle(i); 
        double mt_fitness = landscape.GetFitness(knockout);
        for (int j = 0 ; j < N ; ++j) {
          if (j != i) {
            int ko_pos = j;
            knockout.Toggle(j);
            double ko_fitness = landscape.GetFitness(knockout);
            knockout.Toggle(j);
            kfileout << org_id << "," << mt_pos << "," << ko_pos << "," << mt_fitness << "," << ko_fitness << "," << "\n";
            std::cout << "fitness: " << ko_fitness << std::endl;
          }
        }
        knockout.Toggle(i);
      }
      kfileout.close();
    }
      
  };

  MABE_REGISTER_MODULE(EvalNKRank, "Evaluate bitstrings on an NK fitness lanscape with Rank Epistasis.");
}

#endif
