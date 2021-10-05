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
    std::string mutant_file;
    std::string nk_file;

  public:
    EvalNKRank(mabe::MABE & control,
           const std::string & name="EvalNKRank",
           const std::string & desc="Module to evaluate bitstrings on an NK Fitness Lanscape WITH rank epistasis baked in.",
           size_t _N=100, size_t _K=3, 
           const std::string & _btrait="bits", const std::string & _ftrait="fitness", 
           const std::string & _mfile="mutants.csv", const std::string & _nkfile="nk.csv")
      : Module(control, name, desc)
      , N(_N), K(_K)
      , target_collect(control.GetPopulation(0))
      , bits_trait(_btrait)
      , fitness_trait(_ftrait)
      , mutant_file(_mfile)
      , nk_file(_nkfile)
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
      LinkVar(mutant_file, "mutant_file", "Where should we save the information about mutants?");
      LinkVar(nk_file, "nk_file", "Where should we save the actual NK landscape?");
    }

    void SetupModule() override {
      // Setup the traits.
      AddRequiredTrait<emp::BitVector>(bits_trait);
      AddOwnedTrait<double>(fitness_trait, "NK fitness value", 0.0);

      // Setup the fitness landscape.
      landscape.Config(N, K, control.GetRandom());  // Setup the fitness landscape.
    
      // Output the fitness landscape.
      PrintLandscape(landscape);
    
    }

    void PrintLandscape(NKLandscape nk_landscape) {
      std::ofstream nkFile(nk_file);
      size_t rows = nk_landscape.GetStateCount();
      for (size_t r = 0; r < rows; ++r) {
        for (size_t n = 0; n < N; ++n) {
          nkFile << nk_landscape.GetFitness(n, r) << ",";
        }
        nkFile << "\n";
      }
      nkFile.close();
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
    // 
    void BeforeExit() override {
      // use the existing landscape to evaluate and output the mutants
      std::ofstream mutFile(mutant_file);
      mutFile << "org_ID,pos_REF,pos_MUT,score_REF,score_MUT,\n";
      int org_id = 0;
      const auto & bits = max_bits;
      for (int i = 0; i < N ; ++i) {
        int pos_ref = i;
        auto mutant = bits;
        mutant.Toggle(i); 
        // get fitness of org with single mutation (i)
        double fitness_ref = landscape.GetFitness(mutant);
        for (int j = 0 ; j < N ; ++j) {
          if (j != i) {
            int pos_mut = j;
            mutant.Toggle(j);
            // get fitness of org with dual mutations (i and j)
            double fitness_mut = landscape.GetFitness(mutant);
            mutant.Toggle(j);
            mutFile << org_id << "," << pos_ref << "," << pos_mut << "," << fitness_ref << "," << fitness_mut << "," << "\n";
          }
        }
        mutant.Toggle(i);
      }
      mutFile.close();
    }
      
  };

  MABE_REGISTER_MODULE(EvalNKRank, "Evaluate bitstrings on an NK fitness lanscape with Rank Epistasis.");
}

#endif
