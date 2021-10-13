/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  EvalNKMixed.hpp
 *  @brief MABE Evaluation module for NK Landscapes which hardcodes metrics from rank epistasis and also allows 2 different landscapes. Hopefully soon depreciable. 
 */

#ifndef MABE_EVAL_NK_MIXED_H
#define MABE_EVAL_NK_MIXED_H

#include "../../core/MABE.hpp"
#include "../../core/Module.hpp"
#include "../../tools/NK.hpp"

#include "emp/datastructs/reference_vector.hpp"

namespace mabe {

  class EvalNKMixed : public Module {
  private:
    size_t N;
    size_t K_a;
    size_t K_b;
    NKLandscape landscape_a;
    NKLandscape landscape_b;
    mabe::Collection target_collect;

    std::string nk_type; 
    std::string bits_trait;
    std::string fitness_trait;
    std::string mutant_file;
    std::string nk_file;
    std::string genome_file;

  public:
    EvalNKMixed(mabe::MABE & control,
           const std::string & name="EvalNKMixed",
           const std::string & desc="Module to evaluate bitstrings on a mixed NK Fitness Lanscape WITH rank epistasis baked in.",
           size_t _N=100, size_t _Ka=3, size_t _Kb=3,
           const std::string & _nktype="half",
           const std::string & _btrait="bits", const std::string & _ftrait="fitness", 
           const std::string & _mfile="mutants.csv", const std::string & _nkfile="nk.csv", const std::string & _gfile="ref_genome.csv")
      : Module(control, name, desc)
      , N(_N), K_a(_Ka), K_b(_Kb)
      , target_collect(control.GetPopulation(0))
      , nk_type(_nktype)
      , bits_trait(_btrait)
      , fitness_trait(_ftrait)
      , mutant_file(_mfile)
      , nk_file(_nkfile)
      , genome_file(_gfile)
    {
      SetEvaluateMod(true);
    }
    ~EvalNKMixed() { }

    void SetupConfig() override {
      LinkCollection(target_collect, "target", "Which population(s) should we evaluate?");
      LinkVar(N, "N", "Number of bits required in output");
      LinkVar(K_a, "K_a", "Number of bits used in each gene for landscape a");
      LinkVar(K_b, "K_b", "Number of bits used in each gene for landscape b");
      LinkVar(nk_type,  "nk_type", "Type of mixed landscape [half, merged]");
      LinkVar(bits_trait, "bits_trait", "Which trait stores the bit sequence to evaluate?");
      LinkVar(fitness_trait, "fitness_trait", "Which trait should we store NK fitness in?");
      LinkVar(mutant_file, "mutant_file", "Where should we save the information about mutants?");
      LinkVar(nk_file, "nk_file", "Where should we save the actual NK landscape?");
      LinkVar(genome_file, "genome_file", "Where should we save the maximally performing (i.e. reference) genome (in case we happen to need it later)?");
    }

    const size_t midpt = N/2; 
    emp::BitVector max_bits;
    
    void SetupModule() override {
      // Setup the traits.
      AddRequiredTrait<emp::BitVector>(bits_trait);
      AddOwnedTrait<double>(fitness_trait, "NK fitness value", 0.0);

      // Setup the fitness landscape.
      landscape_a.Config(N, K_a, control.GetRandom());  // Setup the fitness landscape.
      landscape_b.Config(N, K_b, control.GetRandom()); // Setup the 2nd fitness landscape

      // Create the half-landscape
    
      // Output the fitness landscape.
      //PrintLandscape(landscape_a);
      //PrintLandscape(landscape_b);
    
    }

    void PrintLandscape(NKLandscape nk_landscape, const std::string & nk_fname) {
      std::ofstream nkFile(nk_fname);
      size_t rows = nk_landscape.GetStateCount();
      for (size_t r = 0; r < rows; ++r) {
        for (size_t n = 0; n < N; ++n) {
          nkFile << nk_landscape.GetFitness(n, r) << ",";
        }
        nkFile << "\n";
      }
      nkFile.close();
    }

    void PrintGenome(emp::BitVector genome) {
      std::ofstream genomeFile(genome_file);
      genomeFile << genome.ToString();
      genomeFile.close();
    }
    
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

        double fitness = 0;
        if (nk_type == "half") {
          const auto & bits_a = bits.Export(midpt, 0);
          const auto & bits_b = bits.Export(midpt, midpt);
          double fitness_a = landscape_a.GetFitness(bits_a);
          double fitness_b = landscape_b.GetFitness(bits_b);
          fitness = fitness_a + fitness_b;
        } else if (nk_type == "mixed") {
          // to be added
        }
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
      PrintGenome(max_bits);
      // use the existing landscape to evaluate and output the mutants
      std::ofstream mutFile(mutant_file);
      mutFile << "org_ID,pos_REF,pos_MUT,score_REF,score_MUT,\n";
      int org_id = 0;
      const auto & bits = max_bits;
      double fitness_ref;
      double fitness_mut;
      for (size_t i = 0; i < N ; ++i) {
        int pos_ref = i;
        auto genome = bits;
        genome.Toggle(i); 
        // get fitness of org with single mutation (i)
        if (nk_type == "half") {
          auto ref_a = genome.Export(midpt, 0);
          auto ref_b = genome.Export(midpt, midpt);
          double ref_fitness_a = landscape_a.GetFitness(ref_a);
          double ref_fitness_b = landscape_b.GetFitness(ref_b);
          fitness_ref = ref_fitness_a + ref_fitness_b;
        } else if (nk_type == "mixed") {
          // tbd
        }
        for (size_t j = 0 ; j < N ; ++j) {
          if (j != i) {
            int pos_mut = j;
            genome.Toggle(j);
            // get fitness of org with dual mutations (i and j)
            if (nk_type == "half") {
              auto mut_a = genome.Export(midpt, 0);
              auto mut_b = genome.Export(midpt, midpt);
              double mut_fitness_a = landscape_a.GetFitness(mut_a);
              double mut_fitness_b = landscape_b.GetFitness(mut_b);
              double fitness_mut = mut_fitness_a + mut_fitness_b;
        } else if (nk_type == "mixed") {
          // tbd
        }
            genome.Toggle(j);
            mutFile << org_id << "," << pos_ref << "," << pos_mut << "," << fitness_ref << "," << fitness_mut << "," << "\n";
          }
        }
        genome.Toggle(i);
      }
      mutFile.close();
    }
      
  };

  MABE_REGISTER_MODULE(EvalNKMixed, "Evaluate bitstrings on a Mixed NK fitness lanscape with Rank Epistasis.");
}

#endif
