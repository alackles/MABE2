/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  EvalNKVar.hpp
 *  @brief MABE Evaluation module for NK Landscapes which hardcodes metrics from rank epistasis and also allows 2 different landscapes. Hopefully soon depreciable. 
 */

#ifndef MABE_EVAL_NK_VAR_H
#define MABE_EVAL_NK_VAR_H

#include "../../core/MABE.hpp"
#include "../../core/Module.hpp"
#include "../../tools/NK.hpp"

#include "emp/datastructs/reference_vector.hpp"

namespace mabe {

  class EvalNKVar : public Module {
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
    std::string nk_prefix;
    std::string genome_file;

  public:
    EvalNKVar(mabe::MABE & control,
           const std::string & name="EvalNKVar",
           const std::string & desc="Module to evaluate bitstrings on a mixed NK Fitness Lanscape WITH rank epistasis baked in.",
           size_t _N=100, size_t _Ka=3, size_t _Kb=3,
           const std::string & _nktype="half",
           const std::string & _btrait="bits", const std::string & _ftrait="fitness", 
           const std::string & _mfile="mutants.csv", const std::string & _nkprefix="nk", const std::string & _gfile="ref_genome.csv")
      : Module(control, name, desc)
      , N(_N), K_a(_Ka), K_b(_Kb)
      , target_collect(control.GetPopulation(0))
      , nk_type(_nktype)
      , bits_trait(_btrait)
      , fitness_trait(_ftrait)
      , mutant_file(_mfile)
      , nk_prefix(_nkprefix)
      , genome_file(_gfile)
    {
      SetEvaluateMod(true);
    }
    ~EvalNKVar() { }

    // Setup member functions associated with this class.
    static void InitType(emplode::TypeInfo & info) {
      info.AddMemberFunction("EVAL",
                             [](EvalNKVar & mod, Collection list) { return mod.Evaluate(list); },
                             "Use NK landscape to evaluate all orgs in an OrgList.");
      info.AddMemberFunction("RESET",
                             [](EvalNKVar & mod) { mod.landscape_a.Config(mod.N, mod.K_a, mod.control.GetRandom()); mod.landscape_b.Config(mod.N, mod.K_b, mod.control.GetRandom());  return 0; },
                             "Regenerate the NK landscape with current N and K.");
    }
    
    void SetupConfig() override {
      LinkCollection(target_collect, "target", "Which population(s) should we evaluate?");
      LinkVar(N, "N", "Number of bits required in output");
      LinkVar(K_a, "K_a", "Number of bits used in each gene for landscape a");
      LinkVar(K_b, "K_b", "Number of bits used in each gene for landscape b");
      LinkVar(nk_type,  "nk_type", "Type of mixed landscape [half, merged]");
      LinkVar(bits_trait, "bits_trait", "Which trait stores the bit sequence to evaluate?");
      LinkVar(fitness_trait, "fitness_trait", "Which trait should we store NK fitness in?");
      LinkVar(mutant_file, "mutant_file", "Where should we save the information about mutants?");
      LinkVar(nk_prefix, "nk_prefix", "Where should we save NK landscapes (will be suffixed with _[a|b].csv)?");
      LinkVar(genome_file, "genome_file", "Where should we save the maximally performing (i.e. reference) genome (in case we happen to need it later)?");
    }

    void SetupModule() override {
      // Setup the traits.
      AddRequiredTrait<emp::BitVector>(bits_trait);
      AddOwnedTrait<double>(fitness_trait, "NK fitness value", 0.0);

      // Setup the fitness landscape.
      landscape_a.Config(N, K_a, control.GetRandom());  // Setup the fitness landscape.
      landscape_b.Config(N, K_b, control.GetRandom()); // Setup the 2nd fitness landscape

      // Create the half-landscape
    
      // Output the fitness landscape.
      std::string fname_a = nk_prefix + "_a.csv";
      std::string fname_b = nk_prefix + "_b.csv";
      PrintLandscape(landscape_a, fname_a);
      PrintLandscape(landscape_b, fname_b);
    
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

    double NKFitness(emp::BitVector & bitstring) {
      double fitness = 0;
      // Create double-length bitstring to solve wraparound problems
      bitstring.Resize(N*2);
      bitstring |= (bitstring << N);
      if (nk_type == "half") {
        const auto & bits_a = bitstring.Export(N/2, 0); // export first N/2 bits 
        const auto & bits_b = bitstring.Export(N/2, N/2); // export second N/2 bits
        double fitness_a = landscape_a.GetFitness(bits_a); // the chopped up bitstring will be duplicated again but that should be fine
        double fitness_b = landscape_b.GetFitness(bits_b);
        fitness = fitness_a + fitness_b;
      } else if (nk_type == "mixed") {
          for (size_t i = 0; i < N; i++) {
            if (i % 2 == 0) {
              const auto & bits_a = bitstring.Export(K_a+1, i); // export length K+1 bitstring starting at the index of interest
              size_t dec_a = bits_a.GetUInt(0);  // Convert bits representation to decimal representation for table lookup
              fitness += landscape_a.GetFitness(i/2, dec_a); // map evens to landscape A
            } else {
              const auto & bits_b = bitstring.Export(K_b+1, i); // export length K_b+1 bitstring starting at index of interest
              size_t dec_b = bits_b.GetUInt(0); // Convert bits representation to decimal representation for table lookup
              fitness += landscape_b.GetFitness((i-1)/2, dec_b); // map odds to landscape B
            }
          }
      }
      return fitness;
    }
    
    emp::BitVector max_bits;
    
    double Evaluate(const Collection & orgs ) {

      // Loop through the population and evaluate each organism.
      double max_fitness = 0.0;
      emp::Ptr<Organism> max_org = nullptr;
      mabe::Collection alive_orgs( orgs.GetAlive() );
      for (Organism & org : alive_orgs) {
        org.GenerateOutput();
        auto & bits = org.GetTrait<emp::BitVector>(bits_trait);
        if (bits.size() != N) {
          emp::notify::Error("Org returns ", bits.size(), " bits, but ",
                   N, " bits needed for NK landscape.",
                   "\nOrg: ", org.ToString());
        }
        double fitness = NKFitness(bits);
        org.SetTrait<double>(fitness_trait, fitness);

        if (fitness > max_fitness || !max_org) {
          max_fitness = fitness;
          max_org = &org;
          max_bits = bits;
        }
      }
      return max_fitness;
    }


    double Evaluate(Population & pop) {return Evaluate( Collection(pop));}

    double Evaluate(const std::string & in) { return Evaluate( control.ToCollection(in) ); }

    // Rank epistasis analysis on the final population
    // 
    
    void BeforeExit() override {
      PrintGenome(max_bits);
      // use the existing landscape to evaluate and output the mutants
      std::ofstream mutFile(mutant_file);
      mutFile << "org_ID,pos_REF,pos_MUT,score_REF,score_MUT,\n";
      int org_id = 0;
      const auto & bits = max_bits;
      for (size_t i = 0; i < N ; ++i) {
        int pos_ref = i;
        auto genome = bits;
        genome.Toggle(i); 
        // get fitness of org with single mutation (i)
        double fitness_ref = NKFitness(genome);
        for (size_t j = 0 ; j < N ; ++j) {
          if (j != i) {
            int pos_mut = j;
            genome.Toggle(j);
            // get fitness of org with dual mutations (i and j)
            double fitness_mut = NKFitness(genome);
            genome.Toggle(j);
            mutFile << org_id << "," << pos_ref << "," << pos_mut << "," << fitness_ref << "," << fitness_mut << "," << "\n";
          }
        }
        genome.Toggle(i);
      }
      mutFile.close();
    }
      
  };

  MABE_REGISTER_MODULE(EvalNKVar, "Evaluate bitstrings on a Mixed NK fitness lanscape with Rank Epistasis.");
}

#endif
