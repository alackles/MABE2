/**
 *  @note This file is part of MABE, https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019
 *
 *  @file  Organism.h
 *  @brief A wrapper for dealing with organisms is a generic manner.
 *  @note Status: PLANNING
 */

#ifndef MABE_ORGANISM_H
#define MABE_ORGANISM_H

#include "base/assert.h"
#include "data/VarMap.h"
#include "tools/string_utils.h"

#include "OrganismType.h"

namespace mabe {

  class Organism {
  private:
    emp::VarMap var_map;
    emp::Ptr<OrgTypeBase> type_ptr;

  public:
    OrgTypeBase & GetType() { return *type_ptr; }
    const OrgTypeBase & GetType() const { return *type_ptr; }

    bool HasVar(const std::string & name) const { return var_map.Has(name); }
    template <typename T> T & GetVar(const std::string & name) { return var_map.Get<T>(name); }
    template <typename T> const T & GetVar(const std::string & name) const { return var_map.Get<T>(name); }

    template <typename T>
    void SetVar(const std::string & name, const T & value) {
      var_map.Set(name, value);
    }
  };

  concept OrganismWrapper : Organism {
    void MABE_Setup() override { }

    double GetFitness() { return (double) *this; }
    std::string ToString() { return emp::to_string(*this); }
    emp::Ptr<Organism> Clone() { return emp::NewPtr<OrganismWrapper<WRAPPED_T>>(*this); }
    int Mutate() { emp_assert(false, "No default Mutate() available."); return -1; }
  };

}
#endif