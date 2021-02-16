//
// Created by eugene on 16/02/2021.
//

#ifndef ATPIDTASK_PID_MATCHING_FANCYEFFICIENCY_HPP_
#define ATPIDTASK_PID_MATCHING_FANCYEFFICIENCY_HPP_

#include "TEfficiencyHelper.hpp"

#include <TDirectory.h>
#include <TKey.h>

void MakeFancy(TDirectory *dir) {
  for (auto o : *dir->GetListOfKeys()) {
    auto key = (TKey *) o;

    if(TClass::GetClass(key->GetClassName()) == TEfficiency::Class()) {
      auto obj_list = ProjectEfficiency((TEfficiency *) key->ReadObj());

      for (auto o : *obj_list) {
        o->Write();
      }

      delete obj_list;
    }

  }

}

#endif //ATPIDTASK_PID_MATCHING_FANCYEFFICIENCY_HPP_
