//
// Created by eugene on 16/02/2021.
//

#ifndef ATPIDTASK_PID_MATCHING_PLOTEFFICIENCIES_HPP_
#define ATPIDTASK_PID_MATCHING_PLOTEFFICIENCIES_HPP_

#include "TEfficiencyHelper.hpp"

#include <TDirectory.h>
#include <TKey.h>

#include <regex>

inline void MakeFancy(TDirectory *dir) {
  dir->cd();

  TList *key_list = dir->GetListOfKeys();
  for (Int_t ik = 0; ik < key_list->GetEntries(); ++ik) {
    TKey *key = (TKey *) key_list->At(ik);

    if(TClass::GetClass(key->GetClassName()) == TEfficiency::Class()) {
      TList * obj_list = ProjectEfficiency((TEfficiency *) key->ReadObj());

      for (Int_t io = 0; io < obj_list->GetEntries(); ++io) {
        TObject *o = obj_list->At(io);
        o->Write(o->GetName(), TObject::kOverwrite);
      }
    }

  }

}

inline void PlotEfficiencies(const char *filename) {
  TFile f(filename, "UPDATE");
  if (f.IsZombie()) {
    std::cout << "f::IsZombie()" << std::endl;
    return;
  }

  const std::regex re_expr("^efficiency_.*$");

  TList *key_list = f.GetListOfKeys();
  for (Int_t ik = 0; ik < key_list->GetEntries(); ++ik) {
    TKey * key = (TKey *) key_list->At(ik);

    if (TClass::GetClass(key->GetClassName())->InheritsFrom(TDirectory::Class()) &&
        std::regex_match(key->GetName(), re_expr)) {
      std::cout << "Processing '" << key->GetName() << "'..." << std::endl;
      MakeFancy((TDirectory *) key->ReadObj());
    }
  } // keys

  f.Close();


}

#endif //ATPIDTASK_PID_MATCHING_PLOTEFFICIENCIES_HPP_
