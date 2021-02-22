//
// Created by eugene on 16/02/2021.
//

#ifndef ATPIDTASK_PID_MATCHING_PLOTEFFICIENCIES_HPP_
#define ATPIDTASK_PID_MATCHING_PLOTEFFICIENCIES_HPP_

#include "TEfficiencyHelper.hpp"

#include <TDirectory.h>
#include <TKey.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <regex>

inline void ProcessEfficiencyDir(TDirectory *dir) {
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

  TProfile2D * msim_sim_y_pt_prof = dynamic_cast<TProfile2D*>(gDirectory->Get("matched_sim_sim_y_pt_prof"));
  key_list = dir->GetListOfKeys();

  gStyle->SetOptStat(0);
  auto c_ratio = new TCanvas;
  std::string efficiency_ratio_pdf_name{Form("efficiency_ratio_%s.pdf", dir->GetName())};
  c_ratio->SetCanvasSize(800,600);
  c_ratio->SetBatch(true);
  c_ratio->Print((efficiency_ratio_pdf_name + "(").c_str(), "pdf");

  for (Int_t ik = 0; ik < key_list->GetEntries(); ++ik) {
    TKey *key = (TKey *) key_list->At(ik);
    const std::regex re_expr("^matched_sim_sim_centr_y_pt_(\\d+)$");
    std::smatch match_results;
    std::string key_name{key->GetName()};

    if (std::regex_search(key_name, match_results, re_expr)) {
      auto centrality_class = match_results.str(1);
      TProfile2D *ratio_msim_sim_centr_bin_avg = dynamic_cast<TProfile2D*>(key->ReadObj()->
          Clone(Form("ratio_msim_sim_y_pt_centr_%s_total", centrality_class.c_str())));
      ratio_msim_sim_centr_bin_avg->Divide(msim_sim_y_pt_prof);
      ratio_msim_sim_centr_bin_avg->SetTitle(Form("Ratio {%s} / {%s}", ratio_msim_sim_centr_bin_avg->GetTitle(), msim_sim_y_pt_prof->GetTitle()));
      ratio_msim_sim_centr_bin_avg->SetMinimum(0.9);
      ratio_msim_sim_centr_bin_avg->SetMaximum(1.1);
      ratio_msim_sim_centr_bin_avg->Write();

      c_ratio->cd();
      c_ratio->Clear();
      ratio_msim_sim_centr_bin_avg->DrawClone("colz");
      c_ratio->Print(efficiency_ratio_pdf_name.c_str(), "pdf");



      /* lookup vtx_sim_centr_y_pt_%d */
      TProfile2D *residue_vtx_sim_y_pt_centr_bin = dynamic_cast<TProfile2D*>(dir->Get(Form("vtx_sim_centr_y_pt_%s", centrality_class.c_str()))
          ->Clone(Form("residue_vtx_msim_centr_y_pt_%s", centrality_class.c_str())));
      residue_vtx_sim_y_pt_centr_bin->Add(msim_sim_y_pt_prof, -1.);
      residue_vtx_sim_y_pt_centr_bin->SetTitle("N (VtxTracks) - N (Matched SimTracks) / N (SimTracks)");
      residue_vtx_sim_y_pt_centr_bin->SetMinimum(-0.2);
      residue_vtx_sim_y_pt_centr_bin->SetMaximum(0.2);
      residue_vtx_sim_y_pt_centr_bin->Write();
    }
  } // keys

  c_ratio->Print((efficiency_ratio_pdf_name + ")").c_str(), "pdf");
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
      ProcessEfficiencyDir((TDirectory *) key->ReadObj());
    }
  } // keys

  f.Close();


}

#endif //ATPIDTASK_PID_MATCHING_PLOTEFFICIENCIES_HPP_
