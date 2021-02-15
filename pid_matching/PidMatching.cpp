//
// Created by eugene on 09/02/2021.
//

#include "PidMatching.hpp"
#include <boost/program_options.hpp>
#include <AnalysisTree/Matching.hpp>
#include <AnalysisTree/DataHeader.hpp>
#include <pid_new/core/PdgHelper.h>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TH2.h>

TASK_IMPL(PidMatching)

using std::cout;
using std::endl;
using AnalysisTree::Matching;


struct PidMatching::PidEfficiencyQAStruct {
  TDirectory *output_dir{nullptr};

  TH2 *matched_tracks_y_pt{nullptr};
  TH2 *sim_tracks_y_pt{nullptr};

  TEfficiency *efficiency_y_pt_from_hists{nullptr};
  TEfficiency *efficiency_y_pt_filled{nullptr};

};

struct PidMatching::ChargedHadronsEfficiencyStruct {
  TDirectory *output_dir{nullptr};

  TEfficiency *eta_pt_vtx_tracks{nullptr};
  TEfficiency *eta_pt_vtx_tracks_neg{nullptr};
  TEfficiency *eta_pt_vtx_tracks_pos{nullptr};
};

boost::program_options::options_description PidMatching::GetBoostOptions() {
  namespace po = boost::program_options;
  cout << __func__ << endl;
  return {};
}

void PidMatching::PreInit() {
  cout << __func__ << endl;
//  SetInputBranchNames({"VtxTracks2SimTracks"});
//  Matching is stored in different field
}
void PidMatching::Init(std::map<std::string, void *> &map) {
  InitEfficiencies();

  matching_ptr_ = static_cast<Matching *>(map["VtxTracks2SimTracks"]);
  vtx_tracks_ptr = static_cast<AnalysisTree::TrackDetector *>(map["VtxTracks"]);
  sim_track_ptr = static_cast<AnalysisTree::TrackDetector *>(map["SimTracks"]);

  auto vtx_tracks_config = config_->GetBranchConfig("VtxTracks");

  /* input config */
  sim_track_pdg_id_ = VarId("SimTracks/pdg");

  i_dca_x_field_id_ = VarId("VtxTracks/dcax");
  i_dca_y_field_id_ = VarId("VtxTracks/dcay");

  i_nhits_vtpc1_ = VarId("VtxTracks/nhits_vtpc1");
  i_nhits_vtpc2_ = VarId("VtxTracks/nhits_vtpc2");
  i_nhits_mtpc_ = VarId("VtxTracks/nhits_mtpc");
  i_nhits_pot_vtpc1_ = VarId("VtxTracks/nhits_pot_vtpc1");
  i_nhits_pot_vtpc2_ = VarId("VtxTracks/nhits_pot_vtpc2");
  i_nhits_pot_mtpc_ = VarId("VtxTracks/nhits_pot_mtpc");
  i_charge = VarId("VtxTracks/q");

  /* output config */
  NewBranch("MatchedVtxTracks", AnalysisTree::DetType::kParticle);
  y_cm_field_id_ = NewVar<float>("MatchedVtxTracks/y_cm");
  y_field_id_ = NewVar<float>("MatchedVtxTracks/y");

  o_dca_x_field_id_ = NewVar<float>("MatchedVtxTracks/dcax");
  o_dca_y_field_id_ = NewVar<float>("MatchedVtxTracks/dcay");

  o_nhits_total_ = NewVar<int>("MatchedVtxTracks/nhits_total");
  o_nhits_vtpc_ = NewVar<int>("MatchedVtxTracks/nhits_vtpc");
  o_nhits_pot_total_ = NewVar<int>("MatchedVtxTracks/nhits_pot_total");
  o_nhits_ratio_ = NewVar<float>("MatchedVtxTracks/nhits_ratio");

  /* parameters of the matched sim track */
  o_sim_y_cm_ = NewVar<float>("MatchedVtxTracks/sim_y_cm");
  o_sim_pt_ = NewVar<float>("MatchedVtxTracks/sim_pt");
  o_sim_phi_ = NewVar<float>("MatchedVtxTracks/sim_phi");

  matched_particles_config_ = out_config_->GetBranchConfig("MatchedVtxTracks");

  matched_particles_ = new AnalysisTree::Particles;
  out_tree_->Branch("MatchedVtxTracks", &matched_particles_);
}
void PidMatching::Exec() {
  matched_particles_->ClearChannels();

  TLorentzVector sim_momentum;
  TLorentzVector vtx_momentum;

  for (auto &[firstId, secondId] : matching_ptr_->GetMatches()) {
    auto vtx_track = vtx_tracks_ptr->GetChannel(firstId);
    auto sim_track = sim_track_ptr->GetChannel(secondId);

    auto particle = matched_particles_->AddChannel();
    particle->Init(matched_particles_config_);

    auto pdg = sim_track.GetField<int>(sim_track_pdg_id_);
    sim_momentum = sim_track.Get4Momentum(pdg);
    vtx_momentum.SetVectM(vtx_track.GetMomentum3(), PdgHelper::mass(pdg));
    particle->SetMomentum3(vtx_momentum.Vect());
    particle->SetPid(pdg);
    particle->SetMass(PdgHelper::mass(pdg));

    /* y_cm */
    particle->SetField<float>(vtx_momentum.Rapidity(), y_field_id_);
    particle->SetField<float>(vtx_momentum.Rapidity() - data_header_->GetBeamRapidity(), y_cm_field_id_);

    /* dca_x, dca_y */
    particle->SetField<float>(vtx_track.GetField<float>(i_dca_x_field_id_), o_dca_x_field_id_);
    particle->SetField<float>(vtx_track.GetField<float>(i_dca_y_field_id_), o_dca_y_field_id_);

    /* nhits and ratio */
    {
      int nhits_total = vtx_track.GetField<int>(i_nhits_vtpc1_) +
          vtx_track.GetField<int>(i_nhits_vtpc2_) +
          vtx_track.GetField<int>(i_nhits_mtpc_);
      int nhits_vtpc =
          vtx_track.GetField<int>(i_nhits_vtpc1_) +
              vtx_track.GetField<int>(i_nhits_vtpc2_);

      int nhits_pot_total = vtx_track.GetField<int>(i_nhits_pot_vtpc1_) +
          vtx_track.GetField<int>(i_nhits_pot_vtpc2_) +
          vtx_track.GetField<int>(i_nhits_pot_mtpc_);
      particle->SetField(nhits_total, o_nhits_total_);
      particle->SetField<int>(nhits_vtpc, o_nhits_vtpc_);
      particle->SetField(nhits_pot_total, o_nhits_pot_total_);
      particle->SetField(float(nhits_total) / float(nhits_pot_total), o_nhits_ratio_);
    }

    /* matched sim track kinematics */
    particle->SetField<float>(sim_momentum.Rapidity() - data_header_->GetBeamRapidity(), o_sim_y_cm_);
    particle->SetField<float>(sim_momentum.Pt(), o_sim_pt_);
    particle->SetField<float>(sim_momentum.Phi(), o_sim_phi_);

    if (efficiencies.find(pdg) != efficiencies.end()) {
      efficiencies[pdg]->matched_tracks_y_pt->Fill(vtx_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                                   vtx_momentum.Pt());

    }

  } // matched particles

  const auto match_inv = matching_ptr_->GetMatches(true);
  const auto match = matching_ptr_->GetMatches();

  for (size_t i_t = 0; i_t < sim_track_ptr->GetNumberOfChannels(); ++i_t) {
    auto channel = sim_track_ptr->GetChannel(i_t);

    auto pdg = channel.GetField<int>(sim_track_pdg_id_);
    if (efficiencies.find(pdg) != efficiencies.end()) {
      sim_momentum = channel.Get4Momentum(pdg);
      efficiencies[pdg]->sim_tracks_y_pt->Fill(sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                               sim_momentum.Pt());

      /* lookup for match */
      auto has_matched_vtx_track = match_inv.find(i_t) != match_inv.end();
      efficiencies[pdg]->efficiency_y_pt_filled->Fill(has_matched_vtx_track,
                                                      sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                                      sim_momentum.Pt());
    }
  } // sim tracks

  for (size_t i_t = 0; i_t < vtx_tracks_ptr->GetNumberOfChannels(); ++i_t) {
    auto channel = vtx_tracks_ptr->GetChannel(i_t);
    auto momentum = channel.GetMomentum3();

    auto has_matching_sim_track = match.find(i_t) != match.end();
    charged_hadrons_efficiency->eta_pt_vtx_tracks->Fill(has_matching_sim_track, momentum.Eta(), momentum.Pt());

    if (channel.GetField<int>(i_charge) < 0) {
      charged_hadrons_efficiency->eta_pt_vtx_tracks_neg->Fill(has_matching_sim_track, momentum.Eta(), momentum.Pt());
    }

  } // vtx tracks




  cout << "Matched " << matched_particles_->GetNumberOfChannels() << "/" << vtx_tracks_ptr->GetNumberOfChannels()
       << " vertex tracks" << endl;
}
void PidMatching::Finish() {
  cout << __func__ << endl;
  auto cwd = gDirectory;
  for (auto &&[pdg, efficiency] : efficiencies) {
    efficiency->output_dir->cd();
    efficiency->matched_tracks_y_pt->Write();
    efficiency->sim_tracks_y_pt->Write();

    /* clean-up bins */
    for (int i_cell = 0; i_cell < efficiency->matched_tracks_y_pt->GetNcells(); i_cell++) {
      if (efficiency->matched_tracks_y_pt->GetBinContent(i_cell) > efficiency->sim_tracks_y_pt->GetBinContent(i_cell)) {
        efficiency->matched_tracks_y_pt->SetBinContent(i_cell, 0);
        efficiency->sim_tracks_y_pt->SetBinContent(i_cell, 0);
      }
    }

    efficiency->efficiency_y_pt_from_hists = new TEfficiency(*efficiency->matched_tracks_y_pt,
                                                             *efficiency->sim_tracks_y_pt);
    efficiency->efficiency_y_pt_from_hists->SetTitle("N (VtxTracks) / N (SimTracks);#it{y}_{CM};p_{T} (GeV/c)");
    efficiency->efficiency_y_pt_from_hists->Write("efficiency_y_pt_from_hists");
    efficiency->efficiency_y_pt_filled->Write("efficiency_y_pt_filled");

    auto c = new TCanvas;
    c->Divide(2, 1);
    c->cd(1);
    efficiency->efficiency_y_pt_from_hists->SetMarkerSize(0.5);
    efficiency->efficiency_y_pt_from_hists->DrawClone("colz,text");
    c->cd(2);
    efficiency->efficiency_y_pt_filled->SetMarkerSize(0.5);
    efficiency->efficiency_y_pt_filled->DrawClone("colz,text");
    c->Write("efficiency_methods");
  }

  {
    charged_hadrons_efficiency->output_dir->cd();
    charged_hadrons_efficiency->eta_pt_vtx_tracks->Write();
    charged_hadrons_efficiency->eta_pt_vtx_tracks_neg->Write();
    auto c = new TCanvas;
    charged_hadrons_efficiency->eta_pt_vtx_tracks->DrawClone("colz,text");
    c->Write("c_efficiency_all_hadrons");
    charged_hadrons_efficiency->eta_pt_vtx_tracks_neg->DrawClone("colz,text");
    c->Write("c_efficiency_neg");
  }
  cwd->cd();
}
void PidMatching::PostFinish() {
  cout << __func__ << endl;

}

void PidMatching::InitEfficiencies() {
  auto cwd = gDirectory;

  for (int pdg : {211, -211, 2212}) {
    auto qa_struct = new PidEfficiencyQAStruct;
    efficiencies.emplace(pdg, qa_struct);

    qa_struct->output_dir = out_file_->mkdir(Form("efficiency_%d", pdg), "", true);
    qa_struct->output_dir->cd();

    qa_struct->matched_tracks_y_pt = new TH2D("matched_tracks_y_pt",
                                              "Matched tracks;#it{y}_{CM};p_{T} (GeV/c)",
                                              15, -2., 4.,
                                              15, 0., 3.);
    qa_struct->sim_tracks_y_pt = (TH2 *) qa_struct->matched_tracks_y_pt->Clone("sim_tracks_y_pt");

    qa_struct->efficiency_y_pt_filled = new TEfficiency("efficiency_filled",
                                                        "N (Matched SimTracks) / N(SimTracks);#it{y}_{CM};p_{T} (GeV/c)",
                                                        15, -2., 4.,
                                                        15, 0., 3.);

  }

  {
    charged_hadrons_efficiency = new ChargedHadronsEfficiencyStruct;
    charged_hadrons_efficiency->output_dir = out_file_->mkdir("efficiency_charged_hadrons", "", true);
    charged_hadrons_efficiency->eta_pt_vtx_tracks = new TEfficiency("efficiency_all_tracks",
                                                                    "N (VtxTracks) / N (Matched VtxTracks);#eta;p_{T} (GeV/c)",
                                                                    16, 0., 8.,
                                                                    15, 0., 3.);
    charged_hadrons_efficiency->eta_pt_vtx_tracks->SetMarkerSize(0.5);
    charged_hadrons_efficiency->eta_pt_vtx_tracks_neg = (TEfficiency*) charged_hadrons_efficiency->eta_pt_vtx_tracks->Clone("eff_neg");

  }


  // recover gDirectory
  cwd->cd();
}
