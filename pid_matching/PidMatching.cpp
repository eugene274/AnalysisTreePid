//
// Created by eugene on 09/02/2021.
//

#include "PidMatching.hpp"
#include <boost/program_options.hpp>
#include <AnalysisTree/Matching.hpp>
#include <AnalysisTree/DataHeader.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <pid_new/core/PdgHelper.h>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TAxis.h>
#include <TH3.h>

#include "TEfficiencyHelper.hpp"
#include "PlotEfficiencies.hpp"

bool PidMatching::opts_loaded = false;
std::string PidMatching::qa_file_name = "efficiency.root";
bool PidMatching::save_canvases = false;

TASK_IMPL(PidMatching_NoCuts)
TASK_IMPL(PidMatching_StandardCuts)

using std::cout;
using std::endl;
using AnalysisTree::Matching;


struct PidMatching::PidEfficiencyQAStruct {
  TDirectory *output_dir{nullptr};

  TH2 *matched_tracks_y_pt{nullptr};
  TH2 *sim_tracks_y_pt{nullptr};

  TH3 *matched_tracks_centr_y_pt{nullptr};
  TH3 *sim_tracks_centr_y_pt{nullptr};

  TEfficiency *vtx_sim_y_pt{nullptr};
  TEfficiency *matched_sim_sim_y_pt{nullptr};

  TEfficiency *vtx_sim_centr_y_pt{nullptr};
  TEfficiency *matched_sim_sim_centr_y_pt{nullptr};

};

struct PidMatching::ChargedHadronsEfficiencyStruct {
  TDirectory *output_dir{nullptr};

  TEfficiency *eta_pt_vtx_tracks{nullptr};
  TEfficiency *eta_pt_vtx_tracks_neg{nullptr};
  TEfficiency *eta_pt_vtx_tracks_pos{nullptr};
};



boost::program_options::options_description PidMatching::GetBoostOptions() {
  namespace po = boost::program_options;

  if (!opts_loaded) {
    opts_loaded = true;
    po::options_description desc;
    desc.add_options()
        ("save-canvases", po::value(&save_canvases)->default_value(false), "Save canvases")
        ("qa-file-name", po::value(&qa_file_name)->default_value("efficiency_qa.root"));
    return desc;
  }
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
  centrality_ptr = static_cast<AnalysisTree::EventHeader *>(map["Centrality"]);

  auto simTracksConfig = config_->GetBranchConfig("SimTracks");
  simTracksConfig.Print();



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

  i_sim_mother_id = VarId("SimTracks/mother_id");
  i_centrality_id = VarId("Centrality/Centrality_Epsd");

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

  auto centrality = centrality_ptr->GetField<float>(i_centrality_id);

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
      if (CheckVtxTrack(vtx_track) && CheckSimTrack(sim_track)) {
        efficiencies[pdg]->matched_tracks_y_pt->Fill(vtx_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                                     vtx_momentum.Pt());
        efficiencies[pdg]->matched_tracks_centr_y_pt->Fill(
            centrality,
            vtx_momentum.Rapidity() - data_header_->GetBeamRapidity(),
            vtx_momentum.Pt());
      }
    }

  } // matched particles

  const auto match_inv = matching_ptr_->GetMatches(true);
  const auto match = matching_ptr_->GetMatches();

  for (size_t i_t = 0; i_t < sim_track_ptr->GetNumberOfChannels(); ++i_t) {
    auto channel = sim_track_ptr->GetChannel(i_t);

    if (!CheckSimTrack(channel)) continue;

    auto pdg = channel.GetField<int>(sim_track_pdg_id_);
    if (efficiencies.find(pdg) != efficiencies.end()) {
      sim_momentum = channel.Get4Momentum(pdg);
      efficiencies[pdg]->sim_tracks_y_pt->Fill(sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                               sim_momentum.Pt());
      efficiencies[pdg]->sim_tracks_centr_y_pt->Fill(
          centrality,
          sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
          sim_momentum.Pt());

      /* lookup for match */
      auto has_matched_vtx_track = match_inv.find(i_t) != match_inv.end()
          && CheckVtxTrack(vtx_tracks_ptr->GetChannel(match_inv.at(i_t)));
      efficiencies[pdg]->matched_sim_sim_y_pt->Fill(has_matched_vtx_track,
                                                    sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                                    sim_momentum.Pt());
      efficiencies[pdg]->matched_sim_sim_centr_y_pt->Fill(has_matched_vtx_track,
                                                          centrality,
                                                          sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                                          sim_momentum.Pt());
    }
  } // sim tracks

  for (size_t i_t = 0; i_t < vtx_tracks_ptr->GetNumberOfChannels(); ++i_t) {
    auto channel = vtx_tracks_ptr->GetChannel(i_t);
    auto momentum = channel.GetMomentum3();

    if (CheckVtxTrack(channel)) {
      auto has_matching_sim_track = match.find(i_t) != match.end();
      charged_hadrons_efficiency->eta_pt_vtx_tracks->Fill(has_matching_sim_track, momentum.Eta(), momentum.Pt());

      if (channel.GetField<int>(i_charge) < 0) {
        charged_hadrons_efficiency->eta_pt_vtx_tracks_neg->Fill(has_matching_sim_track, momentum.Eta(), momentum.Pt());
      } else if (channel.GetField<int>(i_charge) > 0) {
        charged_hadrons_efficiency->eta_pt_vtx_tracks_pos->Fill(has_matching_sim_track, momentum.Eta(), momentum.Pt());
      }
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

    efficiency->matched_tracks_centr_y_pt->Write();
    efficiency->sim_tracks_centr_y_pt->Write();

    efficiency->matched_sim_sim_centr_y_pt->Write();

    /* clean-up bins */
    for (int i_cell = 0; i_cell < efficiency->matched_tracks_y_pt->GetNcells(); i_cell++) {
      if (efficiency->matched_tracks_y_pt->GetBinContent(i_cell) > efficiency->sim_tracks_y_pt->GetBinContent(i_cell)) {
        efficiency->matched_tracks_y_pt->SetBinContent(i_cell, 0);
        efficiency->sim_tracks_y_pt->SetBinContent(i_cell, 0);
      }
    }

    /* clean-up 3d bins */
    for (int i_cell = 0; i_cell < efficiency->matched_tracks_centr_y_pt->GetNcells(); i_cell++) {
      if (efficiency->matched_tracks_centr_y_pt->GetBinContent(i_cell) > efficiency->sim_tracks_centr_y_pt->GetBinContent(i_cell)) {
        efficiency->matched_tracks_centr_y_pt->SetBinContent(i_cell, 0);
        efficiency->sim_tracks_centr_y_pt->SetBinContent(i_cell, 0);
      }
    }

    efficiency->vtx_sim_y_pt = new TEfficiency(*efficiency->matched_tracks_y_pt,
                                               *efficiency->sim_tracks_y_pt);
    efficiency->vtx_sim_y_pt->SetName("vtx_sim_y_pt");
    efficiency->vtx_sim_y_pt->SetTitle("N (VtxTracks) / N (SimTracks);#it{y}_{CM};p_{T} (GeV/c)");
    efficiency->vtx_sim_y_pt->Write("vtx_sim_y_pt");
    efficiency->matched_sim_sim_y_pt->Write();

    efficiency->vtx_sim_centr_y_pt = new TEfficiency(*efficiency->matched_tracks_centr_y_pt,
                                                     *efficiency->sim_tracks_centr_y_pt);
    efficiency->vtx_sim_centr_y_pt->SetName("vtx_sim_centr_y_pt");
    efficiency->vtx_sim_centr_y_pt->SetTitle("N (VtxTracks) / N (SimTracks);Centrality;#it{y}_{CM};p_{T} (GeV/c)");
    efficiency->vtx_sim_centr_y_pt->Write("vtx_sim_centr_y_pt");

    if (save_canvases) {
      MakeFancy(efficiency->output_dir);
    }

  }

  {
    charged_hadrons_efficiency->output_dir->cd();
    charged_hadrons_efficiency->eta_pt_vtx_tracks->Write();
    charged_hadrons_efficiency->eta_pt_vtx_tracks_neg->Write();
    charged_hadrons_efficiency->eta_pt_vtx_tracks_pos->Write();
    if (save_canvases) {
      MakeFancy(charged_hadrons_efficiency->output_dir);
    }
  }
  cwd->cd();
}
void PidMatching::PostFinish() {
  cout << __func__ << endl;

  if (qa_file_) {
    qa_file_->Close();
    delete qa_file_;
  }

}

void PidMatching::InitEfficiencies() {
  auto cwd = gDirectory;

  qa_file_ = TFile::Open(qa_file_name.c_str(), "RECREATE");

  for (int pdg : {211, -211, 2212}) {
    auto qa_struct = new PidEfficiencyQAStruct;
    efficiencies.emplace(pdg, qa_struct);

    qa_struct->output_dir = qa_file_->mkdir(Form("efficiency_%d", pdg), "", true);
    qa_struct->output_dir->cd();

    qa_struct->matched_tracks_y_pt = new TH2D("matched_tracks_y_pt",
                                              "Matched tracks;#it{y}_{CM};p_{T} (GeV/c)",
                                              30, -2., 4.,
                                              15, 0., 3.);
    qa_struct->sim_tracks_y_pt = (TH2 *) qa_struct->matched_tracks_y_pt->Clone("sim_tracks_y_pt");

    qa_struct->matched_tracks_centr_y_pt = new TH3D("matched_tracks_centr_y_pt",
                                                    "Centrality (%);#it{y}_{CM};p_{T} (GeV/c)",
                                                    20, 0., 100.,
                                                    30,-2., 4.,
                                                    15, 0., 3.);
    qa_struct->sim_tracks_centr_y_pt = (TH3*) qa_struct->matched_tracks_centr_y_pt->Clone("sim_tracks_centr_y_pt");

    qa_struct->matched_sim_sim_y_pt = new TEfficiency("matched_sim_sim_y_pt",
                                                      "N (Matched SimTracks) / N(SimTracks);#it{y}_{CM};p_{T} (GeV/c)",
                                                      30, -2., 4.,
                                                      15, 0., 3.);
    qa_struct->matched_sim_sim_centr_y_pt = new TEfficiency("matched_sim_sim_centr_y_pt",
                                                            "N (Matched SimTracks) / N(SimTracks);Centrality (%);#it{y}_{CM};p_{T} (GeV/c)",
                                                            20,0.,100.,
                                                            30,-2.,4.,
                                                            15,0.,3.);


  }

  {
    charged_hadrons_efficiency = new ChargedHadronsEfficiencyStruct;
    charged_hadrons_efficiency->output_dir = qa_file_->mkdir("efficiency_charged_hadrons", "", true);
    charged_hadrons_efficiency->eta_pt_vtx_tracks = new TEfficiency("efficiency_all_tracks",
                                                                    "N (VtxTracks) / N (Matched VtxTracks);#eta;p_{T} (GeV/c)",
                                                                    16,
                                                                    0.,
                                                                    8.,
                                                                    15,
                                                                    0.,
                                                                    3.);
    charged_hadrons_efficiency->eta_pt_vtx_tracks->SetMarkerSize(0.5);
    charged_hadrons_efficiency->eta_pt_vtx_tracks_neg =
        (TEfficiency *) charged_hadrons_efficiency->eta_pt_vtx_tracks->Clone("eff_neg");
    charged_hadrons_efficiency->eta_pt_vtx_tracks_pos =
        (TEfficiency *) charged_hadrons_efficiency->eta_pt_vtx_tracks->Clone("eff_pos");
  }


  // recover gDirectory
  cwd->cd();
}
