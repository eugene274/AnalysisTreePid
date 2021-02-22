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
void PidMatching::UserInit(std::map<std::string, void *> &map) {
  using AnalysisTree::Types;

  InitEfficiencies();

  matching_ptr_ = static_cast<Matching *>(map["VtxTracks2SimTracks"]);
  vtxt_branch = GetInBranch("VtxTracks");
  simt_branch = GetInBranch("SimTracks");
  vtxt_branch->GetConfig().Print();
  simt_branch->GetConfig().Print();

//  centrality_ptr = static_cast<AnalysisTree::EventHeader *>(map["Centrality"]);


  /* input config */

  /// EVT variable
  evt_centrality = GetVar("Centrality/Centrality_Epsd");

  /// SIM Tracks
  sim_pdg_ = GetVar("SimTracks/pdg");
  sim_mother_id_ = GetVar("SimTracks/mother_id");

  /// VTX Tracks
  vtxt_dca_x_ = GetVar("VtxTracks/dcax");
  vtxt_dca_y_ = GetVar("VtxTracks/dcay");

  vtxt_nhits_vtpc1_ = GetVar("VtxTracks/nhits_vtpc1");
  vtxt_nhits_vtpc2_ = GetVar("VtxTracks/nhits_vtpc2");
  vtxt_nhits_mtpc_ = GetVar("VtxTracks/nhits_mtpc");
  vtxt_nhits_pot_vtpc1_ = GetVar("VtxTracks/nhits_pot_vtpc1");
  vtxt_nhits_pot_vtpc2_ = GetVar("VtxTracks/nhits_pot_vtpc2");
  vtxt_nhits_pot_mtpc_ = GetVar("VtxTracks/nhits_pot_mtpc");
  vtxt_charge = GetVar("VtxTracks/q");

  /// MATCHED TRACKS
  mt_branch = NewBranch("MatchedVtxTracks", AnalysisTree::DetType::kParticle);
  mt_branch->NewVariable("dcax", Types::kFloat);
  mt_branch->NewVariable("dcay", Types::kFloat);
  mt_branch->NewVariable("q", Types::kInteger);
  mt_pid = mt_branch->GetFieldVar("pid"); // internal variable
  mt_mass = mt_branch->GetFieldVar("mass"); // internal variable
  mt_y_cm_ = mt_branch->NewVariable("y_cm", Types::kFloat);
  mt_nhits_total_ = mt_branch->NewVariable("nhits_total", Types::kInteger);
  mt_nhits_vtpc_ = mt_branch->NewVariable("nhits_vtpc", Types::kInteger);
  mt_nhits_pot_total_ = mt_branch->NewVariable("nhits_pot_total", Types::kInteger);
  mt_nhits_ratio_ = mt_branch->NewVariable("nhits_ratio", Types::kFloat);

  /* parameters of the matched sim track */
  mt_sim_y_cm_ = mt_branch->NewVariable("sim_y_cm", Types::kFloat);
  mt_sim_pt_ = mt_branch->NewVariable("sim_pt", Types::kFloat);
  mt_sim_phi_ = mt_branch->NewVariable("sim_phi", Types::kFloat);


  mt_branch->Freeze();
}
void PidMatching::UserExec() {

  using AnalysisTree::Types;
  using AnalysisTree::Particle;
  using AnalysisTree::Track;

  TLorentzVector sim_momentum;
  TLorentzVector vtx_momentum;

  mt_branch->ClearChannels();
  for (auto &&[vtxId, simId] : matching_ptr_->GetMatches()) {
    const auto vtx_track = (*vtxt_branch)[vtxId];
    const auto sim_track = (*simt_branch)[simId];

    auto matched_track = mt_branch->NewChannel();
    matched_track.CopyContents(vtx_track);

    auto pdg = sim_track[sim_pdg_].GetInt();
    sim_momentum.SetVectM(sim_track.DataT<Particle>()->GetMomentum3(), PdgHelper::mass(pdg));
    vtx_momentum.SetVectM(vtx_track.DataT<Track>()->GetMomentum3(), sim_momentum.M());

    matched_track[mt_pid] = pdg;
    matched_track[mt_mass] = float(sim_momentum.M());
    matched_track[mt_y_cm_] = float(vtx_momentum.Rapidity() - data_header_->GetBeamRapidity());
    matched_track[mt_sim_y_cm_] = float(sim_momentum.Rapidity() - data_header_->GetBeamRapidity());
    matched_track[mt_sim_pt_] = float(sim_momentum.Pt());
    matched_track[mt_sim_phi_] = float(sim_momentum.Phi());

    if (efficiencies.find(pdg) != efficiencies.end()) {
      if (CheckVtxTrack(vtx_track) && CheckSimTrack(sim_track)) {
        efficiencies[pdg]->matched_tracks_y_pt->Fill(vtx_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                                     vtx_momentum.Pt());
        efficiencies[pdg]->matched_tracks_centr_y_pt->Fill(
            *evt_centrality,
            vtx_momentum.Rapidity() - data_header_->GetBeamRapidity(),
            vtx_momentum.Pt());
      }
    }

  } // matched particles

  const auto match_inv = matching_ptr_->GetMatches(true);
  const auto match = matching_ptr_->GetMatches();

  for (const auto &sim_track : simt_branch->Loop()) {

    if (!CheckSimTrack(sim_track)) continue;

    auto pdg = sim_track[sim_pdg_].GetInt();
    if (efficiencies.find(pdg) != efficiencies.end()) {
      sim_momentum.SetVectM(sim_track.DataT<Particle>()->GetMomentum3(), PdgHelper::mass(pdg));
      efficiencies[pdg]->sim_tracks_y_pt->Fill(sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                               sim_momentum.Pt());
      efficiencies[pdg]->sim_tracks_centr_y_pt->Fill(
          *evt_centrality,
          sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
          sim_momentum.Pt());

      auto has_matched_vtx_track = match_inv.find(sim_track.GetNChannel()) != match_inv.end()
          && CheckVtxTrack((*vtxt_branch)[match_inv.at(sim_track.GetNChannel())]);
      efficiencies[pdg]->matched_sim_sim_y_pt->Fill(has_matched_vtx_track,
                                                    sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                                    sim_momentum.Pt());
      efficiencies[pdg]->matched_sim_sim_centr_y_pt->Fill(has_matched_vtx_track,
                                                          *evt_centrality,
                                                          sim_momentum.Rapidity() - data_header_->GetBeamRapidity(),
                                                          sim_momentum.Pt());
    }
  } // sim tracks

  for (const auto &vtx_track : vtxt_branch->Loop()) {

    auto momentum = vtx_track.DataT<Track>()->GetMomentum3();
    if (CheckVtxTrack(vtx_track)) {
      auto has_matching_sim_track = match.find(vtx_track.GetNChannel()) != match.end();
      charged_hadrons_efficiency->eta_pt_vtx_tracks->Fill(has_matching_sim_track, momentum.Eta(), momentum.Pt());

      if (vtx_track[vtxt_charge] < 0) {
        charged_hadrons_efficiency->eta_pt_vtx_tracks_neg->Fill(has_matching_sim_track, momentum.Eta(), momentum.Pt());
      } else if (vtx_track[vtxt_charge] > 0) {
        charged_hadrons_efficiency->eta_pt_vtx_tracks_pos->Fill(has_matching_sim_track, momentum.Eta(), momentum.Pt());
      }
    }

  } // vtx tracks


  cout << "Matched " << mt_branch->size() << "/" << vtxt_branch->size()
       << " vertex tracks" << endl;

}
void PidMatching::UserFinish() {
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
      ProcessEfficiencyDir(efficiency->output_dir);
    }

  }

  {
    charged_hadrons_efficiency->output_dir->cd();
    charged_hadrons_efficiency->eta_pt_vtx_tracks->Write();
    charged_hadrons_efficiency->eta_pt_vtx_tracks_neg->Write();
    charged_hadrons_efficiency->eta_pt_vtx_tracks_pos->Write();
    if (save_canvases) {
      ProcessEfficiencyDir(charged_hadrons_efficiency->output_dir);
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
