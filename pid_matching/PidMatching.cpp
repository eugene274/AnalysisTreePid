//
// Created by eugene on 09/02/2021.
//

#include "PidMatching.hpp"
#include <boost/program_options.hpp>
#include <AnalysisTree/Matching.hpp>
#include <AnalysisTree/DataHeader.hpp>
#include <pid_new/core/PdgHelper.h>
#include <TTree.h>
#include <TLorentzVector.h>

TASK_IMPL(PidMatching)

using std::cout;
using std::endl;
using AnalysisTree::Matching;

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
void PidMatching::PostFinish() {
  cout << __func__ << endl;
}
void PidMatching::Init(std::map<std::string, void *> &map) {
  matching_ptr_ = static_cast<Matching *>(map["VtxTracks2SimTracks"]);
  vtx_tracks_ptr = static_cast<AnalysisTree::TrackDetector *>(map["VtxTracks"]);
  sim_track_ptr = static_cast<AnalysisTree::TrackDetector *>(map["SimTracks"]);

  auto vtx_tracks_config = config_->GetBranchConfig("VtxTracks");

  /* input config */
  sim_track_pdg_id_ = config_->GetBranchConfig("SimTracks").GetFieldId("pdg");

  i_dca_x_field_id_ = config_->GetBranchConfig("VtxTracks").GetFieldId("dcax");
  i_dca_y_field_id_ = config_->GetBranchConfig("VtxTracks").GetFieldId("dcay");

  i_nhits_vtpc1_ = vtx_tracks_config.GetFieldId("nhits_vtpc1");
  i_nhits_vtpc2_ = vtx_tracks_config.GetFieldId("nhits_vtpc2");
  i_nhits_mtpc_ = vtx_tracks_config.GetFieldId("nhits_mtpc");
  i_nhits_pot_vtpc1_ = vtx_tracks_config.GetFieldId("nhits_pot_vtpc1");
  i_nhits_pot_vtpc2_ = vtx_tracks_config.GetFieldId("nhits_pot_vtpc2");
  i_nhits_pot_mtpc_ = vtx_tracks_config.GetFieldId("nhits_pot_mtpc");

  /* output config */
  matched_particles_config_.AddField<float>("y_cm");
  matched_particles_config_.AddField<float>("y");
  y_cm_field_id_ = matched_particles_config_.GetFieldId("y_cm");
  y_field_id_ = matched_particles_config_.GetFieldId("y");

  matched_particles_config_.AddField<float>("dcax");
  matched_particles_config_.AddField<float>("dcay");
  o_dca_x_field_id_ = matched_particles_config_.GetFieldId("dcax");
  o_dca_y_field_id_ = matched_particles_config_.GetFieldId("dcay");

  matched_particles_config_.AddField<int>("nhits_total");
  matched_particles_config_.AddField<int>("nhits_vtpc");
  matched_particles_config_.AddField<int>("nhits_pot_total");
  matched_particles_config_.AddField<float>("nhits_ratio");
  o_nhits_total_ = matched_particles_config_.GetFieldId("nhits_total");
  o_nhits_vtpc_ = matched_particles_config_.GetFieldId("nhits_vtpc");
  o_nhits_pot_total_ = matched_particles_config_.GetFieldId("nhits_pot_total");
  o_nhits_ratio_ = matched_particles_config_.GetFieldId("nhits_ratio");

  /* parameters of the matched sim track */
  matched_particles_config_.AddField<float>("sim_y_cm");
  o_sim_y_cm_ = matched_particles_config_.GetFieldId("sim_y_cm");
  matched_particles_config_.AddField<float>("sim_pt");
  o_sim_pt_ = matched_particles_config_.GetFieldId("sim_pt");
  matched_particles_config_.AddField<float>("sim_phi");
  o_sim_phi_ = matched_particles_config_.GetFieldId("sim_phi");

  out_config_->AddBranchConfig(matched_particles_config_);

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
    particle->SetField<float>(o_sim_y_cm_, sim_momentum.Rapidity() - data_header_->GetBeamRapidity());
    particle->SetField<float>(o_sim_pt_, sim_momentum.Pt());
    particle->SetField<float>(o_sim_phi_, sim_momentum.Phi());


  } // particles




  cout << "Matched " << matched_particles_->GetNumberOfChannels() << "/" << vtx_tracks_ptr->GetNumberOfChannels()
       << " vertex tracks" << endl;
}
void PidMatching::Finish() {
  cout << __func__ << endl;
}
