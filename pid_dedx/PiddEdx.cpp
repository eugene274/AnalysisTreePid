//
// Created by eugene on 09/09/2020.
//

#include "PiddEdx.h"
#include "TLorentzVector.h"

#include <AnalysisTree/DataHeader.hpp>
#include <pid_new/core/PdgHelper.h>

TASK_IMPL(PiddEdx)

boost::program_options::options_description PiddEdx::GetBoostOptions() {
  using namespace boost::program_options;

  options_description desc(GetName() + " options");
  desc.add_options()
      ("getter-file", value(&getter_file_)->required(), "Path to ROOT file with getter")
      ("getter-name", value(&getter_name_)->default_value("pid_getter"), "Name of Pid getter")
      ("tracks-branch", value(&tracks_branch_)->default_value("VtxTracks"), "Name of branch with tracks")
      ("dedx-field", value(&dedx_field_name_)->default_value("dedx_total"), "Name of the field with dEdx");
  return desc;
}

void PiddEdx::PreInit() {
  /* Load getter */
  TFile f(getter_file_.c_str(), "read");
  if (!f.IsOpen()) {
    throw std::runtime_error("Unable to open file with getter");
  }
  auto getter_ptr = f.Get(getter_name_.c_str());
  getter_.reset(dynamic_cast<Pid::BaseGetter *>(getter_ptr));
  if (!getter_) {
    throw std::runtime_error("Getter is nullptr");
  }

  SetInputBranchNames({tracks_branch_});
  SetOutputBranchName("RecParticles");
}

void PiddEdx::Init(std::map<std::string, void *> &Map) {
  tracks_ = static_cast<AnalysisTree::TrackDetector *>(Map.at(tracks_branch_));

  /* Input */
  auto tracks_config = config_->GetBranchConfig(tracks_->GetId());
  dedx_field_id_ = tracks_config.GetFieldId(dedx_field_name_);
  charge_field_id_ = tracks_config.GetFieldId("q");
  i_dca_x_field_id_ = tracks_config.GetFieldId("dcax");
  i_dca_y_field_id_ = tracks_config.GetFieldId("dcay");
  i_nhits_vtpc1_ = tracks_config.GetFieldId("nhits_vtpc1");
  i_nhits_vtpc2_ = tracks_config.GetFieldId("nhits_vtpc2");
  i_nhits_mtpc_ = tracks_config.GetFieldId("nhits_mtpc");
  i_nhits_pot_vtpc1_ = tracks_config.GetFieldId("nhits_pot_vtpc1");
  i_nhits_pot_vtpc2_ = tracks_config.GetFieldId("nhits_pot_vtpc2");
  i_nhits_pot_mtpc_ = tracks_config.GetFieldId("nhits_pot_mtpc");

  /* Output */
  rec_particle_config_ = AnalysisTree::BranchConfig(out_branch_, AnalysisTree::DetType::kParticle);
  rec_particle_config_.AddField<float>("y_cm");
  y_cm_field_id_ = rec_particle_config_.GetFieldId("y_cm");

  rec_particle_config_.AddField<float>("dcax");
  rec_particle_config_.AddField<float>("dcay");
  o_dca_x_field_id_ = rec_particle_config_.GetFieldId("dcax");
  o_dca_y_field_id_ = rec_particle_config_.GetFieldId("dcay");

  rec_particle_config_.AddField<int>("nhits_total");
  rec_particle_config_.AddField<int>("nhits_pot_total");
  rec_particle_config_.AddField<float>("nhits_ratio");
  o_nhits_total_ = rec_particle_config_.GetFieldId("nhits_total");
  o_nhits_pot_total_ = rec_particle_config_.GetFieldId("nhits_pot_total");
  o_nhits_ratio_ = rec_particle_config_.GetFieldId("nhits_ratio");

  out_config_->AddBranchConfig(rec_particle_config_);

  rec_particles_ = new AnalysisTree::Particles;
  out_tree_->Branch(out_branch_.c_str(), "AnalysisTree::Particles", &rec_particles_);
}

void PiddEdx::Exec() {

  auto &particle_config = rec_particle_config_;

  rec_particles_->ClearChannels();

  TLorentzVector momentum;

  for (int i_track = 0; i_track < tracks_->GetNumberOfChannels(); ++i_track) {
    auto track = tracks_->GetChannel(i_track);
    auto qp = track.GetP() * track.GetField<int>(charge_field_id_);
    auto dedx = track.GetField<float>(dedx_field_id_);
    auto pid = getter_->GetPid(qp, dedx, 0.9);

    if (pid != -1) {
      auto particle = rec_particles_->AddChannel();
      particle->Init(particle_config);
      particle->SetMomentum3(track.GetMomentum3());
      particle->SetPid(pid);

      /* mass */
      auto mass = PdgHelper::mass(pid);
      particle->SetMass(mass);

      /* y_cm */
      momentum.SetVectM(track.GetMomentum3(), mass);
      particle->SetField<float>(momentum.Rapidity() - data_header_->GetBeamRapidity(), y_cm_field_id_);

      /* dca_x, dca_y */
      particle->SetField<float>(track.GetField<float>(i_dca_x_field_id_), o_dca_x_field_id_);
      particle->SetField<float>(track.GetField<float>(i_dca_y_field_id_), o_dca_y_field_id_);

      /* nhits and ratio */
      {
        int nhits_total = track.GetField<int>(i_nhits_vtpc1_) +
            track.GetField<int>(i_nhits_vtpc2_) +
            track.GetField<int>(i_nhits_mtpc_);
        int nhits_pot_total = track.GetField<int>(i_nhits_pot_vtpc1_) +
            track.GetField<int>(i_nhits_pot_vtpc2_) +
            track.GetField<int>(i_nhits_pot_mtpc_);
        particle->SetField(nhits_total, o_nhits_total_);
        particle->SetField(nhits_pot_total, o_nhits_pot_total_);
        particle->SetField(nhits_total/float(nhits_pot_total), o_nhits_ratio_);
      }

    }

  }

  std::cout << "Identified " << rec_particles_->GetNumberOfChannels() << " particles of " <<
            tracks_->GetNumberOfChannels() << " tracks" << std::endl;
}
