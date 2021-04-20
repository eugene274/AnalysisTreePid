//
// Created by eugene on 09/09/2020.
//

#include "PiddEdx.h"

#include <TKey.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>

#include <AnalysisTree/DataHeader.hpp>

#include <pid_new/core/PdgHelper.h>

#include <regex>
#include <boost/lexical_cast.hpp>

TASK_IMPL(PiddEdx)

struct PiddEdx::Efficiency {

  std::unique_ptr<const TEfficiency> eff{nullptr};

  float Eval(float centrality, float y_cm, float pt) const {
    return eff->GetEfficiency(eff->FindFixBin(centrality, y_cm, pt));
  }
};

boost::program_options::options_description PiddEdx::GetBoostOptions() {
  using namespace boost::program_options;

  options_description desc(GetName() + " options");
  desc.add_options()
      ("getter-file", value(&getter_file_)->required(), "Path to ROOT file with getter")
      ("getter-name", value(&getter_name_)->default_value("pid_getter"), "Name of Pid getter")
      ("tracks-branch", value(&tracks_branch_)->default_value("VtxTracks"), "Name of branch with tracks")
      ("dedx-field", value(&dedx_field_name_)->default_value("dedx_total"), "Name of the field with dEdx")
      ("output-branch", value(&output_branch_name_)->default_value("RecParticles"),
       "Name of the output branch with identified particles")
      ("efficiency-definitions", value(&efficiency_definitions_)->multitoken(), "Efficiency definitions");;
  return desc;
}

void PiddEdx::ProcessBoostVM(const boost::program_options::variables_map &vm) {
  UserTask::ProcessBoostVM(vm);
}

void PiddEdx::PreInit() {
  InitEfficiencyDefinitions();
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
  SetOutputBranchName(output_branch_name_);
}

void PiddEdx::UserInit(std::map<std::string, void *> &Map) {
  tracks_ = static_cast<AnalysisTree::TrackDetector *>(Map.at(tracks_branch_));

  /* Input */
  auto &vtx_tracks_config = config_->GetBranchConfig(tracks_branch_);
  dedx_field_id_ = vtx_tracks_config.GetFieldId(dedx_field_name_);
  charge_field_id_ = vtx_tracks_config.GetFieldId("q");
  i_dca_x_field_id_ = vtx_tracks_config.GetFieldId("dcax");
  i_dca_y_field_id_ = vtx_tracks_config.GetFieldId("dcay");
  i_nhits_vtpc1_ = vtx_tracks_config.GetFieldId("nhits_vtpc1");
  i_nhits_vtpc2_ = vtx_tracks_config.GetFieldId("nhits_vtpc2");
  i_nhits_mtpc_ = vtx_tracks_config.GetFieldId("nhits_mtpc");
  i_nhits_pot_vtpc1_ = vtx_tracks_config.GetFieldId("nhits_pot_vtpc1");
  i_nhits_pot_vtpc2_ = vtx_tracks_config.GetFieldId("nhits_pot_vtpc2");
  i_nhits_pot_mtpc_ = vtx_tracks_config.GetFieldId("nhits_pot_mtpc");
  i_chi2 = vtx_tracks_config.GetFieldId("chi2");
  i_ndf = vtx_tracks_config.GetFieldId("ndf");

  /* Output */
  rec_particle_config_ = AnalysisTree::BranchConfig(out_branch_, AnalysisTree::DetType::kParticle);
  rec_particle_config_.AddField<float>("y_cm");
  rec_particle_config_.AddField<float>("y");
  y_cm_field_id_ = rec_particle_config_.GetFieldId("y_cm");
  y_field_id_ = rec_particle_config_.GetFieldId("y");

  rec_particle_config_.AddField<float>("dcax");
  rec_particle_config_.AddField<float>("dcay");
  o_dca_x_field_id_ = rec_particle_config_.GetFieldId("dcax");
  o_dca_y_field_id_ = rec_particle_config_.GetFieldId("dcay");

  rec_particle_config_.AddField<int>("nhits_total");
  rec_particle_config_.AddField<int>("nhits_vtpc");
  rec_particle_config_.AddField<int>("nhits_pot_total");
  rec_particle_config_.AddField<float>("nhits_ratio");
  o_nhits_total_ = rec_particle_config_.GetFieldId("nhits_total");
  o_nhits_vtpc_ = rec_particle_config_.GetFieldId("nhits_vtpc");
  o_nhits_pot_total_ = rec_particle_config_.GetFieldId("nhits_pot_total");
  o_nhits_ratio_ = rec_particle_config_.GetFieldId("nhits_ratio");

  rec_particle_config_.AddField<float>("chi2_ndf");
  o_chi2_ndf = rec_particle_config_.GetFieldId("chi2_ndf");

  out_config_->AddBranchConfig(rec_particle_config_);

  rec_particles_ = new AnalysisTree::Particles;
  out_tree_->Branch(out_branch_.c_str(), &rec_particles_);
}

void PiddEdx::UserExec() {

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
      particle->SetField<float>(momentum.Rapidity(), y_field_id_);
      particle->SetField<float>(momentum.Rapidity() - data_header_->GetBeamRapidity(), y_cm_field_id_);

      /* dca_x, dca_y */
      particle->SetField<float>(track.GetField<float>(i_dca_x_field_id_), o_dca_x_field_id_);
      particle->SetField<float>(track.GetField<float>(i_dca_y_field_id_), o_dca_y_field_id_);
      particle->SetField<float>(track.GetField<float>(i_chi2)/track.GetField<int>(i_ndf), o_chi2_ndf);
      /* nhits and ratio */
      {
        int nhits_total = track.GetField<int>(i_nhits_vtpc1_) +
            track.GetField<int>(i_nhits_vtpc2_) +
            track.GetField<int>(i_nhits_mtpc_);
        int nhits_vtpc =
            track.GetField<int>(i_nhits_vtpc1_) +
                track.GetField<int>(i_nhits_vtpc2_);

        int nhits_pot_total = track.GetField<int>(i_nhits_pot_vtpc1_) +
            track.GetField<int>(i_nhits_pot_vtpc2_) +
            track.GetField<int>(i_nhits_pot_mtpc_);
        particle->SetField(nhits_total, o_nhits_total_);
        particle->SetField<int>(nhits_vtpc, o_nhits_vtpc_);
        particle->SetField(nhits_pot_total, o_nhits_pot_total_);
        particle->SetField(float(nhits_total) / float(nhits_pot_total), o_nhits_ratio_);
      }

    }

  }

  std::cout << "Identified " << rec_particles_->GetNumberOfChannels() << " particles of " <<
            tracks_->GetNumberOfChannels() << " tracks" << std::endl;
}

void PiddEdx::InitEfficiencyDefinitions() {
  const std::regex tgt_re_expr("^.*tgt:(\\w+).*$");
  const std::regex src_re_expr("^.*src:([^\\s]+).*$");
  const std::regex eff_dir_re_expr("^efficiency_([-\\d]+)$");

  for (auto &eff_def : efficiency_definitions_) {
    std::smatch tgt_match;
    bool tgt_found = std::regex_search(eff_def, tgt_match, tgt_re_expr);
    if (!tgt_found)
      throw std::runtime_error("No 'tgt' entry in the efficiency definition");
    std::string tgt = tgt_match.str(1);

    std::smatch src_match;
    bool src_found = std::regex_search(eff_def, src_match, src_re_expr);
    if (!src_found)
      throw std::runtime_error("No 'src' entry in the efficiency definition");
    std::string src_filename = src_match.str(1);

    /* attempting to reach efficiency src */
    TFile f_efficiency(src_filename.c_str(), "READ");
    if (f_efficiency.IsZombie())
      throw std::runtime_error("Efficiency ROOT file seems to be empty");

    /* Lookup available PID-s */
    for (auto o : *f_efficiency.GetListOfKeys()) {
      std::smatch eff_pid_match;
      std::string obj_name{o->GetName()};
      if (std::regex_search(obj_name, eff_pid_match, eff_dir_re_expr)) {
        auto pid_str = eff_pid_match.str(1);
        auto pid = boost::lexical_cast<int>(pid_str);

        auto key = (TKey *) o;
        if (!TClass::GetClass(key->GetClassName())->InheritsFrom(TDirectory::Class())) {
          throw std::runtime_error("Expected directory for the efficiency");
        }
        auto eff_dir = (TDirectory *) key->ReadObj();

        auto efficiency_matrix = eff_dir->Get<TEfficiency>(efficiency_matrix_name_.c_str());
        if (!efficiency_matrix)
          throw std::runtime_error("Efficiency matrix is not found for " + pid_str);

        /* object is not owned by parent directory */
        efficiency_matrix->SetDirectory(nullptr);

        /* Success, populating structures */
        auto efficiency_struct = std::make_unique<Efficiency>();
        efficiency_struct->eff = std::unique_ptr<TEfficiency>(efficiency_matrix);
        efficiencies_.emplace(pid, std::move(efficiency_struct));
      }
    }

  }

}

