//
// Created by eugene on 09/09/2020.
//

#ifndef ATPIDTASK_PID_DEDX_PIDDEDX_H
#define ATPIDTASK_PID_DEDX_PIDDEDX_H

#include <TFile.h>
#include <TTree.h>

#include <at_task/Task.h>
#include <pid/Getter.h>
#include <AnalysisTree/Detector.hpp>

class PiddEdx : public UserTask {

public:
  void Init(std::map<std::string, void *> &Map) override {
    tracks_ = static_cast<AnalysisTree::TrackDetector *>(Map.at(tracks_branch_));

    auto tracks_config = config_->GetBranchConfig(tracks_->GetId());
    dedx_field_id_ = tracks_config.GetFieldId(dedx_field_name_);

    auto rec_particle_config = AnalysisTree::BranchConfig(out_branch_, AnalysisTree::DetType::kParticle);
    out_config_->AddBranchConfig(rec_particle_config);

    rec_particles_ = new AnalysisTree::Particles;
    out_tree_->Branch(out_branch_.c_str(), "AnalysisTree::Particles", &rec_particles_);
  }
  void Exec() override {

    auto& particle_config = out_config_->GetBranchConfig(out_branch_);

    rec_particles_->ClearChannels();

    for (int i_track = 0; i_track < tracks_->GetNumberOfChannels(); ++i_track) {
      auto track = tracks_->GetChannel(i_track);
      auto abs_momentum = track.GetP();
      auto dedx = track.GetField<float>(dedx_field_id_);
      auto pid = getter_->GetPid(abs_momentum, dedx, 0.9);

      if (pid != -1) {
        auto particle = rec_particles_->AddChannel();
        particle->Init(particle_config);
        particle->SetMomentum3(track.GetMomentum3());
        particle->SetPid(pid);
      }
    }

    std::cout << "Identified " << rec_particles_->GetNumberOfChannels() << " particles of " <<
      tracks_->GetNumberOfChannels() << " tracks" << std::endl;
  }
  void Finish() override {

  }
  boost::program_options::options_description GetBoostOptions() override {
    using namespace boost::program_options;

    options_description desc(GetName() + " options");
    desc.add_options()
        ("getter-file", value(&getter_file_)->required(), "Path to ROOT file with getter")
        ("getter-name", value(&getter_name_)->default_value("pid_getter"), "Name of Pid getter")
        ("tracks-branch", value(&tracks_branch_)->default_value("VtxTracks"), "Name of branch with tracks")
        ("dedx-field", value(&dedx_field_name_)->default_value("dedx_total"), "Name of the field with dEdx")
        ;
    return desc;
  }
  void PreInit() override {
    /* Load getter */
    TFile f(getter_file_.c_str(),"read");
    if (!f.IsOpen()) {
      throw std::runtime_error("Unable to open file with getter");
    }
    auto getter_ptr = f.Get(getter_name_.c_str());
    getter_.reset(dynamic_cast<Pid::BaseGetter*>(getter_ptr));
    if (!getter_) {
      throw std::runtime_error("Getter is nullptr");
    }

    SetInputBranchNames({tracks_branch_});
    SetOutputBranchName("RecParticles");
  }
  void PostFinish() override {
    UserTask::PostFinish();
  }

private:
  std::string getter_file_;
  std::string getter_name_;

  std::string tracks_branch_;
  std::string dedx_field_name_;

  std::shared_ptr<Pid::BaseGetter> getter_;

  AnalysisTree::TrackDetector *tracks_{nullptr};
  int dedx_field_id_{-1};

  AnalysisTree::Particles *rec_particles_{nullptr};


  TASK_DEF(PiddEdx, 0)
};

#endif //ATPIDTASK_PID_DEDX_PIDDEDX_H
