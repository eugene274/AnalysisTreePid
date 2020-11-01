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

class PiddEdx : public UserFillTask {

public:
  void Init(std::map<std::string, void *> &Map) override;
  void Exec() override;
  void Finish() override {

  }
  boost::program_options::options_description GetBoostOptions() override;
  void PreInit() override;
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
  short dedx_field_id_{-1};
  short charge_field_id_{-1};

  AnalysisTree::Particles *rec_particles_{nullptr};

  AnalysisTree::BranchConfig rec_particle_config_;
  short y_cm_field_id_;
  short o_dca_x_field_id_;
  short o_dca_y_field_id_;
  short i_dca_y_field_id_;
  short i_dca_x_field_id_;
  short o_nhits_total_;
  short o_nhits_pot_total_;
  short o_nhits_ratio_;
  short i_nhits_mtpc_;
  short i_nhits_vtpc2_;
  short i_nhits_vtpc1_;
  short i_nhits_pot_mtpc_;
  short i_nhits_pot_vtpc2_;
  short i_nhits_pot_vtpc1_;

TASK_DEF(PiddEdx, 0)
};

#endif //ATPIDTASK_PID_DEDX_PIDDEDX_H
