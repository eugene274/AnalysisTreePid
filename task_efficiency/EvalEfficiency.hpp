//
// Created by eugene on 17/02/2021.
//

#ifndef ATPIDTASK_TASK_EFFICIENCY_EVALEFFICIENCY_HPP_
#define ATPIDTASK_TASK_EFFICIENCY_EVALEFFICIENCY_HPP_

#include <at_task/Task.h>

class EvalEfficiency : public UserFillTask {

 public:
  boost::program_options::options_description GetBoostOptions() override;
  void PreInit() override;
  void PostFinish() override;

  void UserInit(std::map<std::string, void *> &map) override;
  void UserExec() override;
  void UserFinish() override;

 private:
  void LoadEfficiencyFromSource();

  std::string target_branch_name_;
  std::string efficiency_src_file_name_;
  std::string efficiency_field_name_;
  std::string new_branch_name_;

  std::string var_centrality_name_;
  std::string var_pid_name_;
  std::string var_y_cm_name_;
  std::string var_pt_name_;

  ATI2::Branch *rec_particles_branch;

  ATI2::Branch *processed_branch;
  ATI2::Branch *vtx_tracks_branch;

  ATI2::Variable pid_v;
  ATI2::Variable y_cm_v;
  ATI2::Variable pt_v;
 TASK_DEF(EvalEfficiency, 0)
};

#endif //ATPIDTASK_TASK_EFFICIENCY_EVALEFFICIENCY_HPP_
