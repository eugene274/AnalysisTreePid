//
// Created by eugene on 17/02/2021.
//

#include "EvalEfficiency.hpp"

TASK_IMPL(EvalEfficiency)

boost::program_options::options_description EvalEfficiency::GetBoostOptions() {
  using namespace boost::program_options;

  options_description desc(Name());
  desc.add_options()
      ("efficiency-src", value(&efficiency_src_file_name_)->required(),
          "Input file for the efficiency")
      ("target-branch", value(&target_branch_name_)->required(), "Target branch")
      ("efficiency-field-name", value(&efficiency_field_name_)->required(),
          "To which field store the efficiency")
      ("new-branch", value(&new_branch_name_)->default_value(""),
          "Save to the branch with new name")
      ("var-centrality", value(&var_centrality_name_)->required(),
          "Name of centrality variable (outside target branch)")
      ("var-pid", value(&var_pid_name_)->default_value("pid"),
          "Name of the variable with pid (integer)")
      ("var-ycm", value(&var_y_cm_name_)->default_value("y_cm"),
          "Name of variable with y_{CM} (float)")
      ("var-pt", value(&var_pt_name_)->default_value("pT"),
          "Name of variable of with transverse momentum")
      ;
  return desc;
}
void EvalEfficiency::PreInit() {
  UserTask::PreInit();
}
void EvalEfficiency::PostFinish() {
  UserTask::PostFinish();
}
void EvalEfficiency::UserInit(std::map<std::string, void *> &map) {

  /* lookup relevant fields */
  centrality_ = GetVar(var_centrality_name_);
  pdg_ = GetVar(target_branch_name_ + "/" + var_pid_name_);
  ycm_ = GetVar(target_branch_name_ + "/" + var_y_cm_name_);
  pt = GetVar(target_branch_name_ + "/" + var_pt_name_);

  rec_particles = GetBranch(target_branch_name_);

  /* copy input config to the output with new field */
  auto target_branch = config_->GetBranchConfig(target_branch_name_);
  auto output_branch_name = new_branch_name_.empty()? target_branch_name_ : new_branch_name_;
  NewBranch(output_branch_name, target_branch.GetType());


  target_branch.Print();
}
void EvalEfficiency::UserExec() {


  std::cout << "Centrality = " << centrality_.Get<float>() << std::endl;

  for (auto track : rec_particles->Loop()) {
    auto pid = track.Get<float>(pdg_);
    std::cout << pid << std::endl;

  }



}
void EvalEfficiency::UserFinish() {

}
