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
      ("new-branch", value(&new_branch_name_)->default_value("RecParticles"),
          "Save to the branch with new name")
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

  /// INPUT
  vtx_tracks_branch = GetInBranch("VtxTracks");
  rec_particles_branch = GetInBranch(target_branch_name_);

  /// OUTPUT
  processed_branch = NewBranch(new_branch_name_, rec_particles_branch->GetConfig());
  pid_v = processed_branch->GetFieldVar("pid");
  y_cm_v = processed_branch->GetFieldVar(var_y_cm_name_);
  pt_v = processed_branch->GetFieldVar(var_pt_name_);
}
void EvalEfficiency::UserExec() {


  int multiplicity = 0;
  for (auto &vtx_track : vtx_tracks_branch->Loop()) {
    // if good
    ++multiplicity;
  }

  processed_branch->ClearChannels();
  for (auto &rec_particle : rec_particles_branch->Loop()) {
    auto processed_particle = processed_branch->NewChannel();
    processed_particle.CopyContents(rec_particle);

    auto pid = processed_particle[pid_v].GetInt();
    auto y_cm = processed_particle[y_cm_v].GetVal();
    auto pt = processed_particle[pt_v].GetVal();
  }



}
void EvalEfficiency::UserFinish() {

}
