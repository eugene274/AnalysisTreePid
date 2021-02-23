//
// Created by eugene on 17/02/2021.
//

#include "EvalEfficiency.hpp"
#include <TEfficiency.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <regex>
#include <boost/lexical_cast.hpp>

TASK_IMPL(EvalEfficiency)

struct EvalEfficiency::Efficiency {
  TEfficiency *eff_y_pt{nullptr};

  ~Efficiency() {
    delete eff_y_pt;
  }

};

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
      ("weight-name", value(&efficiency_field_name_)->default_value("weight_efficiency"),
          "Name of the variable with efficiency weight")
      ;
  return desc;
}
void EvalEfficiency::PreInit() {
  LoadEfficiencies();
}
void EvalEfficiency::PostFinish() {
  UserTask::PostFinish();
}
void EvalEfficiency::UserInit(std::map<std::string, void *> &map) {

  /// INPUT
  rec_particles_branch = GetInBranch(target_branch_name_);

  /// OUTPUT
  processed_branch = NewBranch(new_branch_name_, rec_particles_branch->GetConfig());
  pid_v = processed_branch->GetFieldVar("pid");
  y_cm_v = processed_branch->GetFieldVar(var_y_cm_name_);
  pt_v = processed_branch->GetFieldVar(var_pt_name_);
  weight_v = processed_branch->NewVariable(efficiency_field_name_, FLOAT);
  processed_branch->Freeze();
}
void EvalEfficiency::UserExec() {


  processed_branch->ClearChannels();
  for (auto &rec_particle : rec_particles_branch->Loop()) {
    auto processed_particle = processed_branch->NewChannel();
    processed_particle.CopyContents(rec_particle);

    auto pid = processed_particle[pid_v].GetInt();
    auto y_cm = processed_particle[y_cm_v].GetVal();
    auto pt = processed_particle[pt_v].GetVal();

    if (efficiencies_.find(pid) != efficiencies_.end()) {
      auto eff_y_pt = efficiencies_[pid]->eff_y_pt;
      auto eff_y_pt_val = eff_y_pt->GetEfficiency(eff_y_pt->FindFixBin(y_cm, pt));
      auto err_lo = eff_y_pt->GetEfficiencyErrorLow(eff_y_pt->FindFixBin(y_cm, pt));
      auto err_hi = eff_y_pt->GetEfficiencyErrorUp(eff_y_pt->FindFixBin(y_cm, pt));
      auto efficiency_eps = (err_hi + err_lo)/eff_y_pt_val;
      auto weight_val = efficiency_eps < efficiency_eps_threshold ? float(1./eff_y_pt_val) : 0.0f;
      processed_particle[weight_v] = weight_val;
    }
  }



}
void EvalEfficiency::UserFinish() {

}
void EvalEfficiency::LoadEfficiencies() {
  TFile eff(efficiency_src_file_name_.c_str(), "READ");
  if (!eff.IsOpen())
    throw std::runtime_error("Unable to open efficiency file");

  const std::regex re_efficiency_dir("^efficiency_(-?\\d+)$");
  std::smatch match_results;
  for (auto key_object : *eff.GetListOfKeys()) {
    std::string obj_name(key_object->GetName());
    if (!std::regex_search(obj_name, match_results, re_efficiency_dir))
      continue;
    std::cout << "Find efficiency dir: " << obj_name << std::endl;

    auto key = (TKey *) key_object;
    auto efficiency_dir = (TDirectory *) key->ReadObj();

    auto efficiency = std::make_shared<Efficiency>();
    efficiency->eff_y_pt = (TEfficiency *) efficiency_dir->Get("matched_sim_sim_y_pt");
    efficiency->eff_y_pt->SetDirectory(nullptr);
    assert(efficiency->eff_y_pt);

    auto pid = boost::lexical_cast<int>(match_results.str(1));

    efficiencies_.emplace(pid, efficiency);
  }

}
