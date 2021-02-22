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



/**
 * @brief Takes dEdx PID from PID Getter
 */
class PiddEdx : public UserFillTask {

public:
 protected:
  bool UseATI2() const override { return false; };
 public:
  void UserInit(std::map<std::string, void *> &Map) override;
  void UserExec() override;
  void UserFinish() override {}
  boost::program_options::options_description GetBoostOptions() override;
  void ProcessBoostVM(const boost::program_options::variables_map &vm) override;
  void PreInit() override;
  void PostFinish() override {
  }

private:
  void InitEfficiencyDefinitions();

  /* SETUP */
  std::string getter_file_;
  std::string getter_name_;

  std::string tracks_branch_;
  std::string dedx_field_name_;

  std::string output_branch_name_;

  /* efficiency */
  std::vector<std::string> efficiency_definitions_;
  std::string efficiency_matrix_name_{"vtx_sim_centr_y_pt"};
  std::string efficiency_centrality_vname_{"Centrality/Centrality_Epsd"};

  struct Efficiency;
  std::map<int, std::unique_ptr<Efficiency>> efficiencies_;




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
  short o_nhits_vtpc_;
  short i_nhits_mtpc_;
  short i_nhits_vtpc2_;
  short i_nhits_vtpc1_;
  short i_nhits_pot_mtpc_;
  short i_nhits_pot_vtpc2_;
  short i_nhits_pot_vtpc1_;

  short y_field_id_;
TASK_DEF(PiddEdx, 0)
};

#endif //ATPIDTASK_PID_DEDX_PIDDEDX_H
