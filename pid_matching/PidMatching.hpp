//
// Created by eugene on 09/02/2021.
//

#ifndef ATPIDTASK_PID_MATCHING_PIDMATCHING_HPP_
#define ATPIDTASK_PID_MATCHING_PIDMATCHING_HPP_

#include <at_task/UserTask.h>
#include <at_task/Task.h>

#include <AnalysisTree/Matching.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Particle.hpp>
#include <AnalysisTree/BranchConfig.hpp>

#include <TEfficiency.h>

class PidMatching : public UserFillTask {

 public:
  boost::program_options::options_description GetBoostOptions() override;
  void PreInit() override;
  void PostFinish() override;

  void Init(std::map<std::string, void *> &map) override;
  void Exec() override;
  void Finish() override;

 private:
  struct PidEfficiencyQAStruct;
  struct ChargedHadronsEfficiencyStruct;

  void InitEfficiencies();

  AnalysisTree::Matching *matching_ptr_{nullptr};
  AnalysisTree::TrackDetector *vtx_tracks_ptr{nullptr};
  AnalysisTree::TrackDetector *sim_track_ptr{nullptr};

  AnalysisTree::BranchConfig matched_particles_config_;
  AnalysisTree::Particles  *matched_particles_{nullptr};


  std::map<int, PidEfficiencyQAStruct *> efficiencies;
  ChargedHadronsEfficiencyStruct *charged_hadrons_efficiency{nullptr};





  TASK_DEF(PidMatching, 0);
  short y_cm_field_id_;
  short y_field_id_;
  short o_dca_x_field_id_;
  short o_dca_y_field_id_;
  short o_nhits_total_;
  short o_nhits_vtpc_;
  short o_nhits_ratio_;
  short o_nhits_pot_total_;
  short sim_track_pdg_id_;
  short i_dca_x_field_id_;
  short i_dca_y_field_id_;
  short i_nhits_vtpc1_;
  short i_nhits_pot_vtpc1_;
  short i_nhits_mtpc_;
  short i_nhits_pot_vtpc2_;
  short i_nhits_pot_mtpc_;
  short i_nhits_vtpc2_;
  short o_sim_y_cm_;
  short o_sim_pt_;
  short o_sim_phi_;
  int i_charge;
};

#endif //ATPIDTASK_PID_MATCHING_PIDMATCHING_HPP_
