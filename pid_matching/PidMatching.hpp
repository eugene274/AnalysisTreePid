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

  /* CONFIG */
  static bool opts_loaded;
  static std::string qa_file_name;
  static bool save_canvases;

  TFile *qa_file_{nullptr};

 protected:
  virtual bool CheckSimTrack(const AnalysisTree::Track &track) const = 0;
  virtual bool CheckVtxTrack(const AnalysisTree::Track &track) const = 0;

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

  short i_sim_mother_id;
};

class PidMatching_NoCuts : public PidMatching {
 protected:
  bool CheckSimTrack(const AnalysisTree::Track &track) const override {
    return track.GetField<int>(i_sim_mother_id) == -1;
  }
  bool CheckVtxTrack(const AnalysisTree::Track &track) const override {
    return true;
  }
 TASK_DEF(PidMatching_NoCuts, 0);
};

class PidMatching_StandardCuts : public PidMatching {
 protected:
  bool CheckSimTrack(const AnalysisTree::Track &track) const override {
    return track.GetField<int>(i_sim_mother_id) == -1;
  }
  bool CheckVtxTrack(const AnalysisTree::Track &track) const override {
    int nhits_total =
        track.GetField<int>(i_nhits_vtpc1_) +
        track.GetField<int>(i_nhits_vtpc2_) +
        track.GetField<int>(i_nhits_mtpc_);
    int nhits_vtpc =
        track.GetField<int>(i_nhits_vtpc1_) +
        track.GetField<int>(i_nhits_vtpc2_);
    int nhits_pot_total =
        track.GetField<int>(i_nhits_pot_vtpc1_) +
        track.GetField<int>(i_nhits_pot_vtpc2_) +
        track.GetField<int>(i_nhits_pot_mtpc_);

    float ratio_nhits_nhits_pot = float(nhits_total) / float(nhits_pot_total);

    auto dca_x = track.GetField<float>(i_dca_x_field_id_);
    auto dca_y = track.GetField<float>(i_dca_y_field_id_);


    return
      nhits_total > 30 && nhits_vtpc > 15 &&
      nhits_pot_total > 0 &&
      ratio_nhits_nhits_pot > 0.55 && ratio_nhits_nhits_pot < 1.1 &&
      abs(dca_x) < 2. && dca_y < 1.;
  }
 TASK_DEF(PidMatching_StandardCuts, 0);
};

#endif //ATPIDTASK_PID_MATCHING_PIDMATCHING_HPP_
