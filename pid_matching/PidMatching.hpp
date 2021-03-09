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
#include <AnalysisTree/EventHeader.hpp>

#include <TEfficiency.h>

class PidMatching : public UserFillTask {

 public:
  boost::program_options::options_description GetBoostOptions() override;
  void PreInit() override;
  void PostFinish() override;

  void UserInit(std::map<std::string, void *> &map) override;
  void UserExec() override;
  void UserFinish() override;

 private:
  struct PidEfficiencyQAStruct;
  struct ChargedHadronsEfficiencyStruct;
  struct ValidateEfficiencyStruct;

  void InitEfficiencies();

  AnalysisTree::Matching *matching_ptr_{nullptr};
  ATI2::Branch *vtxt_branch{nullptr};
  ATI2::Branch *simt_branch{nullptr};

  std::map<int, std::shared_ptr<PidEfficiencyQAStruct>> efficiencies;
  std::map<int, std::shared_ptr<ValidateEfficiencyStruct>> validated_efficiencies;
  ChargedHadronsEfficiencyStruct *charged_hadrons_efficiency{nullptr};

  /* CONFIG */
  static bool opts_loaded;
  static std::string qa_file_name;
  static bool save_canvases;
  static std::string validate_file;

  TFile *qa_file_{nullptr};

 protected:
  virtual bool CheckSimTrack(const ATI2::BranchChannel &sim_track) const = 0;
  virtual bool CheckVtxTrack(const ATI2::BranchChannel &vtx_track) const = 0;

  ATI2::Variable vtxt_dca_x_;
  ATI2::Variable vtxt_dca_y_;
  ATI2::Variable vtxt_nhits_vtpc1_;
  ATI2::Variable vtxt_nhits_pot_vtpc1_;
  ATI2::Variable vtxt_nhits_vtpc2_;
  ATI2::Variable vtxt_nhits_pot_vtpc2_;
  ATI2::Variable vtxt_nhits_mtpc_;
  ATI2::Variable vtxt_charge;
  ATI2::Variable vtxt_nhits_pot_mtpc_;

  ATI2::Branch *mt_branch{nullptr};
  ATI2::Variable mt_y_cm_;
  ATI2::Variable mt_nhits_total_;
  ATI2::Variable mt_nhits_vtpc_;
  ATI2::Variable mt_nhits_ratio_;
  ATI2::Variable mt_nhits_pot_total_;
  ATI2::Variable mt_pid;
  ATI2::Variable mt_mass;
  ATI2::Variable mt_sim_y_cm_;
  ATI2::Variable mt_sim_pt_;
  ATI2::Variable mt_sim_phi_;

  ATI2::Variable sim_mother_id_;
  ATI2::Variable sim_pdg_;

  ATI2::Branch *simtproc_branch;
  ATI2::Variable simtproc_y_cm;
};

class PidMatching_NoCuts : public PidMatching {
 protected:
  bool CheckSimTrack(const ATI2::BranchChannel &sim_track) const override {
    return sim_track[sim_mother_id_].GetInt() == -1;
  }
  bool CheckVtxTrack(const ATI2::BranchChannel &vtx_track) const override {
    return true;
  }
 TASK_DEF(PidMatching_NoCuts, 0);
};

class PidMatching_StandardCuts : public PidMatching {
 protected:
  bool CheckSimTrack(const ATI2::BranchChannel &sim_track) const override {
    return sim_track[sim_mother_id_].GetInt() == -1;
  }
  bool CheckVtxTrack(const ATI2::BranchChannel &vtx_track) const override {
    int nhits_total =
        vtx_track[vtxt_nhits_vtpc1_].GetInt() +
        vtx_track[vtxt_nhits_vtpc2_].GetInt() +
        vtx_track[vtxt_nhits_mtpc_].GetInt();
    int nhits_vtpc =
        vtx_track[vtxt_nhits_vtpc1_].GetInt() +
        vtx_track[vtxt_nhits_vtpc2_].GetInt();
    int nhits_pot_total =
        vtx_track[vtxt_nhits_pot_vtpc1_].GetInt() +
        vtx_track[vtxt_nhits_pot_vtpc2_].GetInt() +
        vtx_track[vtxt_nhits_pot_mtpc_].GetInt();
    float ratio_nhits_nhits_pot = float(nhits_total) / float(nhits_pot_total);

    auto dca_x = vtx_track[vtxt_dca_x_];
    auto dca_y = vtx_track[vtxt_dca_y_];


    return
      nhits_total > 30 && nhits_vtpc > 15 &&
      nhits_pot_total > 0 &&
      ratio_nhits_nhits_pot > 0.55 && ratio_nhits_nhits_pot < 1.1 &&
      abs(dca_x) < 2. && abs(dca_y) < 1.;
  }
 TASK_DEF(PidMatching_StandardCuts, 0);
};

#endif //ATPIDTASK_PID_MATCHING_PIDMATCHING_HPP_
