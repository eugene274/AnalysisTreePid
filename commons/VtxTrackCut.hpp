//
// Created by eugene on 25/02/2021.
//

#ifndef ATPIDTASK_COMMONS_VTXTRACKCUT_HPP_
#define ATPIDTASK_COMMONS_VTXTRACKCUT_HPP_

#include <cassert>
#include <TTree.h>
#include <ati2/ATI2.hpp>

struct VtxTrackCut {
  /* absolute value of DCA */
  const float dcax_max;
  const float dcay_max;

  const int nhits_vtpc_min;
  const int nhits_total_min;

  const float ratio_nhits_nhits_pot_min;
  const float ratio_nhits_nhits_pot_max;

  ATI2::Variable v_dca_x;
  ATI2::Variable v_dca_y;
  ATI2::Variable v_nhits_vtpc1;
  ATI2::Variable v_nhits_vtpc2;
  ATI2::Variable v_nhits_mtpc;
  ATI2::Variable v_nhits_pot_vtpc1;
  ATI2::Variable v_nhits_pot_vtpc2;
  ATI2::Variable v_nhits_pot_mtpc;

  void InitBranch(ATI2::Branch *vtx_branch) {
    assert(vtx_branch);

    std::tie(v_nhits_vtpc1, v_nhits_vtpc2, v_nhits_mtpc,
             v_nhits_pot_vtpc1, v_nhits_pot_vtpc2, v_nhits_pot_mtpc,
             v_dca_x, v_dca_y) = vtx_branch->GetVars(
        "nhits_vtpc1", "nhits_vtpc2", "nhits_mtpc",
        "nhits_pot_vtpc1", "nhits_pot_vtpc2", "nhits_pot_mtpc",
        "dcax", "dcay");
  }

  bool CheckVtxTrack(ATI2::BranchChannel &vtx_track) const {
    int nhits_total =
        vtx_track[v_nhits_vtpc1].GetInt() +
            vtx_track[v_nhits_vtpc2].GetInt() +
            vtx_track[v_nhits_mtpc].GetInt();
    int nhits_vtpc =
        vtx_track[v_nhits_vtpc1].GetInt() +
            vtx_track[v_nhits_vtpc2].GetInt();
    int nhits_pot_total =
        vtx_track[v_nhits_pot_vtpc1].GetInt() +
            vtx_track[v_nhits_pot_vtpc2].GetInt() +
            vtx_track[v_nhits_pot_mtpc].GetInt();
    float ratio_nhits_nhits_pot = float(nhits_total) / float(nhits_pot_total);

    auto dca_x = vtx_track[v_dca_x];
    auto dca_y = vtx_track[v_dca_y];

    return
        nhits_total >= nhits_total_min &&
            nhits_vtpc > nhits_vtpc_min &&
            nhits_pot_total > 0 &&
            ratio_nhits_nhits_pot > ratio_nhits_nhits_pot_min &&
            ratio_nhits_nhits_pot < ratio_nhits_nhits_pot_max &&
            abs(dca_x) < dcax_max &&
            abs(dca_y) < dcay_max;
  }
};

#endif //ATPIDTASK_COMMONS_VTXTRACKCUT_HPP_float ratio_nhits_nhits_pot_min = 0.55;
