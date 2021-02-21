//
// Created by eugene on 16/02/2021.
//

#ifndef ATPIDTASK_PID_MATCHING_TEFFICIENCYHELPER_HPP_
#define ATPIDTASK_PID_MATCHING_TEFFICIENCYHELPER_HPP_

#include <TList.h>
#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TAxis.h>
#include <memory>

inline
TList *
ProjectEfficiency(TEfficiency *eff) {
  auto list = new TList;
  list->SetOwner();

  std::unique_ptr<TH1> total_histo{eff->GetCopyTotalHisto()};
  std::unique_ptr<TH1> passed_histo{eff->GetCopyPassedHisto()};

  if (eff->GetDimension() == 1) {
    /* Convert to TH1* */
  } else if (eff->GetDimension() == 2) {
    /* List of TH1-s */
    auto x_axis = total_histo->GetXaxis();
    auto y_axis = total_histo->GetYaxis();

    /* export as-is */
    auto h2 = new TProfile2D(
        Form("%s_prof", eff->GetName()),
        eff->GetTitle(),
        x_axis->GetNbins(), x_axis->GetXmin(), x_axis->GetXmax(),
        y_axis->GetNbins(), y_axis->GetXmin(), y_axis->GetXmax());
    h2->SetDirectory(nullptr);
    h2->GetXaxis()->SetTitle(x_axis->GetTitle());
    h2->GetYaxis()->SetTitle(y_axis->GetTitle());
    h2->SetMinimum(0.);
    h2->SetMaximum(1.);

    for (int ix = 1; ix < x_axis->GetNbins(); ++ix) {
      for (int iy = 1; iy < y_axis->GetNbins(); ++iy) {
        auto xc = x_axis->GetBinCenter(ix);
        auto yc = y_axis->GetBinCenter(iy);
        bool is_any_passed = passed_histo->GetBinContent(ix, iy) > 0;
        if (is_any_passed)
          h2->Fill(xc, yc, eff->GetEfficiency(eff->GetGlobalBin(ix, iy)));
      }
    }
    list->Add(h2);

    for (int ix = 1; ix < x_axis->GetNbins(); ++ix) {
      auto h1 = new TProfile(
          Form("%s_%d", eff->GetName(), ix),
          Form("%s;  '%s' bin [%f, %f]", eff->GetTitle(), x_axis->GetTitle(),
               x_axis->GetBinLowEdge(ix), x_axis->GetBinUpEdge(ix)),
          y_axis->GetNbins(), y_axis->GetXmin(), y_axis->GetXmax());
      h1->GetXaxis()->SetTitle(y_axis->GetTitle());
      h1->SetDirectory(nullptr);

      for (int iy = 1; iy < y_axis->GetNbins(); ++iy) {
        auto yc = y_axis->GetBinCenter(iy);

        bool is_any_passed = passed_histo->GetBinContent(ix, iy) > 0;
        if (is_any_passed)
          h1->Fill(yc, eff->GetEfficiency(eff->GetGlobalBin(ix, iy)));
      }

      list->Add(h1);
    }

  } else if (eff->GetDimension() == 3) {
    /* List of TH2-s */
    auto x_axis = total_histo->GetXaxis();
    auto y_axis = total_histo->GetYaxis();
    auto z_axis = total_histo->GetZaxis();

    for (int ix = 1; ix < x_axis->GetNbins(); ++ix) {
      auto h2 = new TProfile2D(
          Form("%s_%d", eff->GetName(), ix),
          Form("%s Bin #%d [%f, %f]", eff->GetTitle(), ix, x_axis->GetBinLowEdge(ix), x_axis->GetBinUpEdge(ix)),
          y_axis->GetNbins(), y_axis->GetXmin(), y_axis->GetXmax(),
          z_axis->GetNbins(), z_axis->GetXmin(), z_axis->GetXmax());
      h2->GetXaxis()->SetTitle(y_axis->GetTitle());
      h2->GetYaxis()->SetTitle(z_axis->GetTitle());
      h2->SetDirectory(nullptr);
      h2->SetMinimum(0.);
      h2->SetMaximum(1.);

      for (int iy = 1; iy < y_axis->GetNbins(); ++iy) {
        for (int iz = 1; iz < z_axis->GetNbins(); ++iz) {
          auto yc = y_axis->GetBinCenter(iy);
          auto zc = z_axis->GetBinCenter(iz);
          bool is_any_passed = passed_histo->GetBinContent(ix, iy, iz) > 0;
          if (is_any_passed) {
            h2->Fill(yc, zc, eff->GetEfficiency(eff->GetGlobalBin(ix, iy, iz)));
          }

        }
      }

      list->Add(h2);
    } // ix
  } // ndim == 3

  return list;
}

#endif //ATPIDTASK_PID_MATCHING_TEFFICIENCYHELPER_HPP_
