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

  void Init(std::map<std::string, void *> &map) override;
  void Exec() override;
  void Finish() override;

 private:


 TASK_DEF(EvalEfficiency, 0)
};

#endif //ATPIDTASK_TASK_EFFICIENCY_EVALEFFICIENCY_HPP_
