//
// Created by eugene on 17/02/2021.
//

#include "EvalEfficiency.hpp"

TASK_IMPL(EvalEfficiency)

boost::program_options::options_description EvalEfficiency::GetBoostOptions() {
  return UserTask::GetBoostOptions();
}
void EvalEfficiency::PreInit() {
  UserTask::PreInit();
}
void EvalEfficiency::PostFinish() {
  UserTask::PostFinish();
}
void EvalEfficiency::Init(std::map<std::string, void *> &map) {

}
void EvalEfficiency::Exec() {

}
void EvalEfficiency::Finish() {

}
