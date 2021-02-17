

add_executable(PiddEdx PiddEdx.cpp PiddEdx.h ../task_efficiency/EvalEfficiency.cpp ../task_efficiency/EvalEfficiency.hpp)
target_link_libraries(PiddEdx PUBLIC at_task_main Pid pid_new_core)