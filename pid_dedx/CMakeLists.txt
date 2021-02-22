

add_executable(PiddEdx PiddEdx.cpp PiddEdx.h)
target_link_libraries(PiddEdx PUBLIC at_task_main Pid pid_new_core)