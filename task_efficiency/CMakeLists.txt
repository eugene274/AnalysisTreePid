

add_executable(task_efficiency EvalEfficiency.cpp EvalEfficiency.hpp)
target_link_libraries(task_efficiency PRIVATE at_task_main)