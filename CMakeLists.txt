cmake_minimum_required(VERSION 3.11)
project(ATPidTask)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

FetchContent_Declare(
        ATTaskSkeleton
        GIT_REPOSITORY git@github.com:eugene274/AnalysisTreeTaskSkeleton.git
        GIT_TAG        3577f45a7c3051e283bab18efc3fe6986d2bea81
)

FetchContent_MakeAvailable(ATTaskSkeleton)


set(Pid_BUILD_NEW ON)
FetchContent_Declare(
        Pid
        GIT_REPOSITORY git@git.cbm.gsi.de:pwg-c2f/analysis/pid.git
        GIT_TAG        master
)

FetchContent_MakeAvailable(Pid)


add_subdirectory(commons)
add_subdirectory(pid_dedx)
add_subdirectory(pid_matching)
add_subdirectory(task_efficiency)