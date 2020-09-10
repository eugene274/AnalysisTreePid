cmake_minimum_required(VERSION 3.17)
project(ATPidTask)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

FetchContent_Declare(
        ATTaskSkeleton
        GIT_REPOSITORY git@github.com:eugene274/AnalysisTreeTaskSkeleton.git
        GIT_TAG         0.0.2
)

FetchContent_MakeAvailable(ATTaskSkeleton)

FetchContent_Declare(
        Pid
        GIT_REPOSITORY git@git.cbm.gsi.de:pwg-c2f/analysis/pid.git
        GIT_TAG        master
)

FetchContent_MakeAvailable(Pid)


add_subdirectory(pid_dedx)