
add_library(atpid_commons STATIC VtxTrackCut.cpp VtxTrackCut.hpp)
target_link_libraries(atpid_commons PUBLIC at_task ${ROOT_LIBRARIES})
target_include_directories(atpid_commons PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})