cmake_minimum_required(VERSION 3.13)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(miniproj miniproj.cc)
target_link_libraries(miniproj emdw)

add_executable(miniprojem miniprojem.cc)
target_link_libraries(miniprojem emdw)

add_executable(miniprojemR miniprojemR.cc)
target_link_libraries(miniprojemR emdw)

add_executable(miniprojFull miniprojFull.cc)
target_link_libraries(miniprojFull emdw)

add_executable(miniprojFull2 miniprojFull2.cc)
target_link_libraries(miniprojFull2 emdw)