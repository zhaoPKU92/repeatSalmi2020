cmake_minimum_required(VERSION 3.12)
project(salmi_paper LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(OpenMP REQUIRED)
add_library(salmi_core
  src/numerics.cpp
  src/atmosphere.cpp
  src/beam.cpp
  src/stopping.cpp
  src/radiation.cpp
  src/compton.cpp
  src/redistribution.cpp
  src/redistribution_bins.cpp
  src/checkpoint.cpp
  src/gaunt.cpp
  src/snapshot.cpp
  src/anderson.cpp
)
target_include_directories(salmi_core PUBLIC include)
target_link_libraries(salmi_core PUBLIC OpenMP::OpenMP_CXX)
add_executable(salmi_atmosphere
  experiments/paper_figure3/main.cpp
  experiments/paper_figure3/kernel_library.cpp
  experiments/paper_figure3/hydrostatic.cpp
  experiments/paper_figure3/heating.cpp
  experiments/paper_figure3/opacity.cpp
  experiments/paper_figure3/transfer.cpp
  experiments/paper_figure3/temperature.cpp
  experiments/paper_figure3/balanced_temperature.cpp
  experiments/paper_figure3/coupled_probe.cpp
  experiments/paper_figure3/atlas_numerics.cpp
)
target_include_directories(salmi_atmosphere PRIVATE experiments/paper_figure3)
target_link_libraries(salmi_atmosphere PRIVATE salmi_core)
add_executable(salmi_figure1 app/figure1.cpp)
target_link_libraries(salmi_figure1 PRIVATE salmi_core)
if(MINGW)
  set_target_properties(salmi_atmosphere salmi_figure1 PROPERTIES LINK_FLAGS "-static")
endif()
