#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(AT-TPC)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(avalanche avalanche.cc)
target_link_libraries(avalanche Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file gem/mesh.header gem/mesh.elements gem/mesh.nodes gem/dielectrics.dat gem/gem.result)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

foreach(_file gas_chamber/mesh.header gas_chamber/mesh.elements gas_chamber/mesh.nodes gas_chamber/dielectrics.dat gas_chamber/gas.result)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

