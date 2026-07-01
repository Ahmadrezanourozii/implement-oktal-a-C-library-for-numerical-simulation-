target_sources(oktal PRIVATE ${CMAKE_CURRENT_LIST_DIR}/VtkExport.cpp)
target_link_libraries(oktal PUBLIC HDF5::HDF5)