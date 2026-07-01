add_executable( poisson poisson.cpp )
target_link_libraries( poisson PRIVATE oktal )


add_executable( tgv tgv.cpp )
target_link_libraries( tgv PRIVATE oktal )
