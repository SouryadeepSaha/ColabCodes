
add_library(
        Source
        Grid.h
        Poisson_FD.h

)


set_target_properties(Source PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(Source )