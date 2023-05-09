set(SOURCES Grid.cpp Poisson_FD.cpp Poisson_GMG.cpp ProblemDesc.cpp GMG_solver.cpp)
set(HEADERS Grid.h Poisson_FD.h Poisson_GMG.h ProblemDesc.h GMG_solver.h)

add_library(
        Source
        ${SOURCES}
        ${HEADERS}
)


set_target_properties(Source PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(Source )