list(APPEND ALL_EXE_SRCS
	troueAMR.cpp
	troueAMR_tps.cpp
	troueAMR_p.cpp
	troue.cpp
    contactAMR.cpp
  )


# Include the source directory where mfem.hpp and mfem-performance.hpp are.
include_directories(BEFORE ${PROJECT_BINARY_DIR})

# Add one executable per cpp file
add_mfem_examples(ALL_EXE_SRCS)
