list(APPEND ALL_EXE_SRCS
  cas_point_contact_p.cpp
  cas_point_contact.cpp
  cas_flexion.cpp
  cas_troue.cpp
  cas_troue_p.cpp
  )


# Include the source directory where mfem.hpp and mfem-performance.hpp are.
include_directories(BEFORE ${PROJECT_BINARY_DIR})

# Add one executable per cpp file
add_mfem_examples(ALL_EXE_SRCS)

# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/convergence/carre.msh
  DESTINATION ${PROJECT_BINARY_DIR}/convergence/)
