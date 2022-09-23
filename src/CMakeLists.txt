#Add directories containing the header files for the code
include_directories(include cuda)
set(INC_DIR include)
set(CUDA_INC_DIR cuda)

#set the directory to source code
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
set(CUDA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cuda)
set(PYTHON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/python)
set(MAIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/main)
set(DEV_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dev)

set(SRC_FWI3D ${SRC_DIR}/injector.cpp
                ${SRC_DIR}/misc.cpp
                ${SRC_DIR}/nlsolver.cpp
                ${SRC_DIR}/operator.cpp
                ${SRC_DIR}/param.cpp
                ${SRC_DIR}/spatial_operators.cpp
                ${SRC_DIR}/we_op.cpp
                ${SRC_DIR}/bsplines.cpp
                ${SRC_DIR}/optimization.cpp
	)

if(GPU)
        set(SRC_FWI3D ${SRC_FWI3D} ${CUDA_DIR}/cudaKernels.cu
        ${CUDA_DIR}/cudaWeOp.cpp
        )
else()
#        set(SRC_FWI3D ${SRC_FWI3D} ${SRC_DIR}/dummy.cpp)
endif()

set(HEADERS ${INC_DIR}/injector.hpp
                ${INC_DIR}/IO.hpp
                ${INC_DIR}/misc.hpp
                ${INC_DIR}/nlsolver.hpp
                ${INC_DIR}/operator.hpp
                ${INC_DIR}/optimization.hpp
                ${INC_DIR}/param.hpp
                ${INC_DIR}/spatial_operators.hpp
                ${INC_DIR}/vecReg.hpp
                ${INC_DIR}/we_op.hpp
                ${INC_DIR}/bsplines.hpp
                ${INC_DIR}/mpiWrapper.hpp
        )

set(PYTHON ${PYTHON_DIR}/seppy.py
                ${PYTHON_DIR}/sep2numpy.py)

if (GPU)
        set(HEADERS ${HEADERS} ${CUDA_INC_DIR}/cudaCppWrappers.h
                ${CUDA_INC_DIR}/cudaMisc.h
        )
endif()

#building the project library
add_library(fwi3d ${SRC_FWI3D})
target_link_libraries(fwi3d fftw3f fftw3)
if (MPI_FOUND)
        target_link_libraries(fwi3d ${MPI_CXX_LIBRARIES} )
endif()

#building the main executables
add_executable(WE_MODELING_3D.x ${MAIN_DIR}/we_modeling.cpp)
target_link_libraries(WE_MODELING_3D.x fwi3d sep)

add_executable(FWI3D.x ${MAIN_DIR}/fwi3d.cpp)
target_link_libraries(FWI3D.x fwi3d sep)

add_executable(BSPLINES3D.x ${MAIN_DIR}/bsplines3d.cpp)
target_link_libraries(BSPLINES3D.x fwi3d sep)

#installing libraries and executables
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS fwi3d DESTINATION lib)
install(FILES ${PYTHON} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION bin)
install(TARGETS WE_MODELING_3D.x
                FWI3D.x
                BSPLINES3D.x
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION bin)


# For development purposes
if (DEV)
        add_executable(TEST_FD.x ${DEV_DIR}/test_fd.cpp)
        target_link_libraries(TEST_FD.x fwi3d sep)

        install(TARGETS TEST_FD.x
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION bin)
endif()