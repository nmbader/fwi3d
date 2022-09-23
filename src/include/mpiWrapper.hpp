#pragma once

#include <iostream>
#include "vecReg.hpp"
#ifdef ENABLE_MPI
    #include "mpi.h"
#endif

struct mpiWrapper
{
public:

    mpiWrapper() = delete;

    static void init( int * argc, char *** argv ){

#ifdef ENABLE_MPI
    MPI_Init(argc,argv);
#endif        
    }

    static void setSizeRank( int* size, int* rank ){

#ifdef ENABLE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,size);
    MPI_Comm_rank(MPI_COMM_WORLD,rank);
#else
    *size=1;
    *rank=0;
#endif
    }

    static void finalize(){

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif
    }

    static void barrier(){

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    static void send(const data_t * pin, int n, int recipient, int tag){

#ifdef ENABLE_MPI
        MPI_Datatype mpi_type = MPI_FLOAT;
        if (sizeof(data_t) == 8) mpi_type = MPI_DOUBLE;
        MPI_Send(pin, n, mpi_type, recipient, tag, MPI_COMM_WORLD);
#endif
    }

    static void recv(data_t * pout, int n, int sender, int tag){

#ifdef ENABLE_MPI
        MPI_Datatype mpi_type = MPI_FLOAT;
        if (sizeof(data_t) == 8) mpi_type = MPI_DOUBLE;
        MPI_Recv(pout, n, mpi_type, sender, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
    }

    static void gather(const data_t * pin, data_t * pout, int n){

#ifdef ENABLE_MPI
        MPI_Datatype mpi_type = MPI_FLOAT;
        if (sizeof(data_t) == 8) mpi_type = MPI_DOUBLE;
        MPI_Gather(pin, n, mpi_type, pout, n, mpi_type, 0, MPI_COMM_WORLD);
#else
        memcpy(pout, pin, n*sizeof(data_t));
#endif
    }

    static void allReduceSum(const data_t * pin, data_t * pout, int n){

#ifdef ENABLE_MPI
        MPI_Datatype mpi_type = MPI_FLOAT;
        if (sizeof(data_t) == 8) mpi_type = MPI_DOUBLE;
        if (pin==pout) MPI_Allreduce(MPI_IN_PLACE, pout, n, mpi_type, MPI_SUM, MPI_COMM_WORLD); 
        else MPI_Allreduce(pin, pout, n, mpi_type, MPI_SUM, MPI_COMM_WORLD);
#endif
    }
};