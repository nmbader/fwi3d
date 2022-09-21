#pragma once

#include "vecReg.hpp"
#ifdef ENABLE_MPI
    #include "mpi.h"
#endif

struct MpiWrapper
{
public:

    MpiWrapper() = delete;

    static int init( int * argc, char *** argv ){

#ifdef ENABLE_MPI
    return MPI_Init(argc,argv);
#else
    return 1;
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

    static int send(const data_t * pin, int n, int recipient){

#ifdef ENABLE_MPI
        MPI_Datatype mpi_type = MPI_FLOAT;
        if (sizeof(data_t) == 8) mpi_type = MPI_DOUBLE;
        return MPI_Send(pin, n, mpi_type, recipient, MPI_ANY_TAG, MPI_COMM_WORLD);
#else
        return 1;
#endif
    }

    static int recv(data_t * pout, int n, int sender){

#ifdef ENABLE_MPI
        MPI_Datatype mpi_type = MPI_FLOAT;
        if (sizeof(data_t) == 8) mpi_type = MPI_DOUBLE;
        return MPI_Recv(pout, n, mpi_type, sender, MPI_ANY_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
#else
        return 1;
#endif
    }

    static int gather(const data_t * pin, data_t * pout, int n){

#ifdef ENABLE_MPI
        MPI_Datatype mpi_type = MPI_FLOAT;
        if (sizeof(data_t) == 8) mpi_type = MPI_DOUBLE;
        return MPI_Gather(pin, n, mpi_type, pout, n, mpi_type, 0, MPI_COMM_WORLD);
#else
        memcpy(pout, pin, n*sizeof(data_t));
        return 1;
#endif
    }

    static int allReduceSum(const data_t * pin, data_t * pout, int n){

#ifdef ENABLE_MPI
        MPI_Datatype mpi_type = MPI_FLOAT;
        if (sizeof(data_t) == 8) mpi_type = MPI_DOUBLE;
        if (pin==pout) return MPI_Allreduce(MPI_IN_PLACE, pout, n, mpi_type, MPI_SUM, MPI_COMM_WORLD); 
        else return MPI_Allreduce(pin, pout, n, mpi_type, MPI_SUM, MPI_COMM_WORLD);
#else
        return 1;
#endif
    }
};