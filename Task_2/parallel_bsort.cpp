#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <vector>
#include <random>
#include <cfloat>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <mpi.h>

using namespace std;


double generateRandomDouble(double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}


bool checkSorted(const std::vector<double>& dataArray) {
    bool sorted = true;

    #pragma omp parallel for
    for(int i = 1; i < static_cast<int>(dataArray.size()); i++) {
        if (dataArray[i] < dataArray[i-1]) {
            #pragma omp critical
            {
                sorted = false;
            }
        }
    }

    return sorted;
}


void BatcherMerge(int firstStartInd, int secondStartInd, int distance, int firstPartLen, int secondPartLen, vector<pair<int, int>>& comparators) 
{
    if (firstPartLen * secondPartLen < 1) {
        return;
    }
    if (firstPartLen == 1 && secondPartLen == 1) {
        comparators.push_back(make_pair(firstStartInd, secondStartInd));
        return; 
    }

    BatcherMerge(firstStartInd, secondStartInd, 2 * distance, firstPartLen / 2 + firstPartLen % 2, secondPartLen / 2 + secondPartLen % 2, comparators);
    BatcherMerge(firstStartInd + distance, secondStartInd + distance, 2 * distance, firstPartLen / 2, secondPartLen / 2, comparators);

    for(int i = 1; i < firstPartLen - 1; i += 2) {
        comparators.push_back(make_pair(firstStartInd + distance * i, firstStartInd + distance * (i + 1)));
    }
    int ind = 0;
    if (firstPartLen % 2 == 0) {
        comparators.push_back(make_pair(firstStartInd + distance *(firstPartLen - 1), secondStartInd));
        ind = 1;
    }
    
    for(int i = ind; i < secondPartLen - 1; i += 2) {
        comparators.push_back(make_pair(secondStartInd + distance * i, secondStartInd + distance *(i + 1)));
    }
}


void BatcherSort(int startPos, int distance, int length, vector<pair<int, int>>& comparators) 
{
    if (length < 2) {
        return;
    }

    BatcherSort(startPos, distance, length / 2, comparators);
    BatcherSort(startPos + distance * length / 2, distance, length / 2 + length % 2, comparators);
    BatcherMerge(startPos, startPos + distance * length / 2, distance, length / 2, length / 2 + length % 2, comparators);
}


void Comparator(int proc_ind1, int proc_ind2, int rank, vector<double>& localData) {
    if (rank != proc_ind1 && rank != proc_ind2)
        return;

    MPI_Request sendRequest;
    MPI_Status getStatus;
    int chunkSize = localData.size();
    vector<double> receivedData(chunkSize);
    vector<double> tmp(chunkSize);

    if (rank == proc_ind1) {
        // MPI_Sendrecv(localData.data(), chunkSize, MPI_DOUBLE, proc_ind2, 0, receivedData.data(), chunkSize, MPI_DOUBLE, proc_ind1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Isend(localData.data(), chunkSize, MPI_DOUBLE, proc_ind2, rank, MPI_COMM_WORLD, &sendRequest);
        MPI_Recv(receivedData.data(), chunkSize, MPI_DOUBLE, proc_ind2, proc_ind2, MPI_COMM_WORLD, &getStatus);
        MPI_Wait(&sendRequest, &getStatus);

        // #pragma omp parallel for
        for (int ia = 0, ib = 0, k = 0; k < chunkSize; ++k) {
            // int ia = 0, ib = 0;
            if (localData[ia] < receivedData[ib]) 
                tmp[k] = localData[ia++];
            else
                tmp[k] = receivedData[ib++];
        }

        localData = tmp;
    } else {
        // MPI_Sendrecv(localData.data(), chunkSize, MPI_DOUBLE, proc_ind1, 0, receivedData.data(), chunkSize, MPI_DOUBLE, proc_ind2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Isend(localData.data(), chunkSize, MPI_DOUBLE, proc_ind1, rank, MPI_COMM_WORLD, &sendRequest);
        MPI_Recv(receivedData.data(), chunkSize, MPI_DOUBLE, proc_ind1, proc_ind1, MPI_COMM_WORLD, &getStatus);
        MPI_Wait(&sendRequest, &getStatus);

        // #pragma omp parallel for
        for (int ia = chunkSize - 1, ib = chunkSize - 1, k = chunkSize - 1; k >= 0; --k) {
            // int ia = chunkSize - 1, ib = chunkSize - 1;
            if (localData[ia] > receivedData[ib]) 
                tmp[k] = localData[ia--];
            else 
                tmp[k] = receivedData[ib--];
        }

        localData = tmp;
    }
}


int main(int argc, char** argv) {
    string sizeFlag = "--size=";
    int vectorSize = 1000;

    const int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);

    // cout << num_threads << endl;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg.compare(0, sizeFlag.length(), sizeFlag) == 0) {
            vectorSize = atoi(arg.substr(sizeFlag.length()).c_str());
            break;
        }
    }

    MPI_Init(&argc, &argv);

    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int elementsPerProcess = vectorSize / numProcesses;
    vector<double> localData(elementsPerProcess);

    #pragma omp parallel for
    for (int i = 0; i < elementsPerProcess; i++) {
        localData[i] = generateRandomDouble(0.0, 10000.0);
    }

    vector<pair<int, int>> comparators;
    BatcherSort(0, 1, numProcesses, comparators);

    auto start = MPI_Wtime();
    sort(localData.begin(), localData.end());

    // MPI_Barrier(MPI_COMM_WORLD);
    for (const auto& pair : comparators) {
        Comparator(pair.first, pair.second, rank, localData);
    }
    auto end = MPI_Wtime();
    double delta = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    vector<double> dataArray(vectorSize);
    MPI_Gather(localData.data(), elementsPerProcess, MPI_DOUBLE, dataArray.data(), elementsPerProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double max_time;
    MPI_Reduce(&delta, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = numProcesses * elementsPerProcess; i < vectorSize; i++)
            dataArray[i] = DBL_MAX;

        // cout << "Sorted Array: ";
        // for (const auto& element : dataArray) {
        //     cout << element << " ";
        // }
        // cout << endl;

        // if (checkSorted(dataArray))
        //     cout << "Parallel sort succeded!"<< endl;
        // else
        //     cout << "ERROR: vector is not sorted!" << endl;

        cout << "Sorting time in seconds = " << max_time << endl;
    }

    MPI_Finalize();
    return 0;
}