#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;

void calculateTacts(int first, int second, vector<int>& tacts) {
    // Theoretical count formula
    // return (ceil(log2(2*n)) * (ceil(log2(n))) / 2);

    tacts[first] = tacts[second] = max(tacts[first], tacts[second]) + 1;

    // cout << "pair: " << first << "   " << second << endl;
    // cout << "tact of network: " << tacts[first] << endl;
    // cout << "-----------------" << endl;
}

int getNumTacts(const vector<int>& tacts){
    return *max_element(begin(tacts), end(tacts));
}


void BatcherMerge(int firstStartInd, int secondStartInd, int distance, int firstPartLen, int secondPartLen, vector<pair<int, int>>& comparators, vector<int>& tacts) 
{
    if (firstPartLen * secondPartLen < 1) {
        return;
    }
    if (firstPartLen == 1 && secondPartLen == 1) {
        comparators.push_back(make_pair(firstStartInd, secondStartInd));
        calculateTacts(firstStartInd, secondStartInd, tacts);
        return; 
    }

    BatcherMerge(firstStartInd, secondStartInd, 2 * distance, firstPartLen / 2 + firstPartLen % 2, secondPartLen / 2 + secondPartLen % 2, comparators, tacts);
    BatcherMerge(firstStartInd + distance, secondStartInd + distance, 2 * distance, firstPartLen / 2, secondPartLen / 2, comparators, tacts);

    for(int i = 1; i < firstPartLen - 1; i += 2) {
        comparators.push_back(make_pair(firstStartInd + distance * i, firstStartInd + distance * (i + 1)));
        calculateTacts(firstStartInd + distance * i, firstStartInd + distance * (i + 1), tacts);
    }
    int ind = 0;
    if (firstPartLen % 2 == 0) {
        comparators.push_back(make_pair(firstStartInd + distance *(firstPartLen - 1), secondStartInd));
        calculateTacts(firstStartInd + distance *(firstPartLen - 1), secondStartInd, tacts);
        ind = 1;
    }
    
    for(int i = ind; i < secondPartLen - 1; i += 2) {
        comparators.push_back(make_pair(secondStartInd + distance * i, secondStartInd + distance *(i + 1)));
        calculateTacts(secondStartInd + distance * i, secondStartInd + distance *(i + 1), tacts);
    }
}


void BatcherSort(int startPos, int distance, int length, vector<pair<int, int>>& comparators, vector<int>& tacts) 
{
    if (length < 2) {
        return;
    }

    BatcherSort(startPos, distance, length / 2, comparators, tacts);
    BatcherSort(startPos + distance * length / 2, distance, length / 2 + length % 2, comparators, tacts);
    BatcherMerge(startPos, startPos + distance * length / 2, distance, length / 2, length / 2 + length % 2, comparators, tacts);
}


void writeNetworkScheduleToFile(int n, const vector<pair<int, int>>& comparators, const vector<int>& tacts) {
    ofstream file("schedule.txt");
    if (file.is_open()) {
        file << n << " 0 0" << endl;
        for (const auto& pair : comparators) {
            file << pair.first << " " << pair.second << endl;
        }
        file << comparators.size() << endl << getNumTacts(tacts) << endl;

        file.close();
    }
    else {
        cout << "Cannot open file." << endl;
        exit(1);
    }
}

vector<vector<int>> generateBinaryArrays(int n) {
    vector<vector<int>> binaryArrays;

    for(int i = 1; i < pow(2, n); i++) {
        vector<int> binaryArray;

        int temp = i;
        for(int j = 0; j < n; j++) {
            binaryArray.insert(binaryArray.begin(), temp % 2);
            temp = temp / 2;
        }

        binaryArrays.push_back(binaryArray);
    }

    return binaryArrays;
}

bool checkSorted(const vector<int>& arr) {
    for(int i = 1; i < arr.size(); i++) {
        if (arr[i] < arr[i-1]) {
            return false;
        }
    }
    return true;
}

void testSchedule(int n) {
    ofstream inputArrFile("input_arrays.txt");
    if (!inputArrFile.is_open()) {
        cout << "Cannot open file." << endl;
        exit(1);
    }
    ofstream outputArrFile("output_arrays.txt");
    if (!outputArrFile.is_open()) {
        cout << "Cannot open file." << endl;
        exit(1);
    }
    ofstream comparatorsFile("comparators.txt");
        if (!comparatorsFile.is_open()) {
        cout << "Cannot open file." << endl;
        exit(1);
    }

    for(int i = 1; i < n + 1; i++) {
        cout << "iteration: " << i << endl;
        vector<vector<int>> binaryArrays = generateBinaryArrays(i);

        int k = 1;
        inputArrFile << endl << "-------------------------------( n = " << i << " )-------------------------------" << endl;
        for(const auto& array : binaryArrays) {
            inputArrFile << k << ".  ";
            for(int num : array) {
                inputArrFile << num << " ";
            }
            inputArrFile << endl;
            k++;
        }

        vector<pair<int, int>> comparators;
        vector<int> tacts = vector<int>(i, 0);
        BatcherSort(0, 1, i, comparators, tacts);
        comparatorsFile << endl << "-------------------------------( n = " << i << " )-------------------------------" << endl;
        for (const auto& pair : comparators) {
            comparatorsFile << pair.first << " " << pair.second << endl;
        }
        comparatorsFile << "num of comparators: " << comparators.size() << endl << "num of tacts: " << getNumTacts(tacts) << endl;

        for(auto& array : binaryArrays) {
            for(const auto& pair : comparators) {
                if (array[pair.first] > array[pair.second]) {
                    swap(array[pair.first], array[pair.second]);
                }
            }
        }

        k = 1;
        outputArrFile << endl << "-------------------------------( n = " << i << " )-------------------------------" << endl;
        for(const auto& array : binaryArrays) {
            if (!checkSorted(array)) {
                cout << "Test failed!" << endl;
                inputArrFile.close();
                comparatorsFile.close();
                outputArrFile.close();
                exit(1);
            }
            outputArrFile << k << ".  ";
            for(int num : array) {
                outputArrFile << num << " ";
            }
            outputArrFile << endl;
            k++;
        }
    }

    cout << "All tests succeeded!" << endl;
    inputArrFile.close();
    comparatorsFile.close();
    outputArrFile.close();
}

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        cout << "Invalid number of arguments." << endl;
        return 1;
    }

    int n = atoi(argv[1]);
    if (n < 1 || n > 10000) {
        cout << "Invalid value of n." << endl;
        return 1;
    }

    if (argc == 3) {
        if (strcmp(argv[2], "-t")) {
            cout << "Invalid flag." << endl;
            return 1;
        }
        else {
            testSchedule(n);
        }
    }
    else {
        vector<pair<int, int>> comparators;
        vector<int> tacts = vector<int>(n, 0);

        BatcherSort(0, 1, n, comparators, tacts);

        writeNetworkScheduleToFile(n, comparators, tacts);
    }

    return 0;
}
