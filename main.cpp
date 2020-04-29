#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include<Windows.h>
#include<limits>
#include <algorithm>
#include <iomanip>

using namespace std;

struct Rib {
public:
    int start;
    int end;
    double weight;
};

struct DijkstraAlgorithmResults {
public:
    vector<int> path;
    vector<double> distance;
};

void initializeGraph(int &picks, int &ribs, vector<Rib> &structRibs);

void sortRibs(int &picks, int &ribs, vector<Rib> &structRibs);

vector<double> BelmanFordAlgorithm(int picks, vector<Rib> &structRibs, int startPick);

int findMin(const double *calculationsMatrix, int i, int picks);

int PickExists(int pick, vector<Rib> structRibs, vector<bool> isMarked, vector<bool> isRibUsed);

DijkstraAlgorithmResults DijkstraAlgorithm(const int &picks, int &ribs, vector<Rib> &structRibs, int startPick);

void JonsonAlgorithm(const int &picks, int &ribs, vector<Rib> &structRibs);

void printMatrix(const double *matrix, int picks);

int main() {
    vector<Rib> ribsList;
    int n = 0, m = 0;
    SetConsoleOutputCP(CP_UTF8);
    initializeGraph(n, m, ribsList);
    sortRibs(n, m, ribsList);
    BelmanFordAlgorithm(n, ribsList, 1);
    JonsonAlgorithm(n, m, ribsList);
    return 0;
}

void initializeGraph(int &picks, int &ribs, vector<Rib> &structRibs) {
    Rib myRib{};
    ifstream inFile;
    inFile.open("GraphWithScales.txt");

    if (!inFile.is_open()) cout << "It is not open" << endl;
    inFile >> picks >> ribs;

    for (int i = 0; i < ribs; i++) {
        inFile >> myRib.start >> myRib.end >> myRib.weight;

        structRibs.push_back(myRib);
    }
    inFile.close();
}

void sortRibs(int &picks, int &ribs, vector<Rib> &structRibs) {

    Rib tmp{};
    for (int i = 0; i < ribs - 1; i++) {
        for (int j = 0; j < ribs - 1; j++) {
            if ((structRibs[j].start + structRibs[j].end) > (structRibs[j + 1].start + structRibs[j + 1].end)) {

                tmp = structRibs[j];
                structRibs[j] = structRibs[j + 1];
                structRibs[j + 1] = tmp;

            }
        }
    }
}

vector<double> BelmanFordAlgorithm(int picks, vector<Rib> &structRibs, int startPick, int endP) {
    vector<double> distance(picks);
    for (double &i : distance) {
        i = DBL_MAX;
    }
    distance[startPick - 1] = 0;
    for (int i = 0; i < picks - 1; i++) {
        for (auto &structRib : structRibs) {
            if (distance[structRib.end - 1] > distance[structRib.start - 1] + structRib.weight)
                distance[structRib.end - 1] = distance[structRib.start - 1] + structRib.weight;
        }
    }
    return distance;
}

void JonsonAlgorithm(const int &picks, int &ribs, vector<Rib> &structRibs) {
    vector<Rib> structRibsD = structRibs;

    for (int i = 0; i < picks; i++) {
        Rib tmpRib{};
        tmpRib.start = picks + 1;
        tmpRib.end = i + 1;
        tmpRib.weight = 0;
        structRibsD.push_back(tmpRib);
    }

    vector<double> h = BelmanFordAlgorithm(picks + 1, structRibsD, picks + 1);

    for (auto &structRib : structRibs) {
        structRib.weight = structRib.weight + h[structRib.start - 1] - h[structRib.end - 1];
    }

    auto *distanceMatrix = new double[picks * picks];

    for (int i = 0; i<picks;  i++) {
        vector<double> distance = DijkstraAlgorithm(picks, ribs, structRibs, (i + 1)).distance;
        for (int j = 0; j < picks; j++) {
            *(distanceMatrix + i * picks + j) = distance[j];
        }
    }

    for (int i = 0; i < picks; i++) {
        for (int j = 0; j < picks; j++) {
            *(distanceMatrix + i * picks + j) = *(distanceMatrix + i * picks + j) - h[i] + h[j];
        }
    }
}

DijkstraAlgorithmResults DijkstraAlgorithm(const int &picks, int &ribs, vector<Rib> &structRibs, int startPick) {

    double currentDistance = 0;
    int currentPick = startPick;
    vector<bool> isMarked(picks);
    vector<double> picksDistance(picks);
    auto *calculationsMatrix = new double[picks * picks];
    vector<bool> isRibUsed(structRibs.size());
    vector<int> path;

    for (int i = 0; i < picks; i++) {
        for (int j = 0; j < picks; j++) {
            *(calculationsMatrix + i * picks + j) = DBL_MAX;
        }
    }

    picksDistance[currentPick - 1] = currentDistance;
    isMarked[currentPick - 1] = true;
    *(calculationsMatrix + 0 * picks + (currentPick-1)) = currentDistance;
    path.push_back(currentPick);

    for (int i = 1; i < picks; i++) {
        while (PickExists(currentPick, structRibs, isMarked, isRibUsed)) {

            int ribIndex = PickExists(currentPick, structRibs, isMarked, isRibUsed) - 1;
            if (structRibs[ribIndex].weight < 0) return {};
            isRibUsed[ribIndex] = true;
            int pickIndex = structRibs[ribIndex].end - 1;

            if (*(calculationsMatrix + (i - 1) * picks + pickIndex) > currentDistance + structRibs[ribIndex].weight) {
                *(calculationsMatrix + i * picks + pickIndex) = currentDistance + structRibs[ribIndex].weight;
            } else *(calculationsMatrix + i * picks + pickIndex) = *(calculationsMatrix + (i - 1) * picks + pickIndex);
        }

        int minPickIndex = findMin(calculationsMatrix, i, picks);
        currentPick = minPickIndex + 1;
        path.push_back(currentPick);
        currentDistance = *(calculationsMatrix + i * picks + minPickIndex);
        picksDistance[minPickIndex] = currentDistance;
        isMarked[minPickIndex] = true;
        for(auto && ribUsed : isRibUsed) ribUsed=false;
    }

    DijkstraAlgorithmResults results{};
    results.distance = picksDistance;
    results.path = path;
    return results;
}

int PickExists(int pick, vector<Rib> structRibs, vector<bool> isMarked, vector<bool> isRibUsed) {

    for (int i = 0; i < structRibs.size(); i++) {
        if (structRibs[i].start == pick && isMarked[structRibs[i].end - 1] == false && isRibUsed[i] == false) {

            return i + 1;
        }
    }
    return 0;
}

int findMin(const double *calculationsMatrix, int i, int picks) {
    vector<double> tmpVector(picks);

    for (int j = 0; j < picks; j++) {
        tmpVector[j] = *(calculationsMatrix + i * picks + j);
    }

    auto result = std::min_element(tmpVector.begin(), tmpVector.end());
    return std::distance(tmpVector.begin(), result);
}

void printMatrix(const double *matrix, int picks) {
    for (int i = 0; i < picks; i++) {
        for (int j = 0; j < picks; j++){
            if (*(matrix + i * picks + j) != DBL_MAX) cout << *(matrix + i * picks + j) << setw(5);
            else cout << "∞"<<setw(5);
        }
        cout<< endl;
    }
}