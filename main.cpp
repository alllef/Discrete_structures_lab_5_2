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

struct AlgorithmResults {
public:
    vector<int> path;
    vector<double> distance;
};


void initializeGraph(int &picks, int &ribs, vector<Rib> &structRibs);

void sortRibs(int &picks, int &ribs, vector<Rib> &structRibs);

AlgorithmResults BelmanFordAlgorithm(int picks, vector<Rib> &structRibs, int startPick);

int findMin(const double *calculationsMatrix, int i, int picks);

int PickExists(int pick, vector<Rib> structRibs, vector<bool> isMarked, vector<bool> isRibUsed);

AlgorithmResults DijkstraAlgorithm(const int &picks, int &ribs, vector<Rib> &structRibs, int startPick);

AlgorithmResults JonsonAlgorithm(const int &picks, int &ribs, vector<Rib> &structRibs, int startPick);

void printShortestDistance(const vector<double> &distance, int startPick);

void printShortestPath(const vector<int> &path, int startPick, int endPick);

void makeMenu(int &picks, int &ribs, vector<Rib> &structRibs);

int main() {
    vector<Rib> ribsList;
    int n = 0, m = 0;
    SetConsoleOutputCP(CP_UTF8);
    makeMenu(n, m, ribsList);
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

AlgorithmResults BelmanFordAlgorithm(int picks, vector<Rib> &structRibs, int startPick) {
    vector<double> distance(picks, DBL_MAX);

    vector<int> prev(picks, -1);

    distance[startPick - 1] = 0;
    for (int i = 0; i < picks - 1; i++) {
        for (auto &structRib : structRibs) {
            if (distance[structRib.end - 1] > distance[structRib.start - 1] + structRib.weight) {
                distance[structRib.end - 1] = distance[structRib.start - 1] + structRib.weight;
                prev[structRib.end - 1] = structRib.start;
            }
        }
    }

    for (auto &structRib : structRibs) {
        if (distance[structRib.end - 1] > distance[structRib.start - 1] + structRib.weight) return {};
    }

    AlgorithmResults results;
    results.distance = distance;
    results.path = prev;
    return results;
}

AlgorithmResults JonsonAlgorithm(const int &picks, int &ribs, vector<Rib> &structRibs, int startPick) {
    vector<Rib> structRibsD = structRibs;
    AlgorithmResults results;
    for (int i = 0; i < picks; i++) {
        Rib tmpRib{};
        tmpRib.start = picks + 1;
        tmpRib.end = i + 1;
        tmpRib.weight = 0;
        structRibsD.push_back(tmpRib);
    }

    vector<double> h = BelmanFordAlgorithm(picks + 1, structRibsD, picks + 1).distance;

    for (auto &structRib : structRibs) {
        structRib.weight = structRib.weight + h[structRib.start - 1] - h[structRib.end - 1];
    }

    auto *distanceMatrix = new double[picks * picks];

    for (int i = 0; i < picks; i++) {
        vector<double> distance = DijkstraAlgorithm(picks, ribs, structRibs, (i + 1)).distance;
        if (i + 1 == startPick) results.path = DijkstraAlgorithm(picks, ribs, structRibs, (i + 1)).path;
        for (int j = 0; j < picks; j++) {
            *(distanceMatrix + i * picks + j) = distance[j];
        }
    }

    for (int i = 0; i < picks; i++) {
        for (int j = 0; j < picks; j++) {
            *(distanceMatrix + i * picks + j) = *(distanceMatrix + i * picks + j) - h[i] + h[j];
        }
    }
    vector<double> startPickDistance(picks);

    for (int i = 0; i < picks; i++) {
        startPickDistance[i] = (*(distanceMatrix + (startPick-1) * picks + i));
    }

    results.distance = startPickDistance;
    return results;
}

AlgorithmResults DijkstraAlgorithm(const int &picks, int &ribs, vector<Rib> &structRibs, int startPick) {

    double currentDistance = 0;
    int currentPick = startPick;
    vector<bool> isMarked(picks);
    vector<double> picksDistance(picks);
    auto *calculationsMatrix = new double[picks * picks];
    vector<bool> isRibUsed(structRibs.size());
    vector<int> path(picks, 0);

    for (int i = 0; i < picks; i++) {
        for (int j = 0; j < picks; j++) {
            *(calculationsMatrix + i * picks + j) = DBL_MAX;
        }
    }

    picksDistance[currentPick - 1] = currentDistance;
    isMarked[currentPick - 1] = true;
    *(calculationsMatrix + 0 * picks + (currentPick - 1)) = currentDistance;

    for (int i = 1; i < picks; i++) {
        while (PickExists(currentPick, structRibs, isMarked, isRibUsed)) {

            int ribIndex = PickExists(currentPick, structRibs, isMarked, isRibUsed) - 1;
            if (structRibs[ribIndex].weight < 0) return {};
            isRibUsed[ribIndex] = true;
            int pickIndex = structRibs[ribIndex].end - 1;

            if (*(calculationsMatrix + (i - 1) * picks + pickIndex) > currentDistance + structRibs[ribIndex].weight) {
                *(calculationsMatrix + i * picks + pickIndex) = currentDistance + structRibs[ribIndex].weight;
                path[pickIndex] = currentPick;
            } else *(calculationsMatrix + i * picks + pickIndex) = *(calculationsMatrix + (i - 1) * picks + pickIndex);
        }

        int minPickIndex = findMin(calculationsMatrix, i, picks);
        currentPick = minPickIndex + 1;
        currentDistance = *(calculationsMatrix + i * picks + minPickIndex);
        picksDistance[minPickIndex] = currentDistance;
        isMarked[minPickIndex] = true;
        for (auto &&ribUsed : isRibUsed) ribUsed = false;
    }

    AlgorithmResults results{};
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

void printShortestDistance(const vector<double> &distance, int startPick) {
    cout << "Відстань від " << startPick << " складає\n";
    for (int i = 0; i < distance.size(); i++) {
       if(i!=startPick-1) cout << "До вершини " << i + 1 << ": " << distance[i] << endl;
    }
}

void printShortestPath(const vector<int> &path, int startPick, int endPick) {
    int pickIndex = endPick - 1;
    vector<int> tmpVector;
    cout << "Шлях від " << startPick << " до " << endPick << ": ";
    while (path[pickIndex] != startPick) {
        tmpVector.push_back(path[pickIndex]);
        pickIndex = path[pickIndex] - 1;
    }
    cout << startPick<<setw(2);
    for (int i = tmpVector.size() - 1; i >= 0; i--){
        cout<<tmpVector[i]<<setw(2);
    }
    cout<<endPick;
}

void makeMenu(int &picks, int &ribs, vector<Rib> &structRibs) {

    int startPick;
    int endPick;
    bool isJonson;

    initializeGraph(picks, ribs, structRibs);
    sortRibs(picks, ribs, structRibs);

    cout << " Введіть 0 якщо хочете скористуватись алгоритмом Белмана-Форда та 1 якщо алгоритмом Джонсона";
    cin >> isJonson;
    cout << endl;

    if (isJonson) {
        cout << "Введіть вершину від 1 до " << picks << "від якої хочете почати виконувати алгоритм";

        cin >> startPick;
        AlgorithmResults results = JonsonAlgorithm(picks, ribs, structRibs, startPick);
        if (results.distance.empty()) cout << "В графі є відємні ваги" << endl;
        else {
            printShortestDistance(results.distance, startPick);
            cout << "Введіть маршрут до вершини який ви хочете побачити. Вершини можуть бути від 1 до " << picks;
            cin>>endPick;

            printShortestPath(results.path, startPick, endPick);
        }
    } else {
        cout << "Введіть вершину від 1 до " << picks << "від якої хочете почати виконувати алгоритм";
        cin >> startPick;
        cout << endl;
        AlgorithmResults results = BelmanFordAlgorithm(picks, structRibs, startPick);
        if (results.distance.empty()) cout << "В графі є відємні ваги" << endl;
        else {
            printShortestDistance(results.distance, startPick);
            cout << "Введіть маршрут до вершини який ви хочете побачити. Вершини можуть бути від 1 до " << picks;
            cin>>endPick;
            printShortestPath(results.path, startPick, endPick);
        }
    }
}
