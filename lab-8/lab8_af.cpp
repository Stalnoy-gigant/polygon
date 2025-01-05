#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>

using namespace std;

// Узел графа
struct Node {
    double x, y; // координаты узла
    vector<pair<Node*, double>> neighbors; // соседи узла и веса рёбер

    Node(double x, double y) : x(x), y(y) {}
};

// Граф
class Graph {
private:
    unordered_map<string, Node*> nodes;

    string makeKey(double x, double y) {
        return to_string(x) + "," + to_string(y);
    }

public:
    ~Graph() {
        for (auto& [key, node] : nodes) {
            delete node;
        }
    }

    Node* addNode(double x, double y) {
        string key = makeKey(x, y);
        if (nodes.find(key) == nodes.end()) {
            nodes[key] = new Node(x, y);
        }
        return nodes[key];
    }

    void addEdge(double x1, double y1, double x2, double y2, double weight) {
        Node* node1 = addNode(x1, y1);
        Node* node2 = addNode(x2, y2);
        node1->neighbors.push_back({node2, weight});
    }

    Node* getNode(double x, double y) {
        string key = makeKey(x, y);
        return nodes.find(key) != nodes.end() ? nodes[key] : nullptr;
    }

    Node* find_closet_node(double x, double y){
        double min_x = 999999999, min_y = 9999999999;
        Node* node_f = nullptr;

        for (const auto& node: nodes){
            if (abs(node.second->x - x) <= min_x && abs(node.second->y - y) <= min_y){
                node_f = node.second;
                min_x = abs(node.second->x - x);
                min_y = abs(node.second->y - y);
            }
        }
        return node_f;
    }

    void readFromFile(const string& filename) {
        ifstream fin(filename);
        if (!fin.is_open()) {
            cerr << "Ошибка: не удалось открыть файл " << filename << endl;
            return;
        }

        string line;
        while (getline(fin, line)) {
            stringstream ss(line);
            string nodeData, neighbors;
            getline(ss, nodeData, ':');

            double x1, y1;
            char comma;
            stringstream nodeStream(nodeData);
            nodeStream >> x1 >> comma >> y1;

            while (getline(ss, neighbors, ';')) {
                double x2, y2, weight;
                stringstream neighborStream(neighbors);
                neighborStream >> x2 >> comma >> y2 >> comma >> weight;
                addEdge(x1, y1, x2, y2, weight);
            }
        }

        fin.close();
    }
};

// Поиск в глубину (DFS)
double DFS(Node* start, Node* goal) {
    if (!start || !goal) return -1.0;

    unordered_map<Node*, bool> visited;
    unordered_map<Node*, double> distances;
    vector<Node*> stack;

    stack.push_back(start);
    distances[start] = 0.0;

    while (!stack.empty()) {
        Node* current = stack.back();
        stack.pop_back();

        if (visited[current]) continue;
        visited[current] = true;

        if (current == goal) return distances[current];

        for (auto& [neighbor, weight] : current->neighbors) {
            if (!visited[neighbor]) {
                stack.push_back(neighbor);
                distances[neighbor] = distances[current] + weight;
            }
        }
    }

    return -1.0; // Путь не найден
}

// Поиск в ширину (BFS)
double BFS(Node* start, Node* goal) {
    if (!start || !goal) return -1.0;

    queue<Node*> q;
    unordered_map<Node*, double> distances;
    q.push(start);
    distances[start] = 0.0;

    while (!q.empty()) {
        Node* current = q.front(); q.pop();

        if (current == goal) return distances[current];

        for (auto& [neighbor, weight] : current->neighbors) {
            if (distances.find(neighbor) == distances.end()) {
                q.push(neighbor);
                distances[neighbor] = distances[current] + weight;
            }
        }
    }

    return -1.0; // Путь не найден
}

// Алгоритм Дейкстры
double Dijkstra(Node* start, Node* goal) {
    if (!start || !goal) return -1.0;

    unordered_map<Node*, double> distances;
    priority_queue<pair<double, Node*>, vector<pair<double, Node*>>, greater<>> pq;

    distances[start] = 0.0;
    pq.push({0.0, start});

    while (!pq.empty()) {
        auto [currentDist, current] = pq.top();
        pq.pop();

        if (current == goal) return currentDist;

        for (auto& [neighbor, weight] : current->neighbors) {
            double newDist = currentDist + weight;
            if (distances.find(neighbor) == distances.end() || newDist < distances[neighbor]) {
                distances[neighbor] = newDist;
                pq.push({newDist, neighbor});
            }
        }
    }

    return -1.0; // Путь не найден
}

// Алгоритм A*
double AStar(Node* start, Node* goal) {
    if (!start || !goal) return -1.0;

    auto heuristic = [](Node* a, Node* b) {
        return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2));
    };

    unordered_map<Node*, double> gScore;
    unordered_map<Node*, double> fScore;
    priority_queue<pair<double, Node*>, vector<pair<double, Node*>>, greater<>> pq;

    gScore[start] = 0.0;
    fScore[start] = heuristic(start, goal);
    pq.push({fScore[start], start});

    while (!pq.empty()) {
        auto [currentF, current] = pq.top();
        pq.pop();

        if (current == goal) return gScore[current];

        for (auto& [neighbor, weight] : current->neighbors) {
            double tentativeGScore = gScore[current] + weight;

            if (gScore.find(neighbor) == gScore.end() || tentativeGScore < gScore[neighbor]) {
                gScore[neighbor] = tentativeGScore;
                fScore[neighbor] = tentativeGScore + heuristic(neighbor, goal);
                pq.push({fScore[neighbor], neighbor});
            }
        }
    }

    return -1.0; // Путь не найден
}

void testing() {
    // Простой тест с графом, где все алгоритмы должны найти путь
    Graph graph;
    graph.addEdge(0.0, 0.0, 1.0, 1.0, 1.0);  // edge 1 -> 2
    graph.addEdge(1.0, 1.0, 2.0, 2.0, 1.0);  // edge 2 -> 3
    graph.addEdge(0.0, 0.0, 2.0, 2.0, 2.0);  // direct edge from 1 -> 3

    Node *start = graph.getNode(0.0, 0.0);
    Node *goal = graph.getNode(2.0, 2.0);

    // Test DFS
    auto dfsWeight = DFS(start, goal);
    assert(!(dfsWeight < 0));
    std::cout << "DFS path found with weight: " << dfsWeight << std::endl;

    // Test BFS
    auto bfsWeight = BFS(start, goal);
    assert(!(bfsWeight < 0));
    std::cout << "BFS path found with weight: " << bfsWeight << std::endl;

    // Test Dijkstra
    auto dijkstraWeight = Dijkstra(start, goal);
    assert(!(dijkstraWeight < 0));
    std::cout << "Dijkstra path found with weight: " << dijkstraWeight << std::endl;

    // Test A*
    auto aStarWeight = AStar(start, goal);
    assert(aStarWeight > 0);
    std::cout << "A* path found with weight: " << aStarWeight << std::endl;

    // Тест с несвязанным графом (путь не существует)
    Graph disconnectedGraph;
    disconnectedGraph.addEdge(0.0, 0.0, 1.0, 1.0, 1.0);
    Node *startDisconnected = disconnectedGraph.getNode(0.0, 0.0);
    Node *goalDisconnected = disconnectedGraph.getNode(2.0, 2.0);

    auto dfsWeightDisconnected = DFS(startDisconnected, goalDisconnected);
    assert(dfsWeightDisconnected < 0); // Путь не должен быть найден
    std::cout << "DFS: No path found in disconnected graph." << std::endl;

    auto bfsWeightDisconnected = BFS(startDisconnected, goalDisconnected);
    assert(bfsWeightDisconnected == -1); // Путь не должен быть найден
    std::cout << "BFS: No path found in disconnected graph." << std::endl;

    auto dijkstraWeightDisconnected = Dijkstra(startDisconnected, goalDisconnected);
    assert(dijkstraWeightDisconnected == -1); // Путь не должен быть найден
    std::cout << "Dijkstra: No path found in disconnected graph." << std::endl;

    auto aStarWeightDisconnected = AStar(startDisconnected, goalDisconnected);
    assert(aStarWeightDisconnected == -1); // Путь не должен быть найден
    std::cout << "A*: No path found in disconnected graph." << std::endl;

    // Тест с графом, где несколько путей имеют одинаковую длину
    Graph graphWithMultiplePaths;
    graphWithMultiplePaths.addEdge(0.0, 0.0, 1.0, 1.0, 1.0);
    graphWithMultiplePaths.addEdge(0.0, 0.0, 2.0, 2.0, 1.0);
    graphWithMultiplePaths.addEdge(1.0, 1.0, 2.0, 2.0, 1.0);

    Node *startMultiplePaths = graphWithMultiplePaths.getNode(0.0, 0.0);
    Node *goalMultiplePaths = graphWithMultiplePaths.getNode(2.0, 2.0);

    auto dfsWeightMultiplePaths = DFS(startMultiplePaths, goalMultiplePaths);
    assert(dfsWeightMultiplePaths > 0);
    std::cout << "DFS path (multiple paths) found with weight: " << dfsWeightMultiplePaths << std::endl;

    auto  bfsWeightMultiplePaths = BFS(startMultiplePaths, goalMultiplePaths);
    assert(dfsWeightMultiplePaths > 0);
    std::cout << "BFS path (multiple paths) found with weight: " << bfsWeightMultiplePaths << std::endl;

    auto dijkstraWeightMultiplePaths = Dijkstra(startMultiplePaths, goalMultiplePaths);
    assert(dijkstraWeightMultiplePaths > 0);
    std::cout << "Dijkstra path (multiple paths) found with weight: " << dijkstraWeightMultiplePaths << std::endl;

    auto aStarWeightMultiplePaths = AStar(startMultiplePaths, goalMultiplePaths);
    assert(aStarWeightMultiplePaths > 0);
    std::cout << "A* path (multiple paths) found with weight: " << aStarWeightMultiplePaths << std::endl;
}

int main() {
    Graph graph;

    // Чтение графа из файла
    graph.readFromFile(); // Вставить нормальный путь

    Node* start = graph.find_closet_node(30.372005, 59.934138);
    Node* goal = graph.find_closet_node(30.3132204, 59.9574788);

    // Поиск пути (BFS)
    double distanceBFS = BFS(start, goal);
    cout << "BFS: ";
    if (distanceBFS < 0) {
        cout << "path not find" << endl;
    } else {
        cout << "distance = " << distanceBFS << endl;
    }

    // Поиск пути (DFS)
    double distanceDFS = DFS(start, goal);
    cout << "DFS: ";
    if (distanceDFS < 0) {
        cout << "path not find!" << endl;
    } else {
        cout << "distance = " << distanceDFS << endl;
    }

    // Поиск пути (Dijkstra)
    double distanceDijkstra = Dijkstra(start, goal);
    cout << "Dijkstra: ";
    if (distanceDijkstra < 0) {
        cout << "path not find!" << endl;
    } else {
        cout << "distance = " << distanceDijkstra << endl;
    }

    // Поиск пути (A*)
    double distanceAStar = AStar(start, goal);
    cout << "A*: ";
    if (distanceAStar < 0) {
        cout << "path not find!" << endl;
    } else {
        cout << "distance = " << distanceAStar << endl;
    }

    return 0;
}

