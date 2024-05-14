#include <bits/stdc++.h>
#include <filesystem>

using namespace std;

// Function to calculate Euclidean distance between two points
double euclideanDistance(pair<float, float> p1, pair<float, float> p2) {
    return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
}

// Function to print the minimum spanning tree
void printMinimumSpanningTree(vector<vector<double>>& graph, int parent[]) {
    double MST = 0;

    for (int i = 1; i < graph.size(); i++) {
        MST += graph[i][parent[i]];
    }

    for (int i = 1; i < graph.size(); i++) {
        cout << parent[i] << " - " << i << " \t" << graph[i][parent[i]] << " \n";
    }
}

// Function to compute the minimum spanning tree
void minimumSpanningTree(vector<vector<double>>& graph, const std::string& output_filename) {
    int V = graph.size();
    bool visited[V] = {false};
    double weights[V];
    int parent[V];

    // Initialize weights and parent array
    for (int i = 0; i < V; i++) {
        weights[i] = DBL_MAX;
    }

    weights[0] = 0;
    parent[0] = -1;

    // Find the minimum spanning tree
    for (int i = 0; i < V - 1; i++) {
        int minVertexIndex = -1;
        double minWeight = DBL_MAX;

        // Find the minimum weight edge
        for (int j = 0; j < V; j++) {
            if (!visited[j] && weights[j] < minWeight) {
                minWeight = weights[j];
                minVertexIndex = j;
            }
        }

        if (minVertexIndex == -1) {
            cerr << "Graph is not connected.\n";
            return;
        }

        visited[minVertexIndex] = true;

        // Update weights and parent array
        for (int j = 0; j < V; j++) {
            if (graph[minVertexIndex][j] != 0 && !visited[j] && graph[minVertexIndex][j] < weights[j]) {
                parent[j] = minVertexIndex;
                weights[j] = graph[minVertexIndex][j];
            }
        }
    }

    // Save the minimum spanning tree graph to a text file in the "mstgraphs" directory
    std::string output_path = '<define output folder>' + output_filename;
    ofstream outfile(output_path);

    if (!outfile.is_open()) {
        cerr << "Error opening output file: " << output_path << "\n";
        return;
    }

    // Write the edges and weights to the output file
    for (int i = 1; i < V; i++) {
        outfile << parent[i] << " " << i << " " << graph[i][parent[i]] << " \n";
    }
}

int main() {
    // Specify the directory containing the files
    std::string directory_path = "<define input file>";
    namespace fs = std::filesystem;

    // Iterate through all files in the directory
    for (const auto& entry : fs::directory_iterator(directory_path)) {
        if (entry.path().extension() == ".txt") {
            // Read Gabriel graph from the current file
            std::vector<std::tuple<int, int, double>> edges;
            std::ifstream infile(entry.path());

            if (!infile.is_open()) {
                std::cerr << "Error opening file: " << entry.path() << "\n";
                continue;
            }

            int u, v;
            double weight;
            while (infile >> u >> v >> weight) {
                edges.push_back({u, v, weight});
            }

            int V = 0;
            for (const auto& edge : edges) {
                V = std::max({V, std::get<0>(edge), std::get<1>(edge)});
            }
            V++;  // Increment to account for 0-based indexing

            // Create and populate the graph
            std::vector<std::vector<double>> graph(V, std::vector<double>(V, 0.0));
            for (const auto& edge : edges) {
                int u = std::get<0>(edge);
                int v = std::get<1>(edge);
                double weight = std::get<2>(edge);
                graph[u][v] = graph[v][u] = weight;
            }

            // Run the minimum spanning tree algorithm
            std::string output_filename = entry.path().stem().string() + ".txt";
            minimumSpanningTree(graph, output_filename);
        }
    }

    return 0;
}
