#include <bits/stdc++.h>
#ifdef SHARAELONG
#include "debug.hpp"
#endif
using namespace std;
typedef long long ll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;

// Function to calculate the Euclidean distance between two points
double distance(const pll& p1, const pll& p2) {
    ll dx = p1.first - p2.first;
    ll dy = p1.second - p2.second;
    return sqrt(dx * dx + dy * dy);
}

// Function to calculate the total tour length
double tourLength(const vector<pll>& points, const vector<int>& tour) {
    double totalLength = 0.0;
    int n = tour.size();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        totalLength += distance(points[tour[i]], points[tour[j]]);
    }
    return totalLength;
}

// Function to perform the Chained Lin-Kernighan algorithm
vector<int> chainedLinKernighan(const vector<pll>& points, const chrono::milliseconds& timeLimit) {
    int n = points.size();

    // Initialize the tour to a simple sequential order
    vector<int> tour(n);
    for (int i = 0; i < n; ++i) {
        tour[i] = i;
    }

    double bestLength = numeric_limits<double>::max();
    vector<int> bestTour;

    // Calculate the initial tour length
    double currentLength = tourLength(points, tour);

    // Get the current time
    auto startTime = chrono::high_resolution_clock::now();

    // Main loop of the algorithm
    while (true) {
        bool improvement = false;

        // Try all possible 2-opt exchanges
        for (int i = 0; i < n - 2; ++i) {
            for (int j = i + 2; j < n; ++j) {
                // Calculate the change in tour length
                double delta = 0.0;
                delta += distance(points[tour[i]], points[tour[j]]);
                delta += distance(points[tour[i + 1]], points[tour[(j + 1) % n]]);
                delta -= distance(points[tour[i]], points[tour[i + 1]]);
                delta -= distance(points[tour[j]], points[tour[(j + 1) % n]]);

                // Update the tour length
                double newLength = currentLength + delta;

                // Update the tour if the new length is better
                if (newLength < bestLength) {
                    bestLength = newLength;
                    // Reverse the order of the tour segment between i and j
                    vector<int> newTour = tour;
                    reverse(newTour.begin() + i + 1, newTour.begin() + j + 1);
                    bestTour = newTour;
                    improvement = true;
                }
            }
        }

        // Break the loop if no improvement is made or time limit is exceeded
        auto currentTime = chrono::high_resolution_clock::now();
        auto elapsedTime = chrono::duration_cast<chrono::milliseconds>(currentTime - startTime);
        if (!improvement || elapsedTime > timeLimit) {
            break;
        }

        tour = bestTour;
        currentLength = bestLength;
    }

    return bestTour;
}

// // Function to calculate the Euclidean distance between two points
// double distance(const pll& p1, const pll& p2) {
//     ll dx = p1.first - p2.first;
//     ll dy = p1.second - p2.second;
//     return sqrt(dx * dx + dy * dy);
// }

// Function to perform the K-means clustering algorithm
vector<vector<pll>> kmeans(const vector<pll>& points, int k, int maxIterations) {
    int n = points.size();

    // Randomly initialize the centroids
    vector<pll> centroids(k);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(0, n - 1);
    for (int i = 0; i < k; ++i) {
        centroids[i] = points[dist(gen)];
    }

    vector<int> assignments(n, 0); // Cluster assignments for each point
    vector<int> counts(k, 0); // Number of points assigned to each cluster

    // Main loop of the algorithm
    int iterations = 0;
    while (iterations < maxIterations) {
        // Assign each point to the nearest centroid
        for (int i = 0; i < n; ++i) {
            double minDist = numeric_limits<double>::max();
            int bestCluster = 0;

            for (int j = 0; j < k; ++j) {
                double dist = distance(points[i], centroids[j]);
                if (dist < minDist) {
                    minDist = dist;
                    bestCluster = j;
                }
            }

            assignments[i] = bestCluster;
            counts[bestCluster]++;
        }

        // Update the centroids
        for (int i = 0; i < k; ++i) {
            centroids[i] = pll(0, 0);
        }

        for (int i = 0; i < n; ++i) {
            centroids[assignments[i]].first += points[i].first;
            centroids[assignments[i]].second += points[i].second;
        }

        for (int i = 0; i < k; ++i) {
            if (counts[i] > 0) {
                centroids[i].first /= counts[i];
                centroids[i].second /= counts[i];
            }
        }

        // Clear the counts
        counts.assign(k, 0);

        iterations++;
    }

    // Create the clusters
    vector<vector<pll>> clusters(k);
    for (int i = 0; i < n; ++i) {
        clusters[assignments[i]].push_back(points[i]);
    }

    return clusters;
}

void solve() {
#ifdef SHARAELONG
    // freopen("data/data.in", "r", stdin);
#endif

    int n, k;
    cin >> n >> k;
    n = 8000;
    k = 140;

    map<pll, int> mp;
    vector<pll> points(n);
    for (int i=0; i<n; ++i) {
        int x, y;
        cin >> x >> y;
        points[i] = { x, y };
        mp[points[i]] = i+1;
    }

    const int MAX_COORD = 814000;

    // Set the number of maximum iterations
    int maxIterations = 100;

    // Run the K-means clustering algorithm
    vector<vector<pll>> clusters = kmeans(points, k, maxIterations);

    for (int i=0; i<k; ++i) {
        // Set the time limit as INF
        chrono::milliseconds timeLimit(1000);
        
        // Run the Chained Lin-Kernighan algorithm
        vector<int> tour = chainedLinKernighan(clusters[i], timeLimit);
        
        cout << tour.size() << ' ';
        for (int idx: tour) {
            cout << mp[clusters[i][idx]] << ' ';
        }
        cout << '\n';
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
}
