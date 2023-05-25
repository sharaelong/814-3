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

// Function to perform the 2-opt algorithm
vector<int> twoopt(const vector<pll>& points, const chrono::milliseconds& timeLimit) {
    int n = points.size();

    double bestLength = numeric_limits<double>::max();
    vector<int> bestTour;

    for (int it=0; it<9; ++it) {
        // Initialize the tour to a simple sequential order
        vector<int> tour(n);
        for (int i = 0; i < n; ++i) {
            tour[i] = i;
        }
        
        // Shuffle the vector randomly
        random_device rd;
        mt19937 g(rd());
        shuffle(tour.begin(), tour.end(), g);

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
    }

    return bestTour;
}

void solve() {
#ifdef SHARAELONG
    freopen("data/data.in", "r", stdin);
#endif

    int n, k;
    cin >> n >> k;
    n = 8000;
    k = 140;
    
    vector<pair<pll, int>> points(n);
    for (int i=0; i<n; ++i) {
        int x, y;
        cin >> x >> y;
        points[i] = { { x, y }, i+1 };
    }

    const int MAX_COORD = 814000;

    set<int> s;
    
    int dx = (MAX_COORD + 10) / 10;
    int dy = (MAX_COORD + 14) / 14;
    for (int i=0; i<10; ++i) {
        for (int j=0; j<14; ++j) {
            vector<pll> local_points;
            vector<int> ori_idx;
            for (int l=0; l<n; ++l) {
                auto[p, idx] = points[l];
                auto[x, y] = p;
                if (dx * i <= x && x < dx * (i+1) && dy * j <= y && y < dy * (j+1)) {
                    local_points.push_back(p);
                    ori_idx.push_back(idx);
                }
            }

            // Set the time limit to 32 milliseconds
            chrono::milliseconds timeLimit(32);

            // Run the Chained 2-opt algorithm
            vector<int> tour = twoopt(local_points, timeLimit);
            cout << tour.size() << ' ';
            for (int idx: tour) cout << ori_idx[idx] << ' ';
            cout << '\n';
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
}
