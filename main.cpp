#include <bits/stdc++.h>
#ifdef SHARAELONG
#include "debug.hpp"
#endif
using namespace std;
typedef long long ll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;

// Function to calculate the Euclidean distance between two points
ll calculateDistance(pll p1, pll p2) {
    ll dx = p1.first - p2.first;
    ll dy = p1.second - p2.second;
    return dx * dx + dy * dy;
}

// Function to find the nearest neighbor of a given point
int findNearestNeighbor(int current, vector<pll>& points, vector<bool>& visited) {
    int nearestNeighbor = -1;
    ll minDistance = numeric_limits<ll>::max();

    for (int i = 0; i < points.size(); i++) {
        if (!visited[i] && i != current) {
            ll distance = calculateDistance(points[current], points[i]);
            if (distance < minDistance) {
                minDistance = distance;
                nearestNeighbor = i;
            }
        }
    }

    return nearestNeighbor;
}

// Function to implement the 2-approximation algorithm for TSP
vector<int> tspApproximation(vector<pll>& points) {
    int n = points.size();
    vector<bool> visited(n, false);

    // Start at the first point
    int current = 0;
    visited[current] = true;

    // Store the order of visited points
    vector<int> tour;
    tour.push_back(current);

    // Find nearest neighbors and visit them in order
    for (int i = 0; i < n - 1; i++) {
        int nearestNeighbor = findNearestNeighbor(current, points, visited);
        visited[nearestNeighbor] = true;
        tour.push_back(nearestNeighbor);
        current = nearestNeighbor;
    }

    // Add the distance from the last point back to the starting point
    // tour.push_back(0);
    return tour;
}

void solve() {
    freopen("data.in", "r", stdin);

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
    const int SQRT_K = 11;

    set<int> test;
    int checksum = 0;
    
    int d = (MAX_COORD + SQRT_K) / SQRT_K;
    for (int i=0; i<SQRT_K; ++i) {
        for (int j=0; j<SQRT_K; ++j) {
            vector<pll> local_points;
            vector<int> ori_idx;
            for (int l=(k - SQRT_K * SQRT_K); l<n; ++l) {
                auto[p, idx] = points[l];
                auto[x, y] = p;
                if (d * i <= x && x < d * (i+1) && d * j <= y && y < d * (j+1)) {
                    local_points.push_back(p);
                    ori_idx.push_back(idx);
                }
            }
            vector<int> tour = tspApproximation(local_points);
            cout << local_points.size() << ' ';
            checksum += local_points.size();
            for (int idx: tour) {
                cout << ori_idx[idx] << ' ';
                if (test.find(ori_idx[idx]) != test.end()) {
                    // print(i, j);
                    assert(false);
                }
                test.insert(ori_idx[idx]);
            }
            cout << '\n';
        }
    }
    assert(checksum == 7981);

    for (int i=0; i<k-SQRT_K*SQRT_K; ++i) {
        cout << "1 " << i+1 << '\n';
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
}
