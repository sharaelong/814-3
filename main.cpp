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

// Function to perform the 3-opt algorithm
vector<int> threeOpt(const vector<pll>& points, const chrono::milliseconds& timeLimit) {
    int n = points.size();

    // vector<int> bestTour;

    // Initialize the tour to a simple sequential order
    vector<int> tour(n);
    for (int i = 0; i < n; ++i) {
        tour[i] = i;
    }
        
    // Shuffle the vector randomly
    // random_device rd;
    // mt19937 g(rd());
    // shuffle(tour.begin(), tour.end(), g);

    // Calculate the initial tour length
    // double currentLength = tourLength(points, tour);
    // double bestLength = currentLength;

    // Get the current time
    auto startTime = chrono::high_resolution_clock::now();

    // Main loop of the algorithm
    while (true) {
        bool improvement = false;

        // Try all possible 3-opt exchanges
        for (int i = 0; i < n - 2; ++i) {
            for (int j = i + 2; j < n - 1; ++j) {
                for (int k = j + 2; k < n - 1 + (i > 0); ++k) {
                    const auto& a = points[tour[i]];
                    const auto& b = points[tour[i+1]];
                    const auto& c = points[tour[j]];
                    const auto& d = points[tour[j+1]];
                    const auto& e = points[tour[k]];
                    const auto& f = points[tour[(k+1) % n]];
                    
                    // Calculate the change in tour length
                    double originalLength = 0.0;
                    originalLength += distance(a, b);
                    originalLength += distance(c, d);
                    originalLength += distance(e, f);

                    double c1 = 0.0; // A'BC
                    c1 += distance(a, c);
                    c1 += distance(b, d);
                    c1 += distance(e, f);
                    

                    double c2 = 0.0; // AB'C
                    c2 += distance(a, b);
                    c2 += distance(c, e);
                    c2 += distance(d, f);
                    

                    double c3 = 0.0; // ABC'
                    c3 += distance(a, e);
                    c3 += distance(b, f);
                    c3 += distance(c, d);
                    

                    double c4 = 0.0; // A'B'C
                    c4 += distance(a, c);
                    c4 += distance(b, e);
                    c4 += distance(d, f);
                    

                    double c5 = 0.0; // A'BC'
                    c5 += distance(a, e);
                    c5 += distance(b, d);
                    c5 += distance(c, f);

                    double c6 = 0.0; // AB'C'
                    c6 += distance(a, d);
                    c6 += distance(b, f);
                    c6 += distance(c, e);

                    double c7 = 0.0; // A'B'C'
                    c7 += distance(a, d);
                    c7 += distance(b, e);
                    c7 += distance(c, f);

                    double mn = min({ c1, c2, c3, c4, c5, c6, c7 });
                    if (mn < originalLength) {
                        improvement = true;
                        if (c1 == mn) {
                            reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                        }
                        if (c2 == mn) {
                            reverse(tour.begin() + j + 1, tour.begin() + k + 1);
                        }
                        if (c3 == mn) {
                            reverse(tour.begin() + i + 1, tour.begin() + k + 1); // change body, not segment. It is symmetric.
                        }
                        if (c4 == mn) {
                            reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                            reverse(tour.begin() + j + 1, tour.begin() + k + 1);
                        }
                        if (c5 == mn) {
                            reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                            reverse(tour.begin() + i + 1, tour.begin() + k + 1);
                        }
                        if (c6 == mn) {
                            reverse(tour.begin() + j + 1, tour.begin() + k + 1);
                            reverse(tour.begin() + i + 1, tour.begin() + k + 1);
                        }
                        if (c7 == mn) {
                            reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                            reverse(tour.begin() + j + 1, tour.begin() + k + 1);
                            reverse(tour.begin() + i + 1, tour.begin() + k + 1);
                        }
                    }
                }
            }
        }

        // Break the loop if no improvement is made or time limit is exceeded
        auto currentTime = chrono::high_resolution_clock::now();
        auto elapsedTime = chrono::duration_cast<chrono::milliseconds>(currentTime - startTime);
        if (!improvement || elapsedTime > timeLimit) {
            break;
        }
    }

    return tour;
}

void solve() {
#ifdef SHARAELONG
    // freopen("data/data.in", "r", stdin);
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

            // Run the 3-opt algorithm
            vector<int> tour = threeOpt(local_points, timeLimit);
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
