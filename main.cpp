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

double weightOfCluster(const vector<pair<pll, double>>& points) {
    double ret = 0;
    // for (auto[p, _]: points) {
    //     ret += p.second;
    // }
    ret = 1;
    ret *= sqrt(points.size());
    return ret;
}

// Function to calculate the dot product between two vectors
ll dotProduct(const pll& p1, const pll& p2) {
    return p1.first * p2.first + p1.second * p2.second;
}

// Function to calculate the center of mass of a group of points
pll calculateCenterOfMass(vector<pair<pll, double>>& points) {
    ll sumX = 0, sumY = 0;
    for (const auto& p : points) {
        sumX += p.first.first;
        sumY += p.first.second;
    }
    ll centerX = sumX / points.size();
    ll centerY = sumY / points.size();
    return {centerX, centerY};
}

// Function to find the point in group A mostly oriented in a given direction
int findPointInGroupA(const vector<pair<pll, double>>& groupA, const pll& direction) {
    int bestPointIdx = 0;
    ll maxDotProduct = dotProduct(groupA[0].first, direction);
    for (int i=0; i<groupA.size(); ++i) {
        const auto& p = groupA[i].first;
        ll currentDotProduct = dotProduct(p, direction);
        if (currentDotProduct > maxDotProduct) {
            maxDotProduct = currentDotProduct;
            bestPointIdx = i;
        }
    }
    return bestPointIdx;
}

void clustering(vector<vector<pair<pll, double>>>& groups, const vector<vector<int>>& adj) {
    int k = groups.size();
    for (int it=0; it<500; ++it) {
        vector<double> weights(k);
        for (int i=0; i<k; ++i) {
            weights[i] = weightOfCluster(groups[i]);
            weights[i] *= weights[i];
        }
        // print(weights);

        // Create a discrete distribution based on the weights
        discrete_distribution<int> distribution(weights.begin(), weights.end());

        // Create a random number generator
        random_device rd;
        mt19937 gen(rd());

        // Sample an index randomly based on the weights
        int target = distribution(gen);

        int min_target = -1;
        double min_weight = numeric_limits<double>::max();
        for (int neigh: adj[target]) {
            double weight = weightOfCluster(groups[neigh]);
            if (min_weight > weight) {
                min_weight = weight;
                min_target = neigh;
            }
        }

        pll COMfrom = calculateCenterOfMass(groups[target]);
        pll COMto = calculateCenterOfMass(groups[min_target]);
        pll v = { COMto.first - COMfrom.first, COMto.second - COMto.second };
        int bestPointIdx = findPointInGroupA(groups[target], v);
        auto cpy = groups[target][bestPointIdx];
        groups[target].erase(groups[target].begin() + bestPointIdx);
        groups[min_target].push_back(cpy);
    }
}

void solve() {
#ifdef SHARAELONG
    freopen("data/data1.in", "r", stdin);
#endif

    int n, k;
    cin >> n >> k;
    n = 8000;
    k = 140;
    
    vector<pll> points(n);
    map<pll, int> point_idx;
    for (int i=0; i<n; ++i) {
        int x, y;
        cin >> x >> y;
        points[i] = { x, y };
        point_idx[points[i]] = i+1;
    }

    const int MAX_COORD = 814000;
    int dx = (MAX_COORD + 10) / 10;
    int dy = (MAX_COORD + 14) / 14;

    vector<vector<int>> adj(k);
    for (int i=0; i<10; ++i) {
        for (int j=0; j<14; ++j) {
            int num = 14*i + j;
            if (i < 10 - 1) {
                adj[num].push_back(num + 14);
                adj[num + 14].push_back(num);
            }
            if (j < 14 - 1) {
                adj[num].push_back(num + 1);
                adj[num + 1].push_back(num);
            }
        }
    }
    
    vector<vector<pair<pll, double>>> groups;
    for (int i=0; i<10; ++i) {
        for (int j=0; j<14; ++j) {
            vector<pair<pll, double>> local_points;
            for (int l=0; l<n; ++l) {
                auto[x, y] = points[l];
                if (dx * i <= x && x < dx * (i+1) && dy * j <= y && y < dy * (j+1)) {
                    local_points.push_back({ points[l], 1.0 });
                }
            }
            groups.push_back(local_points);
        }
    }

    // clustering(groups, adj);

    vector<vector<pll>> clusters(k);
    for (int i=0; i<k; ++i) {
        for (auto[p, _]: groups[i]) {
            clusters[i].push_back(p);
        }
    }
    
    for (int i=0; i<k; ++i) {
        // Set the time limit to 32 milliseconds
        chrono::milliseconds timeLimit(32);

        // Run the 3-opt algorithm
        vector<int> tour = threeOpt(clusters[i], timeLimit);
        cout << tour.size() << ' ';
        for (int idx: tour) cout << point_idx[clusters[i][idx]] << ' ';
        cout << '\n';
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
}
