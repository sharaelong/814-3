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

double tourLength(const vector<pll>& points, const vector<int>& tour) {
    double ret = 0;
    int n = points.size();
    for (int i=0; i<n; ++i) {
        ret += distance(points[tour[i]], points[tour[(i+1) % n]]);
    }
    return ret;
}

// Function to perform the 3-opt algorithm
pair<double, vector<int>> threeOpt(const vector<pll>& points) {
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
        // auto currentTime = chrono::high_resolution_clock::now();
        // auto elapsedTime = chrono::duration_cast<chrono::milliseconds>(currentTime - startTime);
        // if (!improvement || elapsedTime > timeLimit) {
        if (!improvement) {
            break;
        }
    }

    return { tourLength(points, tour), tour };
}

// Function to calculate the dot product between two vectors
ll dotProduct(const pll& p1, const pll& p2) {
    return p1.first * p2.first + p1.second * p2.second;
}

// Function to calculate the center of mass of a group of points
pll calculateCenterOfMass(const vector<pll>& points) {
    ll sumX = 0, sumY = 0;
    for (const auto& p : points) {
        sumX += p.first;
        sumY += p.second;
    }
    ll centerX = sumX / points.size();
    ll centerY = sumY / points.size();
    return {centerX, centerY};
}

// Function to find the point in group A mostly oriented in a given direction
int findPointInGroupA(const vector<pll>& groupA, const pll& direction) {
    int bestPointIdx = 0;
    ll maxDotProduct = dotProduct(groupA[0], direction);
    for (int i=0; i<groupA.size(); ++i) {
        const auto& p = groupA[i];
        ll currentDotProduct = dotProduct(p, direction);
        if (currentDotProduct > maxDotProduct) {
            maxDotProduct = currentDotProduct;
            bestPointIdx = i;
        }
    }
    return bestPointIdx;
}

void solve() {
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
    
    vector<vector<pll>> groups;
    for (int i=0; i<10; ++i) {
        for (int j=0; j<14; ++j) {
            vector<pll> local_points;
            for (int l=0; l<n; ++l) {
                auto[x, y] = points[l];
                if (dx * i <= x && x < dx * (i+1) && dy * j <= y && y < dy * (j+1)) {
                    local_points.push_back(points[l]);
                }
            }
            groups.push_back(local_points);
        }
    }

    // clustering(groups, adj);

    vector<double> tourLen(k);
    for (int i=0; i<k; ++i) {
        auto result = threeOpt(groups[i]);
        tourLen[i] = result.first;
    }

    for (int it=1; it<=200; ++it) {
        random_device rd;
        mt19937 gen(rd());

        if (it & 1) {
        vector<double> prob(k);
        for (int i=0; i<k; ++i) {
            prob[i] = exp((tourLen[i]-400000) / 30000);
        }
        discrete_distribution<int> distribution(prob.begin(), prob.end());
        
        // int from = max_element(tourLen.begin(), tourLen.end()) - tourLen.begin();
        int from = distribution(gen);
        
        vector<double> tmp;
        vector<int> tmp_idx;
        for (int there: adj[from]) {
            if (tourLen[there] < tourLen[from]) {
                tmp.push_back(tourLen[there]);
                tmp_idx.push_back(there);
            }
        }
        if (tmp.empty()) continue;

        uniform_int_distribution<int> distribution2(0, tmp.size()-1);
        
        // int to = tmp_idx[min_element(tmp.begin(), tmp.end()) - tmp.begin()];
        int to = tmp_idx[distribution2(gen)];

        pll COMfrom = calculateCenterOfMass(groups[from]);
        pll COMto = calculateCenterOfMass(groups[to]);
        pll v = { COMto.first - COMfrom.first, COMto.second - COMto.second };
        int bestPointIdx = findPointInGroupA(groups[from], v);
        auto cpy = groups[from][bestPointIdx];
        groups[from].erase(groups[from].begin() + bestPointIdx);
        groups[to].push_back(cpy);

        double fromcpy = tourLen[from];
        double tocpy = tourLen[to];
        // print(it & 1, tourLen[from], tourLen[to]);
        tourLen[from] = threeOpt(groups[from]).first;
        tourLen[to] = threeOpt(groups[to]).first;
        if (max(tourLen[from], tourLen[to]) > fromcpy + (30000 / sqrt(it))) {
            tourLen[from] = fromcpy;
            tourLen[to] = tocpy;
            groups[from].push_back(groups[to].back());
            groups[to].pop_back();
        }
        // print(it & 1, tourLen[from], tourLen[to]);
        }
        // } else {
        //     vector<double> prob(k);
        // for (int i=0; i<k; ++i) {
        //     prob[i] = exp((400000-tourLen[i]) / 30000);
        // }
        // discrete_distribution<int> distribution(prob.begin(), prob.end());
        
        // // int from = max_element(tourLen.begin(), tourLen.end()) - tourLen.begin();
        // int from = distribution(gen);
        
        // vector<double> tmp;
        // vector<int> tmp_idx;
        // for (int there: adj[from]) {
        //     if (tourLen[there] > tourLen[from]) {
        //         tmp.push_back(tourLen[there]);
        //         tmp_idx.push_back(there);
        //     }
        // }
        // if (tmp.empty()) continue;

        // uniform_int_distribution<int> distribution2(0, tmp.size()-1);
        
        // // int to = tmp_idx[min_element(tmp.begin(), tmp.end()) - tmp.begin()];
        // int to = tmp_idx[distribution2(gen)];

        // pll COMfrom = calculateCenterOfMass(groups[from]);
        // pll COMto = calculateCenterOfMass(groups[to]);
        // pll v = { COMto.first - COMfrom.first, COMto.second - COMto.second };
        // int bestPointIdx = findPointInGroupA(groups[from], v);
        // auto cpy = groups[from][bestPointIdx];
        // groups[from].erase(groups[from].begin() + bestPointIdx);
        // groups[to].push_back(cpy);

        // double fromcpy = tourLen[from];
        // double tocpy = tourLen[to];
        // // print(it & 1, tourLen[from], tourLen[to]);
        // tourLen[from] = threeOpt(groups[from]).first;
        // tourLen[to] = threeOpt(groups[to]).first;
        // if ((fromcpy < tourLen[from] && tocpy < tourLen[to]) || min(tourLen[from], tourLen[to]) < fromcpy - (30000 / sqrt(it))) {
        //     tourLen[from] = fromcpy;
        //     tourLen[to] = tocpy;
        //     groups[from].push_back(groups[to].back());
        //     groups[to].pop_back();
        // }
        // // print(it & 1, tourLen[from], tourLen[to]);
        // }
    }

    for (int i=0; i<k; ++i) {
        // Set the time limit to 32 milliseconds
        // chrono::milliseconds timeLimit(32);

        // Run the 3-opt algorithm
        auto tmp = threeOpt(groups[i]);
        vector<int> tour = tmp.second;
        cout << tour.size() << ' ';
        for (int idx: tour) cout << point_idx[groups[i][idx]] << ' ';
        cout << '\n';
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
}
