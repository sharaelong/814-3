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

    // Initialize the tour to a simple sequential order
    vector<int> tour(n);
    for (int i = 0; i < n; ++i) {
        tour[i] = i;
    }
        
    // Shuffle the vector randomly
    random_device rd;
    mt19937 g(rd());
    shuffle(tour.begin(), tour.end(), g);

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
        if (!improvement) {
            break;
        }
    }

    return { tourLength(points, tour), tour };
}

vector<double> cut(vector<double> bars) {
    int n = bars.size();
    vector<double> ret(n-1);
    double m = accumulate(bars.begin(), bars.end(), 0.0) / n;
    double accum_m = 0;
    int cnt = 0;
    double used_ratio = 0;
    for (int i=0; i<n; ++i) {
        if (cnt == n-1) break;
        if (accum_m + bars[i] * (1.0 - used_ratio) >= m * (cnt+1)) {
            ret[cnt] = (double)i + used_ratio + (m * (cnt+1) - accum_m) / bars[i];
            used_ratio += (m * (cnt+1) - accum_m) / bars[i];
            accum_m = m * (cnt+1);
            cnt++;
            i--;
        } else {
            accum_m += bars[i] * (1.0 - used_ratio);
            used_ratio = 0;
        }
    }
    return ret;
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
        point_idx[points[i]] = i;
    }

    const int MAX_COORD = 814000;
    int dx = (MAX_COORD + 10) / 10;
    int dy = (MAX_COORD + 14) / 14;

    vector<vector<int>> groups;
    for (int i=0; i<10; ++i) {
        for (int j=0; j<14; ++j) {
            vector<int> local_points;
            for (int l=0; l<n; ++l) {
                auto[x, y] = points[l];
                if (dx * i <= x && x < dx * (i+1) && dy * j <= y && y < dy * (j+1)) {
                    local_points.push_back(l);
                }
            }
            groups.push_back(local_points);
        }
    }

    vector<double> tourLen(k, 1e9);
    for (int i=0; i<k; ++i) {
        // Run the 3-opt algorithm
        vector<pll> convert;
        for (int idx: groups[i]) convert.push_back(points[idx]);
        tourLen[i] = min(tourLen[i], threeOpt(convert).first);
    }
    
    vector<double> w_ysum(14, 0);
    for (int j=0; j<14; ++j) {
        for (int i=0; i<10; ++i) w_ysum[j] += tourLen[14 * i + j];
    }
    vector<double> _ycut = cut(w_ysum);
    vector<int> ycut;
    for (double y: _ycut) ycut.push_back(y * dy);
    // print(ycut);

    for (int i=0; i<10; ++i) {
        vector<vector<int>> local_points(14);
        for (int l=0; l<n; ++l) {
            auto[x, y] = points[l];
            if (dx * i <= x && x < dx * (i+1)) {
                int ynum = lower_bound(ycut.begin(), ycut.end(), y) - ycut.begin();
                local_points[ynum].push_back(l);
            }
        }
        for (int j=0; j<14; ++j) {
            groups[14 * i + j] = local_points[j];
        }
    }
    // print(1);

    for (int i=0; i<k; ++i) {
        // Run the 3-opt algorithm
        vector<pll> convert;
        for (int idx: groups[i]) convert.push_back(points[idx]);
        tourLen[i] = threeOpt(convert).first;
    }
    
    for (int j=0; j<14; ++j) {
        vector<double> w_xsum(10);
        for (int i=0; i<10; ++i) w_xsum[i] = tourLen[14 * i + j];
        vector<double> _xcut = cut(w_xsum);
        vector<int> xcut;
        for (double x: _xcut) xcut.push_back(x * dx);
        // print(j, w_xsum, _xcut, xcut);

        vector<int> targets;
        for (int i=0; i<10; ++i) {
            for (int idx: groups[14 * i + j]) {
                targets.push_back(idx);
            }
        }
        vector<vector<int>> local_points(10);
        for (int l: targets) {
            auto[x, y] = points[l];
            int xnum = lower_bound(xcut.begin(), xcut.end(), x) - xcut.begin();
            local_points[xnum].push_back(l);
        }
        for (int i=0; i<10; ++i) {
            groups[14 * i + j] = local_points[i];
        }
    }
 
    for (int i=0; i<k; ++i) {
        // Run the 3-opt algorithm
        vector<pll> convert;
        for (int idx: groups[i]) convert.push_back(points[idx]);

        auto tmp = threeOpt(convert);

        vector<int> tour = tmp.second;
        // cout << tmp.first << endl;
        cout << tour.size() << ' ';
        for (int idx: tour) cout << groups[i][idx]+1 << ' ';
        cout << '\n';
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
}
