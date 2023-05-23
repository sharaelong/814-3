#include <bits/stdc++.h>
#ifdef SHARAELONG
#include "debug.hpp"
#endif
using namespace std;
typedef long long ll;
typedef pair<int, int> pii;

void solve(int idx) {
    int k = 140;
    vector<vector<int>> tour(k);
    for (int i=0; i<k; ++i) {
        int cnt;
        cin >> cnt;
        for (int j=0; j<cnt; ++j) {
            int x;
            cin >> x;
            tour[i].push_back(x);
        }
    }

    string s1 = "data/data";
    string s2 = to_string(idx);
    string s3 = ".in";
    string s = s1+s2+s3;
    freopen(s.c_str(), "r", stdin);
    
    int n;
    cin >> n >> k;
    vector<pair<ll, ll>> points(n+1);
    for (int i=1; i<=n; ++i) {
        int x, y;
        cin >> x >> y;
        points[i] = { x, y };
    }

    double score = 0; 
    for (auto v: tour) {
        double tmp = 0;
        v.push_back(v[0]);
        for (int i=0; i+1<v.size(); ++i) {
            ll dx = points[v[i]].first - points[v[i+1]].first;
            ll dy = points[v[i]].second - points[v[i+1]].second;
            tmp += sqrt(dx*dx + dy*dy);
        }
        score = max(score, tmp);
    }
    cout << score << ' ';
}

int main(int argc, char **argv) {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    solve(atoi(argv[1]));
}
