#include <bits/stdc++.h>
#ifdef SHARAELONG
#include "debug.hpp"
#endif
using namespace std;
typedef long long ll;
typedef pair<int, int> pii;

int main() {
    // freopen("test.in", "r", stdin);
    // freopen("data.in", "w", stdout);
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<ll> dist(0, (ll)814001 * 814001 - 1);

    set<pair<ll, ll>> s;
    while (s.size() < 8000) {
        ll n = dist(gen);
        ll x = n / 814001;
        ll y = n % 814001;
        s.insert({ x, y });
    }

    cout << "8000 140\n";
    for (auto[x, y]: s) {
        cout << x << ' ' << y << '\n';
    }
}
