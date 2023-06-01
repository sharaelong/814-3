#include <bits/stdc++.h>
#ifdef SHARAELONG
#include "debug.hpp"
#endif
using namespace std;
typedef long long ll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;

bool ge(const ll& a, const ll& b) { return a >= b; }
bool le(const ll& a, const ll& b) { return a <= b; }
bool eq(const ll& a, const ll& b) { return a == b; }
bool gt(const ll& a, const ll& b) { return a > b; }
bool lt(const ll& a, const ll& b) { return a < b; }
int sgn(const ll& a) { return a >= 0 ? a ? 1 : 0 : -1; }

struct pt {
    ll x, y;
    pt() { }
    pt(ll _x, ll _y) : x(_x), y(_y) { }
    pt operator-(const pt& p) const {
        return pt(x - p.x, y - p.y);
    }
    ll cross(const pt& p) const {
        return x * p.y - y * p.x;
    }
    ll cross(const pt& a, const pt& b) const {
        return (a - *this).cross(b - *this);
    }
    ll dot(const pt& p) const {
        return x * p.x + y * p.y;
    }
    ll dot(const pt& a, const pt& b) const {
        return (a - *this).dot(b - *this);
    }
    ll sqrLength() const {
        return this->dot(*this);
    }
    bool operator==(const pt& p) const {
        return eq(x, p.x) && eq(y, p.y);
    }
    bool operator<(const pt& p) const {
        return lt(x, p.x) || lt(y, p.y);
    }
};

const pt inf_pt = pt(1e18, 1e18);

struct QuadEdge {
    pt origin;
    QuadEdge* rot = nullptr;
    QuadEdge* onext = nullptr;
    bool used = false;
    QuadEdge* rev() const {
        return rot->rot;
    }
    QuadEdge* lnext() const {
        return rot->rev()->onext->rot;
    }
    QuadEdge* oprev() const {
        return rot->onext->rot;
    }
    pt dest() const {
        return rev()->origin;
    }
};

QuadEdge* make_edge(pt from, pt to) {
    QuadEdge* e1 = new QuadEdge;
    QuadEdge* e2 = new QuadEdge;
    QuadEdge* e3 = new QuadEdge;
    QuadEdge* e4 = new QuadEdge;
    e1->origin = from;
    e2->origin = to;
    e3->origin = e4->origin = inf_pt;
    e1->rot = e3;
    e2->rot = e4;
    e3->rot = e2;
    e4->rot = e1;
    e1->onext = e1;
    e2->onext = e2;
    e3->onext = e4;
    e4->onext = e3;
    return e1;
}

void splice(QuadEdge* a, QuadEdge* b) {
    swap(a->onext->rot->onext, b->onext->rot->onext);
    swap(a->onext, b->onext);
}

void delete_edge(QuadEdge* e) {
    splice(e, e->oprev());
    splice(e->rev(), e->rev()->oprev());
    delete e->rev()->rot;
    delete e->rev();
    delete e->rot;
    delete e;
}

QuadEdge* connect(QuadEdge* a, QuadEdge* b) {
    QuadEdge* e = make_edge(a->dest(), b->origin);
    splice(e, a->lnext());
    splice(e->rev(), b);
    return e;
}

bool left_of(pt p, QuadEdge* e) {
    return gt(p.cross(e->origin, e->dest()), 0);
}

bool right_of(pt p, QuadEdge* e) {
    return lt(p.cross(e->origin, e->dest()), 0);
}

template <class T>
T det3(T a1, T a2, T a3, T b1, T b2, T b3, T c1, T c2, T c3) {
    return a1 * (b2 * c3 - c2 * b3) - a2 * (b1 * c3 - c1 * b3) +
           a3 * (b1 * c2 - c1 * b2);
}

bool in_circle(pt a, pt b, pt c, pt d) {
    __int128 det = -det3<__int128>(b.x, b.y, b.sqrLength(), c.x, c.y,
                                   c.sqrLength(), d.x, d.y, d.sqrLength());
    det += det3<__int128>(a.x, a.y, a.sqrLength(), c.x, c.y, c.sqrLength(), d.x,
                          d.y, d.sqrLength());
    det -= det3<__int128>(a.x, a.y, a.sqrLength(), b.x, b.y, b.sqrLength(), d.x,
                          d.y, d.sqrLength());
    det += det3<__int128>(a.x, a.y, a.sqrLength(), b.x, b.y, b.sqrLength(), c.x,
                          c.y, c.sqrLength());
    return det > 0;
}

pair<QuadEdge*, QuadEdge*> build_tr(int l, int r, vector<pt>& p) {
    if (r - l + 1 == 2) {
        QuadEdge* res = make_edge(p[l], p[r]);
        return make_pair(res, res->rev());
    }
    if (r - l + 1 == 3) {
        QuadEdge *a = make_edge(p[l], p[l + 1]), *b = make_edge(p[l + 1], p[r]);
        splice(a->rev(), b);
        int sg = sgn(p[l].cross(p[l + 1], p[r]));
        if (sg == 0)
            return make_pair(a, b->rev());
        QuadEdge* c = connect(b, a);
        if (sg == 1)
            return make_pair(a, b->rev());
        else
            return make_pair(c->rev(), c);
    }
    int mid = (l + r) / 2;
    QuadEdge *ldo, *ldi, *rdo, *rdi;
    tie(ldo, ldi) = build_tr(l, mid, p);
    tie(rdi, rdo) = build_tr(mid + 1, r, p);
    while (true) {
        if (left_of(rdi->origin, ldi)) {
            ldi = ldi->lnext();
            continue;
        }
        if (right_of(ldi->origin, rdi)) {
            rdi = rdi->rev()->onext;
            continue;
        }
        break;
    }
    QuadEdge* basel = connect(rdi->rev(), ldi);
    auto valid = [&basel](QuadEdge* e) { return right_of(e->dest(), basel); };
    if (ldi->origin == ldo->origin)
        ldo = basel->rev();
    if (rdi->origin == rdo->origin)
        rdo = basel;
    while (true) {
        QuadEdge* lcand = basel->rev()->onext;
        if (valid(lcand)) {
            while (in_circle(basel->dest(), basel->origin, lcand->dest(),
                             lcand->onext->dest())) {
                QuadEdge* t = lcand->onext;
                delete_edge(lcand);
                lcand = t;
            }
        }
        QuadEdge* rcand = basel->oprev();
        if (valid(rcand)) {
            while (in_circle(basel->dest(), basel->origin, rcand->dest(),
                             rcand->oprev()->dest())) {
                QuadEdge* t = rcand->oprev();
                delete_edge(rcand);
                rcand = t;
            }
        }
        if (!valid(lcand) && !valid(rcand))
            break;
        if (!valid(lcand) ||
            (valid(rcand) && in_circle(lcand->dest(), lcand->origin,
                                       rcand->origin, rcand->dest())))
            basel = connect(rcand, basel->rev());
        else
            basel = connect(basel->rev(), lcand->rev());
    }    return make_pair(ldo, rdo);
}

vector<tuple<pt, pt, pt>> delaunay(vector<pt> p) {
    sort(p.begin(), p.end(), [](const pt& a, const pt& b) {
        return lt(a.x, b.x) || (eq(a.x, b.x) && lt(a.y, b.y));
    });
    auto res = build_tr(0, (int)p.size() - 1, p);
    QuadEdge* e = res.first;
    vector<QuadEdge*> edges = {e};
    while (lt(e->onext->dest().cross(e->dest(), e->origin), 0))
        e = e->onext;
    auto add = [&p, &e, &edges]() {
        QuadEdge* curr = e;
        do {
            curr->used = true;
            p.push_back(curr->origin);
            edges.push_back(curr->rev());
            curr = curr->lnext();
        } while (curr != e);
    };
    add();
    p.clear();
    int kek = 0;
    while (kek < (int)edges.size()) {
        if (!(e = edges[kek++])->used)
            add();
    }
    vector<tuple<pt, pt, pt>> ans;
    for (int i = 0; i < (int)p.size(); i += 3) {
        ans.push_back(make_tuple(p[i], p[i + 1], p[i + 2]));
    }
    return ans;
}

bool isInTriangle(const pt& x, const pt& a, const pt& b, const pt& c) {
    ll sign1 = x.cross(a, b);
    ll sign2 = x.cross(b, c);
    ll sign3 = x.cross(c, a);

    // Check if the signs are all positive or negative
    return (sign1 >= 0 && sign2 >= 0 && sign3 >= 0) || (sign1 <= 0 && sign2 <= 0 && sign3 <= 0);
}

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

struct DisjointSet {
    int n;
    vector<int> parent, rank;
    DisjointSet(int _n) : n(_n) {
        parent.resize(n);
        iota(parent.begin(), parent.end(), 0);
        rank.resize(n, 0);
    }

    int find(int u) {
        return parent[u] = (u == parent[u] ? u : find(parent[u]));
    }

    void merge(int u, int v) {
        u = find(u); v = find(v);
        if (u == v) return;
        if (rank[u] > rank[v]) swap(u, v);
        parent[u] = v;
        if (rank[u] == rank[v]) ++rank[v];
    }
};

const int MAX_N = 8000;

bool adjmat[MAX_N][MAX_N];

bool isConnected(const vector<int>& gp, const vector<vector<int>>& adj) {
    DisjointSet dsu(gp.size());
    int comp_cnt = gp.size();
    for (int i=0; i<gp.size(); ++i) {
        for (int there: adj[gp[i]]) {
            int j = find(gp.begin(), gp.end(), there) - gp.begin();
            if (j < gp.size() && dsu.find(i) != dsu.find(j)) {
                dsu.merge(i, j);
                comp_cnt--;
            }
        }
    }
    return comp_cnt == 1;
}

double mstCost(const vector<int>& gp, const vector<pt>& pts) {
    vector<tuple<ll, int, int>> edges;
    for (int i=0; i<gp.size(); ++i) {
        for (int j=i+1; j<gp.size(); ++j) {
            auto tmp = pts[gp[i]] - pts[gp[j]];
            edges.push_back({ tmp.sqrLength(), i, j });
        }
    }
    
    sort(edges.begin(), edges.end());
    
    ll ret = 0;
    DisjointSet dsu(gp.size());
    for (const auto&[d, i, j]: edges) {
        if (dsu.find(i) != dsu.find(j)) {
            dsu.merge(i, j);
            ret += d;
        }
    }
    
    return sqrt(ret);
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

    vector<pt> pt_wrap(n);
    for (int i=0; i<n; ++i) {
        pt_wrap[i] = { points[i].first, points[i].second };
    }

    auto triangles = delaunay(pt_wrap);

    vector<vector<int>> adj(n);
    for (int i=0; i<triangles.size(); ++i) {
        const auto&[a, b, c] = triangles[i];
        int ai = point_idx[{ a.x, a.y }];
        int bi = point_idx[{ b.x, b.y }];
        int ci = point_idx[{ c.x, c.y }];
        if ((a-b).sqrLength() < 150000000) {
            adj[ai].push_back(bi);
            adj[bi].push_back(ai);
            adjmat[ai][bi] = adjmat[bi][ai] = true;
        }
        if ((b-c).sqrLength() < 150000000) {
            adj[bi].push_back(ci);
            adj[ci].push_back(bi);
            adjmat[bi][ci] = adjmat[ci][bi] = true;
        }
        if ((c-a).sqrLength() < 150000000) {
            adj[ci].push_back(ai);
            adj[ai].push_back(ci);
            adjmat[ci][ai] = adjmat[ai][ci] = true;
        }
    }
    
    vector<double> circle_area(n);
    for (int i=0; i<n; ++i) {
        ll mn = LLONG_MAX;
        for (int j: adj[i]) {
            pt d = pt_wrap[i] - pt_wrap[j];
            mn = min(mn, d.sqrLength());
        }
        circle_area[i] = log2(mn) - 5;
    }

    const int MAX_COORD = 814000;
    int dx = (MAX_COORD + 10) / 10;
    int dy = (MAX_COORD + 14) / 14;

    vector<vector<int>> groups;
    vector<int> belong(n);
    vector<double> weight;
    for (int i=0; i<10; ++i) {
        for (int j=0; j<14; ++j) {
            vector<int> local_points;
            // double w_sum = 0;
            for (int l=0; l<n; ++l) {
                auto[x, y] = points[l];
                if (dx * i <= x && x < dx * (i+1) && dy * j <= y && y < dy * (j+1)) {
                    local_points.push_back(l);
                    belong[l] = groups.size();
                    // w_sum += circle_area[l];
                }
            }
            groups.push_back(local_points);
            weight.push_back(mstCost(local_points, pt_wrap));
            // cout << weight.back() << '\n';
        }
    }
    // print(weight);

    for (int it=2; it<=4000; ++it) {
        random_device rd;
        mt19937 gen(rd());

        vector<int> p(k);
        iota(p.begin(), p.end(), 0);
        sort(p.begin(), p.end(), [&](int a, int b) {
            return weight[a] < weight[b];
        });

        const double ratio = 0.707;
        
        vector<double> prob(k, 1000);
        for (int i=1; i<prob.size(); ++i) {
            prob[i] = prob[i-1] * ratio;
        }
        for (int i=prob.size()-1; i>prob.size()/2; --i) {
            prob[i] = prob[prob.size()-1-i];
        }
        
        discrete_distribution<int> distribution(prob.begin(), prob.end());
        int from = p[distribution(gen)];

        vector<int> neigh_idx;
        for (int here: groups[from]) {
            for (int there: adj[here]) {
                if (belong[there] != from) {
                    neigh_idx.push_back(belong[there]);
                }
            }
        }
        if (neigh_idx.empty()) continue;
        
        sort(neigh_idx.begin(), neigh_idx.end());
        neigh_idx.erase(unique(neigh_idx.begin(), neigh_idx.end()), neigh_idx.end());

        uniform_int_distribution<int> distribution2(0, neigh_idx.size()-1);

        int to = neigh_idx[distribution2(gen)];

        // print(from, weight[from], to, weight[to]);

        bool found = false;
        for (int here: groups[from]) {
            for (int there: adj[here]) {
                if (belong[there] == to) {
                    if (weight[from] < weight[to]) { // if true, why they named from, to? lol
                        groups[to].erase(remove(groups[to].begin(), groups[to].end(), there), groups[to].end());
                        if (!isConnected(groups[to], adj)) {
                            groups[to].push_back(there);
                            continue;
                        }
                        
                        double updated = mstCost(groups[from], pt_wrap);
                        if (updated < max(weight[from], weight[to]) + (50000 / log2(it))) {
                            groups[from].push_back(there);
                            belong[there] = from;
                            weight[to] = updated;
                            weight[from] = mstCost(groups[from], pt_wrap);
                            // print(weight[from], weight[to]);
                            found = true;
                            break;
                        } else {
                            groups[to].push_back(there);
                        }
                    } else if (weight[from] > weight[to]) {
                        groups[from].erase(remove(groups[from].begin(), groups[from].end(), here), groups[from].end());
                        if (!isConnected(groups[from], adj)) {
                            groups[from].push_back(here);
                            continue;
                        }
                        
                        double updated = mstCost(groups[to], pt_wrap);
                        if (updated < max(weight[from], weight[to]) + (50000 / log2(it))) {
                            groups[to].push_back(here);
                            belong[here] = to;
                            weight[from] = updated;
                            weight[to] = mstCost(groups[to], pt_wrap);
                            // print(weight[from], weight[to]);
                            found = true;
                            break;
                        } else {
                            groups[from].push_back(here);
                        }
                    }
                }
            }
            if (found) break;
        }
    }

    for (int i=0; i<k; ++i) {
        // Set the time limit to 32 milliseconds
        // chrono::milliseconds timeLimit(32);

        // Run the 3-opt algorithm
        vector<pll> convert;
        for (int idx: groups[i]) convert.push_back(points[idx]);
        auto tmp = threeOpt(convert);
        vector<int> tour = tmp.second;
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
