#include <ostream>

template<class T, size_t... I>
std::ostream& print_tup(std::ostream& out, const T& tup, std::index_sequence<I...>)
{
    out << "(";
    (..., (out << (I == 0? "" : " ") << std::get<I>(tup)));
    out << ")";
    return out;
}

template <class... T>
std::ostream& operator <<(std::ostream& out, const std::tuple<T...>& tup) {
    print_tup(out, tup, std::make_index_sequence<sizeof...(T)>());
    return out;
}

template <typename T, typename S>
std::ostream& operator <<(std::ostream& out, const std::pair<T, S>& p) {
    out << "(" << p.first << " " << p.second << ")";
    return out;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, const std::vector<T>& v) {
    out << "[";
    for (int i = 0; i < (int)v.size(); i++) {
        out << v[i];
        if (i != (int)v.size()-1) out << " ";
    }
    out << "]";
    return out;
}

template <typename T, size_t size>
std::ostream& operator <<(std::ostream& out, const std::array<T,size>& v) {
    out << '[';
    for (int i = 0; i < (int)v.size(); i++) {
        out << v[i];
        if (i != (int)v.size()-1) out << ' ';
    }
    out << ']';
    return out;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, const std::vector<std::vector<T> >& m) {
    out << std::endl;
    for (int i=0; i<m.size(); ++i) {
        out << m[i] << std::endl;
    }
    return out;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, const std::set<T>& s) {
    out << '{';
    for (T x: s) {
        out << x << ' ';
    }
    out << '}';
    return out;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, const std::multiset<T>& s) {
    out << '{';
    for (T x: s) {
        out << x << ' ';
    }
    out << '}';
    return out;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, std::queue<T> q) {
    out << '[';
    while (!q.empty()) {
        out << q.front();
        q.pop();
        if (!q.empty()) out << ' ';
    }
    out << ']';
    return out;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, std::stack<T> st) {
    std::vector<T> v;
    while (!st.empty()) {
        v.push_back(st.top());
        st.pop();
    }
    
    out << '[';
    while (!v.empty()) {
        out << v.back();
        v.pop_back();
        if (!v.empty()) out << ' ';
    }
    out << ']';
    return out;
}

#define print(...) show(cout, #__VA_ARGS__, __VA_ARGS__)

template<typename H1>
std::ostream& show(std::ostream& out, const char* label, H1&& value) {
    return out << label << " = " << value << std::endl;
}

template<typename H1, typename ...T>
std::ostream& show(std::ostream& out, const char* label, H1&& value, T&&... rest) {
    const char* first_comma = strchr(label, ',');
    const char* left_parenthesis = strchr(label, '(');
    if (left_parenthesis != nullptr && left_parenthesis < first_comma) {
        const char* right_parenthesis = strchr(left_parenthesis, ')');
        assert(right_parenthesis != nullptr);
        const char* pcomma = strchr(right_parenthesis, ',');
        return show(out.write(label, pcomma - label) << " = " << value << ',', pcomma + 1, rest...);
    }
    return show(out.write(label, first_comma - label) << " = " << value << ',', first_comma + 1, rest...);
}
