/*
g++ -static -Wl,--stack=268435456 -O2 -std=c++17 c.cpp -o r.exe
g++ -DEVAL -std=gnu++11 -O2 -pipe -static -s -o r c.cpp
*/

const int base = 1000000000; const int base_digits = 9;
struct bigint {
    vector<int> a; int sign;
 
    bigint() :
        sign(1) {
    }
 
    bigint(long long v) {
        *this = v;
    }
 
    bigint(const string &s) {
        read(s);
    }
 
    void operator=(const bigint &v) {
        sign = v.sign;
        a = v.a;
    }
 
    void operator=(long long v) {
        sign = 1;
        if (v < 0)
            sign = -1, v = -v;
        for (; v > 0; v = v / base)
            a.push_back(v % base);
    }
 
    bigint operator+(const bigint &v) const {
        if (sign == v.sign) {
            bigint res = v;
 
            for (int i = 0, carry = 0; i < (int) max(a.size(), v.a.size()) || carry; ++i) {
                if (i == (int) res.a.size())
                    res.a.push_back(0);
                res.a[i] += carry + (i < (int) a.size() ? a[i] : 0);
                carry = res.a[i] >= base;
                if (carry)
                    res.a[i] -= base;
            }
            return res;
        }
        return *this - (-v);
    }
 
    bigint operator-(const bigint &v) const {
        if (sign == v.sign) {
            if (abs() >= v.abs()) {
                bigint res = *this;
                for (int i = 0, carry = 0; i < (int) v.a.size() || carry; ++i) {
                    res.a[i] -= carry + (i < (int) v.a.size() ? v.a[i] : 0);
                    carry = res.a[i] < 0;
                    if (carry)
                        res.a[i] += base;
                }
                res.trim();
                return res;
            }
            return -(v - *this);
        }
        return *this + (-v);
    }
 
    void operator*=(int v) {
        if (v < 0)
            sign = -sign, v = -v;
        for (int i = 0, carry = 0; i < (int) a.size() || carry; ++i) {
            if (i == (int) a.size())
                a.push_back(0);
            long long cur = a[i] * (long long) v + carry;
            carry = (int) (cur / base);
            a[i] = (int) (cur % base);
            //asm("divl %%ecx" : "=a"(carry), "=d"(a[i]) : "A"(cur), "c"(base));
        }
        trim();
    }
 
    bigint operator*(int v) const {
        bigint res = *this;
        res *= v;
        return res;
    }
 
    friend pair<bigint, bigint> divmod(const bigint &a1, const bigint &b1) {
        int norm = base / (b1.a.back() + 1);
        bigint a = a1.abs() * norm;
        bigint b = b1.abs() * norm;
        bigint q, r;
        q.a.resize(a.a.size());
 
        for (int i = a.a.size() - 1; i >= 0; i--) {
            r *= base;
            r += a.a[i];
            int s1 = r.a.size() <= b.a.size() ? 0 : r.a[b.a.size()];
            int s2 = r.a.size() <= b.a.size() - 1 ? 0 : r.a[b.a.size() - 1];
            int d = ((long long) base * s1 + s2) / b.a.back();
            r -= b * d;
            while (r < 0)
                r += b, --d;
            q.a[i] = d;
        }
 
        q.sign = a1.sign * b1.sign;
        r.sign = a1.sign;
        q.trim();
        r.trim();
        return make_pair(q, r / norm);
    }
 
    bigint operator/(const bigint &v) const {
        return divmod(*this, v).first;
    }
 
    bigint operator%(const bigint &v) const {
        return divmod(*this, v).second;
    }
 
    void operator/=(int v) {
        if (v < 0)
            sign = -sign, v = -v;
        for (int i = (int) a.size() - 1, rem = 0; i >= 0; --i) {
            long long cur = a[i] + rem * (long long) base;
            a[i] = (int) (cur / v);
            rem = (int) (cur % v);
        }
        trim();
    }
 
    bigint operator/(int v) const {
        bigint res = *this;
        res /= v;
        return res;
    }
 
    int operator%(int v) const {
        if (v < 0)
            v = -v;
        int m = 0;
        for (int i = a.size() - 1; i >= 0; --i)
            m = (a[i] + m * (long long) base) % v;
        return m * sign;
    }
 
    void operator+=(const bigint &v) {
        *this = *this + v;
    }
    void operator-=(const bigint &v) {
        *this = *this - v;
    }
    void operator*=(const bigint &v) {
        *this = *this * v;
    }
    void operator/=(const bigint &v) {
        *this = *this / v;
    }
 
    bool operator<(const bigint &v) const {
        if (sign != v.sign)
            return sign < v.sign;
        if (a.size() != v.a.size())
            return a.size() * sign < v.a.size() * v.sign;
        for (int i = a.size() - 1; i >= 0; i--)
            if (a[i] != v.a[i])
                return a[i] * sign < v.a[i] * sign;
        return false;
    }
 
    bool operator>(const bigint &v) const {
        return v < *this;
    }
    bool operator<=(const bigint &v) const {
        return !(v < *this);
    }
    bool operator>=(const bigint &v) const {
        return !(*this < v);
    }
    bool operator==(const bigint &v) const {
        return !(*this < v) && !(v < *this);
    }
    bool operator!=(const bigint &v) const {
        return *this < v || v < *this;
    }
 
    void trim() {
        while (!a.empty() && !a.back())
            a.pop_back();
        if (a.empty())
            sign = 1;
    }
 
    bool isZero() const {
        return a.empty() || (a.size() == 1 && !a[0]);
    }
 
    bigint operator-() const {
        bigint res = *this;
        res.sign = -sign;
        return res;
    }
 
    bigint abs() const {
        bigint res = *this;
        res.sign *= res.sign;
        return res;
    }
 
    long long longValue() const {
        long long res = 0;
        for (int i = a.size() - 1; i >= 0; i--)
            res = res * base + a[i];
        return res * sign;
    }
 
    friend bigint gcd(const bigint &a, const bigint &b) {
        return b.isZero() ? a : gcd(b, a % b);
    }
    friend bigint lcm(const bigint &a, const bigint &b) {
        return a / gcd(a, b) * b;
    }
 
    void read(const string &s) {
        sign = 1;
        a.clear();
        int pos = 0;
        while (pos < (int) s.size() && (s[pos] == '-' || s[pos] == '+')) {
            if (s[pos] == '-')
                sign = -sign;
            ++pos;
        }
        for (int i = s.size() - 1; i >= pos; i -= base_digits) {
            int x = 0;
            for (int j = max(pos, i - base_digits + 1); j <= i; j++)
                x = x * 10 + s[j] - '0';
            a.push_back(x);
        }
        trim();
    }
 
    friend istream& operator>>(istream &stream, bigint &v) {
        string s;
        stream >> s;
        v.read(s);
        return stream;
    }
 
    friend ostream& operator<<(ostream &stream, const bigint &v) {
        if (v.sign == -1)
            stream << '-';
        stream << (v.a.empty() ? 0 : v.a.back());
        for (int i = (int) v.a.size() - 2; i >= 0; --i)
            stream << setw(base_digits) << setfill('0') << v.a[i];
        return stream;
    }
 
    static vector<int> convert_base(const vector<int> &a, int old_digits, int new_digits) {
        vector<long long> p(max(old_digits, new_digits) + 1);
        p[0] = 1;
        for (int i = 1; i < (int) p.size(); i++)
            p[i] = p[i - 1] * 10;
        vector<int> res;
        long long cur = 0;
        int cur_digits = 0;
        for (int i = 0; i < (int) a.size(); i++) {
            cur += a[i] * p[cur_digits];
            cur_digits += old_digits;
            while (cur_digits >= new_digits) {
                res.push_back(int(cur % p[new_digits]));
                cur /= p[new_digits];
                cur_digits -= new_digits;
            }
        }
        res.push_back((int) cur);
        while (!res.empty() && !res.back())
            res.pop_back();
        return res;
    }
 
    typedef vector<long long> vll;
 
    static vll karatsubaMultiply(const vll &a, const vll &b) {
        int n = a.size();
        vll res(n + n);
        if (n <= 32) {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    res[i + j] += a[i] * b[j];
            return res;
        }
 
        int k = n >> 1;
        vll a1(a.begin(), a.begin() + k);
        vll a2(a.begin() + k, a.end());
        vll b1(b.begin(), b.begin() + k);
        vll b2(b.begin() + k, b.end());
 
        vll a1b1 = karatsubaMultiply(a1, b1);
        vll a2b2 = karatsubaMultiply(a2, b2);
 
        for (int i = 0; i < k; i++)
            a2[i] += a1[i];
        for (int i = 0; i < k; i++)
            b2[i] += b1[i];
 
        vll r = karatsubaMultiply(a2, b2);
        for (int i = 0; i < (int) a1b1.size(); i++)
            r[i] -= a1b1[i];
        for (int i = 0; i < (int) a2b2.size(); i++)
            r[i] -= a2b2[i];
 
        for (int i = 0; i < (int) r.size(); i++)
            res[i + k] += r[i];
        for (int i = 0; i < (int) a1b1.size(); i++)
            res[i] += a1b1[i];
        for (int i = 0; i < (int) a2b2.size(); i++)
            res[i + n] += a2b2[i];
        return res;
    }
 
    bigint operator*(const bigint &v) const {
        vector<int> a6 = convert_base(this->a, base_digits, 6);
        vector<int> b6 = convert_base(v.a, base_digits, 6);
        vll a(a6.begin(), a6.end());
        vll b(b6.begin(), b6.end());
        while (a.size() < b.size())
            a.push_back(0);
        while (b.size() < a.size())
            b.push_back(0);
        while (a.size() & (a.size() - 1))
            a.push_back(0), b.push_back(0);
        vll c = karatsubaMultiply(a, b);
        bigint res;
        res.sign = sign * v.sign;
        for (int i = 0, carry = 0; i < (int) c.size(); i++) {
            long long cur = c[i] + carry;
            res.a.push_back((int) (cur % 1000000));
            carry = (int) (cur / 1000000);
        }
        res.a = convert_base(res.a, 6, base_digits);
        res.trim();
        return res;
    }
};

//-----------------------------------------------//
//                  Bitset int                   //
//-----------------------------------------------//
const int nbit = 50;
struct bigh{
    vector<ull> num;
    bigh(){
        num = vector<ull>(nbit, 0);
    }
    ull operator[](int &index){
        return num[index];
    }
    void operator=(bigh &ax){
        num = ax.num;
    }
    void operator=(unsigned int &ax){
        num[0] = ax;
    }
    void operator=(ull &ax){
        num[0] = ax;
    }
    bigh operator+(const bigh &ax) const {
        bigh cx = *this;
        bool carry = false;
        for(int i = 0; i<nbit; i++){
            carry = (cx[i]&ax[i]) || (cx[i]&carry) || (ax[i]&carry);
            cx[i] = ax[i]^cx[i]^carry;
        }
        return cx;
    }
    bigh operator-(const bigh &ax) const {
        bool borrow = false;
        bigh cx = *this;
        for (int i = 0; i < nbit; i++) {
            if (borrow) {
                borrow = !cx[i] || (cx[i] && ax[i]);
                cx[i] = !(cx[i] ^ ax[i]);
                
            } else {
                borrow = !cx[i] && ax[i];
                cx[i] = cx[i] ^ ax[i];    
            }
        }
        return cx;
    }
    bigh operator*(const bigh &bx) const {
        bigh cx, ax = *this;
        for(int i = 0; i<nbit; i++){
            if(ax[i]) cx = cx + (bx << i);
        }
        return cx;
    }
};

//-----------------------------------------------//
//                 ijtoi, itoij                  //
//-----------------------------------------------//
int ijtoi(int i, int j, int mx){
  return (i - 1)*mx + j;
}
pii itoij(int i, int mx){
  return {(i+mx-1)/mx, ((i-1)%mx) + 1};
}


//-----------------------------------------------//
//                    Bitset                     //
//-----------------------------------------------//
const int nbit = 10005;
bitset<nbit> operator+(bitset<nbit> ax, bitset<nbit> bx){
  bitset<nbit> cx;
  bool carry = false;
  for(int i = 0; i<nbit; i++){
    cx[i] = ax[i]^bx[i]^carry;
    carry = (ax[i]&bx[i]) || (ax[i]&carry) || (bx[i]&carry);
  }
  return cx;
}
bitset<nbit> operator-(bitset<nbit> ax, bitset<nbit> bx){
    bool borrow = false;
    bitset<nbit> cx;
    for (int i = 0; i < nbit; i++) {
      if (borrow) {
          cx[i] = !(ax[i] ^ bx[i]);
          borrow = !ax[i] || (ax[i] && bx[i]);
      }
      else {
          cx[i] = ax[i] ^ bx[i];
          borrow = !ax[i] && bx[i];
      }
    }
    return cx;
}
bitset<nbit> operator*(bitset<nbit> ax, bitset<nbit> bx){
  bitset<nbit> cx;
  for(int i = 0; i<nbit; i++){
    if(ax[i]) cx = cx + (bx << i);
  }
  return cx;
}

//-----------------------------------------------//
//                    Matrix                     //
//-----------------------------------------------//
#define ele_maxtrix ll
struct Matrix{
  int n, m;
  vector<vector<ele_maxtrix>> ma;
  Matrix(){}
  Matrix(int _n){
    n = _n; m = _n;
    ele_maxtrix ei;
    ei.set(0,0);
    ma = vector<vector<ele_maxtrix>>(n, vector<ele_maxtrix>(m, ei));
    for(int i = 0; i<n; i++){
      ma[i][i].set(0, 1);
    }
  }
  Matrix(int _n, int _m){
    n = _n; m = _m;
    ele_maxtrix ei;
    ei.set(0,0);
    ma = vector<vector<ele_maxtrix>>(n, vector<ele_maxtrix>(m, ei));
  }
};
Matrix operator*(Matrix ax, Matrix bx){
  if(ax.m != bx.n) return Matrix(0,0);
  Matrix cx = Matrix(ax.n, bx.m);
  for(int i = 0; i<ax.n; i++){
    for(int j = 0; j<bx.m; j++){
      for(int c = 0; c<ax.m; c++){
        cx.ma[i][j] = cx.ma[i][j] + (ax.ma[i][c]*bx.ma[c][j]);
      }
    }
  }
  return cx;
}
void print_matrix(Matrix ax){
  for(int i = 0; i<ax.n; i++){
    for(int j = 0; j<ax.m; j++){
      cout << ax.ma[i][j] << " ";
    }
    cout << "\n";
  }
}

//-----------------------------------------------//
//                  Factorial                    //
//-----------------------------------------------//
struct fac{
    int n; ll mod_fac; vector<ll> f;
    fac(){}
    fac(ll _mod_fac){ n = 0; f = vector<ll>(1); f[0] = 1; mod_fac = _mod_fac; }
    fac(ll _mod_fac, int _n){ n = _n; f = vector<ll>(n + 1, 0); f[0] = 1; mod_fac = _mod_fac; }
    ll get(int index){ if(index <= n) return f[index]; for(int i = n + 1; i<=index; i++){ f.push_back(f.back()*i%mod_fac); n++; } return f[index]; }
} ft;

//-----------------------------------------------//
//                  Math prox                    //
//-----------------------------------------------//
ll pow_hota(ll ax, ll bx){
    if(bx == 0) return 1ll;
    if(bx == 1) return ax%mod;
    ll cx = pow_hota(ax, bx/2);
    cx = (cx*cx)%mod;
    if(bx%2) cx = (cx*ax)%mod;
    return cx;
}

ll Cx(ll _n, ll _k){
    return (ft.get(_n)*pow_hota(ft.get(_k)*ft.get(_n - _k)%mod, mod-2))%mod;
}

//-----------------------------------------------//
//                      DSU                      //
//-----------------------------------------------//
struct hota_DSU{
    int n;
    vector<int> dsu;
    hota_DSU(){}
    hota_DSU(int _n){
        n = _n;
        dsu = vector<int>(n + 1, -1);
    }
    int root(int u){
        if(dsu[u] < 0) return u;
        return dsu[u] = root(dsu[u]);
    }
    bool samer(int u, int v){
        int ru = root(u), rv = root(v);
        return (ru == rv);
    }
    void unce(int u, int v){
        int ru = root(u), rv = root(v);
        if(ru != rv){
            if(dsu[ru] > dsu[rv]) swap(ru, rv);
            dsu[ru] += dsu[rv];
            dsu[rv] = ru;
        }
    }
};


//-----------------------------------------------//
//                  Trie bitmask                 //
//-----------------------------------------------//
struct node_trie{
    int numOfWords, data;
    int val;
    int childs[2];
    bool exit_key(int key){
        return (childs[key] != -1);
    }
    node_trie(){ val = 0; numOfWords = 0; data = 0; childs[0] = -1; childs[1] = -1; }
};
int base = 30;
struct hota_trie{
    int n;
    vector<int> list_empty_nodes;
    vector<node_trie> Trie;
    hota_trie(){
        n = 0;
        list_empty_nodes.clear();
        Trie = vector<node_trie>(1, node_trie());
    }
    hota_trie(int _n){
        n = _n;
        list_empty_nodes.clear();
        for(int i = n; i>=1; i--) list_empty_nodes.push_back(i);
        Trie = vector<node_trie>(n + 1, node_trie());
    }
    int get_empty_node(){
        if(list_empty_nodes.empty()){
            n++;
            Trie.push_back(node_trie());
            return n;
        }
        int id = list_empty_nodes.back(); list_empty_nodes.pop_back();
        return id;
    }
    void add_empty_node(int id){
        list_empty_nodes.push_back(id);
    }
    int get_xor(int key, int i){
        return ((key >> i) & 1);
    }
    bool isEmpty(int id){
        return (Trie[id].numOfWords == 0);
    }
    void insert(int key, int idkey){
        int id = 0;
        vector<int> id_list(0);
        Trie[id].numOfWords++;
        for (int i = 0; i < base; i++){
            id_list.push_back(id);
            int truk = get_xor(key, base - 1 - i);
            if (!Trie[id].exit_key(truk)){
                int fdf = (int)get_empty_node();
                Trie[id].childs[truk] = fdf;
            }
            id = Trie[id].childs[truk];
            Trie[id].numOfWords++;
            Trie[id].data = truk;
        }
        Trie[id].val = idkey;
    }
    void update(vi &id_list){
        int id = id_list.back(); id_list.pop_back();
        bool canDelete = isEmpty(id);
        while(!id_list.empty()){
            int old_id = id;
            char oldKey = Trie[id].data;
            id = id_list.back(); id_list.pop_back();
            if(canDelete){
                Trie[id].childs[oldKey] = -1;
                add_empty_node(old_id);
            }
            Trie[id].numOfWords--;
            canDelete = isEmpty(id);
        }
    }
    void erase(int key){ 
        int id = 0;
        vector<int> id_list(0);
        for (int i = 0; i < base; i++){
            id_list.push_back(id);
            int truk = get_xor(key, base - 1 - i);
            if (!Trie[id].exit_key(truk)) return;
            id = Trie[id].childs[truk];
        }
        if(Trie[id].numOfWords == 0) return;
        id_list.push_back(id);
        Trie[id].numOfWords--;
        update(id_list);
    }
    pii get_val(int p){
        int id = 0;
        int ans = 0;
        for(int i = 0; i<base; i++){
            int bp = get_xor(p, base - 1 - i);
            //Continue
            if(Trie[id].exit_key(bp)){
                id = Trie[id].childs[bp];
            } else if(Trie[id].exit_key(1^bp)){
                ans |= (1 << (base - 1 - i));
                id = Trie[id].childs[1^bp];
            } else break;
        }
        return {Trie[id].val, ans};
    }
} se;

//-----------------------------------------------//
//                  Trie prox                    //
//-----------------------------------------------//
struct node_trie{
    int id, numOfWords, isEndOfWord;
    char data;
    map<char, int> childs;
    bool exit_key(char key){
        return (childs.find(key) != childs.end());
    }
    node_trie(){ numOfWords = 0; isEndOfWord = 0; data = '\0'; childs.clear(); }
};
struct hota_trie{
    int n;
    vector<int> list_empty_nodes;
    vector<node_trie> Trie;
    hota_trie(){
        n = 0;
        list_empty_nodes.clear();
        Trie = vector<node_trie>(1, node_trie());
    }
    hota_trie(int _n){
        n = _n;
        list_empty_nodes.clear();
        for(int i = n; i>=1; i--) list_empty_nodes.push_back(i);
        Trie = vector<node_trie>(n + 1, node_trie());
    }
    int get_empty_node(){
        if(list_empty_nodes.empty()){
            n++;
            Trie.push_back(node_trie());
            return n;
        }
        int id = list_empty_nodes.back(); list_empty_nodes.pop_back();
        return id;
    }
    void add_empty_node(int id){
        list_empty_nodes.push_back(id);
    }
    bool isEmpty(int id){
        return (Trie[id].isEndOfWord == 0 && Trie[id].childs.empty());
    }
    void updateTrie(vector<int> &id_list, int delta){
        int id = id_list.back(); id_list.pop_back();
        bool canDelete = isEmpty(id);
        while(!id_list.empty()){
            int old_id = id;
            char oldKey = Trie[id].data;
            id = id_list.back(); id_list.pop_back();
            if(canDelete){
                Trie[id].childs.erase(oldKey);
                add_empty_node(old_id);
            }
            Trie[id].numOfWords += delta;
            canDelete = isEmpty(id);
        }
    }
    void insert(string key){
        int id = 0;
        vector<int> id_list(0);
        for (int i = 0; i < key.size(); i++){
            id_list.push_back(id);
            if (!Trie[id].exit_key(key[i]))
                Trie[id].childs[key[i]] = get_empty_node();
            id = Trie[id].childs[key[i]];
            Trie[id].data = key[i];
        }
        id_list.push_back(id);
        Trie[id].numOfWords++;
        Trie[id].isEndOfWord++;
        Trie[id].id = Trie[0].numOfWords + 1;
        updateTrie(id_list, 1);
    }
    int exit_nearest(string key){
        int id = 0;
        for (int i = 0; i < key.size(); i++){
            if (!Trie[id].exit_key(key[i]))
                break;
            id = Trie[id].childs[key[i]];
        }
        return Trie[id].isEndOfWord;
    }
    int count_nearest(string key){
        int id = 0;
        for (int i = 0; i < key.size(); i++){
            if (!Trie[id].exit_key(key[i]))
                break;
            id = Trie[id].childs[key[i]];
        }
        return Trie[id].numOfWords;
    }
    int exit(string key){
        int id = 0;
        for (int i = 0; i < key.size(); i++){
            if (!Trie[id].exit_key(key[i]))
                return 0;
            id = Trie[id].childs[key[i]];
        }
        return Trie[id].isEndOfWord;
    }
    int search(string key){
        int id = 0;
        for (int i = 0; i < key.size(); i++){
            if (!Trie[id].exit_key(key[i]))
                return 0;
            id = Trie[id].childs[key[i]];
        }
        return Trie[id].id;
    }
    int next_search(int id, char key){
        if (!Trie[id].exit_key(key)) return -1;
        return Trie[id].childs[key];
    }
    bool exit_id(int id){
        return Trie[id].isEndOfWord;
    }
    int count_words(string key){
        int id = 0;
        for (int i = 0; i < key.size(); i++){
            if (!Trie[id].exit_key(key[i]))
                return 0;
            id = Trie[id].childs[key[i]];
        }
        return Trie[id].numOfWords;
    }
    void erase(string key, bool delete_All = false){ 
        int id = 0;
        vector<int> id_list(0);
        for (int i = 0; i < key.size(); i++){
            id_list.push_back(id);
            if (!Trie[id].exit_key(key[i])) return;
            id = Trie[id].childs[key[i]];
        }
        if(Trie[id].isEndOfWord == 0) return;
        id_list.push_back(id);
        Trie[id].isEndOfWord--;
        Trie[id].numOfWords--;
        int delta = -1;
        if(delete_All){
            delta -= Trie[id].isEndOfWord;
            Trie[id].numOfWords -= Trie[id].isEndOfWord;
            Trie[id].isEndOfWord = 0;
        }
        updateTrie(id_list, delta);
    }
} oniichan;

//-----------------------------------------------//
//                     RMQ                       //
//-----------------------------------------------//
int cha[400005][20], h[100005];
void setupRMQ(){
    for (int j = 1; 1 << j <= n; j++)
        for (int i = 1; i <= n; i++)
            if (cha[i][j - 1] != 0)
                cha[i][j] = cha[cha[i][j - 1]][j - 1];
}
int LCA(int u, int v){
  if(h[u] < h[v]) swap(u, v);
  int logbb = log2(h[u]);
  for(int i = logbb; i>=0; i--){
    if(h[u] - (1 << i) >= h[v]) u = cha[u][i];
  }
  if(u == v) return u;
  for(int i = logbb; i>=0; i--){
    if(cha[u][i] != 0 && cha[u][i] != cha[v][i]){
      u = cha[u][i];
      v = cha[v][i];
    }
  }
  return cha[u][0];
}
void dfs(int u, int pa){
  cha[u][0] = pa;
  for(int v: ed[u]){
    if(v!=pa){
      h[v] = h[u] + 1;
      dfs(v, u);
    }
  }
}

//-----------------------------------------------//
//               Binary search                   //
//-----------------------------------------------//
int bsearch(vi &ax, int bx){
    int l = 0; int r = (int)ax.size() - 1; int ans = -1;
    while(l<=r){
        int mid = (l+r)/2;
        if(ax[mid] <= bx){
            ans = mid;
            l = mid + 1;
        } else r = mid - 1;
    }
    return ans;
}

//-----------------------------------------------//
//               Random a tree                   //
//-----------------------------------------------//
int n;
int nrnow = 0, addg = 137, seed = 94387;
int getr(int l, int r){
    if(nrnow < l) nrnow = l;
    nrnow = ((nrnow - l + 1)*seed%r + addg%r)%r + l;
    return nrnow;
}
void process(){
    srand(time(NULL));
    addg = rand() % (mod - 1000 + 1) + 1000;
    seed = rand() % (mod - 100000 + 1) + 100000;

    int tcase = 1;
    //cin >> tcase;
    for(int ttcase = 1; ttcase <= tcase; ttcase++){
        cin >> n;
        //debug(addg, seed);
        int nin = 1;
        for(int i = 2; i<=n; i++){
            cout << getr(1, nin) << " " << (nin + 1) << "\n";
            nin++;
        }
    }
}


//-----------------------------------------------//
//               Find min max size k             //
//-----------------------------------------------//
deque<int> dq;
while ((!dq.empty()) && getdp(i) <= getdp(dq.back())) dq.pop_back();
dq.push_back(i);
while ((!dq.empty()) && dq.front() <= i - m + 1) dq.pop_front();

//-----------------------------------------------//
//                  IT normal                    //
//-----------------------------------------------//
struct IT{
    int n;
    vector<int> it;
    void init(int _n){
        n = _n;
        it = vector<int>(4*n + 1, 0);
    }
    void update(int i, int l, int r, int u, int val){
        if(l > u || r < u) return;
        if(l >= r){
            it[i] = val;
            return;
        }
        int mid = (l+r)/2;
        update(i*2, l, mid, u, val);
        update(i*2 + 1, mid + 1, r, u, val);
        it[i] = min(it[i*2], it[i*2 + 1]);
    }
    int query(int i, int l, int r, int u, int v){
        if(l > v || r < u) return 1000000007;
        if(l >= u && r <= v) return it[i];
        int mid = (l+r)/2;
        return min(query(i*2, l, mid, u, v), query(i*2 + 1, mid + 1, r, u, v));
    }
    void upd(int u, int val){
        update(1,1,n,u,val);
    }
    int quy(int l, int r){
        return query(1, 1, n, l, r);
    }
} se;

//-----------------------------------------------//
//              IT with lazy update              //
//-----------------------------------------------//
struct sem{
    int l, r, w;
} sems[300005];
struct node{
    int val, lazy;
};
struct IT{
    int n;
    vector<node> it;
    void downnode(int i){
        it[i*2].val += it[i].lazy;
        it[i*2 + 1].val += it[i].lazy;
        it[i*2].lazy += it[i].lazy;
        it[i*2 + 1].lazy += it[i].lazy;
        it[i].lazy = 0;
    }
    void init(int _n){
        n = _n;
        it = vector<node>(4*n + 1, {0, 0});
    }
    void update(int i, int l, int r, int u, int v, int val){
        if(l > v || r < u) return;
        if(l >= u && r <= v){
            it[i].val += val;
            it[i].lazy += val;
            return;
        }
        int mid = (l+r)/2;
        downnode(i);
        update(i*2, l, mid, u, v, val);
        update(i*2 + 1, mid + 1, r, u, v, val);
        it[i].val = min(it[i*2].val, it[i*2 + 1].val);
    }
    int query(int i, int l, int r, int u, int v){
        if(l > v || r < u) return 1000000007;
        if(l >= u && r <= v) return it[i].val;
        int mid = (l+r)/2;
        downnode(i);
        return min(query(i*2, l, mid, u, v), query(i*2 + 1, mid + 1, r, u, v));
    }
    void upd(int l, int r, int val){
        update(1,1,n,l,r,val);
    }
    int quy(int l, int r){
        return query(1, 1, n, l, r);
    }
} se;

//-----------------------------------------------//
//                      BIT                      //
//-----------------------------------------------//
struct BIT{
    int n;
    vi bit;
    void setup(int _n){
        n = _n;
        bit = vi(n + 1, 0);
    }
    void update(int i, int deta){
        for(; i<=n; i+=i&-i){
            bit[i] += deta;
        }
    }
    int query(int i){
        int sum = 0;
        for(; i>0; i-=i&-i){
            sum += bit[i];
        }
        return sum;
    }
    int get(int l, int r){
        return query(r) - query(l-1);
    }
} be;

//-----------------------------------------------//
//                  Magic naruto                 //
//-----------------------------------------------//
#pragma once

//#pragma GCC target("avx2")
#pragma GCC target("avx2,avx512f,avx512vl")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#include <immintrin.h>

using m256 = __m256i;

#define ALIGN __attribute__((aligned(64)))

#define SET(x) _mm256_set1_epi32(x)
#define SET64(x) _mm256_set1_epi64x(x)
#define LOAD(p) _mm256_loadu_si256((__m256i*)(p))
#define STORE(p, A) _mm256_storeu_si256((__m256i*)(p), A)

#define AND(a, b) _mm256_and_si256(a, b)
#define OR(a, b) _mm256_or_si256(a, b)
#define XOR(a, b) _mm256_xor_si256(a, b)

#define ADD(a, b) _mm256_add_epi32(a, b)
#define SUB(a, b) _mm256_sub_epi32(a, b)
#define CMP(a, b) _mm256_cmpgt_epi32(a, b)

#define GETMOD(a, MOD) SUB(a, AND(CMP(a, MOD), MOD))
#define MADD(a, b, MOD) GETMOD(ADD(a, b), MOD)
#define MSUB(a, b, MOD) GETMOD(SUB(ADD(a, MOD), b), MOD)

#define SETLO(a) _mm256_shuffle_epi32(a, 0xA0)
#define SETHI(a) _mm256_shuffle_epi32(a, 0xF9)
#define CAST64(a) AND(a, SET64(0xFFFFFFFF))
#define ADD64(a, b) _mm256_add_epi64(a, b)

//-----------------------------------------------//
//                   Debug magic                 //
//-----------------------------------------------//
void __print(int x) {cout << x;}
void __print(long x) {cout << x;}
void __print(long long x) {cout << x;}
void __print(unsigned x) {cout << x;}
void __print(unsigned long x) {cout << x;}
void __print(unsigned long long x) {cout << x;}
void __print(float x) {cout << x;}
void __print(double x) {cout << x;}
void __print(long double x) {cout << x;}
void __print(char x) {cout << '\'' << x << '\'';}
void __print(const char *x) {cout << '\"' << x << '\"';}
void __print(const string &x) {cout << '\"' << x << '\"';}
void __print(bool x) {cout << (x ? "true" : "false");}
 
template<typename T, typename V>
void __print(const pair<T, V> &x) {cout << '{'; __print(x.first); cout << ','; __print(x.second); cout << '}';}
template<typename T>
void __print(const T &x) {int f = 0; cout << '{'; for (auto &i: x) cout << (f++ ? "," : ""), __print(i); cout << "}";}
void _print() {cout << "]\n";}
template <typename T, typename... V>
void _print(T t, V... v) {__print(t); if (sizeof...(v)) cout << ", "; _print(v...);}
#ifndef ONLINE_JUDGE
#define debug(x...) cout << "["; _print(x)
#else
#define debug(x...)
#endif
//--------------- Debug end --------------

//-----------------------------------------------//
//              Theme                            //
//-----------------------------------------------//
//Hoang Son WIBU lolicon codeforces rate 1681 khong cay
#include <bits/stdc++.h>
#define F first
#define S second
#define times double stime = clock();
#define gettime cout << "\nTime executed: " << (clock() - stime) / CLOCKS_PER_SEC * 1000 << "ms.\n";
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
typedef double dou;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;

const ll mod = 1000000007ll;
const int ooii = 0x7f7f7f7f;
const ll ooll = 0x7f7f7f7f7f7f7f7f;
const double epsilon = 1e-8;
//1000000007ll 998244353ll
bool debug_on = true, debug_test = false;

//1 <= n <= 5*10^4; 1 <= Tmax <= 10^6
void process(){
    int tcase = 1;
    cin >> tcase;
    for(int ttcase = 1; ttcase <= tcase; ttcase++){
    }
}
int online = 0;
void printcopyright();
int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);
    if (online == 0) {
        freopen("in.inp", "r", stdin);
        freopen("out.out", "w", stdout);
    } else if (online == 1) {
        freopen("dddd.inp", "r", stdin);
        freopen("dddd.out", "w", stdout);
    }
    // times
    process();
    printcopyright();
    // gettime
    return 0;
}
void printcopyright() {
    cerr << "\n";
    cerr << "//    _____                   _                       _\n//   / ____|                 | |       "
            "              (_)\n//  | (___   ___  _ __     __| | ___ _ __    ______ _ _\n//   \\___ \\ / _ "
            "\\| '_ \\   / _` |/ _ \\ '_ \\  |_  / _` | |\n//   ____) | (_) | | | | | (_| |  __/ |_) |  / / "
            "(_| | |\n//  |_____/ \\___/|_| |_|  \\__,_|\\___| .__/  /___\\__,_|_|\n//                       "
            "           | |\n//                                  |_|";
    cerr << "\n\n";
}

//-----------------------------------------------//
//                 Hash sha256                   //
//-----------------------------------------------//
#include "D:\c\ctf\challenges.reply.com\2021\weak 1\A1\sha256.h"
#include <cstring>
#include <fstream>
const unsigned int SHA256::sha256_k[64]={0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};void SHA256::transform(const unsigned char*message,unsigned int block_nb) {uint32 w[64];uint32 wv[8];uint32 t1,t2;const unsigned char*sub_block;int i;int j;for(i=0;i<(int)block_nb;i++){sub_block=message+(i<<6);for(j=0;j<16;j++){SHA2_PACK32(&sub_block[j<<2],&w[j]);} for(j=16;j<64;j++){w[j]=SHA256_F4(w[j-2])+w[j-7]+SHA256_F3(w[j-15])+w[j-16];} for(j=0;j<8;j++){wv[j]=m_h[j];} for(j=0;j<64;j++){t1=wv[7]+SHA256_F2(wv[4])+SHA2_CH(wv[4],wv[5],wv[6]) +sha256_k[j]+w[j];t2=SHA256_F1(wv[0])+SHA2_MAJ(wv[0],wv[1],wv[2]);wv[7]=wv[6];wv[6]=wv[5];wv[5]=wv[4];wv[4]=wv[3]+t1;wv[3]=wv[2];wv[2]=wv[1];wv[1]=wv[0];wv[0]=t1+t2;} for(j=0;j<8;j++){m_h[j]+=wv[j];}}} void SHA256::init() {m_h[0]=0x6a09e667;m_h[1]=0xbb67ae85;m_h[2]=0x3c6ef372;m_h[3]=0xa54ff53a;m_h[4]=0x510e527f;m_h[5]=0x9b05688c;m_h[6]=0x1f83d9ab;m_h[7]=0x5be0cd19;m_len=0;m_tot_len=0;} void SHA256::update(const unsigned char*message,unsigned int len) {unsigned int block_nb;unsigned int new_len,rem_len,tmp_len;const unsigned char*shifted_message;tmp_len=SHA224_256_BLOCK_SIZE-m_len;rem_len=len<tmp_len?len:tmp_len;memcpy(&m_block[m_len],message,rem_len);if(m_len+len<SHA224_256_BLOCK_SIZE){m_len+=len;return;} new_len=len-rem_len;block_nb=new_len / SHA224_256_BLOCK_SIZE;shifted_message=message+rem_len;transform(m_block,1);transform(shifted_message,block_nb);rem_len=new_len%SHA224_256_BLOCK_SIZE;memcpy(m_block,&shifted_message[block_nb<<6],rem_len);m_len=rem_len;m_tot_len+=(block_nb+1)<<6;} void SHA256::final(unsigned char*digest) {unsigned int block_nb;unsigned int pm_len;unsigned int len_b;int i;block_nb=(1+((SHA224_256_BLOCK_SIZE-9)<(m_len%SHA224_256_BLOCK_SIZE)));len_b=(m_tot_len+m_len)<<3;pm_len=block_nb<<6;memset(m_block+m_len,0,pm_len-m_len);m_block[m_len]=0x80;SHA2_UNPACK32(len_b,m_block+pm_len-4);transform(m_block,block_nb);for(i=0;i<8;i++){SHA2_UNPACK32(m_h[i],&digest[i<<2]);}} std::string sha256(std::string input) {unsigned char digest[SHA256::DIGEST_SIZE];memset(digest,0,SHA256::DIGEST_SIZE);SHA256 ctx=SHA256();ctx.init();ctx.update((unsigned char*)input.c_str(),input.length());ctx.final(digest);char buf[2*SHA256::DIGEST_SIZE+1];buf[2*SHA256::DIGEST_SIZE]=0;for(int i=0;i<SHA256::DIGEST_SIZE;i++) sprintf(buf+i*2,"%02x",digest[i]);return std::string(buf);}
