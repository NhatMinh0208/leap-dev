#ifndef CPL_TEMPLATE
#define CPL_TEMPLATE
/*
	Normie's Template v2.6
	Changes:
	Added range
*/
// Standard library in one include.
#include <bits/stdc++.h>
using namespace std;
 
// ordered_set library.
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
#define ordered_set(el) tree<el,null_type,less<el>,rb_tree_tag,tree_order_statistics_node_update>
 
// AtCoder library. (Comment out these two lines if you're not submitting in AtCoder.) (Or if you want to use it in other judges, run expander.py first.)
//#include <atcoder/all>
//using namespace atcoder;

//Pragmas (Comment out these three lines if you're submitting in szkopul or USACO.)
#pragma comment(linker, "/stack:200000000")
#pragma GCC optimize("Ofast,unroll-loops,tree-vectorize")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native")
 
//File I/O.
#define FILE_IN "cseq.inp"
#define FILE_OUT "cseq.out"
#define ofile freopen(FILE_IN,"r",stdin);freopen(FILE_OUT,"w",stdout)
 
//Fast I/O.
#define fio ios::sync_with_stdio(0);cin.tie(0)
#define nfio cin.tie(0)
#define endl "\n"
 
//Order checking.
#define ord(a,b,c) ((a>=b)and(b>=c))
 
//min/max redefines, so i dont have to resolve annoying compile errors.
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

// Fast min/max assigns to use with AVX.
// Requires g++ 9.2.0.
// template<typename T>
// __attribute__((always_inline)) void chkmin(T& a, const T& b) {
//     a=(a<b)?a:b;
// }

// template<typename T>
// __attribute__((always_inline)) void chkmax(T& a, const T& b) {
//     a=(a>b)?a:b;
// }
 
//Constants.
#define MOD (ll(998244353))
#define MAX 300001
#define mag 320
const long double PI=3.14159265358979;
 
//Pairs and 3-pairs.
#define p1 first
#define p2 second.first
#define p3 second.second
#define fi first
#define se second
#define pii(element_type) pair<element_type,element_type>
#define piii(element_type) pair<element_type,pii(element_type)>

//Custom begin/end shorthand.
#define rnge(rnge) rnge.begin(),rnge.end()
 
//Quick power of 2.
#define pow2(x) (ll(1)<<x)
 
//Short for-loops.
#define ff(i,__,___) for(int i=__;i<=___;i++)
#define rr(i,__,___) for(int i=__;i>=___;i--)
 
//Typedefs.
#define bi BigInt
typedef long long ll;
typedef long double ld;
typedef short sh;

// Binpow and stuff
ll BOW(ll a, ll x, ll p)
{
	if (!x) return 1;
	ll res=BOW(a,x/2,p);
	res*=res;
	res%=p;
	if (x%2) res*=a;
	return res%p;
}
ll INV(ll a, ll p)
{
	return BOW(a,p-2,p);
}
//---------END-------//
#endif

void raw_hc_sample(int n, int seed, int* h, vector<int>* j,  vector<int>* weight, int* var) {
    srand(seed);
    cout<<"raw run"<<endl;
    cout<<seed<<endl;

    int delta[2048];
    int local_h[2048],local_var[2048];
    vector<int> local_j[2048],local_w[2048];
    ordered_set(int) pos;

    for (int i=0;i<n;i++) {
        local_h[i] = h[i];
        local_j[i] = j[i];
        local_w[i] = weight[i];
        // cout<<"variable "<<i<<endl;
        // cout<<"h = "<<local_h[i]<<endl;
        // cout<<"edges "<<i<<endl;
        for (int k=0;k<local_j[i].size();k++) {
            if (i<local_j[i][k]) {
                cout<<i<<' '<<local_j[i][k]<<' '<<local_w[i][k]<<endl;
            }
        }
    }

    for (int i=0;i<n;i++) {
        local_var[i]=(rand() % 2 * 2) - 1;
    }
    
    for (int i=0;i<n;i++) {
        // local_var[i]=(rand() % 2 * 2) - 1;
        delta[i] = local_h[i] * local_var[i];
        for (int k=0; k<local_j[i].size(); k++) {
            delta[i] += local_var[i] * local_var[local_j[i][k]] * local_w[i][k];
        }
        if (delta[i] > 0) {
            pos.insert(i);
        }
        cout<<"init "<<i<<' '<<local_var[i]<<' '<<delta[i]<<endl;
    }

    while (!pos.empty())
    {
        cout<<"iteration: "<<endl;
        for (auto g : pos) cout<<g<<' ';
        cout<<endl;
        int a = rand() % pos.size();
        // cout<<a<<endl;
        auto it = pos.find_by_order(a);
        int tar = (*it);
        cout<<"select "<<tar<<endl;
        pos.erase(it);
        delta[tar] = -delta[tar];
        local_var[tar] = -local_var[tar];
        for (int i=0;i<local_j[tar].size();i++) {
            delta[local_j[tar][i]] += 2 * (local_var[tar]) * (local_var[local_j[tar][i]]) * (local_w[tar][i]);
            if (delta[local_j[tar][i]] <= 0) {
                pos.erase(local_j[tar][i]);
            }
            else {
                pos.insert(local_j[tar][i]);
            }
        }
        
        for (int i=0;i<n;i++) {
            cout<<"after "<<i<<' '<<local_var[i]<<' '<<delta[i]<<endl;
        }
    }
    

    // for (int i=0;i<n;i++) {
    //     // local_var[i]=(rand() % 2 * 2) - 1;
    //     delta[i] = local_h[i] * local_var[i];
    //     for (int k=0; k<local_j[i].size(); k++) {
    //         delta[i] += local_var[i] * local_var[local_j[i][k]] * local_w[i][k];
    //     }
    //     if (delta[i] > 0) cout<<"wtf";
    //     assert(delta[i] <= 0);
    // }

    for (int i=0;i<n;i++) var[i]=local_var[i];

}

int main() {
    string s;
    int n,opt;
    int h[100],res[100];
    vector<int> j[100];
    vector<int> w[100];
    int i,k,a,b,c;
    cin>>n;
    for (i=0;i<n;i++) {
        cin>>h[i];
    }
    while (true) {
        cin>>a;
        if (a!=-1) {
            cin>>b>>c;
        
        j[a].push_back(b);
        w[a].push_back(c);
        
        j[b].push_back(a);
        w[b].push_back(c);
        }
        else {
            break;
        }
    }
    cin>>opt;
    raw_hc_sample(n,123,h,j,w,res);
    for (i=0;i<n;i++) cout<<i<<' '<<res[i]<<endl;
}