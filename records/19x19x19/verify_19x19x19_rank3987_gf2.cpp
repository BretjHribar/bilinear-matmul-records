#include <fstream>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
using Big=std::vector<uint64_t>;
struct Term { Big u,v,w; std::string us,vs,ws; };
static Big decimal(const std::string&s){Big a(1);for(char ch:s){if(ch<'0'||ch>'9')throw 1;unsigned __int128 carry=ch-'0';for(auto&x:a){carry+=(unsigned __int128)x*10;x=(uint64_t)carry;carry>>=64;}if(carry)a.push_back((uint64_t)carry);}while(a.size()>1&&!a.back())a.pop_back();return a;}
static std::vector<int> bits(const Big&x) { std::vector<int> r;for(size_t k=0;k<x.size();k++){uint64_t z=x[k];while(z){int b=__builtin_ctzll(z);r.push_back(int(k*64+b));z&=z-1;}}return r; }
static bool zero(const Big&x){return x.size()==1&&x[0]==0;}
static int top(const Big&x){auto b=bits(x);return b.empty()?-1:b.back();}
int main(int argc,char**argv){
  if(argc!=2){std::cerr<<"usage: verify_19x19x19_rank3987_gf2 FILE\n";return 2;}
  int n=19,m=19,p=19;
  std::ifstream in(argv[1]); if(!in){std::cerr<<"FAIL: open\n";return 1;}
  std::vector<Term> ts; std::string line; long declared=-1;
  while(std::getline(in,line)){if(line.empty()||line[0]=='#')continue;std::istringstream s(line);std::string a,b,c,d;s>>a;
    if(a=="R")s>>a>>b>>c;else if(!(s>>b)){if(declared<0)declared=std::stol(a);continue;}else s>>c;
    if(a.empty()||b.empty()||c.empty()||(s>>d)){std::cerr<<"FAIL: parse\n";return 1;}try{ts.push_back({decimal(a),decimal(b),decimal(c),a,b,c});}catch(...){std::cerr<<"FAIL: decimal\n";return 1;}
  }
  if(declared>=0&&declared!=(long)ts.size()){std::cerr<<"FAIL: rank header\n";return 1;}
  const int ub=n*m,vb=m*p,wb=n*p;std::unordered_set<uint64_t> r;uint64_t toggles=0;
  std::unordered_set<std::string> seen;
  for(auto&t:ts){if(zero(t.u)||zero(t.v)||zero(t.w)||top(t.u)>=ub||top(t.v)>=vb||top(t.w)>=wb){std::cerr<<"FAIL: factor bounds\n";return 1;}
    std::string key=t.us+" "+t.vs+" "+t.ws;if(!seen.insert(key).second){std::cerr<<"FAIL: duplicate\n";return 1;}
    auto us=bits(t.u),vs=bits(t.v),ws=bits(t.w);toggles+=uint64_t(us.size())*vs.size()*ws.size();
    for(int u:us)for(int v:vs)for(int w:ws){uint64_t q=(uint64_t(u)*vb+v)*wb+w;if(!r.erase(q))r.insert(q);}
  }
  for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int k=0;k<p;k++){uint64_t q=(uint64_t(i*m+j)*vb+j*p+k)*wb+i*p+k;if(!r.erase(q))r.insert(q);}
  std::cout<<(r.empty()?"PASS":"FAIL")<<" dims=<"<<n<<","<<m<<","<<p<<"> rank="<<ts.size()<<" sparse_toggles="<<toggles<<" target_coefficients="<<n*m*p<<" residual="<<r.size()<<"\n";return r.empty()?0:1;
}
