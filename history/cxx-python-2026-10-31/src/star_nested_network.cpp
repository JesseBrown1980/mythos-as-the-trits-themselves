// Mythos blue star as the trits of the nested ternary network.
// The gate network, collision LUTs, streaming and undo are the repository's own
// code (asolaria_ternary_fluid_v1.cpp), included unchanged except the lattice
// is sized to the star: 64 x 64 x 3 (z = slices A, B, C).
// Every cell's 27 channels hold the star seen through the D3Q27 stencil:
//   n_k(x) = star(x + c_k)        -> 27 shifted copies of the star, one per channel.
// Outer network: the 27 channel masses M_k = sum_x n_k(x) are the 27 shell modes.
// Integer only. No floating point anywhere.
#include GATES_INC
#include <cstdio>
static State star_state(const std::vector<std::uint8_t>& star,const array<Velocity,K>& v){
    State f(CELLS*K,0);
    for(int z=0;z<NZ;++z)for(int y=0;y<NY;++y)for(int x=0;x<NX;++x)
      for(int k=0;k<K;++k){
        int sx=mod(x+v[k].x,NX), sy=mod(y+v[k].y,NY), sz=mod(z+v[k].z,NZ);
        f[index(x,y,z)*K+k]=star[(sz*NY+sy)*NX+sx];
      }
    return f;
}
static array<i64,K> modes(const State& f){
    array<i64,K> m{};
    for(int i=0;i<CELLS;++i)for(int k=0;k<K;++k)m[k]+=f[i*K+k];
    return m;
}
static bool same(const Totals&a,const Totals&b){
    return a.m==b.m&&a.px==b.px&&a.py==b.py&&a.pz==b.pz&&a.energy2==b.energy2;
}
// density image rho(x,y,z) and the star-match count: cells whose rest channel equals the star trit
static void dump_rho(const State& f,FILE* fp){
    for(int i=0;i<CELLS;++i){int r=0;for(int k=0;k<K;++k)r+=f[i*K+k];std::uint8_t b=(std::uint8_t)r;fwrite(&b,1,1,fp);}
}
static i64 aligned_channels(const State& f,const std::vector<std::uint8_t>& star){
    // count of (cell,channel) slots equal to the star at that cell: 27*CELLS when every copy sits on the star
    i64 n=0;for(int i=0;i<CELLS;++i)for(int k=0;k<K;++k)n+=(f[i*K+k]==star[i]);return n;
}
int main(){
    try{
    (void)&initial; (void)&simulate; (void)&repo_main; (void)&verify_reverse;
    const auto v=velocities();
    const auto c=make_collisions(v);
    std::vector<std::uint8_t> star(CELLS);
    FILE* in=fopen("star.trits","rb"); if(!in||fread(star.data(),1,CELLS,in)!=(size_t)CELLS)throw std::runtime_error("star.trits");
    fclose(in);
    i64 starsum=0;for(auto t:star)starsum+=t;
    const State start=star_state(star,v);
    const Totals T0=totals(start,v);
    std::cout<<"LATTICE|dims="<<NX<<"x"<<NY<<"x"<<NZ<<"|cells="<<CELLS<<"|trits="<<CELLS*K
             <<"|star_copies=27|collision_gates="<<c.size()<<"|floating_point=0\n";
    std::cout<<"STAR|trit_sum_one_copy="<<starsum<<"|mass_all_copies="<<T0.m<<"|27x_check="<<(T0.m==27*starsum)
             <<"|px="<<T0.px<<"|py="<<T0.py<<"|pz="<<T0.pz<<"|energy2="<<T0.energy2<<"\n";
    FILE* rho=fopen("rho_frames.bin","wb");
    FILE* csv=fopen("outer_modes.csv","w");
    fprintf(csv,"mode,tick");for(int k=0;k<K;++k)fprintf(csv,",M%d",k);fprintf(csv,",Q,aligned\n");
    const int TMAX=192;
    for(int pass=0;pass<2;++pass){
        bool coll=pass==1;
        State cur=start,nxt(CELLS*K,0);
        for(int t=0;t<=TMAX;++t){
            Totals now=totals(cur,v);
            if(!same(now,T0))throw std::runtime_error("conservation violated");
            auto m=modes(cur); i64 Q=0;for(auto x:m)Q+=x;
            i64 al=aligned_channels(cur,star);
            fprintf(csv,"%s,%d",coll?"gates":"free",t);for(auto x:m)fprintf(csv,",%lld",(long long)x);
            fprintf(csv,",%lld,%lld\n",(long long)Q,(long long)al);
            if(t==0||t==1||t==2||t==12||t==72||t==96||t==191||t==192)
                std::cout<<(coll?"GATES":"FREE ")<<"|t="<<t<<"|Q="<<Q<<"|aligned_slots="<<al<<"/"<<(i64)CELLS*K
                         <<"|state==start="<<(cur==start)<<"\n";
            if(coll&&(t%4==0))dump_rho(cur,rho);
            if(t==TMAX)break;
            if(coll)apply_collisions(cur,c,t);
            stream(cur,nxt,v);cur.swap(nxt);
        }
    }
    fclose(rho);fclose(csv);
    // reversibility: 72 forward, 72 back, every one of the 331,776 trits restored
    {
        State cur=start,nxt(CELLS*K,0);
        for(int t=0;t<72;++t){apply_collisions(cur,c,t);stream(cur,nxt,v);cur.swap(nxt);}
        bool moved=cur!=start; i64 diff=0;for(int i=0;i<CELLS*K;++i)diff+=(cur[i]!=start[i]);
        for(int t=71;t>=0;--t){unstream(cur,nxt,v);cur.swap(nxt);undo_collisions(cur,c,t);}
        std::cout<<"REVERSE|steps=72|slots_differing_at_t72="<<diff<<"|moved="<<moved
                 <<"|every_star_trit_restored="<<(cur==start)<<"\n";
        if(cur!=start)throw std::runtime_error("reverse failed");
    }
    // null control: uniform all-ones lattice is a fixed point of gates + streaming
    {
        State u(CELLS*K,1),cur=u,nxt(CELLS*K,0);
        for(int t=0;t<72;++t){apply_collisions(cur,c,t);stream(cur,nxt,v);cur.swap(nxt);}
        std::cout<<"NULL|uniform_state_unchanged_after_72="<<(cur==u)<<"\n";
    }
    std::cout<<"RESULT|mass_momentum_energy_exact_every_tick=PASS\n";
    return 0;
    }catch(const std::exception&e){std::cerr<<"FAIL|"<<e.what()<<"\n";return 1;}
}
