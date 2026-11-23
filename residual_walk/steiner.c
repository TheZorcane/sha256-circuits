// linear-layer optimiser steiner.c -- C port of steiner.py (same algorithm; see that file's docstring).
// Minimum XOR-SLP local search inside a fixed universe of forms given by triples x^y^z = const.
// usage: steiner UNIVERSE.xag n_live TRIPLES.txt OUT.state seed walk_iters rounds [INIT.state]
// Output: OUT.state  ("v c d" per alive XOR node) -- turned into an XAG and re-checked exactly by emit.py.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
static int n_in, n_nodes, nv, base;
static char *kind;            // 'P','A','X'
static int *da, *db;          // current definition
static char *alive; static int *pin;
static int **us; static int *un, *ucap;   // users
static int *toff, *tp, *tq;   // triples CSR: partners of v are (tp[i],tq[i]) for i in toff[v]..toff[v+1]
static uint64_t *pm; static uint32_t *pu; static uint32_t uidc = 1;   // position = (major, uid)
static long nx;
static uint64_t rs[2];
static inline uint64_t rnd(void){ uint64_t a=rs[0],b=rs[1]; rs[0]=b; a^=a<<23; rs[1]=a^b^(a>>17)^(b>>26); return rs[1]+b; }
static inline int lt(int u,int w){ return pm[u]<pm[w] || (pm[u]==pm[w] && pu[u]<pu[w]); }
static void uadd(int u,int n){ if(un[u]==ucap[u]){ ucap[u]=ucap[u]?2*ucap[u]:4; us[u]=realloc(us[u],4*ucap[u]); } us[u][un[u]++]=n; }
static void udel(int u,int n){ for(int i=0;i<un[u];i++) if(us[u][i]==n){ us[u][i]=us[u][--un[u]]; return; } }
// ---- scratch
static int *stk, *Bl, *Fl; static char *mark; // mark: 1=B, 2=F
static int cmp_pos(const void*a,const void*b){ int u=*(const int*)a,w=*(const int*)b; return lt(u,w)?-1:lt(w,u)?1:0; }
static uint64_t *slm; static uint32_t *slu;
// make pos[srcs] < pos[v]; return 0 if v in TFI(srcs)
static int ensure_before(int s0,int s1,int v){
  if(s0==v||s1==v) return 0;
  int nb=0, sp=0; int srcs[2]={s0,s1};
  for(int i=0;i<2;i++){ int c=srcs[i]; if(lt(v,c) && !mark[c]){ mark[c]=1; Bl[nb++]=c; stk[sp++]=c; } }
  if(!nb) return 1;
  int ok=1;
  while(sp){ int u=stk[--sp]; if(kind[u]=='P') continue;
    int d2[2]={da[u],db[u]};
    for(int i=0;i<2;i++){ int w=d2[i]; if(w==v){ ok=0; sp=0; break; }
      if(lt(v,w) && !mark[w]){ mark[w]=1; Bl[nb++]=w; stk[sp++]=w; } } }
  if(!ok){ for(int i=0;i<nb;i++) mark[Bl[i]]=0; return 0; }
  int ub=Bl[0]; for(int i=1;i<nb;i++) if(lt(ub,Bl[i])) ub=Bl[i];
  int nf=0; Fl[nf++]=v; mark[v]=2; stk[sp++]=v;
  while(sp){ int u=stk[--sp];
    for(int i=0;i<un[u];i++){ int w=us[u][i]; if(lt(w,ub) && !mark[w]){ mark[w]=2; Fl[nf++]=w; stk[sp++]=w; } } }
  qsort(Bl,nb,4,cmp_pos); qsort(Fl,nf,4,cmp_pos);
  // merge slots
  int i=0,j=0,k=0;
  while(i<nb||j<nf){ int u; if(j>=nf||(i<nb&&lt(Bl[i],Fl[j]))) u=Bl[i++]; else u=Fl[j++]; slm[k]=pm[u]; slu[k]=pu[u]; k++; }
  k=0; for(i=0;i<nb;i++,k++){ pm[Bl[i]]=slm[k]; pu[Bl[i]]=slu[k]; mark[Bl[i]]=0; }
  for(j=0;j<nf;j++,k++){ pm[Fl[j]]=slm[k]; pu[Fl[j]]=slu[k]; mark[Fl[j]]=0; }
  return 1;
}
// ---- planning
#define MAXP 16
typedef struct { int n,p,q; } PE;
static int DEPTH=2; static double UPHILL=0; static int UPMAX=1;
static inline int inplan(const PE*pl,int n,int w){ for(int i=0;i<n;i++) if(pl[i].n==w) return 1; return 0; }
// plan activation of dead node c given visible plan pl[0..*np); banned nodes in ban[0..nban). returns cost or -1
static int plan_activate(int c,int depth,int*ban,int nban,PE*pl,int*np){
  int best=-1; PE bestp[MAXP]; int bestn=0; int n0=*np;
  for(int t=toff[c];t<toff[c+1];t++){ int p=tp[t],q=tq[t]; int bad=0;
    for(int i=0;i<nban;i++) if(ban[i]==p||ban[i]==q) bad=1;
    if(bad) continue;
    int cost=1; int n1=n0; int ok=1; int w2[2]={p,q};
    for(int i=0;i<2;i++){ int w=w2[i];
      if(alive[w]||inplan(pl,n1,w)) continue;
      if(depth<=1||n1+2>=MAXP){ ok=0; break; }
      ban[nban]=c; int r=plan_activate(w,depth-1,ban,nban+1,pl,&n1);
      if(r<0){ ok=0; break; } cost+=r; }
    if(ok && n1<MAXP && (best<0||cost<best)){ best=cost; bestn=n1-n0; memcpy(bestp,pl+n0,sizeof(PE)*(n1-n0)); bestp[bestn].n=c; bestp[bestn].p=p; bestp[bestn].q=q; bestn++; if(cost==1) break; }
  }
  if(best<0) return -1;
  memcpy(pl+n0,bestp,sizeof(PE)*bestn); *np=n0+bestn; return best;
}
static int *delta, *touched, ntouched; static char *dmark; static int *deadl;
static inline void dl(int u,int x){ if(!delta[u]&&x){ } if(!dmark[u]){ dmark[u]=1; touched[ntouched++]=u; } delta[u]+=x; }
static void setdef(int n,int c,int d){ da[n]=c; db[n]=d; uadd(c,n); uadd(d,n); }
// returns gain (>=0) if applied, -1 otherwise
static int try_move(int v,int c,int d,int accept_zero){
  int oa=da[v], ob=db[v];
  if((oa==c&&ob==d)||(oa==d&&ob==c)) return -1;
  PE pl[MAXP]; int np=0; int ban[8]; int cd[2]={c,d};
  for(int i=0;i<2;i++){ int w=cd[i];
    if(kind[w]=='X' && !alive[w] && !inplan(pl,np,w)){ ban[0]=v; if(plan_activate(w,DEPTH,ban,1,pl,&np)<0) return -1; } }
  int added=np;
  ntouched=0;
  for(int i=0;i<np;i++){ dl(pl[i].p,1); dl(pl[i].q,1); }
  dl(c,1); dl(d,1);
  int nd=0, sp=0;
  dl(oa,-1); stk[sp++]=oa; dl(ob,-1); stk[sp++]=ob;
  while(sp){ int w=stk[--sp];
    if(kind[w]!='X'||!alive[w]||dmark[w]==2||inplan(pl,np,w)) continue;
    if(un[w]+pin[w]+delta[w]==0){ dmark[w]=2; deadl[nd++]=w; dl(da[w],-1); if(dmark[da[w]]!=2) stk[sp++]=da[w]; dl(db[w],-1); if(dmark[db[w]]!=2) stk[sp++]=db[w]; } }
  for(int i=0;i<ntouched;i++){ delta[touched[i]]=0; dmark[touched[i]]=0; }
  for(int i=0;i<nd;i++) dmark[deadl[i]]=0;
  int gain=nd-added;
  if(gain<0){ if(!(accept_zero==2 && gain>=-UPMAX && (rnd()>>11)*(1.0/9007199254740992.0) < UPHILL)) return -1; }
  else if(gain==0&&!accept_zero) return -1;
  // commit
  udel(oa,v); udel(ob,v); setdef(v,c,d);
  for(int i=0;i<np;i++){ int n=pl[i].n; alive[n]=1; setdef(n,pl[i].p,pl[i].q);
    uint64_t m=pm[pl[i].p]>pm[pl[i].q]?pm[pl[i].p]:pm[pl[i].q]; pm[n]=m+1+(rnd()&0xffff); pu[n]=uidc++; }
  int ok=1;
  for(int i=0;i<np&&ok;i++) if(!ensure_before(pl[i].p,pl[i].q,pl[i].n)) ok=0;
  if(ok && !ensure_before(c,d,v)) ok=0;
  for(int i=0;i<np&&ok;i++) if(!(lt(pl[i].p,pl[i].n)&&lt(pl[i].q,pl[i].n))) { if(!ensure_before(pl[i].p,pl[i].q,pl[i].n)) ok=0; }
  if(ok){ if(!(lt(c,v)&&lt(d,v))) ok=0;
    for(int i=0;i<np&&ok;i++){ int n=pl[i].n; if(!(lt(pl[i].p,n)&&lt(pl[i].q,n))) ok=0; for(int k=0;k<un[n]&&ok;k++) if(!lt(n,us[n][k])) ok=0; }
    for(int k=0;k<un[v]&&ok;k++) if(!lt(v,us[v][k])) ok=0; }
  if(!ok){
    for(int i=np-1;i>=0;i--){ int n=pl[i].n; udel(pl[i].p,n); udel(pl[i].q,n); alive[n]=0; }
    udel(c,v); udel(d,v); setdef(v,oa,ob);
    return -1; }
  // garbage collect
  sp=0; stk[sp++]=oa; stk[sp++]=ob; int rd=0;
  while(sp){ int w=stk[--sp]; if(kind[w]!='X'||!alive[w]) continue;
    if(un[w]+pin[w]==0){ alive[w]=0; rd++; udel(da[w],w); udel(db[w],w); stk[sp++]=da[w]; stk[sp++]=db[w]; } }
  nx+=added-rd;
  return rd-added<0?0:rd-added;
}
static void check(void){
  long cnt=0;
  for(int v=base;v<nv;v++) if(alive[v]){
    if(!(alive[da[v]]&&alive[db[v]]&&lt(da[v],v)&&lt(db[v],v))){ fprintf(stderr,"CHECK FAIL at %d\n",v); exit(3); }
    if(kind[v]=='X'){ cnt++; if(un[v]+pin[v]==0){ fprintf(stderr,"CHECK dangling %d\n",v); exit(3);} } }
  if(cnt!=nx){ fprintf(stderr,"CHECK count %ld %ld\n",cnt,nx); exit(3); }
}
static void dump(const char*path){
  char tmp[1024]; snprintf(tmp,sizeof tmp,"%s.tmp",path); FILE*f=fopen(tmp,"w");
  for(int v=base;v<nv;v++) if(alive[v]) fprintf(f,"%d %d %d\n",v,da[v],db[v]);
  fclose(f); rename(tmp,path);
}
typedef struct { int a,b,c; } T3;
static int cmpT(const void*x,const void*y){ return memcmp(x,y,sizeof(T3)); }
int main(int argc,char**argv){
  if(argc<8){ fprintf(stderr,"usage\n"); return 1; }
  FILE*f=fopen(argv[1],"r"); int n_live=atoi(argv[2]); uint64_t seed=strtoull(argv[5],0,10); long iters=atol(argv[6]); int rounds=atoi(argv[7]);
  if(getenv("DEPTH")) DEPTH=atoi(getenv("DEPTH"));
  if(getenv("UPHILL")) UPHILL=atof(getenv("UPHILL"));
  if(getenv("UPMAX")) UPMAX=atoi(getenv("UPMAX"));
  rs[0]=seed*0x9E3779B97F4A7C15ull+12345; rs[1]=0xD1B54A32D192ED03ull^seed; for(int i=0;i<20;i++) rnd();
  char w[64]; int ver; fscanf(f,"%s %d",w,&ver); fscanf(f,"%s %d",w,&n_in); fscanf(f,"%s %d",w,&n_nodes);
  base=n_in+1; nv=base+n_nodes;
  kind=malloc(nv); da=calloc(nv,4); db=calloc(nv,4); alive=calloc(nv,1); pin=calloc(nv,4);
  us=calloc(nv,sizeof(int*)); un=calloc(nv,4); ucap=calloc(nv,4); pm=calloc(nv,8); pu=calloc(nv,4);
  stk=malloc(4*(4*nv+16)); Bl=malloc(4*nv); Fl=malloc(4*nv); mark=calloc(nv,1); slm=malloc(8*nv); slu=malloc(4*nv);
  delta=calloc(nv,4); touched=malloc(4*nv); dmark=calloc(nv,1); deadl=malloc(4*nv);
  size_t tcap=1<<20, nt=0; T3*T=malloc(tcap*sizeof(T3));
  #define ADDT(x,y,z) do{ int q3[3]={x,y,z}; for(int i_=0;i_<3;i_++)for(int j_=i_+1;j_<3;j_++) if(q3[j_]<q3[i_]){int t_=q3[i_];q3[i_]=q3[j_];q3[j_]=t_;} \
     if(q3[0]>0&&q3[0]!=q3[1]&&q3[1]!=q3[2]){ if(nt==tcap){tcap*=2;T=realloc(T,tcap*sizeof(T3));} T[nt].a=q3[0];T[nt].b=q3[1];T[nt].c=q3[2];nt++; } }while(0)
  for(int i=0;i<base;i++){ kind[i]='P'; alive[i]=1; pm[i]=(uint64_t)i<<24; }
  for(int j=0;j<n_nodes;j++){ int v=base+j,a,b; fscanf(f,"%s %d %d",w,&a,&b); kind[v]=w[0]=='A'?'A':'X'; da[v]=a>>1; db[v]=b>>1; pm[v]=(uint64_t)v<<24;
    if(kind[v]=='X') ADDT(v,a>>1,b>>1);
    if(j<n_live||kind[v]=='A'){ alive[v]=1; } }
  int no; fscanf(f,"%s %d",w,&no); for(int i=0;i<no;i++){ int o; fscanf(f,"%d",&o); pin[o>>1]++; }
  fclose(f);
  f=fopen(argv[3],"r"); char line[256];
  while(fgets(line,sizeof line,f)){ int x,y,z; if(line[0]=='T'&&sscanf(line+1,"%d %d %d",&x,&y,&z)==3) ADDT(x,y,z); }
  fclose(f);
  qsort(T,nt,sizeof(T3),cmpT); size_t m=0; for(size_t i=0;i<nt;i++) if(!m||cmpT(&T[i],&T[m-1])) T[m++]=T[i]; nt=m;
  // AND planes (optional): "P N f g h" -> AND N may be redefined with operand pairs (f,g),(f,h),(g,h)
  size_t npl=0, plcap=1<<16; T3*PLN=malloc(plcap*sizeof(T3)); int *plN=malloc(4*plcap);
  if(getenv("PLANES")){ f=fopen(getenv("PLANES"),"r");
    while(fgets(line,sizeof line,f)){ int N,x,y,z; if(line[0]=='P'&&sscanf(line+1,"%d %d %d %d",&N,&x,&y,&z)==4){
      if(npl==plcap){ plcap*=2; PLN=realloc(PLN,plcap*sizeof(T3)); plN=realloc(plN,4*plcap); }
      PLN[npl].a=x; PLN[npl].b=y; PLN[npl].c=z; plN[npl]=N; npl++; } }
    fclose(f); }
  toff=calloc(nv+2,4);
  for(size_t i=0;i<nt;i++){ toff[T[i].a+1]++; toff[T[i].b+1]++; toff[T[i].c+1]++; }
  for(size_t i=0;i<npl;i++) toff[plN[i]+1]+=3;
  for(int v=0;v<=nv;v++) toff[v+1]+=toff[v];
  tp=malloc(4*(3*nt+3*npl+1)); tq=malloc(4*(3*nt+3*npl+1)); int*fill=calloc(nv+1,4);
  for(size_t i=0;i<nt;i++){ int q3[3]={T[i].a,T[i].b,T[i].c}; for(int k=0;k<3;k++){ int v=q3[k]; int o=toff[v]+fill[v]++; tp[o]=q3[(k+1)%3]; tq[o]=q3[(k+2)%3]; } }
  for(size_t i=0;i<npl;i++){ int q3[3]={PLN[i].a,PLN[i].b,PLN[i].c}; int v=plN[i]; for(int k=0;k<3;k++){ int o=toff[v]+fill[v]++; tp[o]=q3[k]; tq[o]=q3[(k+1)%3]; } }
  free(T); free(fill);
  if(argc>8){ // initial state: alive XOR set + definitions; positions by topological sort
    for(int v=base;v<nv;v++) if(kind[v]=='X') alive[v]=0;
    f=fopen(argv[8],"r"); int v,c,d; while(fscanf(f,"%d %d %d",&v,&c,&d)==3){ alive[v]=1; da[v]=c; db[v]=d; } fclose(f);
  }
  for(int v=base;v<nv;v++) if(alive[v]){ uadd(da[v],v); uadd(db[v],v); if(kind[v]=='X') nx++; }
  if(argc>8){ // garbage-collect dangling XOR nodes of the initial state
    long g=0; for(int v=nv-1;v>=base;v--){ int sp=0; stk[sp++]=v;
      while(sp){ int w2=stk[--sp]; if(kind[w2]!='X'||!alive[w2]||un[w2]+pin[w2]) continue;
        alive[w2]=0; nx--; g++; udel(da[w2],w2); udel(db[w2],w2); stk[sp++]=da[w2]; stk[sp++]=db[w2]; } }
    fprintf(stderr,"init GC removed %ld\n",g); }
  if(argc>8){ // Kahn
    int *indeg=calloc(nv,4); int qh=0,qt=0; int*q=malloc(4*nv); uint64_t k=1;
    for(int v=0;v<nv;v++) if(alive[v]){ if(v>=base) indeg[v]=2; else q[qt++]=v; }
    while(qh<qt){ int u=q[qh++]; pm[u]=(k++)<<24; for(int i=0;i<un[u];i++){ int x=us[u][i]; if(--indeg[x]==0) q[qt++]=x; } }
    free(indeg); free(q);
  }
  check();
  int ncand=0; int*cand=malloc(4*nv);
  for(int v=base;v<nv;v++) if((kind[v]=='X'&&toff[v+1]-toff[v]>1)||(kind[v]=='A'&&toff[v+1]>toff[v])) cand[ncand++]=v;
  fprintf(stderr,"universe nodes %d triples %zu start XOR %ld cand %d\n",n_nodes,nt,nx,ncand);
  time_t t0=time(0); long hist0=0,hist1=0; long bestnx=nx;
  for(int rd=0;rd<rounds;rd++){
    for(;;){ long g=0;
      for(int i=0;i<ncand;i++){ int v=cand[i]; if(!alive[v]) continue;
        for(int t=toff[v];t<toff[v+1]&&alive[v];t++){ int r=try_move(v,tp[t],tq[t],0); if(r>0){ g+=r; hist1++; } } }
      if(g==0) break; }
    check(); if(nx<=bestnx){ bestnx=nx; dump(argv[4]); }
    fprintf(stderr,"round %d after sweep XOR %ld best %ld zero-moves %ld gain-moves %ld (%lds)\n",rd,nx,bestnx,hist0,hist1,(long)(time(0)-t0)); fflush(stderr);
    if(rd==rounds-1) break;
    long before=nx;
    int nal=0; static int *al=0; if(!al) al=malloc(4*nv);
    for(long it=0;it<iters;it++){
      if(it%200000==0){ nal=0; for(int i=0;i<ncand;i++) if(alive[cand[i]]) al[nal++]=cand[i]; }
      int v=al[rnd()%nal]; if(!alive[v]) continue;
      int t=toff[v]+rnd()%(toff[v+1]-toff[v]); int r=try_move(v,tp[t],tq[t],UPHILL>0?2:1); if(r==0) hist0++; else if(r>0) hist1++; }
    (void)before;
  }
  return 0;
}
