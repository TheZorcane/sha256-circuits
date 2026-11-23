// linear-layer optimiser deps.c: global low-weight GF(2)-linear dependencies among the signals of an XAG.
// Every signal is an affine form over {1, PIs, AND outputs}; we hash forms (mod constant) with a
// random GF(2)-linear 128-bit hash (seeded), so h(a^b) = h(a)^h(b).  Finds
//   D3: all triples {x,y,z}, x^y^z = 0 (mod const)        [exhaustive over all pairs]
//   D4: all quadruples {x,y,p,q} with x^y^p^q = 0 such that no pair sum is an existing form
//       ("non-trivial": not the sum of two triples)           [exhaustive, bucketed passes]
// usage: deps FILE.xag seed [d4passes(0..256)] > out
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
typedef struct { uint64_t lo, hi; } H;
static uint64_t s[2];
static uint64_t rnd(void){ uint64_t a=s[0],b=s[1]; s[0]=b; a^=a<<23; s[1]=a^b^(a>>17)^(b>>26); return s[1]+b; }
static int n_in, n_nodes, nv; static char *op; static int *fa, *fb; static H *h;
#define TB 20
static int *tab; // open addressing on lo
static inline int look(uint64_t lo, uint64_t hi){
  uint64_t i = (lo*0x9E3779B97F4A7C15ull)>>(64-TB);
  for(;;){ int v=tab[i]; if(v<0) return -1; if(h[v].lo==lo && h[v].hi==hi) return v; i=(i+1)&((1u<<TB)-1);} }
static void ins(int v){ uint64_t i=(h[v].lo*0x9E3779B97F4A7C15ull)>>(64-TB);
  for(;;){ if(tab[i]<0){tab[i]=v;return;} if(h[tab[i]].lo==h[v].lo&&h[tab[i]].hi==h[v].hi){fprintf(stderr,"dup form %d %d\n",tab[i],v);return;} i=(i+1)&((1u<<TB)-1);} }
typedef struct { uint64_t k; int x,y; } P;
static int cmpP(const void*a,const void*b){ uint64_t x=((const P*)a)->k,y=((const P*)b)->k; return x<y?-1:x>y; }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); s[0]=strtoull(argv[2],0,10)*0x9E3779B97F4A7C15ull+1; s[1]=0xD1B54A32D192ED03ull;
  int passes = argc>3?atoi(argv[3]):0; int pstart = argc>4?atoi(argv[4]):0; int skipd3 = argc>5;
  char w[64]; int ver; fscanf(f,"%s %d",w,&ver); fscanf(f,"%s %d",w,&n_in); fscanf(f,"%s %d",w,&n_nodes);
  nv=n_in+1+n_nodes; op=malloc(nv); fa=malloc(4*nv); fb=malloc(4*nv); h=calloc(nv,sizeof(H));
  for(int i=0;i<20;i++) rnd();
  for(int i=1;i<=n_in;i++){ op[i]='I'; h[i].lo=rnd(); h[i].hi=rnd(); }
  op[0]='C';
  for(int j=0;j<n_nodes;j++){ int v=n_in+1+j; int a,b; fscanf(f,"%s %d %d",w,&a,&b); op[v]=w[0]; fa[v]=a>>1; fb[v]=b>>1;
    if(w[0]=='A'){ h[v].lo=rnd(); h[v].hi=rnd(); } else { h[v].lo=h[a>>1].lo^h[b>>1].lo; h[v].hi=h[a>>1].hi^h[b>>1].hi; } }
  tab=malloc(sizeof(int)<<TB); for(int i=0;i<(1<<TB);i++) tab[i]=-1;
  for(int v=1;v<nv;v++) ins(v);
  // D3
  long nt=0, ngate=0;
  if(!skipd3) for(int x=1;x<nv;x++) for(int y=x+1;y<nv;y++){
    int z=look(h[x].lo^h[y].lo,h[x].hi^h[y].hi); if(z>y){
      nt++;
      int g=0; // is it a gate relation?
      if(op[z]=='X' && ((fa[z]==x&&fb[z]==y)||(fa[z]==y&&fb[z]==x))) g=1;
      if(op[y]=='X' && ((fa[y]==x&&fb[y]==z)||(fa[y]==z&&fb[y]==x))) g=1;
      if(g) ngate++; else printf("T %d %d %d\n",x,y,z);
    } }
  fprintf(stderr,"D3: triples %ld, gate relations %ld, extra %ld\n",nt,ngate,nt-ngate);
  if(!passes) return 0;
  // D4
  int *bk[256]; int bn[256]; memset(bn,0,sizeof bn);
  for(int v=1;v<nv;v++) bn[h[v].hi&255]++;
  for(int k=0;k<256;k++){ bk[k]=malloc(4*(bn[k]+1)); bn[k]=0; }
  for(int v=1;v<nv;v++){ int k=h[v].hi&255; bk[k][bn[k]++]=v; }
  size_t cap=(size_t)nv*nv/256/2*12/10+1000000; P*pp=malloc(cap*sizeof(P));
  long nq=0, nnon=0;
  for(int k=pstart;k<passes;k++){
    size_t m=0;
    for(int x=1;x<nv;x++){ int b=(h[x].hi&255)^k; int *L=bk[b]; int c=bn[b];
      for(int i=0;i<c;i++){ int y=L[i]; if(y>x){ if(m>=cap){fprintf(stderr,"cap\n");exit(1);} pp[m].k=h[x].lo^h[y].lo; pp[m].x=x; pp[m].y=y; m++; } } }
    qsort(pp,m,sizeof(P),cmpP);
    for(size_t i=0;i<m;){ size_t j=i+1; while(j<m&&pp[j].k==pp[i].k) j++;
      if(j-i>1){
        for(size_t a=i;a<j;a++) for(size_t b2=a+1;b2<j;b2++){
          int x=pp[a].x,y=pp[a].y,p=pp[b2].x,q=pp[b2].y;
          if((h[x].hi^h[y].hi)!=(h[p].hi^h[q].hi)) continue;
          if(x==p||x==q||y==p||y==q) continue;
          nq++;
          // each quad appears for 3 pairings (in possibly different passes); canonical: x is the min of all four
          if(!(x<p&&x<q) ) { if(!(p<x&&p<y)) continue; int t; t=x;x=p;p=t; t=y;y=q;q=t; }
          // count once: require this pairing to be the one where min's partner is smallest
          if(!(y<p&&y<q)) continue;
          int z1=look(h[x].lo^h[y].lo,h[x].hi^h[y].hi), z2=look(h[x].lo^h[p].lo,h[x].hi^h[p].hi), z3=look(h[x].lo^h[q].lo,h[x].hi^h[q].hi);
          if(z1<0&&z2<0&&z3<0){ nnon++; printf("Q %d %d %d %d\n",x,y,p,q); }
        } }
      i=j; }
    if(k%16==15) fprintf(stderr,"pass %d: pairings %ld nontrivial quads %ld\n",k,nq,nnon);
  }
  fprintf(stderr,"D4: pairings %ld nontrivial quads %ld\n",nq,nnon);
  return 0;
}
