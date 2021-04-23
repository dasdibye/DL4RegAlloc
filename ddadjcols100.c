// KMB 2006 Aug 22
// make gnp_chromatic_number && time nice ./gnp_chromatic_number > gnp_chi.dat
// make gnp_chromatic_number && time nice ./gnp_chromatic_number | graph -Tx -C -g3 -x 0 16 -y 0 7 -m0 -S 16

#include <math.h>
#include "vn_graph.h"
#include <stdlib.h>

//#define GRAPH_MAX_SIZE 50
#define GRAPH_MAX_SIZE 100
#define min(a,b) ((a)<(b)?(a):(b))

void
dump_node_attr_in_list(graph_t g, unsigned int n, FILE* fp)
{
    //fprintf(fp,"dump all edges for nodes = %d and edges = %d\n",n,e);
    int **adjM;
    adjM = (int **)calloc(n,sizeof(int *));
    for (unsigned int i=0; i < n ; i++) 
         adjM[i] = (int *)calloc(n, sizeof(int));

    for ( unsigned int i = 0; i < n; i++ ) {
        for ( unsigned int j = 0; j < n; j++ ) {
           if ( graph_has_edge(g,i,j) && !graph_has_edge(g,j,i)) {
               graph_add_edge(g,j,i);
           }
        }
    }

    unsigned int edges = 0;
    // adjacency of node i stored in bits
    // so for a 64-bit long we can support a graph
    // up to size 64
    // since we will support up to 128 its two longs
    for ( unsigned int i = 0; i < n; i++ ) {
        //adjM[i][i] = 1;

        // For graph of n > 64 we need 2 LONGs
        // We can support up to 128 nodes
        
        // Write the first LONG
        unsigned long long adBits = 0;
        for ( unsigned int j = 0; j < min(64,n); j++ ) {
           if ( graph_has_edge(g,i,j) ) {
               adjM[i][j]=1;
               unsigned long long one = 1;
                if ( j == 0 )
                    adBits |= one;
                else
                    adBits |= ( (unsigned long long )(one << j ) );
                edges++;
            }
        }
        fprintf(fp,"%llu , ",adBits); 

        // Write the second LONG
        
        adBits = 0;
        for ( unsigned int j = min(64,n); j < n; j++ ) {
            if ( graph_has_edge(g,i,j) ) {
                adjM[i][j]=1;
                unsigned long long one = 1;
                if ( j == 0 )
                    adBits |= one;
                else
                    adBits |= ( (unsigned long long )(one << j ) );
                edges++;
            }
        }
        fprintf(fp,"%llu , ",adBits); 
    }

#ifdef PAD
    unsigned int gsz = n;
    unsigned int padsz = GRAPH_MAX_SIZE - gsz ;

    //Pad the remaining sequences  for GRAPH_MAX_SIZE
    // Write 2 zeros for each node
    for ( unsigned int i = 0; i < padsz; i++ ) {
                fprintf(fp,"0 , 0, "); 
    }

    fprintf(stdout,"Total padded size of graph is = %d\n",gsz+padsz);
#endif

    // Now print the colors of all the nodes
    for ( unsigned int i = 0; i < n; i++ ) {
        // color numberis are > 0 
        if (g->c[i] > 0 ) {
            unsigned int m=g->c[i];
            fprintf(fp,"%d , ",m);
        }
    }
    // PAD color numbers
    // colors are padded with 1
    // which is the first color
    for ( unsigned int i = 0; i < padsz; i++ ) {
                fprintf(fp,"1 , "); 
    }
    fprintf(fp,"1"); 
    free(adjM);
}

int main(int argc, char* argv[]) {
  int i,j,k=0,chi,to=0,n;
  double d;
  int p;
  n = atoi(argv[1]);
  graph_t g=graph_new(n);
  fprintf(stdout,"creating graph with %d nodes\n",n);

  // Create all the graphs on n nodes
  // we can have n*(n-1)/2 graphs
  char dotfile[256];
  FILE *outfp = fopen("newadjcols100.csv","a+");

  // generate a max of 100 graphs for each n combo
  // using different prob distributions and using
  // graph_gnp

  // when n is >9 maxgraphs exceeds size 
  unsigned int maxgraphs = 1000;
  //unsigned int maxgraphs = 5;

  if ( n <= 4 ) 
      maxgraphs  = pow(2,(n*(n-1))/2);

  unsigned int gno = 1;
  do {
      //random graph with edge-prob = 0.98
      graph_gnp(g,0.95);
      if ( nedges(g) > 0 && n < 64 ) {
          graph_show(g);
          chi=graph_chromatic_number(g,0);

#ifdef DEBUG
          fprintf(outfp,"***G=<%d,%d>*******\n",n,x); 
#endif
          fprintf(outfp,"%d ,",chi); 
          dump_node_attr_in_list(g,n,outfp);
          fprintf(outfp,"\n"); 

          graph_empty(g);
      }
      gno++;
  } while (gno <= min(5,maxgraphs));

  gno = 1;
  do {
      //random graph with edge-prob = 0.9
      // Create such dense graphs only for small n's
      graph_gnp(g,0.9);
      if ( nedges(g) > 0 && n < 64 ) {
          fprintf(stdout,"***G=<%d>*******\n",n); 
          graph_show(g);
          chi=graph_chromatic_number(g,0);

#ifdef DEBUG
          fprintf(outfp,"***G=<%d,%d>*******\n",n,x); 
#endif
          fprintf(outfp,"%d ,",chi); 
          dump_node_attr_in_list(g,n,outfp);
          fprintf(outfp,"\n"); 

          graph_empty(g);
      }
      gno++;
  } while (gno <= min(10,maxgraphs));

  gno = 1;
  do {
      //random graph with edge-prob = 0.75
      // Create such dense graphs only for small n's
      graph_gnp(g,0.75);
      if ( nedges(g) > 0 && n < 64 ) {
          graph_show(g);
          chi=graph_chromatic_number(g,0);

#ifdef DEBUG
          fprintf(outfp,"***G=<%d,%d>*******\n",n,x); 
#endif
          fprintf(outfp,"%d ,",chi); 
          dump_node_attr_in_list(g,n,outfp);
          fprintf(outfp,"\n"); 

          graph_empty(g);
      }
      gno++;
  } while (gno <= min(25,maxgraphs));

  gno = 1;
  do {
      //random graph with edge-prob = 0.5
      graph_gnp(g,0.5);
      if ( nedges(g) > 0  && n < 80) {
          graph_show(g);
          chi=graph_chromatic_number(g,0);

#ifdef DEBUG
          fprintf(outfp,"***G=<%d,%d>*******\n",n,x); 
#endif
          fprintf(outfp,"%d ,",chi); 
          dump_node_attr_in_list(g,n,outfp);
          fprintf(outfp,"\n"); 

          graph_empty(g);
      }
      gno++;
  } while (gno <= min(25,maxgraphs));

  gno = 1;
  do {
      //random graph with edge-prob = 0.25
      graph_gnp(g,0.25);
      if ( nedges(g) > 0 ) {
          graph_show(g);
          chi=graph_chromatic_number(g,0);

#ifdef DEBUG
          fprintf(outfp,"***G=<%d,%d>*******\n",n,x); 
#endif
          fprintf(outfp,"%d ,",chi); 
          dump_node_attr_in_list(g,n,outfp);
          fprintf(outfp,"\n"); 

          graph_empty(g);
      }
      gno++;
  } while (gno <= min(25,maxgraphs));

  gno = 1;
  do {
      //random graph with edge-prob = 0.05
      graph_gnp(g,0.05);
      if ( nedges(g) > 0 ) {
          graph_show(g);
          chi=graph_chromatic_number(g,0);

#ifdef DEBUG
          fprintf(outfp,"***G=<%d,%d>*******\n",n,x); 
#endif
          fprintf(outfp,"%d ,",chi); 
          dump_node_attr_in_list(g,n,outfp);
          fprintf(outfp,"\n"); 

          graph_empty(g);
      }
      gno++;
  } while (gno <= min(15,maxgraphs));

  graph_clear(g);
  fclose(outfp);
  return 0;
}

