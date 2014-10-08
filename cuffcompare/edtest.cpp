#include "GArgs.h"
#include "GStr.h"

#define USAGE "Usage:\n\
edtest -a '<string_a>' -b '<string_b>' \n\n"

/*
  public class LevenshteinDistance {
    private static int minimum(int a, int b, int c) {
        return Math.min(Math.min(a, b), c);
    }

    public static int computeLevenshteinDistance(CharSequence str1,
            CharSequence str2) {
        int[][] distance = new int[str1.length() + 1][str2.length() + 1];

        for (int i = 0; i <= str1.length(); i++)
            distance[i][0] = i;
        for (int j = 0; j <= str2.length(); j++)
            distance[0][j] = j;

        for (int i = 1; i <= str1.length(); i++)
            for (int j = 1; j <= str2.length(); j++)
                distance[i][j] = minimum(
                        distance[i - 1][j] + 1,
                        distance[i][j - 1] + 1,
                        distance[i - 1][j - 1]
                                + ((str1.charAt(i - 1) == str2.charAt(j - 1)) ? 0
                                        : 1));

        return distance[str1.length()][str2.length()];
    }
}

 */

class IntArray { //two dimensional int array
    int* mem;
    int xsize;
    int ysize;
  public:
   IntArray(int xlen, int ylen) {
     xsize=xlen;
     ysize=ylen;
     GMALLOC(mem, xsize*ysize*sizeof(int));
     }
   ~IntArray() {
     GFREE(mem);
     }
   int& data(int x, int y) {
    return mem[y*xsize+x];
    }
};

class GXConsensus {
 public:
   char* aa;
   int aalen;
   GXConsensus(const char* a) {
     aa=Gstrdup(a);
     aalen=strlen(a);
     }
   ~GXConsensus() {
     GFREE(aa);
     }
};

int aa_diff(GXConsensus* c1, GXConsensus* c2) {
 int diflen=abs(c1->aalen-c2->aalen);
 if (diflen>=5) return diflen;
 IntArray dist(c1->aalen+1, c2->aalen+1);
 for (int i=0;i<=c1->aalen;i++) {
     dist.data(i,0) = i;
     }
 for (int j = 0; j <= c2->aalen; j++) {
     dist.data(0,j) = j;
     }

 for (int i = 1; i <= c1->aalen; i++)
     for (int j = 1; j <= c2->aalen; j++) {
         dist.data(i,j) = GMIN3( dist.data(i-1,j)+1,
             dist.data(i,j-1)+1,
                 dist.data(i-1,j-1)+((c1->aa[i-1] == c2->aa[j-1]) ? 0 : 1) );
         }
 int r=dist.data(c1->aalen,c2->aalen);
 return r;
}

/*
int aa_diff(GXConsensus* c1, GXConsensus* c2) {
 //int diflen=abs(c1->aalen-c2->aalen);
 //if (diflen>=5) return diflen;
 int **distance=NULL;
 GCALLOC(distance,(c1->aalen+1)*sizeof(int*));
 for (int i=0;i<=c1->aalen;i++) {
     GCALLOC(distance[i],(c2->aalen+1)*sizeof(int));
     distance[i][0] = i;
     }
 for (int j = 0; j <= c2->aalen; j++)
     distance[0][j] = j;

 for (int i = 1; i <= c1->aalen; i++)
     for (int j = 1; j <= c2->aalen; j++) {
         distance[i][j] = GMIN3( distance[i-1][j]+1,
                       distance[i][j-1]+1,
                       distance[i-1][j-1]+((c1->aa[i-1] == c2->aa[j-1]) ? 0 : 1) );
         }
 int r=distance[c1->aalen][c2->aalen];
 for (int i=0;i<=c1->aalen;i++)
     GFREE(distance[i]);
 GFREE(distance);
 return r;
}
*/


int main(int argc, char * const argv[]) {
	GArgs args(argc, argv, "Dha:b:");
	int e;
	if ((e=args.isError())>0)
		GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
	if (args.getOpt('h')!=NULL) GError("%s\n", USAGE);
	//bool debug=(args.getOpt('D')!=NULL);
	GXConsensus* ca=NULL;
	GXConsensus* cb=NULL;
    if (args.getOpt('a')!=NULL) {
      ca=new GXConsensus(args.getOpt('a'));
      }
    else GError("%s Required -a option missing!\n", USAGE);
    if (args.getOpt('b')!=NULL) {
      cb=new GXConsensus(args.getOpt('b'));
      }
    else GError("%s Required -b option missing!\n", USAGE);
    int r=aa_diff(ca,cb);
    GMessage("Edit distance is: %d\n",r);
}
