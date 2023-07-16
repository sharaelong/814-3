#pragma GCC optimize("O3")
#include <bits/stdc++.h>
#ifdef SHARAELONG
#include "debug.hpp"
#endif
using namespace std;
typedef long long ll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;

typedef long long GainType;
#ifndef LLONG_MAX
#define LLONG_MAX 9223372036854775807LL
#endif
#ifndef LLONG_MIN
#define LLONG_MIN (-LLONG_MAX - 1LL)
#endif
#define PLUS_INFINITY LLONG_MAX
#define MINUS_INFINITY LLONG_MIN
#define GainFormat "%lld"
#define GainInputFormat "%lld"

#define HashTableSize 65521     /* Largest prime less than USHRT_MAX */
#define MaxLoadFactor 0.75

typedef struct HashTableEntry {
    unsigned Hash;
    GainType Cost;
} HashTableEntry;

typedef struct HashTable {
    HashTableEntry Entry[HashTableSize];
    int Count; /* Number of occupied entries */
} HashTable;

void HashInitialize(HashTable * T) {
    int i;
    for (i = 0; i < HashTableSize; i++) {
        T->Entry[i].Hash = UINT_MAX;
        T->Entry[i].Cost = MINUS_INFINITY;
    }
    T->Count = 0;
}

void HashInsert(HashTable * T, unsigned Hash, GainType Cost)
{
    int i = Hash % HashTableSize;
    if (T->Count >= MaxLoadFactor * HashTableSize) {
        if (Cost > T->Entry[i].Cost)
            return;
    } else {
        int p = Hash % 97 + 1;
        while (T->Entry[i].Cost != MINUS_INFINITY)
            if ((i -= p) < 0)
                i += HashTableSize;
        T->Count++;
    }
    T->Entry[i].Hash = Hash;
    T->Entry[i].Cost = Cost;
}

int HashSearch(HashTable * T, unsigned Hash, GainType Cost)
{
    int i, p;

    i = Hash % HashTableSize;
    p = Hash % 97 + 1;
    while ((T->Entry[i].Hash != Hash || T->Entry[i].Cost != Cost)
           && T->Entry[i].Cost != MINUS_INFINITY)
        if ((i -= p) < 0)
            i += HashTableSize;
    return T->Entry[i].Hash == Hash;
}

typedef struct Node Node;
typedef struct Candidate Candidate;
typedef struct Segment Segment;
typedef struct SSegment SSegment;
typedef struct SwapRecord SwapRecord;
typedef Node *(*MoveFunction) (Node * t1, Node * t2, GainType * G0,
                               GainType * Gain);
typedef int (*CostFunction) (Node * Na, Node * Nb);

/* Genetic.h */
int MaxPopulationSize; /* The maximum size of the population */ 

void AddToPopulation(GainType Cost);
void ApplyCrossover(int i, int j);
void ReplaceIndividualWithTour(int i, GainType Cost);
int ReplacementIndividual(GainType Cost);

// void ERXT();

/* LKH.h */
#define Fixed(a, b) ((a)->FixedTo1 == (b) || (a)->FixedTo2 == (b))
#define FixedOrCommon(a, b) (Fixed(a, b) || IsCommonEdge(a, b))
#define InBestTour(a, b) ((a)->BestSuc == (b) || (b)->BestSuc == (a))
#define InNextBestTour(a, b)                                    \
    ((a)->NextBestSuc == (b) || (b)->NextBestSuc == (a))
#define InInputTour(a, b) ((a)->InputSuc == (b) || (b)->InputSuc == (a))
#define InInitialTour(a, b)                             \
    ((a)->InitialSuc == (b) || (b)->InitialSuc == (a))
#define Near(a, b)                                                      \
    ((a)->BestSuc ? InBestTour(a, b) : (a)->Dad == (b) || (b)->Dad == (a))

#define Link(a, b) { ((a)->Suc = (b))->Pred = (a); }
#define Follow(b, a)                                                    \
    { Link((b)->Pred, (b)->Suc); Link(b, b); Link(b, (a)->Suc); Link(a, b); }
#define Precede(a, b)                                                   \
    { Link((a)->Pred, (a)->Suc); Link(a, a); Link((b)->Pred, a); Link(a, b); }
#define SLink(a, b) { (a)->Suc = (b); (b)->Pred = (a); }

enum Types { TSP, ATSP, SOP, HCP, CVRP, TOUR, HPP };
enum CoordTypes { TWOD_COORDS, THREED_COORDS, NO_COORDS };
enum EdgeWeightTypes { EXPLICIT, EUC_2D, EUC_3D, MAX_2D, MAX_3D, MAN_2D,
                       MAN_3D, CEIL_2D, CEIL_3D, FLOOR_2D, FLOOR_3D,
                       GEO, GEOM, GEO_MEEUS, GEOM_MEEUS, ATT, TOR_2D, TOR_3D,
                       XRAY1, XRAY2, SPECIAL
};
enum EdgeWeightFormats { FUNCTION, FULL_MATRIX, UPPER_ROW, LOWER_ROW,
                         UPPER_DIAG_ROW, LOWER_DIAG_ROW, UPPER_COL, LOWER_COL,
                         UPPER_DIAG_COL, LOWER_DIAG_COL
};
enum CandidateSetTypes { ALPHA, DELAUNAY, NN, POPMUSIC, QUADRANT };
enum InitialTourAlgorithms { BORUVKA, GREEDY, MOORE, NEAREST_NEIGHBOR,
                             QUICK_BORUVKA, SIERPINSKI, WALK
};
enum RecombinationTypes { IPT, GPX2, CLARIST };


/* The Node structure is used to represent nodes (cities) of the problem */

struct Node {
    int Id;     /* Number of the node (1...Dimension) */
    int Loc;    /* Location of the node in the heap 
                   (zero, if the node is not in the heap) */
    int Rank;   /* During the ascent, the priority of the node.
                   Otherwise, the ordinal number of the node in 
                   the tour */
    int V;      /* During the ascent the degree of the node minus 2.
                   Otherwise, the variable is used to mark nodes */
    int LastV;  /* Last value of V during the ascent */
    int Cost;   /* "Best" cost of an edge emanating from the node */
    int NextCost; /* During the ascent, the next best cost of an edge
                     emanating from the node */
    int PredCost, /* The costs of the neighbor edges on the current tour */ 
        SucCost; 
    int SavedCost;
    int Pi;     /* Pi-value of the node */
    int BestPi; /* Currently best pi-value found during the ascent */
    int Beta;   /* Beta-value (used for computing alpha-values) */
    int Subproblem;  /* Number of the subproblem the node is part of */
    int Sons;   /* Number of sons in the minimum spanning tree */
    int *C;     /* A row in the cost matrix */
    Node *Pred, *Suc;  /* Predecessor and successor node in 
                          the two-way list of nodes */
    Node *OldPred, *OldSuc; /* Previous values of Pred and Suc */
    Node *BestSuc,     /* Best and next best successor node in the */
        *NextBestSuc; /* currently best tour */
    Node *Dad;  /* Father of the node in the minimum 1-tree */
    Node *Nearest;     /* Nearest node (used in the greedy heuristics) */
    Node *Next; /* Auxiliary pointer, usually to the next node in a list
                   of nodes (e.g., the list of "active" nodes) */
    Node *Prev; /* Auxiliary pointer, usually to the previous node 
                   in a list of nodes */
    Node *Mark; /* Visited mark */
    Node *FixedTo1,    /* Pointers to the opposite end nodes of fixed edges. */
        *FixedTo2;    /* A maximum of two fixed edges can be incident
                         to a node */
    Node *FixedTo1Saved, /* Saved values of FixedTo1 and FixedTo2 */
        *FixedTo2Saved;
    Node *Head; /* Head of a segment of common edges */
    Node *Tail; /* Tail of a segment of common edges */
    Node *InputSuc;    /* Successor in the INPUT_TOUR file */
    Node *InitialSuc;  /* Successor in the INITIAL_TOUR file */
    Node *SubproblemPred; /* Predecessor in the SUBPROBLEM_TOUR file */
    Node *SubproblemSuc;  /* Successor in the SUBPROBLEM_TOUR file */
    Node *SubBestPred; /* The best predecessor node in a subproblem */
    Node *SubBestSuc;  /* The best successor node in a subproblem */
    Node *Added1, *Added2; /* Pointers to the opposite end nodes
                              of added edges in a submove */
    Node *Deleted1, *Deleted2;  /* Pointers to the opposite end nodes
                                   of deleted edges in a submove */
    Candidate *CandidateSet;    /* Candidate array */
    Candidate *BackboneCandidateSet; /* Backbone candidate array */
    Segment *Parent;   /* Parent segment of a node when the two-level
                          tree representation is used */
    double X, Y, Z;     /* Coordinates of the node */
    double Xc, Yc, Zc;  /* Converted coordinates */
    char Axis;  /* The axis partitioned when the node is part of a KDTree */
    char OldPredExcluded, OldSucExcluded;  /* Booleans used for indicating 
                                              whether one (or both) of the 
                                              adjoining nodes on the old tour 
                                              has been excluded */
};

/* The Candidate structure is used to represent candidate edges */

struct Candidate {
    Node *To;   /* The end node of the edge */
    int Cost;   /* Cost (distance) of the edge */
    int Alpha;  /* Its alpha-value */
};

/* The Segment structure is used to represent the segments in the two-level 
   representation of tours */

struct Segment {
    char Reversed;       /* Reversal bit */
    Node *First, *Last;  /* First and last node in the segment */
    Segment *Pred, *Suc; /* Predecessor and successor in the two-way 
                            list of segments */
    int Rank;   /* Ordinal number of the segment in the list */
    int Size;   /* Number of nodes in the segment */
    SSegment *Parent;    /* The parent super segment */
};

struct SSegment {
    char Reversed;         /* Reversal bit */
    Segment *First, *Last; /* The first and last node in the segment */
    SSegment *Pred, *Suc;  /* The predecessor and successor in the
                              two-way list of super segments */
    int Rank;   /* The ordinal number of the segment in the list */
    int Size;   /* The number of nodes in the segment */
};

/* The SwapRecord structure is used to record 2-opt moves (swaps) */

struct SwapRecord {
    Node *t1, *t2, *t3, *t4;    /* The 4 nodes involved in a 2-opt move */
};


struct LKH {
    /* Sequence.h */
    Node **t;       /* The sequence of nodes to be used in a move */
    Node **T;       /* The currently best t's */
    Node **tSaved;  /* For saving t when using the BacktrackKOptMove function */
    int *p;         /* The permutation corresponding to the sequence in which
                       the t's occur on the tour */
    int *q;         /* The inverse permutation of p */
    int *incl;      /* Array: incl[i] == j, if (t[i], t[j]) is an inclusion edge */
    int *cycle;     /* Array: cycle[i] is cycle number of t[i] */
    GainType *G;    /* For storing the G-values in the BestKOptMove function */
    int K;          /* The value K for the current K-opt move */
    
    int AscentCandidates;
    int BackboneTrials;     /* Number of backbone trials in each run */
    int Backtracking;       /* Specifies whether backtracking is used for 
                               the first move in a sequence of moves */
    GainType BestCost;      
    int *BestTour;  
    GainType BetterCost;    
    int *BetterTour;        /* Table containing the currently best tour 
                               in a run */
    int CacheMask;  
    int *CacheVal;  
    int *CacheSig;  /* Table of the signatures of cached 
                       distances */
    int CandidateFiles;     /* Number of CANDIDATE_FILEs */
    int *CostMatrix;        
    int Dimension;  
    int DimensionSaved;     
    int EdgeFiles;          /* Number of EDGE_FILEs */
    double Excess;  /* Maximum alpha-value allowed for any 
                       candidate edge is set to Excess times the 
                       absolute value of the lower bound of a 
                       solution tour */
    int ExtraCandidates;    /* Number of extra neighbors to be added to 
                               the candidate set of each node */
    Node *FirstActive, *LastActive; /* First and last node in the list 
                                       of "active" nodes */
    Node *FirstNode;        
    Segment *FirstSegment;  /* A pointer to the first segment in the cyclic 
                               list of segments */
    SSegment *FirstSSegment;  /* A pointer to the first super segment in
                                 the cyclic list of segments */
    int Gain23Used; 
    int GainCriterionUsed;  /* Specifies whether L&K's gain criterion is 
                               used */
    double GridSize;        
    int GroupSize;  
    int SGroupSize; 
    int Groups;     
    int SGroups;    
    unsigned Hash;  
    Node **Heap;    
    HashTable *HTable;      
    int InitialPeriod;      
    int InitialStepSize;    
    double InitialTourFraction; /* Fraction of the initial tour to be 
                                   constructed by INITIAL_TOUR_FILE edges */
    char *LastLine; 
    double LowerBound;  
    int Kicks;      /* Specifies the number of K-swap-kicks */
    int KickType;   /* Specifies K for a K-swap-kick */
    int M;          /* The M-value is used when solving an ATSP-
                       instance by transforming it to a STSP-instance */
    int MaxBreadth; /* The maximum number of candidate edges 
                       considered at each level of the search for
                       a move */
    int MaxCandidates;      /* Maximum number of candidate edges to be 
                               associated with each node */
    int MaxMatrixDimension; /* Maximum dimension for an explicit cost
                               matrix */
    int MaxSwaps;   /* Maximum number of swaps made during the 
                       search for a move */
    int MaxTrials;  
    int MoveType;   /* Specifies the sequantial move type to be used 
                       in local search. A value K >= 2 signifies 
                       that a k-opt moves are tried for k <= K */
    Node *NodeSet;  
    int Norm;       /* Measure of a 1-tree's discrepancy from a tour */
    int NonsequentialMoveType; /* Specifies the nonsequential move type to
                                  be used in local search. A value 
                                  L >= 4 signifies that nonsequential
                                  l-opt moves are tried for l <= L */
    GainType Optimum;       /* Known optimal tour length. 
                               If StopAtOptimum is 1, a run will be 
                               terminated as soon as a tour length 
                               becomes equal this value */
    int PatchingA;  /* Specifies the maximum number of alternating
                       cycles to be used for patching disjunct cycles */
    int PatchingC;  /* Specifies the maximum number of disjoint cycles to be 
                       patched (by one or more alternating cycles) */
    int Precision;  /* Internal precision in the representation of 
                       transformed distances */
    int PredSucCostAvailable;  
    int POPMUSIC_InitialTour;  /* Specifies whether the first POPMUSIC tour
                                  is used as initial tour for LK */
    int POPMUSIC_MaxNeighbors; /* Maximum number of nearest neighbors used 
                                  as candidates in iterated 3-opt */
    int POPMUSIC_SampleSize;   
    int POPMUSIC_Solutions;    
    int POPMUSIC_Trials;       /* Maximum trials used for iterated 3-opt */
    unsigned *Rand; 
    int Recombination; /* IPT, GPX2 or CLARIST */
    int RestrictedSearch;      /* Specifies whether the choice of the first 
                                  edge to be broken is restricted */
    short Reversed; /* Boolean used to indicate whether a tour has 
                       been reversed */
    int Run;        
    int Runs;       
    unsigned Seed;  
    double StartTime;       
    int StopAtOptimum;      /* Specifies whether a run will be terminated if 
                               the tour length becomes equal to Optimum */
    int Subgradient;        /* Specifies whether the Pi-values should be 
                               determined by subgradient optimization */
    int SubproblemSize;     
    int SubsequentMoveType; /* Specifies the move type to be used for all 
                               moves following the first move in a sequence 
                               of moves. The value K >= 2 signifies that a 
                               K-opt move is to be used */
    int SubsequentPatching; /* Species whether patching is used for 
                               subsequent moves */
    SwapRecord *SwapStack;  
    int Swaps;      
    double TimeLimit;       
    double TotalTimeLimit;  
    int TraceLevel; /* Specifies the level of detail of the output 
                       given during the solution process. 
                       The value 0 signifies a minimum amount of 
                       output. The higher the value is the more 
                       information is given */
    int Trial;      

    /* The following variables are read by the functions ReadParameters and 
       ReadProblem: */

    char *ParameterFileName, *ProblemFileName, *PiFileName,
        *TourFileName, *OutputTourFileName, *InputTourFileName,
        **CandidateFileName, **EdgeFileName, *InitialTourFileName,
        *SubproblemTourFileName;
    char *Name, *Type, *EdgeWeightType, *EdgeWeightFormat,
        *EdgeDataFormat, *NodeCoordType, *DisplayDataType;
    int CandidateSetSymmetric, CandidateSetType,
        CoordType, DelaunayPartitioning, DelaunayPure,
        ExtraCandidateSetSymmetric, ExtraCandidateSetType,
        InitialTourAlgorithm,
        KarpPartitioning, KCenterPartitioning, KMeansPartitioning,
        MoorePartitioning,
        PatchingAExtended, PatchingARestricted,
        PatchingCExtended, PatchingCRestricted,
        ProblemType,
        RohePartitioning, SierpinskiPartitioning,
        SubproblemBorders, SubproblemsCompressed, WeightType, WeightFormat;

    FILE *ParameterFile, *ProblemFile, *PiFile,
        *TourFile, *InitialTourFile;
    // CostFunction Distance;
    CostFunction c;
    // CostFunction D, C;
    // MoveFunction BestMove, BestSubsequentMove;
    MoveFunction BacktrackMove;
    
    /* LKHmain.c */
    GainType Cost, OldOptimum;
    double Time, LastTime;

    /* Distance.c */
    int Distance(Node * Na, Node * Nb)
    {
        double xd = Na->X - Nb->X, yd = Na->Y - Nb->Y;
        return (int) (sqrt(xd * xd + yd * yd) + 0.5);
    }

    /* C.c */
    int C(Node * Na, Node * Nb)
    {
        return Na->Id < Nb->Id ? Nb->C[Na->Id] : Na->C[Nb->Id];
    }

    int D(Node * Na, Node * Nb)
    {
        return (Na->Id <
                Nb->Id ? Nb->C[Na->Id] : Na->C[Nb->Id]) + Na->Pi + Nb->Pi;
    }

    /* Heap.c */
    int HeapCount;
    int HeapCapacity;

    void HeapMake(int Size) {
        Heap = (Node **) malloc((Size + 1) * sizeof(Node *));
        HeapCapacity = Size;
        HeapCount = 0;
    }
               
    void HeapSiftUp(Node * N) {
        int Loc = N->Loc, Parent = Loc / 2;
        while (Parent && N->Rank < Heap[Parent]->Rank) {
            Heap[Loc] = Heap[Parent];
            Heap[Loc]->Loc = Loc;
            Loc = Parent;
            Parent /= 2;
        }
        Heap[Loc] = N;
        N->Loc = Loc;
    }

    void HeapSiftDown(Node * N) {
        int Loc = N->Loc, Child;
        while (Loc <= HeapCount / 2) {
            Child = 2 * Loc;
            if (Child < HeapCount && Heap[Child + 1]->Rank < Heap[Child]->Rank)
                Child++;
            if (N->Rank <= Heap[Child]->Rank)
                break;
            Heap[Loc] = Heap[Child];
            Heap[Loc]->Loc = Loc;
            Loc = Child;
        }
        Heap[Loc] = N;
        N->Loc = Loc;
    }

    Node *HeapDeleteMin() {
        Node *Remove;
        if (!HeapCount)
            return 0;
        Remove = Heap[1];
        Heap[1] = Heap[HeapCount--];
        Heap[1]->Loc = 1;
        HeapSiftDown(Heap[1]);
        Remove->Loc = 0;
        return Remove;
    }

    void HeapInsert(Node * N) {
        HeapLazyInsert(N);
        HeapSiftUp(N);
    }

    void HeapDelete(Node * N) {
        int Loc = N->Loc;
        if (!Loc)
            return;
        Heap[Loc] = Heap[HeapCount--];
        Heap[Loc]->Loc = Loc;
        if (Heap[Loc]->Rank > N->Rank)
            HeapSiftDown(Heap[Loc]);
        else
            HeapSiftUp(Heap[Loc]);
        N->Loc = 0;
    }

    void HeapLazyInsert(Node * N) {
        assert(HeapCount < HeapCapacity);
        Heap[++HeapCount] = N;
        N->Loc = HeapCount;
    }

    void Heapify() {
        int Loc;
        for (Loc = HeapCount / 2; Loc >= 1; Loc--)
            HeapSiftDown(Heap[Loc]);
    }
    
    /* Segment.h */
#define PRED(a) (Reversed ? (a)->Suc : (a)->Pred)
#define SUC(a) (Reversed ? (a)->Pred : (a)->Suc)
#define BETWEEN(a, b, c) Between(a, b, c)
#define FLIP(a, b, c, d) Flip(a, b, c)

#define Swap1(a1,a2,a3)                         \
    FLIP(a1,a2,a3,0)
#define Swap2(a1,a2,a3, b1,b2,b3)               \
    (Swap1(a1,a2,a3), Swap1(b1,b2,b3))
#define Swap3(a1,a2,a3, b1,b2,b3, c1,c2,c3)             \
    (Swap2(a1,a2,a3, b1,b2,b3), Swap1(c1,c2,c3))
#define Swap4(a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3)           \
    (Swap3(a1,a2,a3, b1,b2,b3, c1,c2,c3), Swap1(d1,d2,d3))
#define Swap5(a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3, e1,e2,e3)         \
    (Swap4(a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3), Swap1(e1,e2,e3))

    /* Make3OptMove.c */
    void
    Make3OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
                 Node * t5, Node * t6, int Case)
    {
        switch (Case) {
        case 1:
        case 2:
            Swap2(t1, t2, t3, t6, t5, t4);
            return;
        case 5:
            Swap3(t1, t2, t4, t6, t5, t4, t6, t2, t3);
            return;
        case 6:
            Swap2(t3, t4, t5, t1, t2, t3);
            return;
        default:
            assert(0);
            return;
        }
    }

    /* Make4OptMove.c */
    void
    Make4OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
                 Node * t5, Node * t6, Node * t7, Node * t8, int Case)
    {
        if (SUC(t1) != t2)
            Reversed ^= 1;
        switch (Case) {
        case 1:
        case 2:
            Swap3(t1, t2, t3, t6, t5, t4, t7, t8, t1);
            return;
        case 3:
        case 4:
            Swap3(t1, t2, t3, t8, t7, t6, t5, t8, t1);
            return;
        case 5:
            if (!BETWEEN(t2, t7, t3))
                Swap3(t5, t6, t7, t2, t1, t4, t1, t4, t5);
            else if (BETWEEN(t2, t7, t6))
                Swap3(t5, t6, t7, t5, t8, t3, t3, t8, t1);
            else
                Swap3(t1, t2, t7, t7, t2, t3, t4, t7, t6);
            return;
        case 6:
            Swap3(t3, t4, t5, t6, t3, t2, t1, t6, t7);
            return;
        case 7:
            Swap3(t6, t5, t8, t2, t1, t4, t8, t5, t4);
            return;
        case 11:
            Swap3(t1, t2, t7, t3, t4, t5, t3, t6, t7);
            return;
        case 12:
            Swap3(t3, t4, t5, t7, t8, t1, t3, t6, t7);
            return;
        case 15:
            Swap3(t3, t4, t5, t3, t6, t7, t8, t3, t2);
            return;
        default:
            assert(0);
            return;
        }
    }

    /* Flip.c */
    void Flip(Node * t1, Node * t2, Node * t3)
    {
        Node *s1, *s2, *t4;
        int R, Temp, Ct2t3, Ct4t1;

        assert(t1->Pred == t2 || t1->Suc == t2);
        if (t3 == t2->Pred || t3 == t2->Suc)
            return;
        t4 = t1->Suc == t2 ? t3->Pred : t3->Suc;
        if (t1->Suc != t2) {
            s1 = t1;
            t1 = t2;
            t2 = s1;
            s1 = t3;
            t3 = t4;
            t4 = s1;
        }
        /* Find the segment with the smallest number of nodes */
        if ((R = t2->Rank - t3->Rank) < 0)
            R += Dimension;
        if (2 * R > Dimension) {
            s1 = t3;
            t3 = t2;
            t2 = s1;
            s1 = t4;
            t4 = t1;
            t1 = s1;
        }
        Ct2t3 = C(t2, t3);
        Ct4t1 = C(t4, t1);
        /* Swap segment (t3 --> t1) */
        R = t1->Rank;
        t1->Suc = 0;
        s2 = t3;
        while ((s1 = s2)) {
            s2 = s1->Suc;
            s1->Suc = s1->Pred;
            s1->Pred = s2;
            s1->Rank = R--;
            Temp = s1->SucCost;
            s1->SucCost = s1->PredCost;
            s1->PredCost = Temp;
        }
        (t3->Suc = t2)->Pred = t3;
        (t4->Suc = t1)->Pred = t4;
        t3->SucCost = t2->PredCost = Ct2t3;
        t1->PredCost = t4->SucCost = Ct4t1;
        SwapStack[Swaps].t1 = t1;
        SwapStack[Swaps].t2 = t2;
        SwapStack[Swaps].t3 = t3;
        SwapStack[Swaps].t4 = t4;
        Swaps++;
        Hash ^= (Rand[t1->Id] * Rand[t2->Id]) ^
            (Rand[t3->Id] * Rand[t4->Id]) ^
            (Rand[t2->Id] * Rand[t3->Id]) ^ (Rand[t4->Id] * Rand[t1->Id]);
    }

    /* Between.c */
    int Between(const Node * ta, const Node * tb, const Node * tc)
    {
        int a, b = tb->Rank, c;

        if (!Reversed) {
            a = ta->Rank;
            c = tc->Rank;
        } else {
            a = tc->Rank;
            c = ta->Rank;
        }
        return a <= c ? b >= a && b <= c : b >= a || b <= c;
    }

    /* Sequence.c */
    Node *tp1;

    void FindPermutation(int k)
    {
        int i, j;

        for (i = j = 1; j <= k; i += 2, j++)
            p[j] = SUC(t[i]) == t[i + 1] ? i : i + 1;
        tp1 = t[p[1]];
        sort(p + 2, p + k + 1, [&](int pa, int pb) {
            return BETWEEN(tp1, t[pa], t[pb]);
        });
        // qsort(p + 2, k - 1, sizeof(int), compare2);
        for (j = 2 * k; j >= 2; j -= 2) {
            p[j - 1] = i = p[j / 2];
            p[j] = i & 1 ? i + 1 : i - 1;
        }
        for (i = 1; i <= 2 * k; i++)
            q[p[i]] = i;
    }

    /*  
     * The FeasibleKOptMove function tests whether the move given by
     * t[1..2k] and incl[1..2k] represents a feasible k-opt move,
     * i.e., making the move on the current tour will result in a tour.
     *   
     * In that case, 1 is returned. Otherwise, 0 is returned. 
     */

    int FeasibleKOptMove(int k)
    {
        int Count, i;

        FindPermutation(k);
        for (Count = 1, i = 2 * k; (i = q[incl[p[i]]] ^ 1); Count++);
        return Count == k;
    }

    /*
     * The Cycles function returns the number of cycles that would appear if 
     * the move given by t[1..2k] and incl[1..2k] was made. 
     * In addition, cycle[i] is assigned the number of the cycle that node t[i] 
     * is a part of (an integer from 1 to Cycles).
     */

    int Cycles(int k)
    {
        int i, j, Count = 0;

        for (i = 1; i <= 2 * k; i++)
            cycle[i] = 0;
        for (i = 1; i <= 2 * k; i++) {
            if (!cycle[p[i]]) {
                Count++;
                j = i;
                do {
                    cycle[p[j]] = Count;
                    j = q[incl[p[j]]];
                    cycle[p[j]] = Count;
                    if ((j ^= 1) > 2 * k)
                        j = 1;
                }
                while (j != i);
            }
        }
        return Count;
    }
    
    int Added(const Node * ta, const Node * tb)
    {
        return ta->Added1 == tb || ta->Added2 == tb;
    }

    /* 
     * The Deleted function is used to test if an edge, (ta,tb), 
     * of the tour has been deleted in the submove under construction.
     */

    int Deleted(const Node * ta, const Node * tb)
    {
        return ta->Deleted1 == tb || ta->Deleted2 == tb;
    }
    
    void MarkAdded(Node * ta, Node * tb)
    {
        if (!ta->Added1)
            ta->Added1 = tb;
        else if (!ta->Added2)
            ta->Added2 = tb;
        if (!tb->Added1)
            tb->Added1 = ta;
        else if (!tb->Added2)
            tb->Added2 = ta;
    }

    /*
     * The MarkDeletedfunction is used to mark an edge, (ta,tb), as deleted
     * in the submove under construction.
     */

    void MarkDeleted(Node * ta, Node * tb)
    {
        if (!ta->Deleted1)
            ta->Deleted1 = tb;
        else if (!ta->Deleted2)
            ta->Deleted2 = tb;
        if (!tb->Deleted1)
            tb->Deleted1 = ta;
        else if (!tb->Deleted2)
            tb->Deleted2 = ta;
    }

    /*
     * The UnmarkAdded function is used mark the edge, (ta,tb), as not
     * added.
     */

    void UnmarkAdded(Node * ta, Node * tb)
    {
        if (ta->Added1 == tb)
            ta->Added1 = 0;
        else if (ta->Added2 == tb)
            ta->Added2 = 0;
        if (tb->Added1 == ta)
            tb->Added1 = 0;
        else if (tb->Added2 == ta)
            tb->Added2 = 0;
    }

    /*
     * The UnmarkDeleted function is used mark the edge, (ta,tb), as not
     * deleted.
     */

    void UnmarkDeleted(Node * ta, Node * tb)
    {
        if (ta->Deleted1 == tb)
            ta->Deleted1 = 0;
        else if (ta->Deleted2 == tb)
            ta->Deleted2 = 0;
        if (tb->Deleted1 == ta)
            tb->Deleted1 = 0;
        else if (tb->Deleted2 == ta)
            tb->Deleted2 = 0;
    }

    /* MakeKOptMove.c */
    void MakeKOptMove(int K)
    {
        int i, j, Best_i = 0, Best_j = 0, BestScore, s;

        FindPermutation(K);
    FindNextReversal:
        /* Find the oriented reversal that has maximal score */
        BestScore = -1;
        for (i = 1; i <= 2 * K - 2; i++) {
            j = q[incl[p[i]]];
            if (j >= i + 2 && (i & 1) == (j & 1) &&
                (s = i & 1 ? Score(i + 1, j, K) :
                 Score(i, j - 1, K)) > BestScore) {
                BestScore = s;
                Best_i = i;
                Best_j = j;
            }
        }
        if (BestScore >= 0) {
            i = Best_i;
            j = Best_j;
            if (i & 1) {
                Swap1(t[p[i + 1]], t[p[i]], t[p[j]]);
                Reverse(i + 1, j);
            } else {
                Swap1(t[p[i - 1]], t[p[i]], t[p[j]]);
                Reverse(i, j - 1);
            }
            goto FindNextReversal;
        }
        /* No more oriented reversals. Cut a simpe hurdle, if any.
         * Note that there can be no super hurdles */
        for (i = 1; i <= 2 * K - 3; i += 2) {
            j = q[incl[p[i]]];
            if (j >= i + 3) {
                Swap1(t[p[i]], t[p[i + 1]], t[p[j]]);
                Reverse(i + 1, j - 1);
                goto FindNextReversal;
            }
        }
    }

    /*
     * The Reverse function reverses the sequence of elements in p[i:j].
     * The inverse permutation q is updated accordingly.
     */

    void Reverse(int i, int j)
    {
        for (; i < j; i++, j--) {
            int pi = p[i];
            q[p[i] = p[j]] = i;
            q[p[j] = pi] = j;
        }
    }

    /*
     * The Score function computes the score of a reversal. The score is the 
     * number of oriented pairs in the resulting reversal.
     */

    int Score(int Left, int Right, int K)
    {
        int Count = 0, i, j;

        Reverse(Left, Right);
        for (i = 1; i <= 2 * K - 2; i++) {
            j = q[incl[p[i]]];
            if (j >= i + 2 && (i & 1) == (j & 1))
                Count++;
        }
        Reverse(Left, Right);
        return Count;
    }

    /* Exclude.c */
    void Exclude(Node * ta, Node * tb)
    {
        if (ta == tb->Pred || ta == tb->Suc)
            return;
        if (ta == tb->OldPred)
            tb->OldPredExcluded = 1;
        else if (ta == tb->OldSuc)
            tb->OldSucExcluded = 1;
        if (tb == ta->OldPred)
            ta->OldPredExcluded = 1;
        else if (tb == ta->OldSuc)
            ta->OldSucExcluded = 1;
    }

    /* Excludable.c */
    int Excludable(Node * ta, Node * tb)
    {
        if (ta == tb->OldPred)
            return !tb->OldPredExcluded;
        if (ta == tb->OldSuc)
            return !tb->OldSucExcluded;
        return 0;
    }

    /* Forbidden.c */
    int Forbidden(const Node * ta, const Node * tb)
    {
        return 0;
    }

    /* IsCommonEdge.c */
    int IsCommonEdge(const Node * ta, const Node * tb)
    {
        return 0;
    }

    /* FixedOrCommonCandidates.c */
    int FixedOrCommonCandidates(Node * N)
    {
        int Count = N->FixedTo2 ? 2 : N->FixedTo1 ? 1 : 0;
        return Count;
    }

    /* IsCandidate.c */
    int IsCandidate(const Node * ta, const Node * tb)
    {
        Candidate *Nta;
        for (Nta = ta->CandidateSet; Nta && Nta->To; Nta++)
            if (Nta->To == tb)
                return 1;
        return 0;
    }

    /* IsPossibleCandidate.c */
    int IsPossibleCandidate(Node * From, Node * To)
    {
        // if (Forbidden(From, To))
        //     return 0;
        // if (InInitialTour(From, To) ||
        //     From->SubproblemSuc == To || To->SubproblemSuc == From ||
        //     FixedOrCommon(From, To))
        //     return 1;
        // if (From->FixedTo2 || To->FixedTo2)
        //     return 0;
        // if (!IsCandidate(From, To) &&
        //     (FixedOrCommonCandidates(From) == 2 ||
        //      FixedOrCommonCandidates(To) == 2))
        //     return 0;
        return 1;
    }

    /* AddCandidate.c */
    int AddCandidate(Node * From, Node * To, int Cost, int Alpha)
    {
        int Count;
        Candidate *NFrom;

        // if (From->Subproblem != FirstNode->Subproblem ||
        //     To->Subproblem != FirstNode->Subproblem ||
        if (Cost == INT_MAX)
            return 0;
        if (From->CandidateSet == 0)
            From->CandidateSet = (Candidate *) calloc(3, sizeof(Candidate));
        if (From == To || To->Subproblem != FirstNode->Subproblem ||
            !IsPossibleCandidate(From, To))
            return 0;
        Count = 0;
        for (NFrom = From->CandidateSet; NFrom->To && NFrom->To != To; NFrom++)
            Count++;
        if (NFrom->To) {
            if (NFrom->Alpha == INT_MAX)
                NFrom->Alpha = Alpha;
            return 0;
        }
        NFrom->Cost = Cost;
        NFrom->Alpha = Alpha;
        NFrom->To = To;
        From->CandidateSet =
            (Candidate *) realloc(From->CandidateSet,
                                  (Count + 2) * sizeof(Candidate));
        From->CandidateSet[Count + 1].To = 0;
        return 1;
    }

    /* AddTourCandidates.c */
    void AddTourCandidates()
    {
        Node *Na, *Nb;
        int i, d, Subproblem = FirstNode->Subproblem;

        /* Add fixed edges */
        Na = FirstNode;
        do {
            if (Na->FixedTo1)
                AddCandidate(Na, Na->FixedTo1, D(Na, Na->FixedTo1), 0);
            if (Na->FixedTo2)
                AddCandidate(Na, Na->FixedTo2, D(Na, Na->FixedTo2), 0);
        }
        while ((Na = Na->Suc) != FirstNode);

        /* Add MERGE_TOUR_FILE edges */
        // removed

        /* Add INITIAL_TOUR_FILE edges */
        Na = FirstNode;
        do {
            Nb = Na->InitialSuc;
            if (!Nb)
                break;
            if (Na->Subproblem == Subproblem && Nb->Subproblem == Subproblem) {
                d = D(Na, Nb);
                AddCandidate(Na, Nb, d, 1);
                AddCandidate(Nb, Na, d, 1);
            }
        }
        while ((Na = Nb) != FirstNode);

        /* Add INPUT_TOUR_FILE edges */
        // removed

        /* Add SUBPROBLEM_TOUR_FILE edges */
        // removed
    }

    /* ResetCandidateSet.c */
    void ResetCandidateSet()
    {
        Node *From;
        Candidate *NFrom, *NN, Temp;

        From = FirstNode;
        /* Loop for all nodes */
        do {
            if (!From->CandidateSet)
                continue;
            /* Reorder the candidate array of From */
            for (NFrom = From->CandidateSet; NFrom->To; NFrom++) {
                Temp = *NFrom;
                for (NN = NFrom - 1;
                     NN >= From->CandidateSet &&
                         (Temp.Alpha < NN->Alpha ||
                          (Temp.Alpha == NN->Alpha && Temp.Cost < NN->Cost)); NN--)
                    *(NN + 1) = *NN;
                *(NN + 1) = Temp;
            }
            NFrom--;
            /* Remove included edges */
            while (NFrom >= From->CandidateSet + 2 && NFrom->Alpha == INT_MAX)
                NFrom--;
            NFrom++;
            NFrom->To = 0;
            /* Remove impossible candidates */
            for (NFrom = From->CandidateSet; NFrom->To; NFrom++) {
                if (!IsPossibleCandidate(From, NFrom->To)) {
                    for (NN = NFrom; NN->To; NN++)
                        *NN = *(NN + 1);
                    NFrom--;
                }
            }
        }
        while ((From = From->Suc) != FirstNode);
    }

    /* TrimCandidateSet.c */
    void TrimCandidateSet(int MaxCandidates)
    {
        Node *From;
        Candidate *NFrom;
        int Count;

        From = FirstNode;
        do {
            Count = 0;
            for (NFrom = From->CandidateSet; NFrom && NFrom->To; NFrom++)
                Count++;
            if (Count > MaxCandidates) {
                From->CandidateSet =
                    (Candidate *) realloc(From->CandidateSet,
                                          (MaxCandidates + 1) * sizeof(Candidate));
                From->CandidateSet[MaxCandidates].To = 0;
            }
        } while ((From = From->Suc) != FirstNode);
    }

    /* SymmetrizeCandidateSet.c */
    void SymmetrizeCandidateSet()
    {
        Node *From, *To;
        Candidate *NFrom;

        From = FirstNode;
        do {
            for (NFrom = From->CandidateSet; NFrom && (To = NFrom->To);
                 NFrom++)
                AddCandidate(To, From, NFrom->Cost, NFrom->Alpha);
        }
        while ((From = From->Suc) != FirstNode);
        ResetCandidateSet();
    }

    /* Random.c */
    #define PRANDMAX INT_MAX

    int a = 0, b = 24, arr[55], initialized = 0;

    unsigned Random()
    {
        int t;

        if (!initialized)
            SRandom(7913);
        if (a-- == 0)
            a = 54;
        if (b-- == 0)
            b = 54;
        if ((t = arr[a] - arr[b]) < 0)
            t += PRANDMAX;
        return (arr[a] = t);
    }

    void SRandom(unsigned Seed)
    {
        int i, ii, last, next;

        Seed %= PRANDMAX;
        arr[0] = last = Seed;
        for (next = i = 1; i < 55; i++) {
            ii = (21 * i) % 55;
            arr[ii] = next;
            if ((next = last - next) < 0)
                next += PRANDMAX;
            last = arr[ii];
        }
        initialized = 1;
        a = 0;
        b = 24;
        for (i = 0; i < 165; i++)
            Random();
    }

    /* SegmentSize.c */
    int SegmentSize(Node * ta, Node * tb)
    {
        int n = !Reversed ? tb->Rank - ta->Rank : ta->Rank - tb->Rank;
        return (n < 0 ? n + Dimension : n) + 1;
    }

    /* Create_POPMUSIC_CandidateSet.c */
#define maxNeighbors POPMUSIC_MaxNeighbors
#define trials POPMUSIC_Trials
#define NB_RES POPMUSIC_Solutions
#define SAMPLE_SIZE POPMUSIC_SampleSize

#define d(a, b) (((a) == (b) ? 0 : D(a, b)              \
                  - (a)->Pi - (b)->Pi) / Precision)
#define less(a, b, dmin) (!c || (c(a, b)                                \
                                 - (a)->Pi - (b)->Pi) / Precision < dmin)

    Node **node;
    Node **node_path;

    void Create_POPMUSIC_CandidateSet(int K)
    {
        int n, i, no_res, setInitialSuc, d;
        int *solution;
        GainType cost, costSum = 0;
        GainType costMin = PLUS_INFINITY, costMax = MINUS_INFINITY;
        Node *N;
        double startTime, entryTime;

        AddTourCandidates();
        // if (MaxCandidates == 0) {
        //     N = FirstNode;
        //     // do {
        //     //     if (!N->CandidateSet)
        //     //         eprintf("MAX_CANDIDATES = 0: No candidates");
        //     // } while ((N = N->Suc) != FirstNode);
        //     return;
        // }

        n = Dimension;
        solution = (int *) malloc((n + 1) * sizeof(int));
        node = (Node **) malloc((n + 1) * sizeof(Node *));
        node_path = (Node **) malloc((n + 1) * sizeof(Node *));

        for (no_res = 1; no_res <= NB_RES; no_res++) {
            n = 0;
            N = FirstNode;
            do {
                node[solution[n] = n] = N;
                n++;
            } while ((N = N->Suc) != FirstNode);
            shuffle(n, solution);
            solution[n] = solution[0];
            node[n] = node[solution[0]];
            build_path(n, solution, SAMPLE_SIZE);
            solution[n] = solution[0];
            node[n] = node[solution[0]];
            cost = length_path(n, solution);
            fast_POPMUSIC(n, solution, SAMPLE_SIZE * SAMPLE_SIZE);
            solution[n] = solution[0];
            node[n] = node[solution[0]];
            cost = length_path(n, solution);
            costSum += cost;
            if (cost > costMax)
                costMax = cost;
            setInitialSuc = 0;
            if (cost < costMin) {
                costMin = cost;
                setInitialSuc = POPMUSIC_InitialTour && !InitialTourFileName;
            }
            for (i = 0; i < n; i++) {
                Node *a = node[solution[i]];
                Node *b = node[solution[i + 1]];
                d = D(a, b);
                AddCandidate(a, b, d, 1);
                AddCandidate(b, a, d, 1);
                if (setInitialSuc)
                    a->InitialSuc = b;
            }
        }
        free(solution);
        free(node);
        free(node_path);
        ResetCandidateSet();
        if (K > 0)
            TrimCandidateSet(K);
        AddTourCandidates();
        if (CandidateSetSymmetric)
            SymmetrizeCandidateSet();
    }

    /************************ Compute the length of a path ************************/
    GainType length_path(int n, int *path)
    {
        Node *a, *b;
        GainType length = 0;
        int i;

        for (i = 1, a = node[path[0]]; i <= n; i++, a = b) {
            b = node[path[i]];
            length += d(a, b);
        }
        return length;
    }

    /*************** Optimize the path between path[0] and path[n] ****************/
    void optimize_path(int n, int *path)
    {
        int i, j;
        int *order, **d;
        GainType length;
        Node *a, *b;

        order = (int *) malloc((n + 1) * sizeof(int));
        d = (int **) malloc((n + 1) * sizeof(int *));
        for (i = 0; i <= n; i++) {
            order[i] = i;
            d[i] = (int *) malloc((n + 1) * sizeof(int));
            node_path[i] = node[path[i]];
        }
        for (i = 0; i < n; i++) {
            a = node_path[i];
            for (j = i + 1; j <= n; j++) {
                b = node_path[j];
                d[i][j] = d[j][i] = 
                    IsPossibleCandidate(a, b) ? d(a, b) : INT_MAX;
            }
        }
        d[0][n] = d[n][0] = 0;
        length = 0;
        for (i = 0; i < n; i++)
            length += d[i][i + 1];
        path_threeOpt(n, d, order, &length);
        for (i = 0; i <= n; i++)
            free(d[i]);
        free(d);
        for (i = 0; i <= n; i++)
            order[i] = path[order[i]];
        for (i = 0; i <= n; i++)
            path[i] = order[i];
        free(order);
    }

    /*********** Build recursively a path between path[0] and path[n] *************/
    void build_path(int n, int *path, int nb_clust)
    {
        int i, j, d, dmin, closest, start, end;
        int *tmp_path, *sample, *assignment, *start_clust, *assigned;

        if (n <= 2)
            return;
        if (n <= nb_clust * nb_clust) {
            optimize_path(n, path);
            return;
        }

        /* Temporary path */
        tmp_path = (int *) malloc((n + 1) * sizeof(int));
        for (i = 0; i <= n; i++)
            tmp_path[i] = path[i];

        /* Let tmp_path[1] be the closest city to path[0] */
        dmin = d(node[tmp_path[1]], node[path[0]]);
        closest = 1;
        for (i = 2; i < n; i++) {
            if (less(node[tmp_path[i]], node[path[0]], dmin) &&
                (d = d(node[tmp_path[i]], node[path[0]])) < dmin) {
                dmin = d;
                closest = i;
            }
        }
        swap(tmp_path + 1, tmp_path + closest);

        /* Let tmp_path[2] be the closest city to path[n] */
        dmin = d(node[tmp_path[2]], node[path[n]]);
        closest = 2;
        for (i = 3; i < n; i++) {
            if (less(node[tmp_path[i]], node[path[n]], dmin) &&
                (d = d(node[tmp_path[i]], node[path[n]])) < dmin) {
                dmin = d;
                closest = i;
            }
        }
        swap(tmp_path + 2, tmp_path + closest);

        /* Choose a sample of nb_clust-2 random cities from tmp_path */
        for (i = 3; i <= nb_clust; i++)
            swap(tmp_path + i, tmp_path + unif(i, n - 1));
        swap(tmp_path + 2, tmp_path + nb_clust);
        sample = (int *) malloc((nb_clust + 2) * sizeof(int));
        for (i = 0; i <= nb_clust; i++)
            sample[i] = tmp_path[i];
        sample[nb_clust + 1] = path[n];
        optimize_path(nb_clust + 1, sample);

        /* Assign each city of path to the closest of the sample */
        assignment = (int *) malloc((n + 1) * sizeof(int));
        for (i = 1; i < n; i++) {
            dmin = INT_MAX;
            for (j = 1; j <= nb_clust; j++) {
                if (path[i] == sample[j]) {
                    closest = j;
                    break;
                }
                if (less(node[path[i]], node[sample[j]], dmin) &&
                    (d = d(node[path[i]], node[sample[j]])) < dmin) {
                    dmin = d;
                    closest = j;
                }
            }
            assignment[i] = closest;
        }

        /* Build clusters: ith cluster has (start_clust[i+1]-start_clust[i])
           cities */
        start_clust = (int *) calloc(nb_clust + 1, sizeof(int));
        for (i = 1; i < n; i++)
            start_clust[assignment[i]]++;
        for (i = 1; i <= nb_clust; i++)
            start_clust[i] += start_clust[i - 1];

        /* Clusters are stored in tmp_path in the order given by sample */
        assigned = (int *) calloc(nb_clust + 1, sizeof(int));
        for (i = 1; i < n; i++)
            tmp_path[start_clust[assignment[i] - 1] +
                     assigned[assignment[i]]++] = path[i];

        /* Reorder original path */
        for (i = 1; i < n; i++)
            path[i] = tmp_path[i - 1];
        free(tmp_path);
        free(sample);
        free(assigned);
        free(assignment);
        for (i = 0; i < nb_clust; i++) {
            start = start_clust[i];
            end = start_clust[i + 1] + 1;
            /* Recursively optimize sub-path corresponding to each cluster */
            build_path(end - start, path + start, nb_clust);
        }
        free(start_clust);
    }

    void reverse(int *path, int i, int j)
    {
        while (i < j)
            swap(path + i++, path + j--);
    }

    void circular_right_shift(int n, int *path, int positions)
    {
        reverse(path, 0, positions - 1);
        reverse(path, positions, n - 1);
        reverse(path, 0, n - 1);
    }

    /* Fast POPMUSIC: Optimize independent subpaths of R cities */
    void fast_POPMUSIC(int n, int *path, int R)
    {
        int scan, i;

        if (R > n)
            R = n;
        /* Optimize subpaths of R cities with R/2 overlap; 2 scans */
        for (scan = 1; scan <= 2; scan++) {
            if (scan == 2) {
                circular_right_shift(n, path, R / 2);
                path[n] = path[0];
                node[n] = node[path[0]];
            }
            for (i = 0; i < n / R; i++)
                optimize_path(R, path + R * i);
            if (n % R != 0) /* Optimize last portion of the path */
                optimize_path(R, path + n - R);
        }
    }

    /* Iterated 3-opt code */

    int n;
    int **dist;
    int *tour, *pos;
    int **neighbor;
    int *neighbors;
    int reversed;
    char *dontLook;
    GainType tourLength;

    void path_threeOpt(int N, int **D, int *best_sol,
                       GainType * best_cost)
    {
        int i, j, trial;
        int *bestTour;
        GainType bestTourLength;

        n = N + 1;
        tour = (int *) malloc(n * sizeof(int));
        pos = (int *) malloc(n * sizeof(int));
        bestTour = (int *) malloc(n * sizeof(int));
        for (i = 0; i < n; i++)
            pos[bestTour[i] = tour[i] = best_sol[i]] = i;
        dist = D;
        createNeighbors();
        dontLook = (char *) calloc(n, sizeof(char));
        bestTourLength = tourLength = *best_cost;
        if (POPMUSIC_Trials == 0)
            trials = n;
        for (trial = 1; trial <= trials; trial++) {
            threeOpt();
            if (tourLength < bestTourLength) {
                for (i = 0; i < n; i++)
                    bestTour[i] = tour[i];
                bestTourLength = tourLength;
            } else {
                for (i = 0; i < n; i++) {
                    pos[tour[i] = bestTour[i]] = i;
                    tourLength = bestTourLength;
                }
            }
            if (n <= 5 || trial == trials)
                break;
            doubleBridgeKick();
        }
        *best_cost = bestTourLength;
        for (i = 0; i < n; i++)
            pos[tour[i] = bestTour[i]] = i;
        reversed = next(0) == N;
        for (i = 0, j = 0; j < n; i = NEXT(i), j++)
            best_sol[j] = i;
        free(tour);
        free(pos);
        free(bestTour);
        free(neighbors);
        for (i = 0; i < n; i++)
            free(neighbor[i]);
        free(neighbor);
        free(dontLook);
    }

    int unif(int low, int high)
    {
        return low + Random() % (high - low + 1);
    }

    void shuffle(int n, int *path)
    {
        int i;
        for (i = 1; i < n; i++)
            swap(path + i, path + unif(0, i));
    }

    void swap(int *a, int *b)
    {
        int tmp = *a;
        *a = *b;
        *b = tmp;
    }

    int fixed(int a, int b)
    {
        return (a == 0 && b == n - 1) || (a == n - 1 && b == 0) ||
            FixedOrCommon(node_path[a], node_path[b]);
    }

    int prev(int v)
    {
        return tour[pos[v] > 0 ? pos[v] - 1 : n - 1];
    }

    int next(int v)
    {
        return tour[pos[v] < n - 1 ? pos[v] + 1 : 0];
    }

    int between(int v1, int v2, int v3)
    {
        int a = pos[v1], b = pos[v2], c = pos[v3];
        return a <= c ? b >= a && b <= c : b <= c || b >= a;
    }

    int PREV(int v)
    {
        return reversed ? next(v) : prev(v);
    }

    int NEXT(int v)
    {
        return reversed ? prev(v) : next(v);
    }

    int _BETWEEN(int v1, int v2, int v3)
    {
        return !reversed ? between(v1, v2, v3) : between(v3, v2, v1);
    }

    void flip(int from, int to)
    {
        int i, j, size, tmp;

        if (from == to)
            return;
        if (reversed) {
            tmp = from;
            from = to;
            to = tmp;
        }
        i = pos[from], j = pos[to];
        size = j - i;
        if (size < 0)
            size += n;
        if (size >= n / 2) {
            tmp = i;
            i = ++j < n ? j : 0;
            j = --tmp >= 0 ? tmp : n - 1;
        }
        while (i != j) {
            tmp = tour[i];
            pos[tour[i] = tour[j]] = i;
            pos[tour[j] = tmp] = j;
            if (++i == n)
                i = 0;
            if (i != j && --j < 0)
                j = n - 1;
        }
    }

    void threeOpt()
    {
        int improved = 1, a, b, c, d, e, f, xa, xc, xe, i, j;
        GainType g0, g1, g2, g3, gain;

        while (improved) {
            improved = 0;
            for (b = 0; b < n; b++) {
                if (dontLook[b])
                    continue;
                dontLook[b] = 1;
                for (xa = 1; xa <= 2; xa++, reversed = !reversed) {
                    a = PREV(b);
                    if (fixed(a, b))
                        continue;
                    g0 = dist[a][b];
                    for (i = 0; i < neighbors[b]; i++) {
                        c = neighbor[b][i];
                        if (c == prev(b) || c == next(b))
                            continue;
                        g1 = g0 - dist[b][c];
                        if (g1 <= 0)
                            break;
                        for (xc = 1; xc <= 2; xc++) {
                            d = xc == 1 ? PREV(c) : NEXT(c);
                            if (d == a || fixed(c, d))
                                continue;
                            g2 = g1 + dist[c][d];
                            if (xc == 1) {
                                gain = g2 - dist[d][a];
                                if (gain > 0) {
                                    flip(b, d);
                                    tourLength -= gain;
                                    dontLook[a] = dontLook[b] = 0;
                                    dontLook[c] = dontLook[d] = 0;
                                    improved = 1;
                                    i = neighbors[b];
                                    break;
                                }
                            }
                            for (j = 0; j < neighbors[d]; j++) {
                                e = neighbor[d][j];
                                if (e == prev(d) || e == next(d) ||
                                    (xc == 2 && !_BETWEEN(b, e, c)))
                                    continue;
                                g3 = g2 - dist[d][e];
                                if (g3 <= 0)
                                    break;
                                for (xe = 1; xe <= xc; xe++) {
                                    if (xc == 1)
                                        f = _BETWEEN(b, e,
                                                     c) ? NEXT(e) : PREV(e);
                                    else
                                        f = xe == 1 ? PREV(e) : NEXT(e);
                                    if (f == a || fixed(e, f))
                                        continue;
                                    gain = g3 + dist[e][f] - dist[f][a];
                                    if (gain > 0) {
                                        if (xc == 1) {
                                            flip(b, d);
                                            if (f == PREV(e))
                                                flip(e, a);
                                            else
                                                flip(a, e);
                                        } else if (xe == 1) {
                                            flip(e, c);
                                            if (b == NEXT(a))
                                                flip(b, f);
                                            else
                                                flip(f, b);
                                        } else {
                                            flip(d, a);
                                            if (f == NEXT(e))
                                                flip(f, c);
                                            else
                                                flip(c, f);
                                            if (b == NEXT(d))
                                                flip(b, e);
                                            else
                                                flip(e, b);
                                        }
                                        tourLength -= gain;
                                        dontLook[a] = dontLook[b] = 0;
                                        dontLook[c] = dontLook[d] = 0;
                                        dontLook[e] = dontLook[f] = 0;
                                        improved = 1;
                                        goto Next_b;
                                    }
                                }
                            }
                        }
                    }
                Next_b:;
                }
            }
        }
    }

    void createNeighbors()
    {
        int i, j, k, d;

        neighbor = (int **) malloc(n * sizeof(int *));
        for (i = 0; i < n; i++)
            neighbor[i] = (int *) malloc((maxNeighbors + 1) * sizeof(int));
        neighbors = (int *) calloc(n, sizeof(int));
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (i == j)
                    continue;
                d = dist[i][j];
                k = neighbors[i] <
                    maxNeighbors ? neighbors[i]++ : maxNeighbors;
                while (k > 0 && d < dist[i][neighbor[i][k - 1]]) {
                    neighbor[i][k] = neighbor[i][k - 1];
                    k--;
                }
                neighbor[i][k] = j;
            }
        }
    }

    int legal(int r, int i, int t[4])
    {
        while (--i >= 0)
            if (r == t[i])
                return 0;
        return !fixed(r, next(r));
    }

    int select_t(int i, int t[4])
    {
        int r, r0;
        r = r0 = unif(0, n - 1);
        while (!legal(r, i, t)) {
            if (++r == n)
                r = 0;
            if (r == r0)
                return -1;
        }
        return r;
    }

    void doubleBridgeKick()
    {
        int t[4], a, b, c, d, e, f, g, h, i;

        reversed = 0;
        for (i = 0; i <= 3; i++) {
            t[i] = select_t(i, t);
            if (t[i] < 0)
                return;
        }
        if (pos[t[0]] > pos[t[1]])
            swap(t, t + 1);
        if (pos[t[2]] > pos[t[3]])
            swap(t + 2, t + 3);
        if (pos[t[0]] > pos[t[2]])
            swap(t, t + 2);
        if (pos[t[1]] > pos[t[3]])
            swap(t + 1, t + 3);
        if (pos[t[1]] > pos[t[2]])
            swap(t + 1, t + 2);
        a = t[0];
        b = next(a);
        c = t[2];
        d = next(c);
        e = t[1];
        f = next(e);
        g = t[3];
        h = next(g);
        flip(b, c);
        if (f == next(e))
            flip(f, h);
        else
            flip(h, f);
        if (b == next(d))
            flip(c, d);
        else
            flip(d, c);
        dontLook[a] = dontLook[b] = dontLook[c] = dontLook[d] = 0;
        dontLook[e] = dontLook[f] = dontLook[g] = dontLook[h] = 0;
        tourLength -= dist[a][b] - dist[b][c] +
            dist[c][d] - dist[d][a] +
            dist[e][f] - dist[f][g] + dist[g][h] - dist[h][e];
    }


    /* PatchCycles */
    int CurrentCycle, Patchwork = 0, RecLevel = 0;
#define MaxPatchwork Dimension

    /*
     * The PatchCycles function tries to find a gainful move by patching the
     * cycles that would occur if the move represented by t[1..2k] and incl[1..2k]
     * was made using one one or more alternating cycles.
     * The alternating cycles are put in continuation of t, starting at 2k+1.
     */

    GainType PatchCycles(int k, GainType Gain)
    {
        Node *s1, *s2, *sStart, *sStop;
        GainType NewGain;
        int M, i;

        FindPermutation(k);
        M = Cycles(k);
        if (M == 1 && Gain > 0) {
            MakeKOptMove(k);
            return Gain;
        }
        if (M == 1 || M > PatchingC || k + M > NonsequentialMoveType)
            return 0;
        if (RecLevel == 0)
            Patchwork = 0;
        CurrentCycle = ShortestCycle(M, k);
        for (i = 0; i < k; i++) {
            if (cycle[p[2 * i]] != CurrentCycle)
                continue;
            sStart = t[p[2 * i]];
            sStop = t[p[2 * i + 1]];
            for (s1 = sStart; s1 != sStop; s1 = s2) {
                s2 = SUC(s1);
                if (FixedOrCommon(s1, s2))
                    continue;
                if (++Patchwork > MaxPatchwork)
                    return 0;
                t[2 * k + 1] = s1;
                t[2 * k + 2] = s2;
                MarkDeleted(s1, s2);
                /* Find a set of gainful alternating cycles */
                NewGain = PatchCyclesRec(k, 2, M, Gain + C(s1, s2));
                UnmarkDeleted(s1, s2);
                if (NewGain > 0)
                    return NewGain;
            }
        }
        return 0;
    }

    GainType PatchCyclesRec(int k, int m, int M, GainType G0)
    {
        Node *s1, *s2, *s3, *s4, *s5, *s6, *S3 = 0, *S4 = 0;
        Candidate *Ns2, *Ns4;
        GainType G1, G2, G3, G4, Gain, CloseUpGain,
            BestCloseUpGain = PatchingAExtended ? MINUS_INFINITY : 0;
        int X4, X6;
        int i, NewCycle, *cycleSaved = 0, *pSaved = 0;
        int Breadth2 = 0, Breadth4;

        s1 = t[2 * k + 1];
        s2 = t[i = 2 * (k + m) - 2];
        incl[incl[i] = i + 1] = i;

        /* Choose (s2,s3) as a candidate edge emanating from s2 */
        for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
            if (s3 == s2->Pred || s3 == s2->Suc || Added(s2, s3) ||
                (NewCycle = Cycle(s3, k)) == CurrentCycle)
                continue;
            if (++Breadth2 > MaxBreadth)
                break;
            MarkAdded(s2, s3);
            t[2 * (k + m) - 1] = s3;
            G1 = G0 - Ns2->Cost;
            /* Choose s4 as one of s3's two neighbors on the tour */
            for (X4 = 1; X4 <= 2; X4++) {
                s4 = X4 == 1 ? s3->Pred : s3->Suc;
                if (FixedOrCommon(s3, s4) || Deleted(s3, s4))
                    continue;
                MarkDeleted(s3, s4);
                t[2 * (k + m)] = s4;
                G2 = G1 + C(s3, s4);
                if (M > 2) {
                    if (!cycleSaved) {
                        cycleSaved = (int *) malloc(2 * k * sizeof(int));
                        memcpy(cycleSaved, cycle + 1, 2 * k * sizeof(int));
                    }
                    for (i = 1; i <= 2 * k; i++)
                        if (cycle[i] == NewCycle)
                            cycle[i] = CurrentCycle;
                    /* Extend the current alternating path */
                    if ((Gain = PatchCyclesRec(k, m + 1, M - 1, G2)) > 0) {
                        UnmarkAdded(s2, s3);
                        UnmarkDeleted(s3, s4);
                        goto End_PatchCyclesRec;
                    }
                    memcpy(cycle + 1, cycleSaved, 2 * k * sizeof(int));
                    if (PatchingA >= 2 && Patchwork < MaxPatchwork &&
                        k + M < NonsequentialMoveType &&
                        !Forbidden(s4, s1) &&
                        (!PatchingARestricted || IsCandidate(s4, s1))) {
                        GainType Bound = BestCloseUpGain >= 0 ||
                            IsCandidate(s4, s1) ? BestCloseUpGain : 0;
                        if ((!c || G2 - c(s4, s1) > Bound) &&
                            (CloseUpGain = G2 - C(s4, s1)) > Bound) {
                            S3 = s3;
                            S4 = s4;
                            BestCloseUpGain = CloseUpGain;
                        }
                    }
                } else if (!Forbidden(s4, s1) && (!c || G2 - c(s4, s1) > 0)
                           && (Gain = G2 - C(s4, s1)) > 0) {
                    incl[incl[2 * k + 1] = 2 * (k + m)] = 2 * k + 1;
                    MakeKOptMove(k + m);
                    UnmarkAdded(s2, s3);
                    UnmarkDeleted(s3, s4);
                    goto End_PatchCyclesRec;
                }
                UnmarkDeleted(s3, s4);
            }
            UnmarkAdded(s2, s3);
        }
        if (M == 2 && !PatchingCRestricted) {
            /* Try to patch the two cycles by a sequential 3-opt move */
            incl[incl[2 * (k + m)] = 2 * (k + m) + 1] = 2 * (k + m);
            incl[incl[2 * k + 1] = 2 * (k + m) + 2] = 2 * k + 1;
            Breadth2 = 0;
            /* Choose (s2,s3) as a candidate edge emanating from s2 */
            for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
                if (s3 == s2->Pred || s3 == s2->Suc || Added(s2, s3))
                    continue;
                if (++Breadth2 > MaxBreadth)
                    break;
                t[2 * (k + m) - 1] = s3;
                G1 = G0 - Ns2->Cost;
                NewCycle = Cycle(s3, k);
                /* Choose s4 as one of s3's two neighbors on the tour */
                for (X4 = 1; X4 <= 2; X4++) {
                    s4 = X4 == 1 ? s3->Pred : s3->Suc;
                    if (FixedOrCommon(s3, s4) || Deleted(s3, s4))
                        continue;
                    t[2 * (k + m)] = s4;
                    G2 = G1 + C(s3, s4);
                    Breadth4 = 0;
                    /* Choose (s4,s5) as a candidate edge emanating from s4 */
                    for (Ns4 = s4->CandidateSet; (s5 = Ns4->To); Ns4++) {
                        if (s5 == s4->Pred || s5 == s4->Suc || s5 == s1 ||
                            Added(s4, s5) ||
                            (NewCycle == CurrentCycle &&
                             Cycle(s5, k) == CurrentCycle))
                            continue;
                        if (++Breadth4 > MaxBreadth)
                            break;
                        G3 = G2 - Ns4->Cost;
                        /* Choose s6 as one of s5's two neighbors on the tour */
                        for (X6 = 1; X6 <= 2; X6++) {
                            s6 = X6 == 1 ? s5->Pred : s5->Suc;
                            if (s6 == s1 || Forbidden(s6, s1)
                                || FixedOrCommon(s5, s6)
                                || Deleted(s5, s6)
                                || Added(s6, s1))
                                continue;
                            G4 = G3 + C(s5, s6);
                            if ((!c || G4 - c(s6, s1) > 0) &&
                                (Gain = G4 - C(s6, s1)) > 0) {
                                if (!pSaved) {
                                    pSaved = (int *) malloc(2 * k * sizeof(int));
                                    memcpy(pSaved, p + 1, 2 * k * sizeof(int));
                                }
                                t[2 * (k + m) + 1] = s5;
                                t[2 * (k + m) + 2] = s6;
                                if (FeasibleKOptMove(k + m + 1)) {
                                    MakeKOptMove(k + m + 1);
                                    goto End_PatchCyclesRec;
                                }
                                memcpy(p + 1, pSaved, 2 * k * sizeof(int));
                                for (i = 1; i <= 2 * k; i++)
                                    q[p[i]] = i;
                            }
                        }
                    }
                }
            }
        }
        Gain = 0;
        if (S4) {
            int OldCycle = CurrentCycle;
            if (!pSaved) {
                pSaved = (int *) malloc(2 * k * sizeof(int));
                memcpy(pSaved, p + 1, 2 * k * sizeof(int));
            }
            t[2 * (k + m) - 1] = S3;
            t[2 * (k + m)] = S4;
            incl[incl[2 * k + 1] = 2 * (k + m)] = 2 * k + 1;
            /* Find a new alternating cycle */
            PatchingA--;
            RecLevel++;
            MarkAdded(s2, S3);
            MarkDeleted(S3, S4);
            MarkAdded(S4, s1);
            Gain = PatchCycles(k + m, BestCloseUpGain);
            UnmarkAdded(s2, S3);
            UnmarkDeleted(S3, S4);
            UnmarkAdded(S4, s1);
            RecLevel--;
            PatchingA++;
            if (Gain <= 0) {
                memcpy(cycle + 1, cycleSaved, 2 * k * sizeof(int));
                memcpy(p + 1, pSaved, 2 * k * sizeof(int));
                for (i = 1; i <= 2 * k; i++)
                    q[p[i]] = i;
                CurrentCycle = OldCycle;
            }
        }

    End_PatchCyclesRec:
        free(cycleSaved);
        free(pSaved);
        return Gain;
    }

    /*
     * The Cycle function returns the number of the cycle containing
     * a given node, N.
     *
     * Time complexity: O(log k).
     */

    int Cycle(Node * N, int k)
    {
        /* Binary search */
        int Low = 1, High = k;
        while (Low < High) {
            int Mid = (Low + High) / 2;
            if (BETWEEN(t[p[2 * Low]], N, t[p[2 * Mid + 1]]))
                High = Mid;
            else
                Low = Mid + 1;
        }
        return cycle[p[2 * Low]];
    }

    /*
     * The ShortestCycle function returns the number of the cycle with
     * the smallest number of nodes. Note however that if the two-level
     * list is used, the number of nodes of each cycle is only approximate
     * (for efficiency reasons).
     *
     * Time complexity: O(k + M), where M = Cycles(k).
     *
     * The function may only be called after a call of the Cycles function.
     */

    int ShortestCycle(int M, int k)
    {
        int i, Cycle, MinCycle = 0;
        int *Size, MinSize = INT_MAX;

        Size = (int *) calloc(1 + M, sizeof(int));
        p[0] = p[2 * k];
        for (i = 0; i < 2 * k; i += 2)
            Size[cycle[p[i]]] += SegmentSize(t[p[i]], t[p[i + 1]]);
        for (Cycle = 1; Cycle <= M; Cycle++) {
            if (Size[Cycle] < MinSize) {
                MinSize = Size[Cycle];
                MinCycle = Cycle;
            }
        }
        free(Size);
        return MinCycle;
    }

    /* RestoreTour.c */
    void RestoreTour()
    {
        Node *t1, *t2, *t3, *t4;

        /* Loop as long as the stack is not empty */
        while (Swaps > 0) {
            /* Undo topmost 2-opt move */
            Swaps--;
            t1 = SwapStack[Swaps].t1;
            t2 = SwapStack[Swaps].t2;
            t3 = SwapStack[Swaps].t3;
            t4 = SwapStack[Swaps].t4;
            Swap1(t3, t2, t1);
            Swaps--;
            /* Make edges (t1,t2) and (t2,t3) excludable again */
            t1->OldPredExcluded = t1->OldSucExcluded = 0;
            t2->OldPredExcluded = t2->OldSucExcluded = 0;
            t3->OldPredExcluded = t3->OldSucExcluded = 0;
            t4->OldPredExcluded = t4->OldSucExcluded = 0;
        }
    }


    /* BestKOptMove.c */
    GainType BestG2;

    Node *BestSubsequentMove(Node * t1, Node * t2, GainType * G0, GainType * Gain)
    {
        K = Swaps == 0 ? MoveType : SubsequentMoveType;
        *Gain = 0;
        t[1] = t1;
        t[2] = t2;
        T[2 * K] = 0;
        BestG2 = MINUS_INFINITY;

        /* 
         * Determine (T[3],T[4], ..., T[2K]) = (t[3],t[4], ..., t[2K])
         * such that
         *
         *     G[2 * K] = *G0 - C(t[2],T[3]) + C(T[3],T[4])
         *                    - C(T[4],T[5]) + C(T[5],T[6])
         *                      ...
         *                    - C(T[2K-3],T[2K-2]) + C(T[2K-1],T[2K])
         *
         * is maximum, and (T[2K-1],T[2K]) has not previously been included.
         * If during this process a legal move with *Gain > 0 is found, then 
         * make the move and exit BestKOptMove immediately.
         */

        MarkDeleted(t1, t2);
        *Gain = BestKOptMoveRec(2, *G0);
        UnmarkDeleted(t1, t2);

        if (*Gain <= 0 && T[2 * K]) {
            int i;
            memcpy(t + 1, T + 1, 2 * K * sizeof(Node *));
            for (i = 2; i < 2 * K; i += 2)
                incl[incl[i] = i + 1] = i;
            incl[incl[1] = 2 * K] = 1;
            MakeKOptMove(K);
            for (i = 1; i < 2 * K; i += 2)
                Exclude(T[i], T[i + 1]);
            *G0 = BestG2;
            return T[2 * K];
        }
        return 0;
    }

    
    Node *BestMove(Node * t1, Node * t2, GainType * G0, GainType * Gain)
    {
        K = Swaps == 0 ? MoveType : SubsequentMoveType;
        *Gain = 0;
        t[1] = t1;
        t[2] = t2;
        T[2 * K] = 0;
        BestG2 = MINUS_INFINITY;

        /* 
         * Determine (T[3],T[4], ..., T[2K]) = (t[3],t[4], ..., t[2K])
         * such that
         *
         *     G[2 * K] = *G0 - C(t[2],T[3]) + C(T[3],T[4])
         *                    - C(T[4],T[5]) + C(T[5],T[6])
         *                      ...
         *                    - C(T[2K-3],T[2K-2]) + C(T[2K-1],T[2K])
         *
         * is maximum, and (T[2K-1],T[2K]) has not previously been included.
         * If during this process a legal move with *Gain > 0 is found, then 
         * make the move and exit BestKOptMove immediately.
         */

        MarkDeleted(t1, t2);
        *Gain = BestKOptMoveRec(2, *G0);
        UnmarkDeleted(t1, t2);

        if (*Gain <= 0 && T[2 * K]) {
            int i;
            memcpy(t + 1, T + 1, 2 * K * sizeof(Node *));
            for (i = 2; i < 2 * K; i += 2)
                incl[incl[i] = i + 1] = i;
            incl[incl[1] = 2 * K] = 1;
            MakeKOptMove(K);
            for (i = 1; i < 2 * K; i += 2)
                Exclude(T[i], T[i + 1]);
            *G0 = BestG2;
            return T[2 * K];
        }
        return 0;
    }

    // Node 
    
    Node *BestKOptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain)
    {
        K = Swaps == 0 ? MoveType : SubsequentMoveType;
        *Gain = 0;
        t[1] = t1;
        t[2] = t2;
        T[2 * K] = 0;
        BestG2 = MINUS_INFINITY;

        /* 
         * Determine (T[3],T[4], ..., T[2K]) = (t[3],t[4], ..., t[2K])
         * such that
         *
         *     G[2 * K] = *G0 - C(t[2],T[3]) + C(T[3],T[4])
         *                    - C(T[4],T[5]) + C(T[5],T[6])
         *                      ...
         *                    - C(T[2K-3],T[2K-2]) + C(T[2K-1],T[2K])
         *
         * is maximum, and (T[2K-1],T[2K]) has not previously been included.
         * If during this process a legal move with *Gain > 0 is found, then 
         * make the move and exit BestKOptMove immediately.
         */

        MarkDeleted(t1, t2);
        *Gain = BestKOptMoveRec(2, *G0);
        UnmarkDeleted(t1, t2);

        if (*Gain <= 0 && T[2 * K]) {
            int i;
            memcpy(t + 1, T + 1, 2 * K * sizeof(Node *));
            for (i = 2; i < 2 * K; i += 2)
                incl[incl[i] = i + 1] = i;
            incl[incl[1] = 2 * K] = 1;
            MakeKOptMove(K);
            for (i = 1; i < 2 * K; i += 2)
                Exclude(T[i], T[i + 1]);
            *G0 = BestG2;
            return T[2 * K];
        }
        return 0;
    }

    GainType BestKOptMoveRec(int k, GainType G0)
    {
        Candidate *Nt2;
        Node *t1, *t2, *t3, *t4;
        GainType G1, G2, G3, Gain;
        int X4, i;
        int Breadth2 = 0;

        t1 = t[1];
        t2 = t[i = 2 * k - 2];
        incl[incl[i] = i + 1] = i;
        incl[incl[1] = i + 2] = 1;
        /* Choose (t2,t3) as a candidate edge emanating from t2 */
        for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
            if (t3 == t2->Pred || t3 == t2->Suc ||
                ((G1 = G0 - Nt2->Cost) <= 0 && GainCriterionUsed &&
                 ProblemType != HCP && ProblemType != HPP) ||
                Added(t2, t3))
                continue;
            if (++Breadth2 > MaxBreadth)
                break;
            MarkAdded(t2, t3);
            t[2 * k - 1] = t3;
            G[2 * k - 2] = G1 + t3->Pi;
            /* Choose t4 as one of t3's two neighbors on the tour */
            for (X4 = 1; X4 <= 2; X4++) {
                t4 = X4 == 1 ? PRED(t3) : SUC(t3);
                if (FixedOrCommon(t3, t4) || Deleted(t3, t4))
                    continue;
                t[2 * k] = t4;
                G2 = G1 + C(t3, t4);
                G3 = MINUS_INFINITY;
                if (t4 != t1 && !Forbidden(t4, t1) && !Added(t4, t1) &&
                    (!c || G2 - c(t4, t1) > 0) &&
                    (G3 = G2 - C(t4, t1)) > 0 && FeasibleKOptMove(k)) {
                    UnmarkAdded(t2, t3);
                    MakeKOptMove(k);
                    return G3;
                }
                if (Backtracking && !Excludable(t3, t4))
                    continue;
                MarkDeleted(t3, t4);
                G[2 * k - 1] = G2 - t4->Pi;
                if (k < K) {
                    if ((Gain = BestKOptMoveRec(k + 1, G2)) > 0) {
                        UnmarkAdded(t2, t3);
                        UnmarkDeleted(t3, t4);
                        return Gain;
                    }
                    incl[incl[1] = 2 * k] = 1;
                }
                if (t4 != t1 && !Forbidden(t4, t1) &&
                    k + 1 < NonsequentialMoveType &&
                    PatchingC >= 2 && PatchingA >= 1 &&
                    (Swaps == 0 || SubsequentPatching)) {
                    if (G3 == MINUS_INFINITY)
                        G3 = G2 - C(t4, t1);
                    if ((PatchingCRestricted ? G3 > 0 && IsCandidate(t4, t1) :
                         PatchingCExtended ? G3 > 0
                         || IsCandidate(t4, t1) : G3 > 0)
                        && (Gain = PatchCycles(k, G3)) > 0) {
                        UnmarkAdded(t2, t3);
                        UnmarkDeleted(t3, t4);
                        return Gain;
                    }
                }
                UnmarkDeleted(t3, t4);
                if (k == K && t4 != t1 && t3 != t1 && G3 <= 0 &&
                    !Added(t4, t1) &&
                    (!GainCriterionUsed || G2 - Precision >= t4->Cost)) {
                    if (!Backtracking || Swaps > 0) {
                        if ((G2 > BestG2 ||
                             (G2 == BestG2 && !Near(t3, t4) &&
                              Near(T[2 * K - 1], T[2 * K]))) &&
                            Swaps < MaxSwaps &&
                            Excludable(t3, t4) && !InInputTour(t3, t4)) {
                            if (RestrictedSearch && K > 2 &&
                                ProblemType != HCP && ProblemType != HPP) {
                                /* Ignore the move if the gain does not vary */
                                G[0] = G[2 * K - 2];
                                G[1] = G[2 * K - 1];
                                for (i = 2 * K - 3; i >= 2; i--)
                                    if (G[i] != G[i % 2])
                                        break;
                                if (i < 2)
                                    continue;
                            }
                            if (FeasibleKOptMove(K)) {
                                BestG2 = G2;
                                memcpy(T + 1, t + 1, 2 * K * sizeof(Node *));
                            }
                        }
                    } else if (MaxSwaps > 0 && FeasibleKOptMove(K)) {
                        Node *SUCt1 = SUC(t1);
                        MakeKOptMove(K);
                        for (i = 1; i < 2 * k; i += 2) {
                            Exclude(t[i], t[i + 1]);
                            UnmarkDeleted(t[i], t[i + 1]);
                        }
                        for (i = 2; i < 2 * k; i += 2)
                            UnmarkAdded(t[i], t[i + 1]);
                        memcpy(tSaved + 1, t + 1, 2 * k * sizeof(Node *));
                        while ((t4 = BestSubsequentMove(t1, t4, &G2, &Gain)));
                        if (Gain > 0) {
                            UnmarkAdded(t2, t3);
                            return Gain;
                        }
                        RestoreTour();
                        K = k;
                        memcpy(t + 1, tSaved + 1, 2 * K * sizeof(Node *));
                        for (i = 1; i < 2 * K - 2; i += 2)
                            MarkDeleted(t[i], t[i + 1]);
                        for (i = 2; i < 2 * K; i += 2)
                            MarkAdded(t[i], t[i + 1]);
                        for (i = 2; i < 2 * K; i += 2)
                            incl[incl[i] = i + 1] = i;
                        incl[incl[1] = 2 * K] = 1;
                        if (SUCt1 != SUC(t1))
                            Reversed ^= 1;
                        T[2 * K] = 0;
                    }
                }
            }
            UnmarkAdded(t2, t3);
            if (t3 == t1)
                continue;
            /* Try to delete an added edge, (_,t3) or (t3,_) */
            for (i = 2 * k - 4; i >= 2; i--) {
                if (t3 == t[i]) {
                    t4 = t[i ^ 1];
                    if (t4 == t1 || Forbidden(t4, t1) || FixedOrCommon(t3, t4)
                        || Added(t4, t1))
                        continue;
                    G2 = G1 + C(t3, t4);
                    if ((!c || G2 - c(t4, t1) > 0)
                        && (Gain = G2 - C(t4, t1)) > 0) {
                        incl[incl[i ^ 1] = 1] = i ^ 1;
                        incl[incl[i] = 2 * k - 2] = i;
                        if (FeasibleKOptMove(k - 1)) {
                            MakeKOptMove(k - 1);
                            return Gain;
                        }
                        incl[incl[i ^ 1] = i] = i ^ 1;
                    }
                }
            }
            incl[1] = 2 * k;
            incl[2 * k - 2] = 2 * k - 1;
        }
        return 0;
    }

    /* MinimumSpanningTree */
    void MinimumSpanningTree(int Sparse)
    {
        Node *Blue;         /* Points to the last node included in the tree */
        Node *NextBlue = 0; /* Points to the provisional next node to be included */
        Node *N;
        Candidate *NBlue;
        int d;

        Blue = N = FirstNode;
        Blue->Dad = 0;      /* The root of the tree has no father */
        if (Sparse && Blue->CandidateSet) {
            /* The graph is sparse */
            /* Insert all nodes in the heap */
            Blue->Loc = 0;  /* A blue node is not in the heap */
            while ((N = N->Suc) != FirstNode) {
                N->Dad = Blue;
                N->Cost = N->Rank = INT_MAX;
                HeapLazyInsert(N);
            }
            /* Update all neighbors to the blue node */
            for (NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++) {
                if (FixedOrCommon(Blue, N)) {
                    N->Dad = Blue;
                    N->Cost = NBlue->Cost + Blue->Pi + N->Pi;
                    N->Rank = INT_MIN;
                    HeapSiftUp(N);
                } else if (!Blue->FixedTo2 && !N->FixedTo2) {
                    N->Dad = Blue;
                    N->Cost = N->Rank = NBlue->Cost + Blue->Pi + N->Pi;
                    HeapSiftUp(N);
                }
            }
            /* Loop as long as there are more nodes to include in the tree */
            while ((NextBlue = HeapDeleteMin())) {
                Follow(NextBlue, Blue);
                Blue = NextBlue;
                /* Update all neighbors to the blue node */
                for (NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++) {
                    if (!N->Loc)
                        continue;
                    if (FixedOrCommon(Blue, N)) {
                        N->Dad = Blue;
                        N->Cost = NBlue->Cost + Blue->Pi + N->Pi;
                        N->Rank = INT_MIN;
                        HeapSiftUp(N);
                    } else if (!Blue->FixedTo2 && !N->FixedTo2 &&
                               (d =
                                NBlue->Cost + Blue->Pi + N->Pi) < N->Cost) {
                        N->Dad = Blue;
                        N->Cost = N->Rank = d;
                        HeapSiftUp(N);
                    }
                }
            }
        } else {
            /* The graph is dense */
            while ((N = N->Suc) != FirstNode)
                N->Cost = INT_MAX;
            /* Loop as long as there a more nodes to include in the tree */
            while ((N = Blue->Suc) != FirstNode) {
                int Min = INT_MAX;
                /* Update all non-blue nodes (the successors of Blue in the list) */
                do {
                    if (FixedOrCommon(Blue, N)) {
                        N->Dad = Blue;
                        N->Cost = D(Blue, N);
                        NextBlue = N;
                        Min = INT_MIN;
                    } else {
                        if (!Blue->FixedTo2 && !N->FixedTo2 &&
                            !Forbidden(Blue, N) &&
                            (!c || c(Blue, N) < N->Cost) &&
                            (d = D(Blue, N)) < N->Cost) {
                            N->Cost = d;
                            N->Dad = Blue;
                        }
                        if (N->Cost < Min) {
                            Min = N->Cost;
                            NextBlue = N;
                        }
                    }
                }
                while ((N = N->Suc) != FirstNode);
                Follow(NextBlue, Blue);
                Blue = NextBlue;
            }
        }
    }

    /* Connect.c */
    void Connect(Node * N1, int Max, int Sparse)
    {
        Node *N;
        Candidate *NN1;
        int d;

        N1->Next = 0;
        N1->NextCost = INT_MAX;
        if (!Sparse || N1->CandidateSet == 0 ||
            N1->CandidateSet[0].To == 0 || N1->CandidateSet[1].To == 0) {
            /* Find the requested edge in a dense graph */
            N = FirstNode;
            do {
                if (N == N1 || N == N1->Dad || N1 == N->Dad)
                    continue;
                if (FixedOrCommon(N1, N)) {
                    N1->NextCost = D(N1, N);
                    N1->Next = N;
                    return;
                }
                if (!N1->FixedTo2 && !N->FixedTo2 &&
                    !Forbidden(N1, N) &&
                    (!c || c(N1, N) < N1->NextCost) &&
                    (d = D(N1, N)) < N1->NextCost) {
                    N1->NextCost = d;
                    if (d <= Max)
                        return;
                    N1->Next = N;
                }
            }
            while ((N = N->Suc) != FirstNode);
        } else {
            /* Find the requested edge in a sparse graph */
            for (NN1 = N1->CandidateSet; (N = NN1->To); NN1++) {
                if (N == N1->Dad || N1 == N->Dad)
                    continue;
                if (FixedOrCommon(N1, N)) {
                    N1->NextCost = NN1->Cost + N1->Pi + N->Pi;
                    N1->Next = N;
                    return;
                }
                if (!N1->FixedTo2 && !N->FixedTo2 &&
                    !Forbidden(N1, N) &&
                    (d = NN1->Cost + N1->Pi + N->Pi) < N1->NextCost) {
                    N1->NextCost = d;
                    if (d <= Max)
                        return;
                    N1->Next = N;
                }
            }
        }
    }

    
    /* Minimum1TreeCost.c */
    GainType Minimum1TreeCost(int Sparse)
    {
        Node *N, *N1 = 0;
        GainType Sum = 0;
        int Max = INT_MIN;

        MinimumSpanningTree(Sparse);
        N = FirstNode;
        do {
            N->V = -2;
            Sum += N->Pi;
        }
        while ((N = N->Suc) != FirstNode);
        Sum *= -2;
        while ((N = N->Suc) != FirstNode) {
            N->V++;
            N->Dad->V++;
            Sum += N->Cost;
            N->Next = 0;
        }
        FirstNode->Dad = FirstNode->Suc;
        FirstNode->Cost = FirstNode->Suc->Cost;
        do {
            if (N->V == -1) {
                Connect(N, Max, Sparse);
                if (N->NextCost > Max && N->Next) {
                    N1 = N;
                    Max = N->NextCost;
                }
            }
        }
        while ((N = N->Suc) != FirstNode);
        assert(N1);
        N1->Next->V++;
        N1->V++;
        Sum += N1->NextCost;
        Norm = 0;
        do
            Norm += N->V * N->V;
        while ((N = N->Suc) != FirstNode);
        if (N1 == FirstNode)
            N1->Suc->Dad = 0;
        else {
            FirstNode->Dad = 0;
            Precede(N1, FirstNode);
            FirstNode = N1;
        }
        if (Norm == 0) {
            for (N = FirstNode->Dad; N; N1 = N, N = N->Dad)
                Follow(N, N1);
            for (N = FirstNode->Suc; N != FirstNode; N = N->Suc) {
                N->Dad = N->Pred;
                N->Cost = D(N, N->Dad);
            }
            FirstNode->Suc->Dad = 0;
        }
        return Sum;
    }

    /* FreeStructures.c */
    void FreeCandidateSets()
    {
        Node *N = FirstNode;
        if (!N)
            return;
        do {
            N->CandidateSet = 0;
            N->BackboneCandidateSet = 0;
        }
        while ((N = N->Suc) != FirstNode);
    }

    /* GenerateCandidates.c */
    void GenerateCandidates(int MaxCandidates, GainType MaxAlpha,
                            int Symmetric)
    {
        Node *From, *To;
        Candidate *NFrom, *NN;
        int a, d, Count;

        if (MaxAlpha < 0 || MaxAlpha > INT_MAX)
            MaxAlpha = INT_MAX;
        /* Initialize CandidateSet for each node */
        FreeCandidateSets();
        From = FirstNode;
        do
            From->Mark = 0;
        while ((From = From->Suc) != FirstNode);

        if (MaxCandidates > 0) {
            do {
                From->CandidateSet =
                    (Candidate *) malloc((MaxCandidates + 1) *
                                         sizeof(Candidate));
                From->CandidateSet[0].To = 0;
            }
            while ((From = From->Suc) != FirstNode);
        } else {
            AddTourCandidates();
            // do {
            //     if (!From->CandidateSet)
            //         eprintf("MAX_CANDIDATES = 0: No candidates");
            // } while ((From = From->Suc) != FirstNode);
            return;
        }

        /* Loop for each node, From */
        do {
            NFrom = From->CandidateSet;
            if (From != FirstNode) {
                From->Beta = INT_MIN;
                for (To = From; To->Dad != 0; To = To->Dad) {
                    To->Dad->Beta =
                        !FixedOrCommon(To, To->Dad) ?
                        max(To->Beta, To->Cost) : To->Beta;
                    To->Dad->Mark = From;
                }
            }
            Count = 0;
            /* Loop for each node, To */
            To = FirstNode;
            do {
                if (To == From)
                    continue;
                d = c && !FixedOrCommon(From, To) ? c(From, To) : D(From, To);
                if (From == FirstNode)
                    a = To == From->Dad ? 0 : d - From->NextCost;
                else if (To == FirstNode)
                    a = From == To->Dad ? 0 : d - To->NextCost;
                else {
                    if (To->Mark != From)
                        To->Beta =
                            !FixedOrCommon(To, To->Dad) ?
                            max(To->Dad->Beta, To->Cost) : To->Dad->Beta;
                    a = d - To->Beta;
                }
                if (FixedOrCommon(From, To))
                    a = INT_MIN;
                else {
                    if (From->FixedTo2 || To->FixedTo2 || Forbidden(From, To))
                        continue;
                    if (InInputTour(From, To)) {
                        a = 0;
                        if (c)
                            d = D(From, To);
                    } else if (c) {
                        if (a > MaxAlpha ||
                            (Count == MaxCandidates &&
                             (a > (NFrom - 1)->Alpha ||
                              (a == (NFrom - 1)->Alpha
                               && d >= (NFrom - 1)->Cost))))
                            continue;
                        if (To == From->Dad) {
                            d = From->Cost;
                            a = 0;
                        } else if (From == To->Dad) {
                            d = To->Cost;
                            a = 0;
                        } else {
                            a -= d;
                            a += (d = D(From, To));
                        }
                    }
                }
                if (a <= MaxAlpha && IsPossibleCandidate(From, To)) {
                    /* Insert new candidate edge in From->CandidateSet */
                    NN = NFrom;
                    while (--NN >= From->CandidateSet) {
                        if (a > NN->Alpha || (a == NN->Alpha && d >= NN->Cost))
                            break;
                        *(NN + 1) = *NN;
                    }
                    NN++;
                    NN->To = To;
                    NN->Cost = d;
                    NN->Alpha = a;
                    if (Count < MaxCandidates) {
                        Count++;
                        NFrom++;
                    }
                    NFrom->To = 0;
                }
            }
            while ((To = To->Suc) != FirstNode);
        }
        while ((From = From->Suc) != FirstNode);

        AddTourCandidates();
        if (Symmetric)
            SymmetrizeCandidateSet();
    }

    /* Ascent.c */
    GainType Ascent()
    {
        Node *t;
        GainType BestW, W, W0, Alpha, MaxAlpha;
        int T, Period, P, InitialPhase, BestNorm;

    Start:
        /* Initialize Pi and BestPi */
        t = FirstNode;
        do
            t->Pi = t->BestPi = 0;
        while ((t = t->Suc) != FirstNode);

        /* Compute the cost of a minimum 1-tree */
        W = Minimum1TreeCost(CandidateSetType == DELAUNAY ||
                             CandidateSetType == POPMUSIC ||
                             MaxCandidates == 0);

        /* Return this cost 
           if either
           (1) subgradient optimization is not wanted, or
           (2) the norm of the tree (its deviation from a tour) is zero
           (in that case the true optimum has been found).
        */
        if (!Subgradient || !Norm)
            return W;

        if (MaxCandidates > 0) {
            /* Generate symmetric candididate sets for all nodes */
            MaxAlpha = INT_MAX;
            if (Optimum != MINUS_INFINITY
                && (Alpha = Optimum * Precision - W) >= 0)
                MaxAlpha = Alpha;
            GenerateCandidates(AscentCandidates, MaxAlpha, 1);
        }

        /* Set LastV of every node to V (the node's degree in the 1-tree) */
        t = FirstNode;
        do
            t->LastV = t->V;
        while ((t = t->Suc) != FirstNode);

        BestW = W0 = W;
        BestNorm = Norm;
        InitialPhase = 1;
        /* Perform subradient optimization with decreasing period length 
           and decreasing step size */
        for (Period = InitialPeriod, T = InitialStepSize * Precision;
             Period > 0 && T > 0 && Norm != 0; Period /= 2, T /= 2) {
            /* Period and step size are halved at each iteration */
            for (P = 1; T && P <= Period && Norm != 0; P++) {
                /* Adjust the Pi-values */
                t = FirstNode;
                do {
                    if (t->V != 0) {
                        t->Pi += T * (7 * t->V + 3 * t->LastV) / 10;
                        if (t->Pi > INT_MAX / 10)
                            t->Pi = INT_MAX / 10;
                        else if (t->Pi < INT_MIN / 10)
                            t->Pi = INT_MIN / 10;
                    }
                    t->LastV = t->V;
                }
                while ((t = t->Suc) != FirstNode);
                /* Compute a minimum 1-tree in the sparse graph */
                W = Minimum1TreeCost(1);
                /* Test if an improvement has been found */
                if (W > BestW || (W == BestW && Norm < BestNorm)) {
                    /* If the lower bound becomes greater than twice its
                       initial value it is taken as a sign that the graph might be
                       too sparse */
                    if (W - W0 > (W0 >= 0 ? W0 : -W0) && AscentCandidates > 0
                        && AscentCandidates < Dimension) {
                        W = Minimum1TreeCost(CandidateSetType == DELAUNAY ||
                                             CandidateSetType == POPMUSIC ||
                                             MaxCandidates == 0);
                        if (W < W0) {
                            /* Double the number of candidate edges 
                               and start all over again */
                            if ((AscentCandidates *= 2) > Dimension)
                                AscentCandidates = Dimension;
                            goto Start;
                        }
                        W0 = W;
                    }
                    BestW = W;
                    BestNorm = Norm;
                    /* Update the BestPi-values */
                    t = FirstNode;
                    do
                        t->BestPi = t->Pi;
                    while ((t = t->Suc) != FirstNode);
                    /* If in the initial phase, the step size is doubled */
                    if (InitialPhase && T * sqrt((double) Norm) > 0)
                        T *= 2;
                    /* If the improvement was found at the last iteration of the 
                       current period, then double the period */
                    if (CandidateSetType != DELAUNAY &&
                        CandidateSetType != POPMUSIC &&
                        P == Period && (Period *= 2) > InitialPeriod)
                        Period = InitialPeriod;
                } else {
                    if (InitialPhase && P > Period / 2) {
                        /* Conclude the initial phase */
                        InitialPhase = 0;
                        P = 0;
                        T = 3 * T / 4;
                    }
                }
            }
        }

        t = FirstNode;
        do {
            t->Pi = t->BestPi;
            t->BestPi = 0;
        } while ((t = t->Suc) != FirstNode);

        /* Compute a minimum 1-tree */
        W = BestW = Minimum1TreeCost(CandidateSetType == DELAUNAY ||
                                     CandidateSetType == POPMUSIC ||
                                     MaxCandidates == 0);

        if (MaxCandidates > 0 && CandidateSetType != POPMUSIC) {
            FreeCandidateSets();
        } else {
            Candidate *Nt;
            t = FirstNode;
            do {
                for (Nt = t->CandidateSet; Nt && Nt->To; Nt++)
                    Nt->Cost += t->Pi + Nt->To->Pi;
            }
            while ((t = t->Suc) != FirstNode);
        }
        return W;
    }

    /* KSwapKick.c */
    static int compare(const void *Na, const void *Nb)
    {
        return (*(Node **) Na)->Rank - (*(Node **) Nb)->Rank;
    }

    void KSwapKick(int K)
    {
        Node **s, *N;
        int Count, i;

        s = (Node **) malloc(K * sizeof(Node *));
        Count = 0;
        N = FirstNode;
        do {
            N->Rank = ++Count;
            N->V = 0;
        } while ((N = N->Suc) != FirstNode);
        N = s[0] = RandomNode();
        if (!N)
            goto End_KSwapKick;
        N->V = 1;
        for (i = 1; i < K; i++) {
            N = s[i] = RandomNode();
            if (!N)
                K = i;
            else
                N->V = 1;
        }
        if (K < 4)
            goto End_KSwapKick;
        qsort(s, K, sizeof(Node *), compare);
        for (i = 0; i < K; i++)
            s[i]->OldSuc = s[i]->Suc;
        for (i = 0; i < K; i++)
            Link(s[(i + 2) % K], s[i]->OldSuc);
    End_KSwapKick:
        free(s);
    }

    Node *RandomNode()
    {
        Node *N;
        int Count;

        if (Dimension == DimensionSaved)
            N = &NodeSet[1 + Random() % Dimension];
        else {
            N = FirstNode;
            for (Count = Random() % Dimension; Count > 0; Count--)
                N = N->Suc;
        }
        Count = 0;
        while ((N->V || FixedOrCommon(N, N->Suc)) && Count < Dimension) {
            N = N->Suc;
            Count++;
        }
        return Count < Dimension ? N : 0;
    }

    /* RecordBetterTour.c */
    void RecordBetterTour()
    {
        Node *N = FirstNode, *Stop = N;

        if (ProblemType != ATSP) {
            int i = 1;
            do
                BetterTour[i++] = N->Id;
            while ((N = N->Suc) != Stop);
        }
        BetterTour[0] = BetterTour[DimensionSaved];
        N = FirstNode;
        do {
            N->NextBestSuc = N->BestSuc;
            N->BestSuc = N->Suc;
        }
        while ((N = N->Suc) != FirstNode);
    }

    /* RecordBestTour.c */
    void RecordBestTour()
    {
        for (int i = 0; i <= DimensionSaved; i++)
            BestTour[i] = BetterTour[i];
    }

    /* StoreTour.c */
    void StoreTour()
    {
        Node *t, *u;
        Candidate *Nt;
        int i;

        while (Swaps > 0) {
            Swaps--;
            for (i = 1; i <= 4; i++) {
                t = i == 1 ? SwapStack[Swaps].t1 :
                    i == 2 ? SwapStack[Swaps].t2 :
                    i == 3 ? SwapStack[Swaps].t3 : SwapStack[Swaps].t4;
                Activate(t);
                t->OldPred = t->Pred;
                t->OldSuc = t->Suc;
                t->OldPredExcluded = t->OldSucExcluded = 0;
                t->Cost = INT_MAX;
                for (Nt = t->CandidateSet; (u = Nt->To); Nt++)
                    if (u != t->Pred && u != t->Suc && Nt->Cost < t->Cost)
                        t->Cost = Nt->Cost;
            }
        }
    }


    /* ChooseInitialTour.c */
    void ChooseInitialTour()
    {
        Node *N, *NextN, *FirstAlternative, *Last;
        Candidate *NN;
        int Alternatives, Count, i;

        if (KickType > 0 && Kicks > 0 && Trial > 1) {
            for (Last = FirstNode; (N = Last->BestSuc) != FirstNode; Last = N)
                Follow(N, Last);
            for (i = 1; i <= Kicks; i++)
                KSwapKick(KickType);
            return;
        }

    Start:
        /* Mark all nodes as "not chosen" by setting their V field to zero */
        N = FirstNode;
        do
            N->V = 0;
        while ((N = N->Suc) != FirstNode);
        Count = 0;

        /* Choose FirstNode without two incident fixed or common candidate edges */
        do {
            if (FixedOrCommonCandidates(N) < 2)
                break;
        }
        while ((N = N->Suc) != FirstNode);
        FirstNode = N;

        /* Move nodes with two incident fixed or common candidate edges in
           front of FirstNode */
        for (Last = FirstNode->Pred; N != Last; N = NextN) {
            NextN = N->Suc;
            if (FixedOrCommonCandidates(N) == 2)
                Follow(N, Last);
        }

        /* Mark FirstNode as chosen */
        FirstNode->V = 1;
        N = FirstNode;

        /* Loop as long as not all nodes have been chosen */
        while (N->Suc != FirstNode) {
            FirstAlternative = 0;
            Alternatives = 0;
            Count++;

            /* Case A */
            for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
                if (!NextN->V && Fixed(N, NextN)) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }

            // removed, mergetourfiles

            if (Alternatives == 0 && FirstNode->InitialSuc && Trial == 1 &&
                Count <= InitialTourFraction * Dimension) {
                /* Case B */
                for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
                    if (!NextN->V && InInitialTour(N, NextN)) {
                        Alternatives++;
                        NextN->Next = FirstAlternative;
                        FirstAlternative = NextN;
                    }
                }
            }
            if (Alternatives == 0 && Trial > 1 &&
                ProblemType != HCP && ProblemType != HPP) {
                /* Case C */
                for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
                    if (!NextN->V && FixedOrCommonCandidates(NextN) < 2 &&
                        NN->Alpha == 0 && (InBestTour(N, NextN) ||
                                           InNextBestTour(N, NextN))) {
                        Alternatives++;
                        NextN->Next = FirstAlternative;
                        FirstAlternative = NextN;
                    }
                }
            }
            if (Alternatives == 0) {
                /* Case D */
                for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
                    if (!NextN->V && FixedOrCommonCandidates(NextN) < 2) {
                        Alternatives++;
                        NextN->Next = FirstAlternative;
                        FirstAlternative = NextN;
                    }
                }
            }
            if (Alternatives == 0) {
                /* Case E (actually not really a random choice) */
                NextN = N->Suc;
                while ((FixedOrCommonCandidates(NextN) == 2 ||
                        Forbidden(N, NextN)) && NextN->Suc != FirstNode)
                    NextN = NextN->Suc;
                if (FixedOrCommonCandidates(NextN) == 2 || Forbidden(N, NextN)) {
                    FirstNode = N;
                    goto Start;
                }
            } else {
                NextN = FirstAlternative;
                if (Alternatives > 1) {
                    /* Select NextN at random among the alternatives */
                    i = Random() % Alternatives;
                    while (i--)
                        NextN = NextN->Next;
                }
            }
            /* Include NextN as the successor of N */
            Follow(NextN, N);
            N = NextN;
            N->V = 1;
        }
        if (Forbidden(N, N->Suc)) {
            FirstNode = N;
            goto Start;
        }
        if (MaxTrials == 0) {
            GainType Cost = 0;
            N = FirstNode;
            do
                Cost += C(N, N->Suc) - N->Pi - N->Suc->Pi;
            while ((N = N->Suc) != FirstNode);
            Cost /= Precision;
            if (Cost < BetterCost) {
                BetterCost = Cost;
                RecordBetterTour();
            }
        }
    }

    /* Activate.c */
    void Activate(Node * N)
    {
        if (N->Next != 0)
            return;
        if (FirstActive == 0)
            FirstActive = LastActive = N;
        else
            LastActive = LastActive->Next = N;
        LastActive->Next = FirstActive;
    }

    /* RemoveFirstActive.c */
    Node *RemoveFirstActive()
    {
        Node *N = FirstActive;
        if (FirstActive == LastActive)
            FirstActive = LastActive = 0;
        else
            LastActive->Next = FirstActive = FirstActive->Next;
        if (N)
            N->Next = 0;
        return N;
    }

    /* BridgeGain.c */
    GainType
    BridgeGain(Node * s1, Node * s2, Node * s3, Node * s4,
               Node * s5, Node * s6, Node * s7, Node * s8, int Case6,
               GainType G)
    {
        Node *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *u2 = 0, *u3 = 0;
        Candidate *Nt2, *Nt4, *Nt6;
        GainType G0, G1, G2, G3, G4, G5, G6, Gain;
        int X4;
        int Breadth2, Breadth4, Breadth6;

        /* From the original tour select a segment (u2 --> u3) which contains 
           as few nodes as possible */
        switch (Case6) {
        case 3:
            if (2 * SegmentSize(s5, s4) <= Dimension) {
                u2 = s5;
                u3 = s4;
            } else {
                u2 = s3;
                u3 = s6;
            }
            break;
        case 4:
            if (2 * SegmentSize(s2, s5) <= Dimension) {
                u2 = s2;
                u3 = s5;
            } else {
                u2 = s6;
                u3 = s1;
            }
            break;
        case 0:
        case 7:
            if (2 * SegmentSize(s2, s3) <= Dimension) {
                u2 = s2;
                u3 = s3;
            } else {
                u2 = s4;
                u3 = s1;
            }
        }

        /* Choose t1 between u2 and u3 */
        for (t1 = u2; t1 != u3; t1 = t2) {
            /* Choose t2 as the successor of t1 */
            t2 = SUC(t1);
            if ((t1 == s1 && t2 == s2) ||
                (t1 == s2 && t2 == s1) ||
                (t1 == s3 && t2 == s4) ||
                (t1 == s4 && t2 == s3) ||
                (t1 == s5 && t2 == s6) ||
                (t1 == s6 && t2 == s5) ||
                (t1 == s7 && t2 == s8) ||
                (t1 == s8 && t2 == s7) || FixedOrCommon(t1, t2))
                continue;
            G0 = G + C(t1, t2);
            /* Choose (t2,t3) as a candidate edge emanating from t2. 
               t3 must not be between u2 and u3 */
            Breadth2 = 0;
            for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
                if (t3 == t2->Pred || t3 == t2->Suc || BETWEEN(u2, t3, u3))
                    continue;
                G1 = G0 - Nt2->Cost;
                if (++Breadth2 > MaxBreadth)
                    break;
                /* Choose t4 as one of t3's two neighbors on the tour */
                for (X4 = 1; X4 <= 2; X4++) {
                    t4 = X4 == 1 ? SUC(t3) : PRED(t3);
                    if (t4 == t2 ||
                        (t3 == s1 && t4 == s2) ||
                        (t3 == s2 && t4 == s1) ||
                        (t3 == s3 && t4 == s4) ||
                        (t3 == s4 && t4 == s3) ||
                        (t3 == s5 && t4 == s6) ||
                        (t3 == s6 && t4 == s5) ||
                        (t3 == s7 && t4 == s8) ||
                        (t3 == s8 && t4 == s7) || FixedOrCommon(t3, t4))
                        continue;
                    G2 = G1 + C(t3, t4);
                    /* Test if an improvement can be obtained */
                    if (!Forbidden(t4, t1) &&
                        (!c || G2 - c(t4, t1) > 0) &&
                        (t4 != s2 || t1 != s3) &&
                        (t4 != s3 || t1 != s2) &&
                        (t4 != s4 || t1 != s5) &&
                        (t4 != s5 || t1 != s4) &&
                        (t4 != s6 || t1 != s7) &&
                        (t4 != s7 || t1 != s6) &&
                        (t4 != s8 || t1 != s1) &&
                        (t4 != s1 || t1 != s8) &&
                        (s8 ||
                         ((t4 != s6 || t1 != s1) &&
                          (t4 != s1 || t1 != s6))) &&
                        (Gain = G2 - C(t4, t1)) > 0) {
                        switch (Case6) {
                        case 0:
                            if (X4 == 1)
                                Swap3(s1, s2, s4, t3, t4, t1, s1, s3, s2);
                            else
                                Swap2(t1, t2, t3, s1, s2, s3);
                            return Gain;
                        case 3:
                            if ((X4 == 1) ==
                                (!BETWEEN(s2, t1, s6) && !BETWEEN(s2, t3, s6)))
                                Swap3(s1, s2, s3, t1, t2, t3, s5, s6, s1);
                            else
                                Swap4(s1, s2, s3, t1, t2, t4, s5, s6, s1, t2,
                                      t4, t1);
                            if (s8)
                                Swap1(s7, s8, s1);
                            return Gain;
                        case 4:
                            if ((X4 == 1) ==
                                (!BETWEEN(s3, t1, s5) && !BETWEEN(s3, t3, s5)))
                                Swap3(s1, s2, s3, t1, t2, t3, s5, s6, s1);
                            else
                                Swap4(s1, s2, s3, t1, t2, t4, s5, s6, s1, t2,
                                      t4, t1);
                            if (s8)
                                Swap1(s7, s8, s1);
                            return Gain;
                        case 7:
                            if ((X4 == 1) ==
                                (!BETWEEN(s4, t1, s6) && !BETWEEN(s4, t3, s6)))
                                Swap3(s5, s6, s1, t1, t2, t3, s3, s4, s5);
                            else
                                Swap4(s5, s6, s1, t1, t2, t4, s3, s4, s5, t2,
                                      t4, t1);
                            if (s8)
                                Swap1(s7, s8, s1);
                            return Gain;
                        }
                    }
                    /* If BridgeGain has been called with a nonfeasible 2-opt move,
                       then try to find a 3-opt or 4-opt move which, when composed 
                       with the 2-opt move, results in an improvement of the tour */
                    if (Case6 != 0)
                        continue;
                    Breadth4 = 0;
                    /* Choose (t4,t5) as a candidate edge emanating from t4 */
                    for (Nt4 = t4->CandidateSet; (t5 = Nt4->To); Nt4++) {
                        if (t5 == t4->Pred || t5 == t4->Suc || t5 == t1
                            || t5 == t2)
                            continue;
                        /* Choose t6 as one of t5's two neighbors on the tour.
                           Only one choice! */
                        t6 = X4 == 1
                            || BETWEEN(u2, t5, u3) ? PRED(t5) : SUC(t5);
                        if ((t5 == s1 && t6 == s2) || (t5 == s2 && t6 == s1)
                            || (t5 == s3 && t6 == s4) || (t5 == s4 && t6 == s3)
                            || FixedOrCommon(t5, t6))
                            continue;
                        G3 = G2 - Nt4->Cost;
                        G4 = G3 + C(t5, t6);
                        if (!Forbidden(t6, t1) &&
                            (!c || G4 - c(t6, t1) > 0) &&
                            (Gain = G4 - C(t6, t1)) > 0) {
                            if (X4 == 1)
                                Swap4(s1, s2, s4, t3, t4, t1, s1, s3, s2, t5,
                                      t6, t1);
                            else
                                Swap3(t1, t2, t3, s1, s2, s3, t5, t6, t1);
                            return Gain;
                        }
                        if (++Breadth4 > MaxBreadth)
                            break;
                        Breadth6 = 0;
                        /* Choose (t7,t8) as a candidate edge emanating from t7.
                           Only one choice! */
                        for (Nt6 = t6->CandidateSet; (t7 = Nt6->To); Nt6++) {
                            if (t7 == t6->Pred || t7 == t6->Suc)
                                continue;
                            /* Choose t8 as one of t7's two neighbors on the tour.
                               Only one choice! */
                            if (X4 == 1)
                                t8 = (BETWEEN(u2, t5, t1) ? BETWEEN(t5, t7, t1)
                                      : BETWEEN(t2, t5, u3) ? BETWEEN(u2, t7,
                                                                      t1)
                                      || BETWEEN(t5, t7, u3) : BETWEEN(SUC(u3),
                                                                       t5,
                                                                       t3) ?
                                      BETWEEN(u2, t7, u3)
                                      || BETWEEN(t5, t7, t3) : !BETWEEN(t4, t7,
                                                                        t6)) ?
                                    PRED(t7) : SUC(t7);
                            else
                                t8 = (BETWEEN(u2, t5, t1) ?
                                      !BETWEEN(u2, t7, t6)
                                      && !BETWEEN(t2, t7, u3) : BETWEEN(t2, t5,
                                                                        u3) ?
                                      !BETWEEN(t2, t7, t6) : BETWEEN(SUC(u3),
                                                                     t5,
                                                                     t4) ?
                                      !BETWEEN(SUC(u3), t7, t5)
                                      && !BETWEEN(t3, t7,
                                                  PRED(u2)) : !BETWEEN(t3, t7,
                                                                       t5)) ?
                                    PRED(t7) : SUC(t7);
                            if (t8 == t1
                                || (t7 == t1 && t8 == t2) || (t7 == t3
                                                              && t8 == t4)
                                || (t7 == t4 && t8 == t3) || (t7 == s1
                                                              && t8 == s2)
                                || (t7 == s2 && t8 == s1) || (t7 == s3
                                                              && t8 == s4)
                                || (t7 == s4 && t8 == s3))
                                continue;
                            if (FixedOrCommon(t7, t8) || Forbidden(t8, t1))
                                continue;
                            G5 = G4 - Nt6->Cost;
                            G6 = G5 + C(t7, t8);
                            /* Test if an improvement can be achieved */
                            if ((!c || G6 - c(t8, t1) > 0) &&
                                (Gain = G6 - C(t8, t1)) > 0) {
                                if (X4 == 1)
                                    Swap4(s1, s2, s4, t3, t4, t1, s1, s3, s2,
                                          t5, t6, t1);
                                else
                                    Swap3(t1, t2, t3, s1, s2, s3, t5, t6, t1);
                                Swap1(t7, t8, t1);
                                return Gain;
                            }
                            if (++Breadth6 > MaxBreadth)
                                break;
                        }
                    }
                }
            }
        }
        /* No improvement has been found */
        return 0;
    }

    /* Gain23.c */
    GainType Gain23()
    {
        static Node *s1 = 0;
        static short OldReversed = 0;
        Node *s2, *s3, *s4, *s5, *s6 = 0, *s7, *s8 = 0, *s1Stop;
        Candidate *Ns2, *Ns4, *Ns6;
        GainType G0, G1, G2, G3, G4, G5, G6, Gain, Gain6;
        int X2, X4, X6, X8, Case6 = 0, Case8 = 0;
        int Breadth2, Breadth4, Breadth6;
    
        if (!s1 || s1->Subproblem != FirstNode->Subproblem)
            s1 = FirstNode;
        s1Stop = s1;
        for (X2 = 1; X2 <= 2; X2++) {
            Reversed = X2 == 1 ? OldReversed : (OldReversed ^= 1);
            do {
                s2 = SUC(s1);
                if (FixedOrCommon(s1, s2))
                    continue;
                G0 = C(s1, s2);
                Breadth2 = 0;
                /* Choose (s2,s3) as a candidate edge emanating from s2 */
                for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
                    if (s3 == s2->Pred || s3 == s2->Suc)
                        continue;
                    if (++Breadth2 > MaxBreadth)
                        break;
                    G1 = G0 - Ns2->Cost;
                    for (X4 = 1; X4 <= 2; X4++) {
                        s4 = X4 == 1 ? SUC(s3) : PRED(s3);
                        if (FixedOrCommon(s3, s4))
                            continue;
                        G2 = G1 + C(s3, s4);
                        /* Try any gainful nonfeasible 2-opt move
                           followed by a 2-, 3- or 4-opt move */
                        if (X4 == 1 && s4 != s1 && !Forbidden(s4, s1) &&
                            2 * SegmentSize(s2, s3) <= Dimension &&
                            (!c || G2 - c(s4, s1) > 0) &&
                            (G3 = G2 - C(s4, s1)) > 0 &&
                            (Gain = BridgeGain(s1, s2, s3, s4, 0, 0, 0, 0, 0,
                                               G3)) > 0)
                            return Gain;
                        if (X4 == 2 &&
                            !Forbidden(s4, s1) &&
                            (!c || G2 - c(s4, s1) > 0) &&
                            (Gain = G2 - C(s4, s1)) > 0) {
                            Swap1(s1, s2, s3);
                            return Gain;
                        }
                        if (G2 - s4->Cost <= 0)
                            continue;
                        Breadth4 = 0;
                        /* Try any gainful nonfeasible 3- or 4-opt move
                           folllowed by a 2-opt move */
                        /* Choose (s4,s5) as a candidate edge emanating from s4 */
                        for (Ns4 = s4->CandidateSet; (s5 = Ns4->To); Ns4++) {
                            if (s5 == s4->Pred || s5 == s4->Suc ||
                                (G3 = G2 - Ns4->Cost) <= 0)
                                continue;
                            if (++Breadth4 > MaxBreadth)
                                break;
                            /* Choose s6 as one of s5's two neighbors on the tour */
                            for (X6 = 1; X6 <= 2; X6++) {
                                if (X4 == 2) {
                                    if (X6 == 1) {
                                        Case6 = 1 + !BETWEEN(s2, s5, s4);
                                        s6 = Case6 == 1 ? SUC(s5) : PRED(s5);
                                    } else {
                                        s6 = s6 ==
                                            s5->Pred ? s5->Suc : s5->Pred;
                                        if (s5 == s1 || s6 == s1)
                                            continue;
                                        Case6 += 2;
                                    }
                                } else if (BETWEEN(s2, s5, s3)) {
                                    Case6 = 4 + X6;
                                    s6 = X6 == 1 ? SUC(s5) : PRED(s5);
                                    if (s6 == s1)
                                        continue;
                                } else {
                                    if (X6 == 2)
                                        break;
                                    Case6 = 7;
                                    s6 = PRED(s5);
                                }
                                if (FixedOrCommon(s5, s6))
                                    continue;
                                G4 = G3 + C(s5, s6);
                                Gain6 = 0;
                                if (!Forbidden(s6, s1) &&
                                    (!c || G4 - c(s6, s1) > 0) &&
                                    (Gain6 = G4 - C(s6, s1)) > 0) {
                                    if (Case6 <= 2 || Case6 == 5 || Case6 == 6) {
                                        Make3OptMove(s1, s2, s3, s4, s5, s6,
                                                     Case6);
                                        return Gain6;
                                    }
                                    if ((Gain =
                                         BridgeGain(s1, s2, s3, s4, s5, s6, 0,
                                                    0, Case6, Gain6)) > 0)
                                        return Gain;
                                }
                                Breadth6 = 0;
                                /* Choose (s6,s7) as a candidate edge
                                   emanating from s6 */
                                for (Ns6 = s6->CandidateSet; (s7 = Ns6->To);
                                     Ns6++) {
                                    if (s7 == s6->Pred || s7 == s6->Suc
                                        || (s6 == s2 && s7 == s3) || (s6 == s3
                                                                      && s7 ==
                                                                      s2)
                                        || (G5 = G4 - Ns6->Cost) <= 0)
                                        continue;
                                    if (++Breadth6 > MaxBreadth)
                                        break;
                                    /* Choose s8 as one of s7's two neighbors
                                       on the tour */
                                    for (X8 = 1; X8 <= 2; X8++) {
                                        if (X8 == 1) {
                                            Case8 = Case6;
                                            switch (Case6) {
                                            case 1:
                                                s8 = BETWEEN(s2, s7,
                                                             s5) ? SUC(s7) :
                                                    PRED(s7);
                                                break;
                                            case 2:
                                                s8 = BETWEEN(s3, s7,
                                                             s6) ? SUC(s7) :
                                                    PRED(s7);
                                                break;
                                            case 3:
                                                if (BETWEEN(s5, s7, s4))
                                                    s8 = SUC(s7);
                                                else {
                                                    s8 = BETWEEN(s3, s7,
                                                                 s1) ? PRED(s7)
                                                        : SUC(s7);
                                                    Case8 = 17;
                                                }
                                                break;
                                            case 4:
                                                if (BETWEEN(s2, s7, s5))
                                                    s8 = BETWEEN(s2, s7,
                                                                 s4) ? SUC(s7)
                                                        : PRED(s7);
                                                else {
                                                    s8 = PRED(s7);
                                                    Case8 = 18;
                                                }
                                                break;
                                            case 5:
                                                s8 = PRED(s7);
                                                break;
                                            case 6:
                                                s8 = BETWEEN(s2, s7,
                                                             s3) ? SUC(s7) :
                                                    PRED(s7);
                                                break;
                                            case 7:
                                                if (BETWEEN(s2, s7, s3))
                                                    s8 = SUC(s7);
                                                else {
                                                    s8 = BETWEEN(s5, s7,
                                                                 s1) ? PRED(s7)
                                                        : SUC(s7);
                                                    Case8 = 19;
                                                }
                                            }
                                        } else {
                                            if (Case8 >= 17 ||
                                                (Case6 != 3 && Case6 != 4
                                                 && Case6 != 7))
                                                break;
                                            s8 = s8 ==
                                                s7->Pred ? s7->Suc : s7->Pred;
                                            Case8 += 8;
                                        }
                                        if (s8 == s1 ||
                                            (s7 == s1 && s8 == s2) ||
                                            (s7 == s3 && s8 == s4) ||
                                            (s7 == s4 && s8 == s3))
                                            continue;
                                        if (FixedOrCommon(s7, s8)
                                            || Forbidden(s8, s1))
                                            continue;
                                        G6 = G5 + C(s7, s8);
                                        if ((!c || G6 - c(s8, s1) > 0) &&
                                            (Gain = G6 - C(s8, s1)) > 0) {
                                            if (Case8 <= 15) {
                                                Make4OptMove(s1, s2, s3, s4,
                                                             s5, s6, s7, s8,
                                                             Case8);
                                                return Gain;
                                            }
                                            if (Gain > Gain6 &&
                                                (Gain =
                                                 BridgeGain(s1, s2, s3, s4, s5,
                                                            s6, s7, s8, Case6,
                                                            Gain)) > 0)
                                                return Gain;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            while ((s1 = s2) != s1Stop);
        }
        return 0;
    }

    /* LinKernighan.c */
    GainType LinKernighan()
    {
        Node *t1, *t2, *SUCt1;
        GainType Gain, G0, Cost;
        int X2, i, it = 0;
        Candidate *Nt1;
        Segment *S;
        SSegment *SS;

        Reversed = 0;
        S = FirstSegment;
        i = 0;
        do {
            S->Size = 0;
            S->Rank = ++i;
            S->Reversed = 0;
            S->First = S->Last = 0;
        }
        while ((S = S->Suc) != FirstSegment);
        SS = FirstSSegment;
        i = 0;
        do {
            SS->Size = 0;
            SS->Rank = ++i;
            SS->Reversed = 0;
            SS->First = SS->Last = 0;
        }
        while ((SS = SS->Suc) != FirstSSegment);

        FirstActive = LastActive = 0;
        Swaps = 0;

        /* Compute the cost of the initial tour, Cost.
           Compute the corresponding hash value, Hash.
           Initialize the segment list.
           Make all nodes "active" (so that they can be used as t1). */
        Cost = 0;
        Hash = 0;
        i = 0;
        t1 = FirstNode;
        do {
            t2 = t1->OldSuc = t1->Suc;
            t1->OldPred = t1->Pred;
            t1->Rank = ++i;
            Cost += (t1->SucCost = t2->PredCost = C(t1, t2)) - t1->Pi - t2->Pi;
            Hash ^= Rand[t1->Id] * Rand[t2->Id];
            t1->Cost = INT_MAX;
            for (Nt1 = t1->CandidateSet; (t2 = Nt1->To); Nt1++)
                if (t2 != t1->Pred && t2 != t1->Suc && Nt1->Cost < t1->Cost)
                    t1->Cost = Nt1->Cost;
            t1->Parent = S;
            S->Size++;
            if (S->Size == 1)
                S->First = t1;
            S->Last = t1;
            if (SS->Size == 0)
                SS->First = S;
            S->Parent = SS;
            SS->Last = S;
            if (S->Size == GroupSize) {
                S = S->Suc;
                SS->Size++;
                if (SS->Size == SGroupSize)
                    SS = SS->Suc;
            }
            t1->OldPredExcluded = t1->OldSucExcluded = 0;
            t1->Next = 0;
            if (KickType == 0 || Kicks == 0 || Trial == 1 ||
                !InBestTour(t1, t1->Pred) || !InBestTour(t1, t1->Suc))
                Activate(t1);
        }
        while ((t1 = t1->Suc) != FirstNode);
        if (S->Size < GroupSize)
            SS->Size++;
        Cost /= Precision;
        PredSucCostAvailable = 1;

        /* Loop as long as improvements are found */
        do {
            /* Choose t1 as the first "active" node */
            while ((t1 = RemoveFirstActive())) {
                SUCt1 = SUC(t1);
                /* Choose t2 as one of t1's two neighbors on the tour */
                for (X2 = 1; X2 <= 2; X2++) {
                    t2 = X2 == 1 ? PRED(t1) : SUCt1;
                    if (FixedOrCommon(t1, t2) ||
                        (RestrictedSearch && Near(t1, t2) &&
                         (Trial == 1 ||
                          (Trial > BackboneTrials &&
                           (KickType == 0 || Kicks == 0)))))
                        continue;
                    G0 = C(t1, t2);
                    Gain = 0;
                    /* Try to find a tour-improving chain of moves */
                    do
                        t2 = Swaps == 0 ? BestMove(t1, t2, &G0, &Gain) :
                            BestSubsequentMove(t1, t2, &G0, &Gain);
                    while (t2);
                    if (Gain > 0) {
                        /* An improvement has been found */
                        assert(Gain % Precision == 0);
                        Cost -= Gain / Precision;
                        StoreTour();
                        if (HashSearch(HTable, Hash, Cost))
                            goto End_LinKernighan;
                        /* Make t1 "active" again */
                        Activate(t1);
                        break;
                    }
                    RestoreTour();
                    if (Dimension != DimensionSaved && SUC(t1) != SUCt1)
                        Reversed ^= 1;
                }
            }
            if (HashSearch(HTable, Hash, Cost))
                goto End_LinKernighan;
            HashInsert(HTable, Hash, Cost);
            /* Try to find improvements using non-sequential 4/5-opt moves */
            Gain = 0;
            if (Gain23Used && (Gain = Gain23()) > 0) {
                /* An improvement has been found */
                assert(Gain % Precision == 0);
                Cost -= Gain / Precision;
                StoreTour();
                if (HashSearch(HTable, Hash, Cost))
                    goto End_LinKernighan;
            }
        }
        while (Gain > 0);

    End_LinKernighan:
        PredSucCostAvailable = 0;
        // NormalizeNodeList();
        {
            Node *t1, *t2;

            t1 = FirstNode;
            do {
                t2 = SUC(t1);
                t1->Pred = PRED(t1);
                t1->Suc = t2;
                t1->Parent = 0;
            }
            while ((t1 = t2) != FirstNode);
            Reversed = 0;
        }

        // NormalizeSegmentList();
        {
            Segment *s1, *s2;

            s1 = FirstSegment;
            do {
                if (!s1->Parent->Reversed)
                    s2 = s1->Suc;
                else {
                    s2 = s1->Pred;
                    s1->Pred = s1->Suc;
                    s1->Suc = s2;
                }
            }
            while ((s1 = s2) != FirstSegment);
        }

        return Cost;
    }

    /* IsBackboneCandidate.c */
    int IsBackboneCandidate(const Node * ta, const Node * tb)
    {
        Candidate *Nta;

        for (Nta = ta->BackboneCandidateSet; Nta && Nta->To; Nta++)
            if (Nta->To == tb)
                return 1;
        return 0;
    }

    /* AdjustCandidateSet.c */
    void AdjustCandidateSet()
    {
        Candidate *NFrom, *NN, Temp;
        Node *From = FirstNode, *To;

        /* Extend and reorder candidate sets */
        do {
            if (!From->CandidateSet)
                From->CandidateSet = (Candidate *) calloc(3, sizeof(Candidate));
            /* Extend */
            for (To = From->Pred; To; To = To == From->Pred ? From->Suc : 0) {
                int Count = 0;
                if ((ProblemType == HCP || ProblemType == HPP) &&
                    !IsBackboneCandidate(From, To))
                    continue;
                for (NFrom = From->CandidateSet; NFrom->To && NFrom->To != To;
                     NFrom++)
                    Count++;
                if (!NFrom->To) {
                    /* Add new candidate edge */
                    NFrom->Cost = C(From, To);
                    NFrom->To = To;
                    NFrom->Alpha = INT_MAX;
                    From->CandidateSet =
                        (Candidate *) realloc(From->CandidateSet,
                                              (Count + 2) * sizeof(Candidate));
                    From->CandidateSet[Count + 1].To = 0;
                }
            }
            /* Reorder */
            for (NFrom = From->CandidateSet + 1; (To = NFrom->To); NFrom++)
                if (InBestTour(From, To) && InNextBestTour(From, To)) {
                    /* Move the edge to the start of the candidate table */
                    Temp = *NFrom;
                    for (NN = NFrom - 1; NN >= From->CandidateSet; NN--)
                        *(NN + 1) = *NN;
                    *(NN + 1) = Temp;
                }
        }
        while ((From = From->Suc) != FirstNode);
    }


    /* FindTour.c */
    GainType OrdinalTourCost;

    void SwapCandidateSets()
    {
        Node *t = FirstNode;
        do {
            Candidate *Temp = t->CandidateSet;
            t->CandidateSet = t->BackboneCandidateSet;
            t->BackboneCandidateSet = Temp;
        } while ((t = t->Suc) != FirstNode);
    }

    GainType FindTour()
    {
        GainType Cost;
        Node *t;
        int i;

        t = FirstNode;
        do
            t->OldPred = t->OldSuc = t->NextBestSuc = t->BestSuc = 0;
        while ((t = t->Suc) != FirstNode);
        if (Run == 1 && Dimension == DimensionSaved) {
            OrdinalTourCost = 0;
            for (i = 1; i < Dimension; i++)
                OrdinalTourCost += C(&NodeSet[i], &NodeSet[i + 1])
                    - NodeSet[i].Pi - NodeSet[i + 1].Pi;
            OrdinalTourCost += C(&NodeSet[Dimension], &NodeSet[1])
                - NodeSet[Dimension].Pi - NodeSet[1].Pi;
            OrdinalTourCost /= Precision;
        }
        BetterCost = PLUS_INFINITY;
        // if (MaxTrials > 0)
        HashInitialize(HTable);
        // else {
        //     Trial = 1;
        //     ChooseInitialTour();
        // }

        for (Trial = 1; Trial <= MaxTrials; Trial++) {
            /* Choose FirstNode at random */
            if (Dimension == DimensionSaved)
                FirstNode = &NodeSet[1 + Random() % Dimension];
            else
                for (i = Random() % Dimension; i > 0; i--)
                    FirstNode = FirstNode->Suc;
            ChooseInitialTour();
            Cost = LinKernighan();
            if (Cost < BetterCost) {
                BetterCost = Cost;
                RecordBetterTour();
                if (StopAtOptimum && BetterCost == Optimum)
                    break;
                AdjustCandidateSet();
                HashInitialize(HTable);
                HashInsert(HTable, Hash, Cost);
            }
        }
        t = FirstNode;
        if (Norm == 0 || MaxTrials == 0) {
            do
                t = t->BestSuc = t->Suc;
            while (t != FirstNode);
        }
        Hash = 0;
        do {
            (t->Suc = t->BestSuc)->Pred = t;
            Hash ^= Rand[t->Id] * Rand[t->Suc->Id];
        } while ((t = t->BestSuc) != FirstNode);
        if (Trial > MaxTrials)
            Trial = MaxTrials;
        ResetCandidateSet();
        return BetterCost;
    }

    LKH(int __n, const vector<pll>& points) : Dimension(__n) {
        ProblemFileName = PiFileName = InputTourFileName =
            OutputTourFileName = TourFileName = 0;
        CandidateFiles = 0;
        AscentCandidates = 50;
        BackboneTrials = 0;
        Backtracking = 0;
        CandidateSetSymmetric = 0;
        CandidateSetType = ALPHA;
        // Crossover = ERXT;
        DelaunayPartitioning = 0;
        DelaunayPure = 0;
        Excess = -1;
        ExtraCandidates = 0;
        ExtraCandidateSetSymmetric = 0;
        ExtraCandidateSetType = QUADRANT;
        // Gain23Used = 1;
        Gain23Used = 0;
        GainCriterionUsed = 1;
        GridSize = 1000000.0;
        InitialPeriod = -1;
        InitialStepSize = 0;
        InitialTourAlgorithm = WALK;
        InitialTourFraction = 1.0;
        KarpPartitioning = 0;
        KCenterPartitioning = 0;
        KMeansPartitioning = 0;
        Kicks = 1;
        KickType = 0;
        MaxBreadth = INT_MAX;
        // MaxCandidates = 5;
        MaxCandidates = 4;
        MaxPopulationSize = 0;
        MaxSwaps = -1;
        MaxTrials = -1;
        MoorePartitioning = 0;
        MoveType = 5;
        NonsequentialMoveType = -1;
        Optimum = MINUS_INFINITY;
        PatchingA = 1;
        PatchingC = 0;
        // Optimum = 12345;
        // PatchingA = 2;
        // PatchingC = 3;
        PatchingAExtended = 0;
        PatchingARestricted = 0;
        PatchingCExtended = 0;
        PatchingCRestricted = 0;
        Precision = 100;
        POPMUSIC_InitialTour = 0;
        POPMUSIC_MaxNeighbors = 5;
        POPMUSIC_SampleSize = 10;
        POPMUSIC_Solutions = 50;
        POPMUSIC_Trials = 1;
        Recombination = IPT;
        RestrictedSearch = 1;
        RohePartitioning = 0;
        // Runs = 10;
        Runs = 1;
        Seed = 1;
        SierpinskiPartitioning = 0;
        StopAtOptimum = 1;
        Subgradient = 1;
        SubproblemBorders = 0;
        SubproblemsCompressed = 0;
        SubproblemSize = 0;
        SubsequentMoveType = 0;
        SubsequentPatching = 1;
        TimeLimit = DBL_MAX;
        TotalTimeLimit = DBL_MAX;
        TraceLevel = 1;

        MaxMatrixDimension = 20000;

        /* ReadProblem() */

        CostMatrix = 0;
        BestTour = 0;
        BetterTour = 0;
        SwapStack = 0;
        HTable = 0;
        Rand = 0;
        CacheSig = 0;
        CacheVal = 0;
        Name = 0;
        Heap = 0;
        t = 0;
        T = 0;
        tSaved = 0;
        p = 0;
        q = 0;
        incl = 0;
        cycle = 0;
        G = 0;

        FirstNode = 0;
        WeightType = WeightFormat = ProblemType = -1;
        CoordType = NO_COORDS;
        // Name = Copy("Unnamed");
        Type = EdgeWeightType = EdgeWeightFormat = 0;
        EdgeDataFormat = NodeCoordType = DisplayDataType = 0;
        // Distance = 0;
        GridSize = 1000000.0;
        // C = 0;
        c = 0;
        ProblemType = TSP;
        DimensionSaved = Dimension;
        WeightType = EUC_2D;
        // Distance = Distance_EUC_2D;
        // c = c_EUC_2D;
        CoordType = TWOD_COORDS;

        /* Read coord section */
        {
            Node *N;
            int Id, i;
            if (!FirstNode) {
                Node *Prev = 0, *N = 0;
                int i;

                NodeSet = (Node *) calloc(Dimension + 1, sizeof(Node));
                for (i = 1; i <= Dimension; i++, Prev = N) {
                    N = &NodeSet[i];
                    if (i == 1)
                        FirstNode = N;
                    else
                        Link(Prev, N);
                    N->Id = i;
                }
                Link(N, FirstNode);
            }
            N = FirstNode;
            do
                N->V = 0;
            while ((N = N->Suc) != FirstNode);

            for (i = 1; i <= Dimension; i++) {
                N = &NodeSet[i];
                N->V = 1;
                N->X = points[i-1].first;
                N->Y = points[i-1].second;
            }
            N = FirstNode;
            do
                if (!N->V && N->Id <= Dimension)
                    break;
            while ((N = N->Suc) != FirstNode);
        }
        Swaps = 0;

        /* Adjust parameters */
        if (Seed == 0)
            Seed = (unsigned) time(0);
        if (Precision == 0)
            Precision = 100;
        if (InitialStepSize == 0)
            InitialStepSize = 1;
        if (MaxSwaps < 0)
            MaxSwaps = Dimension;
        if (KickType > Dimension / 2)
            KickType = Dimension / 2;
        if (Runs == 0)
            // Runs = 10;
            Runs = 1;
        if (MaxCandidates > Dimension - 1)
            MaxCandidates = Dimension - 1;
        if (ExtraCandidates > Dimension - 1)
            ExtraCandidates = Dimension - 1;
        if (SubproblemSize >= Dimension)
            SubproblemSize = Dimension;
        else if (SubproblemSize == 0) {
            if (AscentCandidates > Dimension - 1)
                AscentCandidates = Dimension - 1;
            if (InitialPeriod < 0) {
                InitialPeriod = Dimension / 2;
                if (InitialPeriod < 100)
                    InitialPeriod = 100;
            }
            if (Excess < 0)
                Excess = 1.0 / Dimension;
            if (MaxTrials == -1)
                // MaxTrials = Dimension;
                MaxTrials = 1;
            HeapMake(Dimension);
        }
        if (POPMUSIC_MaxNeighbors > Dimension - 1)
            POPMUSIC_MaxNeighbors = Dimension - 1;
        if (POPMUSIC_SampleSize > Dimension)
            POPMUSIC_SampleSize = Dimension;

        {
            Node *Ni, *Nj;
            CostMatrix = (int *) calloc((size_t) Dimension * (Dimension - 1) / 2,
                                        sizeof(int));
            Ni = FirstNode->Suc;
            do {
                Ni->C =
                    &CostMatrix[(size_t) (Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
                if (ProblemType != HPP || Ni->Id < Dimension)
                    for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                        Ni->C[Nj->Id] = Fixed(Ni, Nj) ? 0 : Distance(Ni, Nj);
                else
                    for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                        Ni->C[Nj->Id] = 0;
            }
            while ((Ni = Ni->Suc) != FirstNode);
            WeightType = EXPLICIT;
            c = 0;
        }
        if (Precision > 1 && (WeightType == EXPLICIT || ProblemType == ATSP)) {
            int j, n = ProblemType == ATSP ? Dimension / 2 : Dimension;
            for (int i = 2; i <= n; i++) {
                Node *N = &NodeSet[i];
                // for (j = 1; j < i; j++)
                //     if (N->C[j] * Precision / Precision != N->C[j])
                //         eprintf("PRECISION (= %d) is too large", Precision);
            }
        }
        // C = WeightType == EXPLICIT ? C_EXPLICIT : C_FUNCTION;
        // D = WeightType == EXPLICIT ? D_EXPLICIT : D_FUNCTION;
        if (SubsequentMoveType == 0)
            SubsequentMoveType = MoveType;
        K = MoveType >= SubsequentMoveType
            || !SubsequentPatching ? MoveType : SubsequentMoveType;
        if (PatchingC > K)
            PatchingC = K;
        if (PatchingA > 1 && PatchingA >= PatchingC)
            PatchingA = PatchingC > 2 ? PatchingC - 1 : 1;
        if (NonsequentialMoveType == -1 ||
            NonsequentialMoveType > K + PatchingC + PatchingA - 1)
            NonsequentialMoveType = K + PatchingC + PatchingA - 1;
        if (PatchingC >= 1) {
            // We copy this body
            // BestMove = BestSubsequentMove = BestKOptMove;
        }
        
        /* AllocateStructures() */
        {
            int i, K;

            HeapMake(Dimension);
            BestTour = (int *) calloc(1 + Dimension, sizeof(int));
            BetterTour = (int *) calloc(1 + Dimension, sizeof(int));
            HTable = (HashTable *) malloc(sizeof(HashTable));
            HashInitialize((HashTable *) HTable);
            SRandom(Seed);
            Rand = (unsigned *) malloc((Dimension + 1) * sizeof(unsigned));
            for (i = 1; i <= Dimension; i++)
                Rand[i] = Random();
            SRandom(Seed);
            if (WeightType != EXPLICIT) {
                for (i = 0; (1 << i) < (Dimension << 1); i++);
                i = 1 << i;
                CacheSig = (int *) calloc(i, sizeof(int));
                CacheVal = (int *) calloc(i, sizeof(int));
                CacheMask = i - 1;
            }

            // AllocateSegments();
            {
                Segment *S = 0, *SPrev;
                SSegment *SS = 0, *SSPrev;
                int i;
                
                FirstSegment = 0;
                FirstSSegment = 0;
                GroupSize = Dimension;
                Groups = 0;
                for (i = Dimension, SPrev = 0; i > 0; i -= GroupSize, SPrev = S) {
                    S = (Segment *) malloc(sizeof(Segment));
                    S->Rank = ++Groups;
                    if (!SPrev)
                        FirstSegment = S;
                    else
                        SLink(SPrev, S);
                }
                SLink(S, FirstSegment);
                SGroupSize = Dimension;
                SGroups = 0;
                for (i = Groups, SSPrev = 0; i > 0; i -= SGroupSize, SSPrev = SS) {
                    SS = (SSegment *) malloc(sizeof(SSegment));
                    SS->Rank = ++SGroups;
                    if (!SSPrev)
                        FirstSSegment = SS;
                    else
                        SLink(SSPrev, SS);
                }
                SLink(SS, FirstSSegment);
            }

            K = MoveType;
            if (SubsequentMoveType > K)
                K = SubsequentMoveType;
            T = (Node **) malloc((1 + 2 * K) * sizeof(Node *));
            G = (GainType *) malloc(2 * K * sizeof(GainType));
            t = (Node **) malloc(6 * K * sizeof(Node *));
            tSaved = (Node **) malloc((1 + 2 * K) * sizeof(Node *));
            p = (int *) malloc(6 * K * sizeof(int));
            q = (int *) malloc(6 * K * sizeof(int));
            incl = (int *) malloc(6 * K * sizeof(int));
            cycle = (int *) malloc(6 * K * sizeof(int));
            SwapStack =
                (SwapRecord *) malloc((MaxSwaps + 6 * K) * sizeof(SwapRecord));
        }

        /* CreateCandidateSet() */
        {
            GainType Cost, MaxAlpha, A;
            Node *Na;
            int CandidatesRead = 0, i;

            Norm = 9999;
            // if (C == C_EXPLICIT) {
            {
                Na = FirstNode;
                do {
                    for (i = 1; i < Na->Id; i++)
                        Na->C[i] *= Precision;
                }
                while ((Na = Na->Suc) != FirstNode);
            }

            // if (!ReadPenalties()) {
            {
                Na = FirstNode;
                do
                    Na->Pi = 0;
                while ((Na = Na->Suc) != FirstNode);
                CandidatesRead = 0;
                Cost = Ascent();
                if (Subgradient && SubproblemSize == 0) {
                    PiFile = 0;
                }
            }
            
            LowerBound = (double) Cost / Precision;

            MaxAlpha = (GainType) fabs(Excess * Cost);
            if ((A = Optimum * Precision - Cost) > 0 && A < MaxAlpha)
                MaxAlpha = A;
            GenerateCandidates(MaxCandidates, MaxAlpha, CandidateSetSymmetric);

        End_CreateCandidateSet:
            ResetCandidateSet();
            if (MaxTrials > 0 ||
                (InitialTourAlgorithm != SIERPINSKI &&
                 InitialTourAlgorithm != MOORE)) {
                Na = FirstNode;
                do {
                    // if (!Na->CandidateSet || !Na->CandidateSet[0].To) {
                    //     if (MaxCandidates == 0)
                    //         eprintf
                    //             ("MAX_CANDIDATES = 0: Node %d has no candidates",
                    //              Na->Id);
                    //     else
                    //         eprintf("Node %d has no candidates", Na->Id);
                    // }
                }
                while ((Na = Na->Suc) != FirstNode);
            }
            // if (C == C_EXPLICIT) {
            {
                Na = FirstNode;
                do
                    for (i = 1; i < Na->Id; i++)
                        Na->C[i] += Na->Pi + NodeSet[i].Pi;
                while ((Na = Na->Suc) != FirstNode);
            }
        }
        
        BestCost = PLUS_INFINITY;

        for (Run = 1; Run <= Runs; Run++) {
            Cost = FindTour(); // LKH heuristic
            if (Cost < BestCost) {
                BestCost = Cost;
                RecordBetterTour();
                RecordBestTour();
            }
            OldOptimum = Optimum;
            Optimum = min(Optimum, Cost);
            if (StopAtOptimum && Cost == OldOptimum && MaxPopulationSize >= 1) {
                Runs = Run;
                break;
            }
            SRandom(++Seed);
        }
    }

    pair<double, vector<int>> return_result() {
        vector<int> ret(DimensionSaved);
        for (int i=0; i<DimensionSaved; ++i) {
            ret[i] = BestTour[i+1];
        }
        return { BestCost, ret };
    }
};

vector<double> cut(vector<double> bars) {
    int n = bars.size();
    vector<double> ret(n-1);
    double m = accumulate(bars.begin(), bars.end(), 0.0) / n;
    double accum_m = 0;
    int cnt = 0;
    double used_ratio = 0;
    for (int i=0; i<n; ++i) {
        if (cnt == n-1) break;
        if (accum_m + bars[i] * (1.0 - used_ratio) >= m * (cnt+1)) {
            ret[cnt] = (double)i + used_ratio + (m * (cnt+1) - accum_m) / bars[i];
            used_ratio += (m * (cnt+1) - accum_m) / bars[i];
            accum_m = m * (cnt+1);
            cnt++;
            i--;
        } else {
            accum_m += bars[i] * (1.0 - used_ratio);
            used_ratio = 0;
        }
    }
    return ret;
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

    const int MAX_COORD = 814000;
    int dx = (MAX_COORD + 10) / 10;
    int dy = (MAX_COORD + 14) / 14;

    vector<vector<int>> groups;
    for (int i=0; i<10; ++i) {
        for (int j=0; j<14; ++j) {
            vector<int> local_points;
            for (int l=0; l<n; ++l) {
                auto[x, y] = points[l];
                if (dx * i <= x && x < dx * (i+1) && dy * j <= y && y < dy * (j+1)) {
                    local_points.push_back(l);
                }
            }
            groups.push_back(local_points);
        }
    }

    vector<double> tourLen(k, 1e9);
    for (int i=0; i<k; ++i) {
        // Run the 3-opt algorithm
        vector<pll> convert;
        for (int idx: groups[i]) convert.push_back(points[idx]);
        if (i == 1) {
            // for (int j=0; j<convert.size(); ++j) {
            //     cerr << j+1 << ' ' << convert[j].first << ' ' << convert[j].second << '\n';
            // }

        }
        LKH lkh(convert.size(), convert);
        tourLen[i] = lkh.return_result().first;
        // if (i % 14 == 0) print(i);
        // tourLen[i] = min(tourLen[i], threeOpt(convert).first);
    }
    // print(tourLen);
    
    vector<double> w_ysum(14, 0);
    for (int j=0; j<14; ++j) {
        for (int i=0; i<10; ++i) w_ysum[j] += tourLen[14 * i + j];
    }
    vector<double> _ycut = cut(w_ysum);
    vector<int> ycut;
    for (double y: _ycut) ycut.push_back(y * dy);
    // print(ycut);

    for (int i=0; i<10; ++i) {
        vector<vector<int>> local_points(14);
        for (int l=0; l<n; ++l) {
            auto[x, y] = points[l];
            if (dx * i <= x && x < dx * (i+1)) {
                int ynum = lower_bound(ycut.begin(), ycut.end(), y) - ycut.begin();
                local_points[ynum].push_back(l);
            }
        }
        for (int j=0; j<14; ++j) {
            groups[14 * i + j] = local_points[j];
        }
    }
    // print(1);

    for (int i=0; i<k; ++i) {
        // Run the 3-opt algorithm
        vector<pll> convert;
        for (int idx: groups[i]) convert.push_back(points[idx]);
        LKH lkh(convert.size(), convert);
        tourLen[i] = lkh.return_result().first;
        // if (i % 14 == 0) print(i);
        // tourLen[i] = threeOpt(convert).first;
    }
    // print(tourLen);
    
    for (int j=0; j<14; ++j) {
        vector<double> w_xsum(10);
        for (int i=0; i<10; ++i) w_xsum[i] = tourLen[14 * i + j];
        vector<double> _xcut = cut(w_xsum);
        vector<int> xcut;
        for (double x: _xcut) xcut.push_back(x * dx);
        // print(j, w_xsum, _xcut, xcut);

        vector<int> targets;
        for (int i=0; i<10; ++i) {
            for (int idx: groups[14 * i + j]) {
                targets.push_back(idx);
            }
        }
        vector<vector<int>> local_points(10);
        for (int l: targets) {
            auto[x, y] = points[l];
            int xnum = lower_bound(xcut.begin(), xcut.end(), x) - xcut.begin();
            local_points[xnum].push_back(l);
        }
        for (int i=0; i<10; ++i) {
            groups[14 * i + j] = local_points[i];
        }
    }

    vector<double> _tourLen(k);
    for (int i=0; i<k; ++i) {
        // Run the 3-opt algorithm
        vector<pll> convert;
        for (int idx: groups[i]) convert.push_back(points[idx]);

        LKH lkh(convert.size(), convert);
        // auto tmp = threeOpt(convert);
        _tourLen[i] = lkh.return_result().first;

        vector<int> tour = lkh.return_result().second;
        // cout << tmp.first << endl;
        cout << tour.size() << ' ';
        for (int idx: tour) cout << groups[i][idx-1]+1 << ' ';
        cout << '\n';
    }

    // for (int j=13; j>=0; --j) {
    //     for (int i=0; i<10; ++i) {
    //         cerr << _tourLen[14*i+j] << ' ';
    //     }
    //     cerr << '\n';
    // }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
}
