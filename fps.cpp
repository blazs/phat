#include "fps.h"

///////////////////////
// Farthest-point-sampling 
// This naive algorithm does O(l^2 * n) distance function evaluations: n+2n+3n+...+ln
// Note that we usually set l = log n, which gives O(n (log n)^2)
void TFpsUtil::GetLandmarks(const TVec<TIntFltKdV>& S, const int& l, TVec<TInt>& L) {
    EAssertR(l > 0, "Number of landmark points must be positive.");
    EAssertR(S.Len() >= l, "Number of landmark points exceeds the sample size");
    
    TRnd Rnd;
    // Pick the initial landmark uniformly at random 
    L.Add(Rnd.GetUniDevInt(0, S.Len()));
    while (L.Len() < l) { L.Add(NextLandmark(S, L)); }
}

// Given a vector of sparse vectors S and a set of landmark indices L, it 
// returns the index of the point x in S that maximizes distance between
// x and the closest landmark: argmax_x (min_l d(x,l)) 
int TFpsUtil::NextLandmark(const TVec<TIntFltKdV>& S, const TVec<TInt>& L) {
    double MxMnDst = 0.0, CrrDst = 0.0;
    int MxMnIdx = 0;
    for (int PtIdx = 0; PtIdx < S.Len(); ++PtIdx) {
        CrrDst = dist(S.GetVal(PtIdx), S.GetVal(L.GetVal(0)));
        for (auto LmIt = L.BegI(); LmIt != L.EndI(); ++LmIt) {
            CrrDst = TMath::Mn(CrrDst, dist(S.GetVal(PtIdx), S.GetVal(*LmIt)));
        }
        if (CrrDst > MxMnDst) { MxMnIdx = PtIdx; }
    }
    return MxMnIdx;
}

// Projects the set S into R^l, mapping each x from S into
// (d(x,p1), d(x,p2), ..., d(x,pl)),  a vector of distances 
// of x to each landmark 
// Does O(l * n) distance function evaluations; for l = log n, this is O(n log n) 
void TFpsUtil::Project(const TVec<TIntFltKdV>& S, const TVec<TInt>& L, TVec<TFltV>& ProjS) {
    for (int PtIdx = 0; PtIdx < S.Len(); ++PtIdx) {
        TFltV CrrV;
        for (int LmIdx = 0; LmIdx < L.Len(); ++LmIdx) {
            const double CrrDst = dist(S.GetVal(PtIdx), S.GetVal(LmIdx));
            CrrV.Add(CrrDst);
        }
        ProjS.Add(CrrV);
        CrrV.Clr();
    }
}

// Given q, find p that maximizes d(p,q) 
int TFpsUtil::GetFarthestPoint(const TVec<TIntFltKdV>& S, const int& PtIdx) {
    int FstIdx = 0;
    double Dst = TFpsUtil::dist(S[PtIdx], S[FstIdx]);
    for (int Idx = 0; Idx < S.Len(); ++Idx) {
        double TmpDst = TFpsUtil::dist(S[Idx], S[PtIdx]);
        if (TmpDst > Dst) {
            Dst = TmpDst; FstIdx = Idx;
        }
    }
    return FstIdx;
}

///////////////////////
// Landmarks-By-Dimension
// Picks the set of the first l points from S, so that each pairs of points
// differs in a least one dimension 
void TDimUtil::GetLandmarks(const TVec<TIntFltKdV>& S, const int& l, TVec<TInt>& L) {
    for (int PtIdx = 0; PtIdx < S.Len(); ++PtIdx) {
        TFltV CrrV;
        bool IsOk = (L.Len() == 0);
        for (int LmIdx = 0; LmIdx < L.Len(); ++LmIdx) {
            // Check whether they "point" in different "directions" 
            const auto V1 = S.GetVal(PtIdx);
            const auto V2 = S.GetVal(L.GetVal(LmIdx));
            auto It1 = V1.BegI();
            auto It2 = V2.BegI();
            while (It1 != V1.EndI() && It2 != V2.EndI() && It1->Key == It2->Key)
                { ++It1; ++It2; }
            IsOk = !(It1 == V1.EndI() && It2 == V2.EndI());
        }
        // printf("%d", IsOk);
        if (IsOk) { L.Add(PtIdx); }
        if (L.Len() >= l) { break; } // OK, enough landmarks 
    }
}

///////////////////////
// Random-Sampling
// Picks a random l-subset of [n], using the algorithm from Knuth's TAOCP 
// Volume 2, section 3.4.2, p. 142. (XXX: Consider reservoir sampling?)
void TRndSampleUtil::GetLandmarks(const int& n, const int& l, TVec<TInt>& L) {
    EAssertR(n >= l, "Can't select l-subset of items from a n-set if n < l.");
    TRnd Rnd(0);
    // Idx ... total number of records we've seen so far 
    int m = 0; // Number of items selected so far 
    for (int Idx = 0; Idx < n && m < l; ++Idx) {
        const double U = Rnd.GetUniDev(); // Random 0 <= U <= 1
        if ((n-Idx)*U < l-m) {
            ++m;
            L.Add(Idx);
        }
    }
    EAssertR(m == l, "Oops: l != m");
}

///////////////////////
// Sparse-Random-Projection
// Multiplies d-times-n matrix S with a random n-times-k matrix A=(a_ij),
// where k << d and a_ij ~ (-1, 0, 1) with probabilities (1/6 2/3 1/6) 
// See "Database-friendly random projections: Johnson-Lindenstrauss with
// binary coins" by Dimitris Achlioptas
void TRandProjUtil::Project(const TVec<TIntFltKdV>& S, TVec<TFltV>& ProjS, const int& Dim, TIntV& ZerosV) {
    TVec<TIntFltKdV> R; // Random matrix 
    TRnd Rnd(0);
    // TVec<TIntFltKdV> St;
    // TLinAlg::Transpose(S, St);
    for (auto It = S.BegI(); It != S.EndI(); ++It) {
        TFltV TmpV;
        int ZerosN = 0;
        for (int r = 0; r < Dim; ++r) {
            double Res = 0.0;
            for (auto VecIt = It->BegI(); VecIt != It->EndI(); ++VecIt) {
                // Pr[X = 1] = 1/3; Pr[X = 0] = 2/3;
                int TmpRnd = (Rnd.GetUniDevInt(0, 2) == 1); 
                if (TmpRnd != 0) {
                    // Pr[sgn = -] = 1/2; Pr[sgn = +] = 1/2;
                    TmpRnd = Rnd.GetUniDevInt(0, 1) == 1 ? TmpRnd : -TmpRnd;
                    TmpRnd *= TMath::Sqrt(3);
                    Res += TmpRnd*VecIt->Dat;
                }
            }
            TmpV.Add(Res);
            if (Res == 0.0) { ++ZerosN; }
        }
        ZerosV.Add(ZerosN);
        // printf("ZerosN = %d\n", ZerosN);
        ProjS.Add(TmpV);
    }
}

///////////////////////
// KMeans-Sampling
// template<class V, class LA, class M>
// typedef TKMeans<TIntFltKdV, TLinAlg, TSparseColMatrix> TSparseKMeans;
void TKMeansSampleUtil::GetLandmarks(const TVec<TIntFltKdV>& S, const int& l, TVec<TInt>& L) {
    EFailR("Need to reimplement.");
    
    TIntV IdxV;
    printf("Selecting a random sample of the data set...\n");
    TRndSampleUtil::GetLandmarks(S.Len(), S.Len()/100, IdxV);
    TVec<TIntFltKdV> Sub;
    for (auto It = IdxV.BegI(); It != IdxV.EndI(); ++It) {
        Sub.Add(S.GetVal(It->Val));
    }
    
    printf("Converting to sparse column matrix...\n");
    TSparseColMatrix SpColMatrix(Sub);
    printf("SpColMatrix: %d x %d\n", SpColMatrix.RowN, SpColMatrix.ColN);
    printf("KMeans on %d data points...\n", Sub.Len());
    TSparseKMeans SpKMeans(&SpColMatrix, l, 10000, TStdNotify::New());
    printf("Initializing k-means...\n");
    SpKMeans.Init();
    printf("Running k-means...\n");
    SpKMeans.Apply();
    
    const TIntV& AssignmentV = SpKMeans.GetAssignment();
    for (auto It = AssignmentV.BegI(); It != AssignmentV.EndI(); ++It) {
        printf("%d\n", It->Val);
    }
    
}

