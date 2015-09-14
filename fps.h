#ifndef FPS_H
#define FPS_H

#include <base.h>
#include <mine.h>

// typedef TIntFltKdV TSpVec;

///////////////////////
// Farthest-Point-Sampling 
class TFpsUtil {
public:
    // Distance (Euclidean) between sparse vectors 
    static double dist(const TIntFltKdV& x, const TIntFltKdV& y) {
        // AssertR(1.0-TLinAlg::Norm(x) <= 1e-4 && 1.0-TLinAlg::Norm(y) <= 1e-4, "Not L2-normalized.");
        TIntFltKdV z; z.Reserve(x.Len()+y.Len()); // z = x - y 
        TLinAlg::LinComb(1.0, x, -1.0, y, z);
        return TLinAlg::Norm(z);
    }
    // Distance (Euclidean) between "ordinary" vectors 
    static double dist(const TFltV& x, const TFltV& y) {
        // printf("%d :: %d\n", x.Len(), y.Len());
        TFltV z(x); // z.Reserve(x.Len()); // z = x - y 
        // printf("z: %d\n", z.Len());
        AssertR(x.Len() == y.Len(), "x.len() == y.len()");
        TLinAlg::LinComb(1.0, x, -1.0, y, z);
        return TLinAlg::Norm(z); // XXX: return TLinAlg::DotProduct(x, y);
    }
    static double dist_cos(const TIntFltKdV& x, const TIntFltKdV& y) {
        const double cosDist = 1.0-TLinAlg::DotProduct(x, y);
        AssertR(cosDist <= 1.0 && cosDist >= -1e-9, "Oopst, cosine distance larger than 1 or smaller than 0.");
        return TMath::Log2(1.0+cosDist);
    }
    static void GetLandmarks(const TVec<TIntFltKdV>& S, const int& l, TVec<TInt>& L);
    static int NextLandmark(const TVec<TIntFltKdV>& S, const TVec<TInt>& L);
    static void Project(const TVec<TIntFltKdV>& S, const TVec<TInt>& L, TVec<TFltV>& ProjS);
    static int GetFarthestPoint(const TVec<TIntFltKdV>& S, const int& PtIdx);
};

///////////////////////
// Landmarks-By-Dimension
class TDimUtil {
public:
    static void GetLandmarks(const TVec<TIntFltKdV>& S, const int& l, TVec<TInt>& L);
};

///////////////////////
// Random-Sampling
class TRndSampleUtil {
public:
    static void GetLandmarks(const TVec<TIntFltKdV>& S, const int& l, TVec<TInt>& L) {
        GetLandmarks(S.Len(), l, L);
    }
    static void GetLandmarks(const int& n, const int& l, TVec<TInt>& L);
};

///////////////////////
// Sparse-Random-Projection
class TRandProjUtil {
public:
    static void Project(const TVec<TIntFltKdV>& S, TVec<TFltV>& ProjS, const int& Dim, TIntV& ZerosV);
};

///////////////////////
// K-Means-Sampling
// TODO: Find a scalable implementation of k-means clustering for sparse data 
class TKMeansSampleUtil {
public:
    static void GetLandmarks(const TVec<TIntFltKdV>& S, const int& l, TVec<TInt>& L);
};

#endif

