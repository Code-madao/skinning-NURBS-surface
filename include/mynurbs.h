#pragma once

#include "tinynurbs/tinynurbs.h"
#include "glm/glm.hpp"
#include "MathUtils.h"

#include <vector>
#include <set>
#include <algorithm>

#define MyEpsilon 1e-6f
#define MaxDistance 1e6

using namespace std;

// using vec3f = glm::vec3;
using vec3f = glm::f64vec3;

namespace myNurbs
{
    // return index of first knot less than u
    inline size_t getKnotSpanIndex(unsigned degree, const vector<double> &knots, double u)
    {
        int n = knots.size() - degree - 2;
        assert(n >= 0);
        if (knots[n + 1] - u < MyEpsilon)
            return n;
        if (u - knots[degree] < MyEpsilon)
            return degree;

        int low = degree, high = n + 1;
        int mid;
        while (low < high)
        {
            mid = (low + high) / 2;
            if (knots[mid] >= u)
                high = mid;
            else
                low = mid + 1;
        }
        return low - 1;
    }

    inline vector<double> getBsplineBasis(size_t spanIndex, unsigned degree, const vector<double> &knots, double u)
    {
        // return tinynurbs::bsplineBasis(degree, spanIndex, knots, u);
        vector<double> basisFunc(degree + 1, 0.);
        vector<double> leftDiff(degree + 1, 0.);
        vector<double> rightDiff(degree + 1, 0.);

        // non-recursively compute triangle
        basisFunc[0] = 1.0;
        for (size_t i = 1; i <= degree; ++i)
        {
            leftDiff[i] = u - knots[spanIndex + 1 - i];
            rightDiff[i] = knots[spanIndex + i] - u;

            double saved = 0;

            for (size_t j = 0; j < i; ++j)
            {
                double temp = basisFunc[j] / (rightDiff[j + 1] + leftDiff[i - j]);
                basisFunc[j] = saved + rightDiff[j + 1] * temp;
                saved = leftDiff[i - j] * temp;
            }
            basisFunc[i] = saved;
        }
        return basisFunc;
    }

    inline vector<glm::dvec4> getWeightCurveCtrPnts(const tinynurbs::RationalCurve3d& curve)
    {
        vector<glm::dvec4> res;
        for (size_t i = 0; i < curve.control_points.size(); i++)
        {
            glm::dvec4 temp = glm::dvec4(curve.control_points[i]*curve.weights[i], curve.weights[i]);
            res.emplace_back(temp);
        }
        return res;
    }
    
    inline void decomposeWeightCurveCtrPnts(const vector<glm::dvec4>& wCtrPnts, vector<vec3f>& ctrPnts, vector<double>& weight)
    {
        ctrPnts.resize(wCtrPnts.size());
        weight.resize(wCtrPnts.size());

        for (size_t i = 0; i < wCtrPnts.size(); i++)
        {
            auto temp = wCtrPnts[i];
            weight[i]=temp.w;
            temp/=temp.w;
            ctrPnts[i] = vec3f(temp);
        }
    }
    
    template <typename T>
    T myCurvePointTemplate(unsigned degree, const vector<double> &knots, const vector<T> &ctrPoints, double u)
    {
        size_t spanIndex = getKnotSpanIndex(degree, knots, u);
        vector<double> N = getBsplineBasis(spanIndex, degree, knots, u);

        T point(T(0));
        for (size_t i = 0; i <= degree; ++i)
        {
            point += ctrPoints[spanIndex - degree + i] * N[i];
        }
        return point;
    }

    // NURBS curve
    inline vec3f myRationalCurvePoint(const tinynurbs::RationalCurve3d &crv, double u)
    {
        // return tinynurbs::curvePoint(crv, u);
        vector<glm::dvec4> ctrPntsWithWeights;
        for (size_t i = 0; i < crv.control_points.size(); ++i)
        {
            ctrPntsWithWeights.push_back(glm::dvec4(crv.control_points[i] * crv.weights[i], crv.weights[i])); // MAYBEBUG
        }

        glm::dvec4 point = myCurvePointTemplate<glm::dvec4>(crv.degree, crv.knots, ctrPntsWithWeights, u);
        // glm::dvec4 point = tinynurbs::internal::curvePoint<4,double>(crv.degree, crv.knots, ctrPntsWithWeights, u);

        return vec3f(point / point.w); // MAYBEBUG
    }

    inline vector<vector<double>> getBasisFuncDerivatives(size_t spanIndex, unsigned degree, int derivativeNum,
                                                          const vector<double> &knots, double u)
    {
        vector<vector<double>> basisFuncDer(derivativeNum + 1, vector<double>(degree + 1));

        // tinynurbs::array2<double> temp= tinynurbs::bsplineDerBasis(degree, spanIndex, knots, u, derivativeNum);
        // for (size_t i = 0; i < temp.rows(); ++i)
        //{
        //    for (size_t j = 0; j < temp.cols(); ++j)
        //    {
        //       basisFuncDer[i][j] = temp(i, j);
        //    }
        // }
        // return basisFuncDer;

        vector<vector<double>> ndu(degree + 1, vector<double>(degree + 1));
        ndu[0][0] = 1.f;

        vector<double> left(degree + 1);
        vector<double> right(degree + 1);

        for (size_t i = 1; i <= degree; i++)
        {
            left[i] = u - knots[spanIndex + 1 - i];
            right[i] = knots[spanIndex + i] - u;

            double temp = 0.f, saved = 0.f;
            for (size_t j = 0; j < i; j++)
            {
                ndu[i][j] = right[j + 1] + left[i - j];
                temp = ndu[j][i - 1] / ndu[i][j];

                ndu[j][i] = saved + right[j + 1] * temp;
                saved = left[i - j] * temp;
            }
            ndu[i][i] = saved;
        }

        for (size_t i = 0; i <= degree; i++)
            basisFuncDer[0][i] = ndu[i][degree];

        vector<vector<double>> a(2, vector<double>(degree + 1));
        for (size_t i = 0; i <= degree; i++)
        {
            int s1 = 0;
            int s2 = 1;
            a[0][0] = 1.0;
            // MAYBEBUG
            for (int k = 1; k <= derivativeNum; k++)
            {
                double d = 0.0;
                int rk = i - k;
                int pk = degree - k;

                if (i >= k)
                {
                    a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                    d = a[s2][0] * ndu[rk][pk];
                }

                int j1 = 0;
                int j2 = 0;

                if (rk >= -1)
                    j1 = 1;
                else
                    j1 = -rk;

                if (i - 1 <= pk)
                    j2 = k - 1;
                else
                    j2 = degree - i;

                for (int j = j1; j <= j2; j++)
                {
                    a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                    d += a[s2][j] * ndu[rk + j][pk];
                }
                if (i <= pk)
                {
                    a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][i];
                    d += a[s2][k] * ndu[i][pk];
                }
                basisFuncDer[k][i] = d;

                swap(s1, s2);
            }
        }
        int temp = degree;
        for (size_t i = 1; i <= derivativeNum; i++)
        {
            for (size_t j = 0; j <= degree; j++)
            {
                basisFuncDer[i][j] *= temp;
            }
            temp *= degree - i;
        }
        return basisFuncDer;
    }

    // The NURBS book page93
    template <typename T>
    vector<T> myCurveDerivativeTemplate(unsigned degree, const vector<double> &knots,
                                        const vector<T> &ctrPoints, int derivativeNum, double u)
    {
        vector<T> curveDerivatives(derivativeNum + 1, T(0));

        size_t spanIndex = getKnotSpanIndex(degree, knots, u);
        vector<vector<double>> basisFuncDer = getBasisFuncDerivatives(spanIndex, degree, derivativeNum, knots, u);

        int D = min(derivativeNum, (int)degree);
        for (size_t i = 0; i <= D; i++)
        {
            for (size_t j = 0; j <= degree; j++)
            {
                curveDerivatives[i] += basisFuncDer[i][j] * ctrPoints[spanIndex - degree + j];
            }
        }
        return curveDerivatives;
    }

    inline int getFactorial(int num)
    {
        if (num == 0)
            return 1;
        else
            return num * getFactorial(num - 1);
    }

    // compute Cnk
    inline double binomial(int n, int k)
    {
        return getFactorial(n) / (getFactorial(k) * getFactorial(n - k));
    }

    inline vector<vec3f> myRationalCurveDerivative(const tinynurbs::RationalCurve3d &crv, int derivativeNum, double u)
    {
        vector<vec3f> result;
        vector<glm::dvec4> ctrPntsWithWeights;
        for (size_t i = 0; i < crv.control_points.size(); ++i)
            ctrPntsWithWeights.push_back(glm::dvec4(crv.control_points[i] * crv.weights[i], crv.weights[i])); // MAYBEBUG

        vector<glm::dvec4> tempDer = myCurveDerivativeTemplate<glm::dvec4>(crv.degree, crv.knots, ctrPntsWithWeights, derivativeNum, u);

        vector<vec3f> tempDerTrunc;
        vector<double> tempDerWeight;
        for (const auto &ele : tempDer)
        {
            tempDerTrunc.push_back(vec3f(ele)); // MAYBEBUG
            tempDerWeight.push_back(ele.w);
        }

        for (size_t i = 0; i <= derivativeNum; i++)
        {
            vec3f t = tempDerTrunc[i];
            for (size_t j = 1; j <= i; j++)
            {
                t -= binomial(i, j) * tempDerWeight[j] * result[i - j];
            }
            result.push_back(t / tempDerWeight[0]);
        }
        return result;
    }

    // NURBS Book page103
    template <typename T>
    T mySurfPointTemplate(unsigned uDegree, const vector<double> &uKnots, double u,
                          unsigned vDegree, const vector<double> &vKnots, double v,
                          const tinynurbs::array2<T> &ctrPoints)
    {
        size_t uSpanIndex = getKnotSpanIndex(uDegree, uKnots, u);
        size_t vSpanIndex = getKnotSpanIndex(vDegree, vKnots, v);
        vector<double> uN = getBsplineBasis(uSpanIndex, uDegree, uKnots, u);
        vector<double> vN = getBsplineBasis(vSpanIndex, vDegree, vKnots, v);

        T result(T(0));
        size_t uInd = uSpanIndex - uDegree;
        size_t vInd = vSpanIndex - vDegree;
        for (size_t i = 0; i <= vDegree; i++)
        {
            T temp(T(0));
            for (size_t j = 0; j <= uDegree; j++)
            {
                temp += uN[j] * ctrPoints(uInd + j, vInd + i);
            }
            result += vN[i] * temp;
        }
        return result;
    }

    inline vec3f myRationalSurfPoint(const tinynurbs::RationalSurface3d &surf, double u, double v)
    {
        size_t rows = surf.control_points.rows();
        size_t cols = surf.control_points.cols();

        tinynurbs::array2<glm::dvec4> ctrPointsWithWeight;
        ctrPointsWithWeight.resize(rows, cols);

        for (size_t i = 0; i < surf.control_points.rows(); i++)
        {
            for (size_t j = 0; j < surf.control_points.cols(); j++)
            {
                ctrPointsWithWeight(i, j) = (glm::dvec4(surf.control_points(i, j) * surf.weights(i, j), surf.weights(i, j)));
            }
        }

        glm::dvec4 result = mySurfPointTemplate<glm::dvec4>(surf.degree_u, surf.knots_u, u,
                                                            surf.degree_v, surf.knots_v, v,
                                                            ctrPointsWithWeight);
        return vec3f(result / result.w); // MAYBEBUG
    }

    template <typename T>
    vector<vector<T>> mySurfDerivativeTemplate(unsigned uDegree, const vector<double> &uKnots, double u,
                                               unsigned vDegree, const vector<double> &vKnots, double v,
                                               const tinynurbs::array2<T> &ctrPoints, int derivativeNum)
    {
        vector<vector<T>> result(derivativeNum + 1, vector<T>(derivativeNum + 1, T(0)));

        size_t uSpanIndex = getKnotSpanIndex(uDegree, uKnots, u);
        size_t vSpanIndex = getKnotSpanIndex(vDegree, vKnots, v);
        vector<vector<double>> uDer = getBasisFuncDerivatives(uSpanIndex, uDegree, derivativeNum, uKnots, u);
        vector<vector<double>> vDer = getBasisFuncDerivatives(vSpanIndex, vDegree, derivativeNum, vKnots, v);
        // auto uDer = tinynurbs::bsplineDerBasis(uDegree, uSpanIndex, uKnots, u, derivativeNum);
        // auto vDer = tinynurbs::bsplineDerBasis(vDegree, vSpanIndex, vKnots, v, derivativeNum);

        int uD = min(derivativeNum, (int)uDegree);
        int vD = min(derivativeNum, (int)vDegree);

        vector<T> temp(vDegree + 1);
        for (size_t k = 0; k <= uD; k++)
        {
            for (size_t i = 0; i <= vDegree; i++)
            {
                temp[i] = T(0);
                for (size_t j = 0; j <= uDegree; j++)
                {
                    temp[i] += uDer[k][j] * ctrPoints(uSpanIndex - uDegree + j, vSpanIndex - vDegree + i);
                    // temp[i] += uDer(k,i) * ctrPoints(uSpanIndex - uDegree + j, vSpanIndex - vDegree + i);
                }
            }
            int kk = min(derivativeNum - (int)k, vD);
            for (size_t i = 0; i <= kk; i++)
            {
                for (size_t j = 0; j <= vDegree; j++)
                {
                    result[k][i] += vDer[i][j] * temp[j];
                    // result[k][i] += vDer(i,j) * temp[j];
                }
            }
        }
        return result;
    }

    inline tinynurbs::array2<vec3f> myRationalSurfDerivative(const tinynurbs::RationalSurface3d &surf, int derivativeNum, double u, double v)
    {
        tinynurbs::array2<glm::dvec4> ctrPntsWithWeights;
        ctrPntsWithWeights.resize(surf.control_points.rows(), surf.control_points.cols());
        for (size_t i = 0; i < surf.control_points.rows(); i++)
            for (size_t j = 0; j < surf.control_points.cols(); j++)
                ctrPntsWithWeights(i, j) = glm::dvec4(surf.control_points(i, j) * surf.weights(i, j), surf.weights(i, j));

        vector<vector<glm::dvec4>> tempDer = mySurfDerivativeTemplate<glm::dvec4>(surf.degree_u, surf.knots_u, u,
                                                                                  surf.degree_v, surf.knots_v, v,
                                                                                  ctrPntsWithWeights, derivativeNum);
        // tinynurbs::array2<glm::dvec4> tempDer = tinynurbs::internal::surfaceDerivatives(surf.degree_u, surf.degree_v, surf.knots_u, surf.knots_v, ctrPntsWithWeights, derivativeNum, u, v);
        tinynurbs::array2<vec3f> result;
        vector<vector<vec3f>> tempDerTrunc;
        vector<vector<double>> tempDerWeight;
        result.resize(derivativeNum + 1, derivativeNum + 1);
        tempDerTrunc.resize(derivativeNum + 1, vector<vec3f>(derivativeNum + 1));
        tempDerWeight.resize(derivativeNum + 1, vector<double>(derivativeNum + 1));
        for (size_t i = 0; i < tempDer.size(); i++)
        // for (size_t i = 0; i < tempDer.rows(); i++)
        {
            for (size_t j = 0; j < tempDer[0].size(); j++)
            // for (size_t j = 0; j < tempDer.cols(); j++)
            {
                tempDerTrunc[i][j] = tempDer[i][j];
                tempDerWeight[i][j] = tempDer[i][j].w;
                // tempDerTrunc[i][j] = tempDer(i,j);
                // tempDerWeight[i][j] = tempDer(i, j).w;
            }
        }

        for (size_t k = 0; k <= derivativeNum; k++)
        {
            for (size_t i = 0; i <= derivativeNum - k; i++)
            {
                vec3f t = tempDerTrunc[k][i];
                for (size_t j = 1; j <= i; j++)
                {
                    t -= binomial(i, j) * tempDerWeight[0][j] * result(k, i - j);
                }
                for (size_t j = 1; j <= k; j++)
                {
                    t -= binomial(k, j) * tempDerWeight[j][0] * result(k - j, i);

                    vec3f tVec(0.f);
                    for (size_t m = 1; m <= i; m++)
                    {
                        tVec += binomial(i, m) * tempDerWeight[j][m] * result(k - j, i - m);
                    }
                    t -= double(binomial(k, j)) * tVec;
                }
                result(k, i) = t / tempDerWeight[0][0];
            }
        }
        return result;
    }

    inline vec3f myRationalSurfNormal(const tinynurbs::RationalSurface3d &surf, double u, double v)
    {
        // return tinynurbs::surfaceNormal(surf, u, v);

        tinynurbs::array2<vec3f> der = myRationalSurfDerivative(surf, 1, u, v);
        // tinynurbs::array2<vec3f> der = tinynurbs::surfaceDerivatives(surf, 1, u, v);

        vec3f result = glm::cross(der(0, 1), der(1, 0));
        double length = glm::length(result);
        if (length > MyEpsilon)
        {
            result /= length;
        }
        return result;
    }

    inline double getChordLength(const vector<vector<glm::dvec4>> &ctrPnts, int index)
    {
        double totalLength = 0;
        for (size_t i = 1; i < ctrPnts.size(); i++)
            totalLength += glm::distance(ctrPnts[i][index], ctrPnts[i - 1][index]);
        return totalLength;
    }

    //Algo A5.9 in the NURBS book 2nd p206
    inline tinynurbs::RationalCurve3d elevateCurveDegree(const tinynurbs::RationalCurve3d &curve, int times)
    {
        tinynurbs::RationalCurve3d result;
        int degree = curve.degree;
        vector<double> knotVector = curve.knots;
        auto controlPoints = getWeightCurveCtrPnts(curve);

        int n = controlPoints.size() - 1;
        int m = n + degree + 1;
        int ph = degree + times;
        int ph2 = floor(ph / 2);

        // compute Bezier degree elevation coefficients
        vector<vector<double>> bezalfs(degree + times + 1, vector<double>(degree + 1));
        bezalfs[0][0] = bezalfs[ph][degree] = 1.0;

        for (int i = 1; i <= ph2; i++)
        {
            double inv = 1.0 / MathUtils::Binomial(ph, i);
            int mpi = min(degree, i);

            for (int j = max(0, i - times); j <= mpi; j++)
            {
                bezalfs[i][j] = inv * MathUtils::Binomial(degree, j) * MathUtils::Binomial(times, i - j);
            }
        }

        for (int i = ph2 + 1; i <= ph - 1; i++)
        {
            int mpi = min(degree, i);
            for (int j = max(0, i - times); j <= mpi; j++)
            {
                bezalfs[i][j] = bezalfs[ph - i][degree - j];
            }
        }

        int mh = ph;
        int kind = ph + 1;
        int r = -1;
        int a = degree;
        int b = degree + 1;
        int cind = 1;
        double ua = knotVector[0];

        int moresize = n * times * 2;
        vector<glm::dvec4> updatedControlPoints(moresize, glm::dvec4(MaxDistance, MaxDistance, MaxDistance, 1));
        updatedControlPoints[0] = controlPoints[0];

        vector<double> updatedKnotVector(moresize + ph + 1, MaxDistance);
        for (int i = 0; i <= ph; i++)
        {
            updatedKnotVector[i] = ua;
        }

        // initialize first Bezier seg
        vector<glm::dvec4> bpts(degree + 1);
        for (int i = 0; i <= degree; i++)
        {
            bpts[i] = controlPoints[i];
        }

        vector<glm::dvec4> nextbpts(degree - 1);

        while (b < m)
        {
            int i = b;
            while (b < m && MathUtils::IsAlmostEqualTo(knotVector[b], knotVector[b + 1]))
            {
                b = b + 1;
            }
            int mul = b - i + 1;
            mh += mul + times;
            double ub = knotVector[b];

            int oldr = r;
            r = degree - mul;

            int lbz = oldr > 0 ? floor((oldr + 2) / 2) : 1;
            int rbz = r > 0 ? floor(ph - (r + 1) / 2) : ph;
            //insert knot to get Bezier seg
            if (r > 0)
            {
                double numer = ub - ua;
                vector<double> alfs(degree - 1);
                for (int k = degree; k > mul; k--)
                {
                    alfs[k - mul - 1] = numer / (knotVector[a + k] - ua);
                }
                for (int j = 1; j <= r; j++)
                {
                    int save = r - j;
                    int s = mul + j;

                    for (int k = degree; k >= s; k--)
                    {
                        bpts[k] = alfs[k - s] * bpts[k] + (1.0 - alfs[k - s]) * bpts[k - 1];
                    }
                    nextbpts[save] = bpts[degree];
                }
            }

            vector<glm::dvec4> ebpts(degree + times + 1);
            for (int i = lbz; i <= ph; i++)
            {
                ebpts[i] = glm::dvec4(0);
                int mpi = min(degree, i);
                for (int j = max(0, i - times); j <= mpi; j++)
                {
                    ebpts[i] += bezalfs[i][j] * bpts[j];
                }
            }
            //
            if (oldr > 1)
            {
                int first = kind - 2;
                int last = kind;
                double den = ub - ua;
                double bet = (ub - updatedKnotVector[kind - 1]) / den;
                //knot remove
                for (int tr = 1; tr < oldr; tr++)
                {
                    int i = first;
                    int j = last;
                    int kj = j - kind + 1;

                    while (j - i > tr)
                    {
                        if (i < cind)
                        {
                            double alf = (ub - updatedKnotVector[i]) / (ua - updatedKnotVector[i]);
                            updatedControlPoints[i] = alf * updatedControlPoints[i] + (1.0 - alf) * updatedControlPoints[i - 1];
                        }

                        if (j >= lbz)
                        {
                            if (j - tr <= kind - ph + oldr)
                            {
                                double gam = (ub - updatedKnotVector[j - tr]) / den;
                                ebpts[kj] = gam * ebpts[kj] + (1.0 - gam) * ebpts[kj + 1];
                            }
                            else
                            {
                                ebpts[kj] = bet * ebpts[kj] + (1.0 - bet) * ebpts[kj + 1];
                            }
                        }

                        i = i + 1;
                        j = j - 1;
                        kj = kj - 1;
                    }

                    first -= 1;
                    last += 1;
                }
            }

            if (a != degree)
            {
                for (int i = 0; i < ph - oldr; i++)
                {
                    updatedKnotVector[kind] = ua;
                    kind++;
                }
            }

            for (int j = lbz; j <= rbz; j++)
            {
                updatedControlPoints[cind] = ebpts[j];
                cind++;
            }

            if (b < m)
            {
                for (int j = 0; j < r; j++)
                {
                    bpts[j] = nextbpts[j];
                }
                for (int j = r; j <= degree; j++)
                {
                    bpts[j] = controlPoints[b - degree + j];
                }

                a = b;
                b = b + 1;
                ua = ub;
            }
            else
            {
                for (int i = 0; i <= ph; i++)
                {
                    updatedKnotVector[kind + i] = ub;
                }
            }
        }

        for (int i = updatedControlPoints.size() - 1; i > 0; i--)
        {
            if (MathUtils::IsAlmostEqualTo(updatedControlPoints[i][0], MaxDistance) &&
                MathUtils::IsAlmostEqualTo(updatedControlPoints[i][1], MaxDistance) &&
                MathUtils::IsAlmostEqualTo(updatedControlPoints[i][2], MaxDistance))
            {
                updatedControlPoints.pop_back();
                continue;
            }
            break;
        }
        for (int i = updatedKnotVector.size() - 1; i > 0; i--)
        {
            if (MathUtils::IsAlmostEqualTo(updatedKnotVector[i], MaxDistance))
            {
                updatedKnotVector.pop_back();
                continue;
            }
            break;
        }
        result.degree = ph;
        result.knots = updatedKnotVector;
        decomposeWeightCurveCtrPnts(updatedControlPoints, result.control_points, result.weights);

        return result;
    }

    //Algo A5.4 in the NURBS book 2nd p164
    inline void refineKnotsCurve(const tinynurbs::RationalCurve3d& curve, const vector<double>& knotInserted, tinynurbs::RationalCurve3d& result)
    {
        assert(!knotInserted.empty());
        printf("refining knots...\n");
        int degree = curve.degree;
        std::vector<double> knotVector = curve.knots;
        std::vector<glm::dvec4> controlPoints = getWeightCurveCtrPnts(curve);

        int n = controlPoints.size() - 1;
        int m = n + degree + 1;
        int r = knotInserted.size() - 1;

        int a = getKnotSpanIndex(degree, knotVector, knotInserted[0]);
        int b = getKnotSpanIndex(degree, knotVector, knotInserted[r]) + 1;

        std::vector<double> insertedKnotVector(m + r + 2);
        for (int j = 0; j <= a; j++)
        {
            insertedKnotVector[j] = knotVector[j];
        }
        for (int j = b + degree; j <= m; j++)
        {
            insertedKnotVector[j + r + 1] = knotVector[j];
        }

        std::vector<glm::dvec4> updatedControlPoints(n + r + 2);
        for (int j = 0; j <= a - degree; j++)
        {
            updatedControlPoints[j] = controlPoints[j];
        }
        for (int j = b - 1; j <= n; j++)
        {
            updatedControlPoints[j + r + 1] = controlPoints[j];
        }

        int i = b + degree - 1;
        int k = b + degree + r;
        for (int j = r; j >= 0; j--)
        {
            while (knotInserted[j] <= knotVector[i] && i > a)
            {
                updatedControlPoints[k - degree - 1] = controlPoints[i - degree - 1];
                insertedKnotVector[k] = knotVector[i];
                k = k - 1;
                i = i - 1;
            }

            updatedControlPoints[k - degree - 1] = updatedControlPoints[k - degree];
            for (int l = 1; l <= degree; l++)
            {
                int ind = k - degree + l;
                double alpha = insertedKnotVector[k + l] - knotInserted[j];
                if (MathUtils::IsAlmostEqualTo(abs(alpha), 0.0))
                {
                    updatedControlPoints[ind - 1] = updatedControlPoints[ind];
                }
                else
                {
                    alpha = alpha / (insertedKnotVector[k + l] - knotVector[i - degree + l]);
                    updatedControlPoints[ind - 1] = alpha * updatedControlPoints[ind - 1] + (1.0 - alpha) * updatedControlPoints[ind];
                }
            }
            insertedKnotVector[k] = knotInserted[j];
            k = k - 1;
        }
        result.degree = degree;
        result.knots = insertedKnotVector;
        decomposeWeightCurveCtrPnts(updatedControlPoints,result.control_points,result.weights);
    }

    // inline vector<vector<double>> solveLinearSystem(const vector<vector<double>> &matrixA, const vector<vector<double>> &right)
    // {
    //     vector<vector<double>> res;

    //     return res;
    // }

    template <typename T>
    void globalInterpolation(const vector<T> &targetPnts, size_t degree, const vector<double> &params, const vector<double> &knots, tinynurbs::RationalCurve3d &resCurve)
    {
        assert(targetPnts.size() == params.size());

        size_t n = targetPnts.size() - 1;
        vector<vector<double>> matrixA(n + 1, vector<double>(n + 1, 0));
        for (size_t i = 1; i < n; i++)
        // for (size_t i = 0; i <= n; i++)
        {
            int spanIndex = getKnotSpanIndex(degree, knots, params[i]);
            vector<double> basis = getBsplineBasis(spanIndex, degree, knots, params[i]);
            for (size_t j = 0; j <= degree; j++)
            {
                matrixA[i][spanIndex - degree + j] = basis[j];
            }
        }
        matrixA[0][0] = 1.;
        matrixA[n][n] = 1.;

        int dim = T(0).length();
        assert(dim == 4);
        vector<vector<double>> right(n + 1, vector<double>(dim));
        for (size_t i = 0; i < n + 1; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                right[i][j] = targetPnts[i][j];
            }
        }

        // MathUtils aMathUtils;
        // auto resMatrix = aMathUtils.SolveLinearSystem(matrixA, right);
        auto resMatrix = MathUtils::SolveLinearSystem(matrixA, right);

        resCurve.control_points.resize(n + 1, vec3f(0));
        resCurve.weights.resize(n + 1, 1);
        for (size_t i = 0; i < resMatrix.size(); i++)
        {
            T temp(0);
            for (size_t j = 0; j < dim; j++)
            {
                temp[j] = resMatrix[i][j];
            }

            if (dim == 4)
            { // MAYBEBUG
                resCurve.weights[i] = temp.w;
                temp /= temp.w;
            }
            resCurve.control_points[i] = vec3f(temp);
        }
    }

    //inline vector<double> mergeSameSizeKnots(const tinynurbs::RationalCurve3d& a, const tinynurbs::RationalCurve3d& b)
    inline vector<double> mergeKnots(const vector<double>& a, const vector<double>& b)
    {
       //vector<double> res(b.knots);
       //int abegin = a.degree + 1;
       //int aend = a.knots.size() - 1 - a.degree - 1;
       //int bbegin = b.degree + 1;
       //int bend = b.knots.size() - 1 - b.degree - 1;
       //for (size_t ia = abegin, ib = bbegin; ia <= aend && ib <= bend; ++ia)
       //{
       //   if (MathUtils::IsAlmostEqualTo(a.knots[ia], b.knots[ib]))
       //      ++ib;
       //   else
       //      res.push_back(a.knots[ia]);
       //}

       //sort(res.begin(), res.end());
       vector<double> res;
       int ia = 0, ib = 0;
       for (; ia < a.size() && ib < b.size();)
       {
          if (MathUtils::IsAlmostEqualTo(a[ia], b[ib]))
          {
             res.push_back(a[ia]);
             ++ia;
             ++ib;
          }
          else
          {
             if (a[ia] < b[ib])
             {
                res.push_back(a[ia]);
                ++ia;
             }
             if (a[ia] > b[ib])
             {
                res.push_back(b[ib]);
                ++ib;
             }
          }
       }

       for (size_t i = ia; i < a.size(); ++i)
          res.push_back(a[i]);

       for (size_t i = ib; i < b.size(); ++i)
          res.push_back(b[i]);

       //sort(res.begin(), res.end());
       return res;
    }

    inline vector<double> getUniformKnots(const vector<tinynurbs::RationalCurve3d>& sections)
    {
       //auto setComp = [](const double& a, const double& b)->bool
       //{
       //   if (fabs(a - b) < 1e-4)
       //      return false;
       //   else
       //      return a < b;
       //};
       ////set<double, decltype(setComp)> container(setComp);
       //multiset<double, decltype(setComp)> container(setComp);

       //for (const auto& ele : sections)
       //{
       //   size_t begin = ele.degree + 1;
       //   size_t end = ele.knots.size() - 1 - ele.degree - 1;
       //   for (size_t i = begin; i <= end; i++)
       //   {
       //      container.insert(ele.knots[i]);
       //   }
       //}
       //printf("container.size: %d\n", container.size());

       //int maxDegree = sections[0].degree;
       //vector<double> uniformKnotsU(maxDegree + 1, 0);
       //for (const auto& ele : container)
       //{
       //   uniformKnotsU.push_back(ele);
       //}
       //for (size_t i = 0; i < maxDegree + 1; i++)
       //{
       //   uniformKnotsU.push_back(1.);
       //}
       vector<double> uniformKnotsU = sections[0].knots;
       for (size_t i = 1; i < sections.size(); ++i)
       {
          uniformKnotsU = mergeKnots(uniformKnotsU, sections[i].knots);
       }
       //printf("uniformKnotsU.size: %d\n", uniformKnotsU.size());
       //for(const auto& ele: uniformKnotsU)
       //   printf("%.3lf ",ele);
       //printf("\n");
       return uniformKnotsU;
    }

    inline vector<double> getKnotsInserted(const vector<double>& uniformKnots, const vector<double>& origin)
    {
       //for (size_t i = 0, ii = 0; i < uniformKnotsU.size() && ii < ele.knots.size(); i++)
       //{
       //   if (MathUtils::IsAlmostEqualTo(uniformKnotsU[i], ele.knots[ii]))
       //      ++ii;
       //   else
       //      knotsInserted.push_back(uniformKnotsU[i]);
       //}
       vector<double> res;
       int ui = 0, oi = 0;
       for (; ui < uniformKnots.size() && oi < origin.size();)
       {
          if (MathUtils::IsAlmostEqualTo(uniformKnots[ui], origin[oi]))
          {
             ++ui;
             ++oi;
          }
          else
          {
             if (uniformKnots[ui] < origin[oi])
             {
                res.push_back(uniformKnots[ui]);
                ++ui;
             }
             if (uniformKnots[ui] > origin[oi])
             {
                ++oi;
             }
          }
       }
       for (; ui < uniformKnots.size(); ++ui)
       {
          res.push_back(uniformKnots[ui]);
       }
       return res;
    }

    inline void createLoftSurface(const vector<tinynurbs::RationalCurve3d> &inputSections, tinynurbs::RationalSurface3d &result)
    {
        typedef glm::dvec4 vec4d;
        vector<tinynurbs::RationalCurve3d> sections = inputSections;
        size_t sectionSize = sections.size();
        int k = sectionSize - 1;

        int maxDegree = 0;
        for (const auto &it : sections)
        {
            if (it.degree > maxDegree)
                maxDegree = it.degree;
        }
        result.degree_u = maxDegree;
        result.degree_v = maxDegree > sectionSize ? sectionSize : maxDegree;//degree_v should be <= sectionSize

        /******preprocess: make sure uniform degree and knots******/
        ///1. uniform degree
        for (auto &ele : sections)
        {
            if (ele.degree < maxDegree)
            {
                printf("elevating curve degree...\n");
                auto curveRes = elevateCurveDegree(ele, maxDegree - ele.degree);
                ele = curveRes;
                assert(ele.degree == curveRes.degree && ele.knots == curveRes.knots);
                assert(ele.degree == maxDegree);
            }
        }
        ///2. get uniform knots after elevating degree
        //MAYBEBUG
        auto uniformKnotsU = getUniformKnots(sections);

         //result.knots_u = sections[0].knots;
        result.knots_u = uniformKnotsU;

        ///3. refine each section's knots
        for (auto& ele: sections)
        {
            if (ele.knots != uniformKnotsU)
            {
                assert(uniformKnotsU.size() >= ele.knots.size());
                
                vector<double> knotsInserted;
                knotsInserted = getKnotsInserted(uniformKnotsU, ele.knots);//MAYBEBUG

                tinynurbs::RationalCurve3d resCurve;
                refineKnotsCurve(ele, knotsInserted, resCurve);
                ele = resCurve;
            }
        }
        /******end of preprocess******/

        ///4. build skinning surface

        printf("4.0\n");
        vector<vector<vec4d>> wCtrPointsVector2;
        for (const auto &ele : sections)
        {
            assert(ele.degree == maxDegree);
            assert(ele.knots == uniformKnotsU);
            auto tempWCtrPoints = getWeightCurveCtrPnts(ele);
            // for (size_t i = 0; i < ele.control_points.size(); i++)
            // {
            //     tempWCtrPoints.emplace_back(vec4d(ele.control_points[i] * ele.weights[i], ele.weights[i]));
            // }
            wCtrPointsVector2.emplace_back(tempWCtrPoints); // MAYBEBUG
        }

        ///4.1 compute param vi of ith curve
        printf("4.1\n");
        vector<double> vParaVector(sectionSize);
        vParaVector[0] = 0;
        vParaVector[k] = 1;
        size_t ctrPntSize = wCtrPointsVector2[0].size();
        size_t n = ctrPntSize - 1;

        vector<double> di(ctrPntSize, 0);
        for (size_t i = 0; i < ctrPntSize; ++i) // get di
        {
            di[i] = getChordLength(wCtrPointsVector2, i);
        }
        for (size_t i = 1; i <= k - 1; i++)
        {
            double sumOfChordLengthProportion = 0;
            for (size_t j = 0; j <= n; j++)
            {
                auto tChordLength = glm::distance(wCtrPointsVector2[i][j], wCtrPointsVector2[i - 1][j]);
                sumOfChordLengthProportion += tChordLength / di[j];
            }
            // equation 10.8 in The NURBS book 2nd.
            // vParaVector[i] = vParaVector[i - 1] + 1. / (ctrPntSize + 1) * sumOfChordLengthProportion;
            vParaVector[i] = vParaVector[i - 1] + 1. / (n + 1) * sumOfChordLengthProportion;
        }

        ///4.2 average knot v
        printf("4.2\n");
        vector<double> knotsV(sectionSize + result.degree_v + 1, 0);
        for (size_t i = knotsV.size() - result.degree_v - 1; i < knotsV.size(); i++)
        {
            knotsV[i] = 1;
        }
        // equation 9.8 in the NURBS book 2nd
        // for (size_t i = result.degree_v+1; i <= ctrPntSize; i++)
        //{
        //    double sumOfVParam=0;
        //    for (size_t j = i-result.degree_v; j <= i-1; j++)
        //    {
        //        sumOfVParam+=vParaVector[j];
        //    }
        //    knotsV[i]=sumOfVParam/result.degree_v;
        //}
        int vCtrPntsNum = vParaVector.size();
        for (size_t i = 1; i <= vCtrPntsNum - 1 - result.degree_v; ++i)
        {
            double sunOfVParam = 0;
            for (size_t j = i; j <= i + result.degree_v - 1; ++j)
            {
                sunOfVParam += vParaVector[j];
            }
            knotsV[i + result.degree_v] = sunOfVParam / result.degree_v;
        }
        result.knots_v = knotsV;

        ///4.3 interpolation
        printf("4.3\n");
        result.control_points.resize(ctrPntSize, vCtrPntsNum);
        result.weights.resize(ctrPntSize, vCtrPntsNum);

        for (size_t i = 0; i < ctrPntSize; i++)
        {  
            printf("global interpolating...\n");
            // vector<vec3f> dataPnts;
            // for (size_t j = 0; j < sectionSize; j++)
            //     dataPnts.emplace_back(sections[j].control_points[i]);
            vector<vec4d> dataPnts;
            for (size_t j = 0; j < sectionSize; j++)
                dataPnts.emplace_back(wCtrPointsVector2[j][i]);

            tinynurbs::RationalCurve3d resCurve;
            globalInterpolation<vec4d>(dataPnts, result.degree_v, vParaVector, result.knots_v, resCurve);

            for (size_t j = 0; j < resCurve.control_points.size(); j++)
            {
                result.control_points(i, j) = resCurve.control_points[j];
                result.weights(i, j) = resCurve.weights[j];
            }
        }
    }

} // end of namespace myNurbs
