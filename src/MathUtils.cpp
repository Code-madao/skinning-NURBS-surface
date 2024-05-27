﻿#include "../include/MathUtils.h"
#include <limits>

void GetCoFactor(const std::vector<std::vector<double>> &matrix, std::vector<std::vector<double>> &temp, int row, int column, int n)
{
    int i = 0, j = 0;
    for (int r = 0; r < n; r++)
    {
        for (int col = 0; col < n; col++)
        {
            if (r != row && col != column)
            {
                temp[i][j++] = matrix[r][col];
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

void GetAdjointMatrix(const std::vector<std::vector<double>> &matrix, std::vector<std::vector<double>> &adjointMatrix)
{
    int n = matrix.size();
    if (n == 1)
    {
        adjointMatrix[0][0] = 1;
        return;
    }

    int sign = 1;
    std::vector<std::vector<double>> temp(n, std::vector<double>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            GetCoFactor(matrix, temp, i, j, n);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            adjointMatrix[j][i] = (sign) * (MathUtils::GetDeterminant(temp, n - 1));
        }
    }
}

bool MathUtils::IsAlmostEqualTo(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    double eps = (abs(value1) + abs(value2) + 10) * tolerance;
    double delta = value1 - value2;
    return (-eps < delta) && (eps > delta);
}

bool MathUtils::IsGreaterThan(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return value1 > value2 && !IsAlmostEqualTo(value1, value2, tolerance);
}

bool MathUtils::IsGreaterThanOrEqual(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return (value1 - value2 > tolerance) || IsAlmostEqualTo(value1, value2, tolerance);
}

bool MathUtils::IsLessThan(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return value1 < value2 && !IsAlmostEqualTo(value1, value2, tolerance);
}

bool MathUtils::IsLessThanOrEqual(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return value1 < value2 || IsAlmostEqualTo(value1, value2, tolerance);
}

bool MathUtils::IsInfinite(double value)
{
    constexpr double maxValue = (std::numeric_limits<double>::max)();
    double minValue = -maxValue;
    return !(minValue <= value && value <= maxValue);
}

bool MathUtils::IsNaN(double value)
{
    return value != value;
}

double MathUtils::RadiansToAngle(double radians)
{
    return radians * 180.0 / Pi;
}

double MathUtils::AngleToRadians(double angle)
{
    return angle * Pi / 180.0;
}

int MathUtils::Factorial(int number)
{
    if (number == 0)
        return 1;
    else
        return number * Factorial(number - 1);
}

double MathUtils::Binomial(int number, int i)
{
    // MAYBEBUG:
    return Factorial(number) / (Factorial(i) * Factorial(number - i));
    // return Factorial(number) / (Factorial(i) * Factorial(number - 1));
}

double MathUtils::ComputerCubicEquationsWithOneVariable(double cubic, double quadratic, double linear, double constant)
{
    double result;
    double initial = 0.001;
    result = initial - ((cubic * pow(initial, 3) + quadratic * pow(initial, 2) + linear * initial + constant) / (3 * cubic * pow(initial, 2) + 2 * quadratic * initial + linear));
    while (MathUtils::IsGreaterThan(abs(result - initial), DoubleEpsilon))
    {
        initial = result;
        result = initial - ((cubic * pow(initial, 3) + quadratic * pow(initial, 2) + linear * initial + constant) / (3 * cubic * pow(initial, 2) + 2 * quadratic * initial + linear));
    }
    return result;
}

std::vector<std::vector<double>> MathUtils::MatrixMultiply(const std::vector<std::vector<double>> &left, const std::vector<std::vector<double>> &right)
{
    int m = left.size();
    int n = left[0].size();
    int p = right[0].size();

    std::vector<std::vector<double>> result(m, std::vector<double>(p, 0.0));
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < p; j++)
        {
            for (int k = 0; k < n; k++)
            {
                result[i][j] += left[i][k] * right[k][j];
            }
        }
    }
    return result;
}

std::vector<std::vector<double>> MathUtils::MakeDiagonal(int size)
{
    std::vector<std::vector<double>> result(size);
    for (int i = 0; i < size; i++)
    {
        result[i].resize(size);
        for (int j = 0; j < size; j++)
        {
            if (i == j)
            {
                result[i][j] = 1;
            }
            else
            {
                result[i][j] = 0;
            }
        }
    }
    return result;
}

std::vector<std::vector<double>> MathUtils::CreateMatrix(int row, int column)
{
    std::vector<std::vector<double>> result;
    for (int i = 0; i < row; i++)
    {
        std::vector<double> v(column, 0);
        result.emplace_back(v);
    }
    return result;
}

bool MathUtils::IsSquareMatrix(const std::vector<std::vector<double>> &matrix)
{
    int row = matrix.size();
    int column = matrix[0].size();
    return row == column;
}

double MathUtils::GetDeterminant(const std::vector<std::vector<double>> &matrix, int dimension)
{
    if (!IsSquareMatrix(matrix))
    {
        return 0.0;
    }

    std::vector<std::vector<double>> temp = matrix;
    double result = 0.0;
    if (dimension == 1)
    {
        return matrix[0][0];
    }

    int sign = 1;
    for (int f = 0; f < dimension; f++)
    {
        GetCoFactor(matrix, temp, 0, f, dimension);
        result += sign * matrix[0][f] * GetDeterminant(temp, dimension - 1);
        sign = -sign;
    }

    return result;
}

bool MathUtils::MakeInverse(const std::vector<std::vector<double>> &matrix, std::vector<std::vector<double>> &inverse)
{
    if (!IsSquareMatrix(matrix))
    {
        return false;
    }

    int n = matrix.size();
    double det = GetDeterminant(matrix, n);
    if (IsAlmostEqualTo(det, 0.0))
    {
        if (n == 1 && !IsAlmostEqualTo(matrix[0][0], 0.0))
        {
            inverse = {{1.0 / (double)matrix[0][0]}};
            return true;
        }
        else
        {
            return false;
        }
    }

    std::vector<std::vector<double>> lower;
    std::vector<std::vector<double>> upper;
    if (LUDecomposition(matrix, lower, upper))
    {
        std::vector<std::vector<double>> inverseLower(n, std::vector<double>(n));
        std::vector<std::vector<double>> inverseUpper(n, std::vector<double>(n));

        for (int i = 0; i < n; i++)
        {
            inverseUpper[i][i] = 1 / upper[i][i];
            for (int k = i - 1; k >= 0; k--)
            {
                double s = 0;
                for (int j = k + 1; j <= i; j++)
                {
                    s = s + upper[k][j] * inverseUpper[j][i];
                }
                if (IsAlmostEqualTo(abs(s), 0.0))
                {
                    inverseUpper[k][i] = 0.0;
                }
                else
                {
                    inverseUpper[k][i] = -s / upper[k][k];
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            inverseLower[i][i] = 1;
            for (int k = i + 1; k < n; k++)
            {
                for (int j = i; j <= k - 1; j++)
                {
                    double temp = inverseLower[k][i] - lower[k][j] * inverseLower[j][i];
                    if (IsAlmostEqualTo(temp, 0.0))
                    {
                        inverseLower[k][i] = 0.0;
                    }
                    else
                    {
                        inverseLower[k][i] = temp;
                    }
                }
            }
        }

        inverse = MatrixMultiply(inverseUpper, inverseLower);
        return true;
    }
    else
    {
        bool rs = false;
        inverse.resize(n);
        for (int k = 0; k < n; k++)
        {
            inverse[k].resize(n);
        }
        for (int k = 0; k < n; k++)
        {
            std::vector<double> b(n, 0.0);
            b[k] = 1;
            std::vector<double> pivot;
            if (!LUPDecomposition(matrix, lower, upper, pivot))
            {
                rs = false;
                break;
            }
            std::vector<double> x(n);
            std::vector<double> y(n);

            for (int i = 0; i < n; i++)
            {
                y[i] = b[static_cast<int>(pivot[i])];
                for (int j = 0; j < i; j++)
                {
                    y[i] = y[i] - lower[i][j] * y[j];
                }
            }
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = y[i];
                for (int j = n - 1; j > i; j--)
                {
                    x[i] = x[i] - upper[i][j] * x[j];
                }
                x[i] /= upper[i][i];
            }

            for (int i = 0; i < n; i++)
            {
                inverse[i][k] = x[i];
            }
        }
        if (rs)
        {
            return true;
        }
        else
        {
            std::vector<std::vector<double>> adjoint(n, std::vector<double>(n));
            GetAdjointMatrix(matrix, adjoint);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    inverse[i][j] = adjoint[i][j] / float(det);
                }
            }
            return true;
        }
    }
}

bool MathUtils::LUDecomposition(const std::vector<std::vector<double>> &matrix, std::vector<std::vector<double>> &lowerTriMatrix, std::vector<std::vector<double>> &upperTriMatrix)
{
    if (!IsSquareMatrix(matrix))
    {
        return false;
    }

    int n = matrix.size();
    lowerTriMatrix.resize(n);
    upperTriMatrix.resize(n);
    for (int i = 0; i < n; i++)
    {
        lowerTriMatrix[i].resize(n);
        upperTriMatrix[i].resize(n);
        for (int j = 0; j < n; j++)
        {
            upperTriMatrix[i][j] = 0;
            if (i == j)
            {
                lowerTriMatrix[i][j] = 1.0;
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = i; j < n; j++)
        {
            for (int k = 0; k <= i - 1; k++)
            {
                sum += lowerTriMatrix[i][k] * upperTriMatrix[k][j];
            }
            double temp = matrix[i][j] - sum;
            if (IsAlmostEqualTo(temp, 0.0))
            {
                if (i == j)
                    return false;
                upperTriMatrix[i][j] = 0.0;
            }
            else
            {
                upperTriMatrix[i][j] = temp;
            }
            sum = 0.0;
        }

        for (int x = i + 1; x < n; x++)
        {
            for (int k = 0; k <= i - 1; k++)
            {
                sum += lowerTriMatrix[x][k] * upperTriMatrix[k][i];
            }
            double temp = matrix[x][i] - sum;
            if (IsAlmostEqualTo(temp, 0.0))
            {
                lowerTriMatrix[x][i] = 0.0;
            }
            else
            {
                lowerTriMatrix[x][i] = temp / upperTriMatrix[i][i];
            }
            sum = 0.0;
        }
    }
    return true;
}

bool MathUtils::LUPDecomposition(const std::vector<std::vector<double>> &matrix, std::vector<std::vector<double>> &lowerTriMatrix, std::vector<std::vector<double>> &upperTriMatrix, std::vector<double> &pivot)
{
    if (!IsSquareMatrix(matrix))
    {
        return false;
    }

    int n = matrix.size();
    std::vector<std::vector<double>> copy(n, std::vector<double>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            copy[i][j] = matrix[i][j];
            lowerTriMatrix[i][j] = 0.0;
            upperTriMatrix[i][j] = 0.0;
        }
    }
    pivot.resize(n);
    for (int i = 0; i < n; i++)
    {
        pivot[i] = i;
    }

    int row = 0;
    for (int i = 0; i < n - 1; i++)
    {
        double p = 0.0;
        for (int j = i; j < n; j++)
        {
            if (IsGreaterThan(abs(copy[j][i]), p))
            {
                p = abs(copy[j][i]);
                row = j;
            }
        }
        if (IsAlmostEqualTo(p, 0.0))
        {
            return false;
        }

        int tmp = pivot[i];
        pivot[i] = pivot[row];
        pivot[row] = tmp;

        double tmp2 = 0.0;
        for (int j = 0; j < n; j++)
        {
            tmp2 = copy[i][j];
            copy[i][j] = copy[row][j];
            copy[row][j] = tmp2;
        }

        double u = copy[i][i];
        double l = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            l = copy[j][i] / u;
            copy[j][i] = l;
            for (int k = i + 1; k < n; k++)
            {
                copy[j][k] = copy[j][k] - copy[i][k] * l;
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            if (i != j)
            {
                lowerTriMatrix[i][j] = copy[i][j];
            }
            else
            {
                lowerTriMatrix[i][j] = 1;
            }
        }
        for (int k = i; k < n; k++)
        {
            upperTriMatrix[i][k] = copy[i][k];
        }
    }
    return true;
}

std::vector<std::vector<double>> MathUtils::SolveLinearSystem(const std::vector<std::vector<double>> &matrix, const std::vector<std::vector<double>> &right)
{
    std::vector<std::vector<double>> result;
    std::vector<std::vector<double>> inverse;
    bool canInverse = MakeInverse(matrix, inverse);
    if (canInverse)
    {
        result = MatrixMultiply(inverse, right);
    }
    return result;
}
