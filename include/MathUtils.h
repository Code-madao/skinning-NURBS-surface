#pragma once

#include <vector>

const double DefaultTolerance = 0.1;
const double DoubleEpsilon = 1E-6;
const double DistanceEpsilon = 1E-4;
const double MaxDistance = 1E9;
const double Pi = 3.14159265358979323846;

namespace MathUtils
{
	 bool IsAlmostEqualTo(double value1, double value2, double tolerance = DoubleEpsilon);

	 bool IsGreaterThan(double value1, double value2, double tolerance = DoubleEpsilon);

	 bool IsGreaterThanOrEqual(double value1, double value2, double tolerance = DoubleEpsilon);

	 bool IsLessThan(double value1, double value2, double tolerance = DoubleEpsilon);

	 bool IsLessThanOrEqual(double value1, double value2, double tolerance = DoubleEpsilon);

	 bool IsInfinite(double value);

	 bool IsNaN(double value);

	 double RadiansToAngle(double radians);

	 double AngleToRadians(double angle);

	 int Factorial(int number);

	 double Binomial(int number, int i);

	/// <summary>
	/// The NURBS Book 2nd Edition Page445
	/// Equation 9.102.
	/// </summary>
	double ComputerCubicEquationsWithOneVariable(double cubic, double quadratic, double linear, double constant);

	template <typename T>
	 void Transpose(const std::vector<std::vector<T>> &matrix, std::vector<std::vector<T>> &transposed)
	{
		std::vector<T> temp;

		for (int i = 0; i < matrix[0].size(); i++)
		{
			for (int j = 0; j < matrix.size(); j++)
			{
				temp.emplace_back(matrix[j][i]);
			}
			transposed.emplace_back(temp);
			temp.erase(temp.begin(), temp.end());
		}
	}

	template <typename T>
	 std::vector<T> GetColumn(const std::vector<std::vector<T>> &matrix, int columnIndex)
	{
		int size = matrix.size();
		std::vector<T> result(size);
		for (int i = 0; i < size; i++)
		{
			result[i] = matrix[i][columnIndex];
		}
		return result;
	}

	 std::vector<std::vector<double>> MatrixMultiply(const std::vector<std::vector<double>> &left, const std::vector<std::vector<double>> &right);

	 std::vector<std::vector<double>> MakeDiagonal(int size);

	 std::vector<std::vector<double>> CreateMatrix(int row, int column);

	 bool IsSquareMatrix(const std::vector<std::vector<double>> &matrix);

	 double GetDeterminant(const std::vector<std::vector<double>> &matrix, int dimension);

	 bool MakeInverse(const std::vector<std::vector<double>> &matrix, std::vector<std::vector<double>> &inverse);

	 bool LUDecomposition(const std::vector<std::vector<double>> &matrix, std::vector<std::vector<double>> &lowerTriMatrix, std::vector<std::vector<double>> &upperTriMatrix);

	 bool LUPDecomposition(const std::vector<std::vector<double>> &matrix, std::vector<std::vector<double>> &lowerTriMatrix, std::vector<std::vector<double>> &upperTriMatrix, std::vector<double> &pivot);

	// matrix * result = right.
	 std::vector<std::vector<double>> SolveLinearSystem(const std::vector<std::vector<double>> &matrix, const std::vector<std::vector<double>> &right);
};
