#pragma once

#include <vector>
#include <fstream>
#include <array>
#include <Eigen/Dense>

namespace Discregrid
{

class DiscreteGrid2D
{
public:

	using CoefficientVector = Eigen::Matrix<double, 12, 1>;
	using ContinuousFunction = std::function<double(Eigen::Vector2d const&)>;
	using MultiIndex = std::array<unsigned int, 2>;
	using Predicate = std::function<bool(Eigen::Vector2d const&, double)>;
	using SamplePredicate = std::function<bool(Eigen::Vector2d const&)>;

	DiscreteGrid2D() = default;
	DiscreteGrid2D(Eigen::AlignedBox2d const& domain, std::array<unsigned int, 2> const& resolution)
		: m_domain(domain), m_resolution(resolution), m_n_fields(0u)
	{
		auto n = Eigen::Matrix<unsigned int, 2, 1>::Map(resolution.data());
		m_cell_size = domain.diagonal().cwiseQuotient(n.cast<double>());
		m_inv_cell_size = m_cell_size.cwiseInverse();
		m_n_cells = n.prod();
	}
	virtual ~DiscreteGrid2D() = default;

	virtual void save(std::string const& filename) const = 0;
	virtual void load(std::string const& filename) = 0;

	virtual unsigned int addFunction(ContinuousFunction const& func, bool verbose = false,
		SamplePredicate const& pred = nullptr) = 0;

	double interpolate(Eigen::Vector2d const& xi, Eigen::Vector2d* gradient = nullptr) const
	{
		return interpolate(0u, xi, gradient);
	}

	virtual double interpolate(unsigned int field_id, Eigen::Vector2d const& xi,
		Eigen::Vector2d* gradient = nullptr) const = 0;

	/**
	 * @brief Determines the shape functions for the discretization with ID field_id at point xi.
	 * 
	 * @param field_id Discretization ID
	 * @param x Location where the shape functions should be determined
	 * @param cell cell of x
	 * @param c0 vector required for the interpolation
	 * @param N	shape functions for the cell of x
	 * @param dN (Optional) derivatives of the shape functions, required to compute the gradient
	 * @return Success of the function.
	 */
	virtual bool determineShapeFunctions(unsigned int field_id, Eigen::Vector2d const &x,
		std::array<unsigned int, 12> &cell, Eigen::Vector2d &c0, Eigen::Matrix<double, 12, 1> &N,
		Eigen::Matrix<double, 12, 2> *dN = nullptr) const = 0;

	/**
	 * @brief Evaluates the given discretization with ID field_id at point xi.
	 * 
	 * @param field_id Discretization ID
	 * @param xi Location where the discrete function is evaluated
	 * @param cell cell of xi
	 * @param c0 vector required for the interpolation
	 * @param N	shape functions for the cell of xi
	 * @param gradient (Optional) if a pointer to a vector is passed the gradient of the discrete function will be evaluated
	 * @param dN (Optional) derivatives of the shape functions, required to compute the gradient
	 * @return double Results of the evaluation of the discrete function at point xi
	 */
	virtual double interpolate(unsigned int field_id, Eigen::Vector2d const& xi, const std::array<unsigned int, 12> &cell, const Eigen::Vector2d &c0, const Eigen::Matrix<double, 12, 1> &N,
		Eigen::Vector2d* gradient = nullptr, Eigen::Matrix<double, 12, 2> *dN = nullptr) const = 0;

	virtual void reduceField(unsigned int field_id, Predicate pred) {}


	MultiIndex singleToMultiIndex(unsigned int i) const;
	unsigned int multiToSingleIndex(MultiIndex const& ij) const;

	Eigen::AlignedBox2d subdomain(MultiIndex const& ij) const;
	Eigen::AlignedBox2d subdomain(unsigned int l) const;

	Eigen::AlignedBox2d const& domain() const { return m_domain; }
	std::array<unsigned int, 2> const& resolution() const { return m_resolution; }
	Eigen::Vector2d const& cellSize() const { return m_cell_size;}
	Eigen::Vector2d const& invCellSize() const { return m_inv_cell_size;}

protected:


	Eigen::AlignedBox2d m_domain;
	std::array<unsigned int, 2> m_resolution;
	Eigen::Vector2d m_cell_size;
	Eigen::Vector2d m_inv_cell_size;
	std::size_t m_n_cells;
	std::size_t m_n_fields;
};
}