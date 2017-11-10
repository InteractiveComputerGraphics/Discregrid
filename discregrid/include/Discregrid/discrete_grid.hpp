#pragma once

#include <vector>
#include <fstream>
#include <array>
#include <Eigen/Dense>

namespace Discregrid
{

class DiscreteGrid
{
public:

	using CoefficientVector = Eigen::Matrix<double, 32, 1>;
	using ContinuousFunction = std::function<double(Eigen::Vector3d const&)>;
	using MultiIndex = std::array<unsigned int, 3>;
	using Predicate = std::function<bool(Eigen::Vector3d const&, double)>;
	using SamplePredicate = std::function<bool(Eigen::Vector3d const&)>;

	DiscreteGrid() = default;
	DiscreteGrid(Eigen::AlignedBox3d const& domain, std::array<unsigned int, 3> const& resolution)
		: m_domain(domain), m_resolution(resolution), m_n_fields(0u)
	{
		auto n = Eigen::Matrix<unsigned int, 3, 1>::Map(resolution.data());
		m_cell_size = domain.diagonal().cwiseQuotient(n.cast<double>());
		m_inv_cell_size = m_cell_size.cwiseInverse();
		m_n_cells = n.prod();
	}
	virtual ~DiscreteGrid() = default;

	virtual void save(std::string const& filename) const = 0;
	virtual void load(std::string const& filename) = 0;

	virtual unsigned int addFunction(ContinuousFunction const& func, bool verbose = false,
		SamplePredicate const& pred = nullptr) = 0;

	double interpolate(Eigen::Vector3d const& xi, Eigen::Vector3d* gradient = nullptr) const
	{
		return interpolate(0u, xi, gradient);
	}

	virtual double interpolate(unsigned int field_id, Eigen::Vector3d const& xi,
		Eigen::Vector3d* gradient = nullptr) const = 0;

	virtual void reduceField(unsigned int field_id, Predicate pred) {}


	MultiIndex singleToMultiIndex(unsigned int i) const;
	unsigned int multiToSingleIndex(MultiIndex const& ijk) const;

	Eigen::AlignedBox3d subdomain(MultiIndex const& ijk) const;
	Eigen::AlignedBox3d subdomain(unsigned int l) const;

	Eigen::AlignedBox3d const& domain() const { return m_domain; }
	std::array<unsigned int, 3> const& resolution() const { return m_resolution; };
	Eigen::Vector3d const& cellSize() const { return m_cell_size;}
	Eigen::Vector3d const& invCellSize() const { return m_inv_cell_size;}

protected:

	Eigen::AlignedBox3d m_domain;
	std::array<unsigned int, 3> m_resolution;
	Eigen::Vector3d m_cell_size;
	Eigen::Vector3d m_inv_cell_size;
	std::size_t m_n_cells;
	std::size_t m_n_fields;

};
}