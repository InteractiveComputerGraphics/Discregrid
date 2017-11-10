#pragma once

#include "discrete_grid.hpp"

namespace Discregrid
{

class CubicLagrangeDiscreteGrid : public DiscreteGrid
{
public:

	CubicLagrangeDiscreteGrid(std::string const& filename);
	CubicLagrangeDiscreteGrid(Eigen::AlignedBox3d const& domain,
		std::array<unsigned int, 3> const& resolution);

	void save(std::string const& filename) const override;
	void load(std::string const& filename) override;

	unsigned int addFunction(ContinuousFunction const& func, bool verbose = false,
		SamplePredicate const& pred = nullptr) override;


	std::size_t nCells() const { return m_n_cells; };
	double interpolate(unsigned int field_id, Eigen::Vector3d const& xi,
		Eigen::Vector3d* gradient = nullptr) const override;

	void reduceField(unsigned int field_id, Predicate pred) override;

	void forEachCell(unsigned int field_id,
		std::function<void(unsigned int, Eigen::AlignedBox3d const&, unsigned int)> const& cb) const;

private:

	Eigen::Vector3d indexToNodePosition(unsigned int l) const;


private:

	std::vector<std::vector<double>> m_nodes;
	std::vector<std::vector<std::array<unsigned int, 32>>> m_cells;
	std::vector<std::vector<unsigned int>> m_cell_map;
};

}