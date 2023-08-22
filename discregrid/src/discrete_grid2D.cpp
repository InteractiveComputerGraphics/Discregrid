#include <discrete_grid2D.hpp>

using namespace Eigen;

namespace Discregrid
{


DiscreteGrid2D::MultiIndex
DiscreteGrid2D::singleToMultiIndex(unsigned int l) const
{
	const auto n0 = m_resolution[0];
	const auto j = l / n0;
	const auto i = l % n0;
	return {{i, j}};
}

unsigned int
DiscreteGrid2D::multiToSingleIndex(MultiIndex const & ij) const
{
	return ij[1] * m_resolution[0] + ij[0];
}

AlignedBox2d
DiscreteGrid2D::subdomain(MultiIndex const& ij) const
{
	auto origin = m_domain.min() + Map<Matrix<unsigned int, 2, 1> const>(
		ij.data()).cast<double>().cwiseProduct(m_cell_size);
	return { origin, origin + m_cell_size};
}

AlignedBox2d
DiscreteGrid2D::subdomain(unsigned int l) const
{
	return subdomain(singleToMultiIndex(l));
}


}