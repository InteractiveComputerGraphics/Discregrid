#include "cubic_lagrange_discrete_grid2D.hpp"

#include "data/z_sort_table.hpp"
#include <utility/serialize.hpp>
#include "utility/spinlock.hpp"

#include <iostream>
#include <iomanip>
#include <atomic>
#include <numeric>
#include <set>
#include <chrono>
#include <future>

using namespace Eigen;

namespace Discregrid
{

namespace
{

// See Zienkiewicz The Finite Element Method its Basis and Fundamentals 2013 page 159 (2D) and 168 (3D)
Matrix<double, 12, 1>
shape_function_(Vector2d const &xi, Matrix<double, 12, 2> *gradient = nullptr)
{
	auto res = Matrix<double, 12, 1>{};

	const auto x = xi[0];
	const auto y = xi[1];

	const auto x2 = x * x;
	const auto y2 = y * y;
	
	const auto _1mx = 1.0 - x;
	const auto _1my = 1.0 - y;
	
	const auto _1px = 1.0 + x;
	const auto _1py = 1.0 + y;
	
	const auto _1m3x = 1.0 - 3.0 * x;
	const auto _1m3y = 1.0 - 3.0 * y;
	
	const auto _1p3x = 1.0 + 3.0 * x;
	const auto _1p3y = 1.0 + 3.0 * y;
	
	const auto _1mxt1my = _1mx * _1my;
	const auto _1mxt1py = _1mx * _1py;
	const auto _1pxt1my = _1px * _1my;
	const auto _1pxt1py = _1px * _1py;
	
	const auto _1mx2 = 1.0 - x2;
	const auto _1my2 = 1.0 - y2;

	// Corner nodes.
	auto factor = (1.0 / 32.0) * (9.0 * (x2 + y2) - 10.0);
	res[0] = factor * _1mxt1my;
	res[1] = factor * _1pxt1my;
	res[2] = factor * _1mxt1py;
	res[3] = factor * _1pxt1py;

	// Edge nodes.
	factor = (9.0 / 32.0) * _1mx2;
	auto factort1m3x = factor * _1m3x;	// Note that x is multiplied with xi (+-1/3) so 9 * xi * x becomes 1/3 * x
	auto factort1p3x = factor * _1p3x;
	res[4] = factort1m3x * _1my;
	res[5] = factort1p3x * _1my;
	res[6] = factort1m3x * _1py;
	res[7] = factort1p3x * _1py;

	factor = (9.0 / 32.0) * _1my2;
	auto factort1m3y = factor * _1m3y;
	auto factort1p3y = factor * _1p3y;
	res[8]  = factort1m3y * _1mx;
	res[9]  = factort1p3y * _1mx;
	res[10] = factort1m3y * _1px;
	res[11] = factort1p3y * _1px;

	if (gradient)
	{
		auto &dN = *gradient;

		const auto _27tx2p9ty2m10 = 27 * x2 + 9 * y2 - 10;
		const auto _27ty2p9tx2m10 = 27 * y2 + 9 * x2 - 10;

		const auto _18x = 18 * x;
		const auto _18y = 18 * y;

		dN(0, 0) = _1my * (_18x - _27tx2p9ty2m10);
		dN(0, 1) = _1mx * (_18y - _27ty2p9tx2m10);
		dN(1, 0) = _1my * (_18x + _27tx2p9ty2m10);
		dN(1, 1) = _1px * (_18y - _27ty2p9tx2m10);
		dN(2, 0) = _1py * (_18x - _27tx2p9ty2m10);
		dN(2, 1) = _1mx * (_18y + _27ty2p9tx2m10);
		dN(3, 0) = _1py * (_18x + _27tx2p9ty2m10);
		dN(3, 1) = _1px * (_18y + _27ty2p9tx2m10);

		dN.topRows(4) /= 32.0;

		const auto _9x2 = 9 * x2;
		const auto _2x = 2 * x;
		const auto _3x = 3 * x;
		const auto _3x3 = 3 * x2 * x;
		const auto _9x2m2xm3 = _9x2 - _2x - 3;
		const auto _3x3mx2m3xp1 = _3x3 - x2 - _3x + 1;
		const auto _9x2p2xm3 = _9x2 + _2x - 3;
		const auto _3x3px2m3xm1 = _3x3 + x2 - _3x - 1;
		
		dN(4, 0) = _1my * _9x2m2xm3;
		dN(4, 1) = -_3x3mx2m3xp1;
		dN(5, 0) = _1my * -_9x2p2xm3;
		dN(5, 1) = _3x3px2m3xm1;
		dN(6, 0) = _1py * _9x2m2xm3;
		dN(6, 1) = _3x3mx2m3xp1;
		dN(7, 0) = _1py * -_9x2p2xm3;
		dN(7, 1) = -_3x3px2m3xm1;
		
		const auto _9y2 = 9 * y2;
		const auto _2y = 2 * y;
		const auto _3y = 3 * y;
		const auto _3y3 = 3 * y2 * y;
		const auto _3y3my2m3yp1 = _3y3 - y2 - _3y + 1;
		const auto _9y2m2ym3 = _9y2 - _2y - 3;
		const auto _3y3py2m3ym1 = _3y3 + y2 - _3y - 1;
		const auto _9y2p2ym3 = _9y2 + _2y - 3;

		dN(8, 0) = -_3y3my2m3yp1;
		dN(8, 1) = _1mx * _9y2m2ym3;
		dN(9, 0) = _3y3py2m3ym1;
		dN(9, 1) = _1mx * -_9y2p2ym3;
		dN(10, 0) = _3y3my2m3yp1;
		dN(10, 1) = _1px * _9y2m2ym3;
		dN(11, 0) = -_3y3py2m3ym1;
		dN(11, 1) = _1px * -_9y2p2ym3;

		dN.bottomRows(12u - 4u) *= 9.0 / 32.0;
	}

	return res;
}

// Determines Morten value according to z-curve.
inline uint64_t
zValue(Vector2d const &x, double invCellSize)
{
	std::array<int, 2> key;
	for (unsigned int i(0); i < 2; ++i)
	{
		if (x[i] >= 0.0)
			key[i] = static_cast<int>(invCellSize * x[i]);
		else
			key[i] = static_cast<int>(invCellSize * x[i]) - 1;
	}

	std::array<unsigned int, 2> p = {
		static_cast<unsigned int>(static_cast<int64_t>(key[0]) - (std::numeric_limits<int>::lowest() + 1)),
		static_cast<unsigned int>(static_cast<int64_t>(key[1]) - (std::numeric_limits<int>::lowest() + 1))};

	return morton_lut(p);
}
} // namespace

	Vector2d CubicLagrangeDiscreteGrid2D::indexToNodePosition(unsigned int l) const
{
	auto x = Vector2d{};

	auto n = Matrix<unsigned int, 2, 1>::Map(m_resolution.data());

	// Amount of vertices and edges
	auto nv = (n[0] + 1) * (n[1] + 1);
	auto ne_x = n[0] * (n[1] + 1);

	auto ij = Matrix<unsigned int, 2, 1>{};
	// 0 -> nv: l is a cell's corner vertex
	if (l < nv)
	{
		ij(1) = l / (n[0] + 1);
		ij(0) = l % (n[0] + 1);

		x = m_domain.min() + m_cell_size.cwiseProduct(ij.cast<double>());
	}
	// nv -> nv + 2 * ne_x: l is a non-corner vertex in x direction
	else if (l < nv + 2 * ne_x)
	{
		l -= nv;
		// Identify cell index, then split into x and y dimensions.
		auto e_ind = l / 2;
		ij(1) = e_ind / n[0];
		ij(0) = e_ind % n[0];
		// Calculate cell position
		x = m_domain.min() + m_cell_size.cwiseProduct(ij.cast<double>());
		// Add sub-cell node position to cell position: Either 1/3 or 2/3 of cell size in x direction
		x(0) += (1.0 + static_cast<double>(l % 2)) / 3.0 * m_cell_size[0];
	}
	// nv + 2 * ne_x -> nv + 2 * (ne_x + ne_y): l is  non-corner vertex in y direction
	else
	{
		l -= nv + 2 * ne_x;
		// Identify cell index, then split into x and y dimensions.
		auto e_ind = l / 2;
		ij(1) = e_ind % n[1];
		ij(0) = e_ind / n[1];

		// Calculate cell position
		x = m_domain.min() + m_cell_size.cwiseProduct(ij.cast<double>());
		// Add sub-cell node position to cell position: Either 1/3 or 2/3 of cell size in y direction
		x(1) += (1.0 + static_cast<double>(l % 2)) / 3.0 * m_cell_size[1];
	}

	return x;
}

CubicLagrangeDiscreteGrid2D::CubicLagrangeDiscreteGrid2D(std::string const &filename)
{
	load(filename);
}

CubicLagrangeDiscreteGrid2D::CubicLagrangeDiscreteGrid2D(AlignedBox2d const &domain, std::array<unsigned int, 2> const &resolution) : DiscreteGrid2D(domain, resolution)
{
}

void CubicLagrangeDiscreteGrid2D::save(std::string const &filename) const
{
	auto out = std::ofstream(filename, std::ios::binary);
	serialize::write(*out.rdbuf(), m_domain);
	serialize::write(*out.rdbuf(), m_resolution);
	serialize::write(*out.rdbuf(), m_cell_size);
	serialize::write(*out.rdbuf(), m_inv_cell_size);
	serialize::write(*out.rdbuf(), m_n_cells);
	serialize::write(*out.rdbuf(), m_n_fields);

	serialize::write(*out.rdbuf(), m_nodes.size());
	for (auto const &nodes : m_nodes)
	{
		serialize::write(*out.rdbuf(), nodes.size());
		for (auto const &node : nodes)
		{
			serialize::write(*out.rdbuf(), node);
		}
	}

	serialize::write(*out.rdbuf(), m_cells.size());
	for (auto const &cells : m_cells)
	{
		serialize::write(*out.rdbuf(), cells.size());
		for (auto const &cell : cells)
		{
			serialize::write(*out.rdbuf(), cell);
		}
	}

	serialize::write(*out.rdbuf(), m_cell_map.size());
	for (auto const &maps : m_cell_map)
	{
		serialize::write(*out.rdbuf(), maps.size());
		for (auto const &map : maps)
		{
			serialize::write(*out.rdbuf(), map);
		}
	}

	out.close();
}

void CubicLagrangeDiscreteGrid2D::load(std::string const &filename)
{
	auto in = std::ifstream(filename, std::ios::binary);

	if (!in.good())
	{
		std::cerr << "ERROR: Discrete grid can not be loaded. Input file does not exist!" << std::endl;
		return;
	}

	serialize::read(*in.rdbuf(), m_domain);
	serialize::read(*in.rdbuf(), m_resolution);
	serialize::read(*in.rdbuf(), m_cell_size);
	serialize::read(*in.rdbuf(), m_inv_cell_size);
	serialize::read(*in.rdbuf(), m_n_cells);
	serialize::read(*in.rdbuf(), m_n_fields);

	auto n_nodes = std::size_t{};
	serialize::read(*in.rdbuf(), n_nodes);
	m_nodes.resize(n_nodes);
	for (auto &nodes : m_nodes)
	{
		serialize::read(*in.rdbuf(), n_nodes);
		nodes.resize(n_nodes);
		for (auto &node : nodes)
		{
			serialize::read(*in.rdbuf(), node);
		}
	}

	auto n_cells = std::size_t{};
	serialize::read(*in.rdbuf(), n_cells);
	m_cells.resize(n_cells);
	for (auto &cells : m_cells)
	{
		serialize::read(*in.rdbuf(), n_cells);
		cells.resize(n_cells);
		for (auto &cell : cells)
		{
			serialize::read(*in.rdbuf(), cell);
		}
	}

	auto n_cell_maps = std::size_t{};
	serialize::read(*in.rdbuf(), n_cell_maps);
	m_cell_map.resize(n_cell_maps);
	for (auto &cell_maps : m_cell_map)
	{
		serialize::read(*in.rdbuf(), n_cell_maps);
		cell_maps.resize(n_cell_maps);
		for (auto &cell_map : cell_maps)
		{
			serialize::read(*in.rdbuf(), cell_map);
		}
	}

	in.close();
}

unsigned int CubicLagrangeDiscreteGrid2D::addFunction(ContinuousFunction const &func, bool verbose, SamplePredicate const &pred)
{
	using namespace std::chrono;

	const auto t0_construction = high_resolution_clock::now();

	auto n = Matrix<unsigned int, 2, 1>::Map(m_resolution.data());

	// Amount of vertices and edges
	const auto nv = (n[0] + 1) * (n[1] + 1);
	const auto ne_x = n[0] * (n[1] + 1);
	const auto ne_y = (n[0] + 1) * n[1];
	const auto ne = ne_x + ne_y;
	const auto n_nodes = nv + 2 * ne;

	m_nodes.emplace_back();
	auto &coeffs = m_nodes.back();
	coeffs.resize(n_nodes);

	std::atomic_uint counter(0u);
	SpinLock mutex;
	auto t0 = high_resolution_clock::now();

	// Evaluate function at every node
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static) nowait
		for (int l = 0; l < static_cast<int>(n_nodes); ++l)
		{
			auto x = indexToNodePosition(l);
			auto &c = coeffs[l];

			if (!pred || pred(x))
				c = func(x);
			else
				c = std::numeric_limits<double>::max();

			if (verbose && (++counter == n_nodes || duration_cast<milliseconds>(high_resolution_clock::now() - t0).count() > 1000u))
			{
				std::async(std::launch::async, [&]
				{
					mutex.lock();
					t0 = high_resolution_clock::now();
					std::cout << "\r"
							  << "Construction " << std::setw(20)
							  << 100.0 * static_cast<double>(counter) / static_cast<double>(n_nodes) << "%";
					mutex.unlock();
				});
			}
		}
	}

	m_cells.emplace_back();
	auto &cells = m_cells.back();
	cells.resize(m_n_cells);
	const auto nx = n[0];
	const auto ny = n[1];
	for (auto l = 0u; l < m_n_cells; ++l)
	{
		// Build translation table for easy querying without having to recalculate the index positions
		const auto j = l / nx;
		const auto i = l % nx;

		auto &cell = cells[l];
		cell[0] = (nx + 1) * j + i;
		cell[1] = cell[0] + 1;
		cell[2] = (nx + 1) * (j + 1) + i;
		cell[3] = cell[2] + 1;

		auto offset = nv;
		cell[4] = offset + 2 * (nx * j + i);
		cell[5] = cell[4] + 1;
		cell[6] = offset + 2 * (nx * (j + 1) + i);
		cell[7] = cell[6] + 1;

		offset += 2 * ne_x;
		cell[8] = offset + 2 * (ny * i + j);
		cell[9] = cell[8] + 1;
		cell[10] = offset + 2 * (ny * (i + 1) + j);
		cell[11] = cell[10] + 1;
	}

	m_cell_map.emplace_back();
	auto &cell_map = m_cell_map.back();
	cell_map.resize(m_n_cells);
	std::iota(cell_map.begin(), cell_map.end(), 0u);

	if (verbose)
	{
		std::cout << "\rConstruction took " << std::setw(15) << static_cast<double>(duration_cast<milliseconds>(high_resolution_clock::now() - t0_construction).count()) / 1000.0 << "s" << std::endl;
	}

	return static_cast<unsigned int>(m_n_fields++);
}

bool CubicLagrangeDiscreteGrid2D::determineShapeFunctions(unsigned int field_id, Eigen::Vector2d const &x, std::array<unsigned int, 12> &cell, Eigen::Vector2d &c0, Eigen::Matrix<double, 12, 1> &N, Eigen::Matrix<double, 12, 2> *dN) const
{
	if (!m_domain.contains(x))
		return false;

	auto mi = (x - m_domain.min()).cwiseProduct(m_inv_cell_size).cast<unsigned int>().eval();
	if (mi[0] >= m_resolution[0])
		mi[0] = m_resolution[0] - 1;
	if (mi[1] >= m_resolution[1])
		mi[1] = m_resolution[1] - 1;
	auto i = multiToSingleIndex({ { mi(0), mi(1) } });
	auto i_ = m_cell_map[field_id][i];
	if (i_ == std::numeric_limits<unsigned int>::max())
		return false;

	auto sd = subdomain(i);
	i = i_;

	// Mapping to isoparametric space from world (actually local?) space
	const auto denom = (sd.max() - sd.min()).eval();
	c0 = Vector2d::Constant(2.0).cwiseQuotient(denom).eval();
	const auto c1 = (sd.max() + sd.min()).cwiseQuotient(denom).eval();
	const auto xi = (c0.cwiseProduct(x) - c1).eval();

	cell = m_cells[field_id][i];
	N = shape_function_(xi, dN);
	return true;
}

double CubicLagrangeDiscreteGrid2D::interpolate(unsigned int field_id, Eigen::Vector2d const& xi, const std::array<unsigned int, 12> &cell, const Eigen::Vector2d &c0, const Eigen::Matrix<double, 12, 1> &N, Eigen::Vector2d* gradient, Eigen::Matrix<double, 12, 2> *dN) const
{
	if (!gradient)
	{
		auto phi = 0.0;
		for (auto j = 0u; j < 12u; ++j)
		{
			auto v = cell[j];
			auto c = m_nodes[field_id][v];
			if (c == std::numeric_limits<double>::max())
			{
				return std::numeric_limits<double>::max();
			}
			phi += c * N[j];
		}

		return phi;
	}

	auto phi = 0.0;
	gradient->setZero();
	for (auto j = 0u; j < 12u; ++j)
	{
		auto v = cell[j];
		auto c = m_nodes[field_id][v];
		if (c == std::numeric_limits<double>::max())
		{
			gradient->setZero();
			return std::numeric_limits<double>::max();
		}
		phi += c * N[j];
		(*gradient)(0) += c * (*dN)(j, 0);
		(*gradient)(1) += c * (*dN)(j, 1);
	}
	gradient->array() *= c0.array();

	return phi;
}

double CubicLagrangeDiscreteGrid2D::interpolate(unsigned int field_id, Vector2d const &x, Vector2d *gradient) const
{
	if (!m_domain.contains(x))
		return std::numeric_limits<double>::max();

	auto mi = (x - m_domain.min()).cwiseProduct(m_inv_cell_size).cast<unsigned int>().eval();
	if (mi[0] >= m_resolution[0])
		mi[0] = m_resolution[0] - 1;
	if (mi[1] >= m_resolution[1])
		mi[1] = m_resolution[1] - 1;
	auto i = multiToSingleIndex({{mi(0), mi(1)}});
	auto i_ = m_cell_map[field_id][i];
	if (i_ == std::numeric_limits<unsigned int>::max())
		return std::numeric_limits<double>::max();

	auto sd = subdomain(i);
	i = i_;
	auto d = sd.diagonal().eval();
	
	auto denom = (sd.max() - sd.min()).eval();
	auto c0 = Vector2d::Constant(2.0).cwiseQuotient(denom).eval();
	auto c1 = (sd.max() + sd.min()).cwiseQuotient(denom).eval();
	auto xi = (c0.cwiseProduct(x) - c1).eval();

	auto const &cell = m_cells[field_id][i];
	if (!gradient)
	{
		auto phi = 0.0;
		auto N = shape_function_(xi, nullptr);
		for (auto j = 0u; j < 12u; ++j)
		{
			auto v = cell[j];
			auto c = m_nodes[field_id][v];
			if (c == std::numeric_limits<double>::max())
			{
				return std::numeric_limits<double>::max();
			}
			phi += c * N[j];
		}

		return phi;
	}

	auto dN = Matrix<double, 12, 2>{};
	auto N = shape_function_(xi, &dN);
	
	auto phi = 0.0;
	gradient->setZero();
	for (auto j = 0u; j < 12u; ++j)
	{
		auto v = cell[j];
		auto c = m_nodes[field_id][v];
		if (c == std::numeric_limits<double>::max())
		{
			gradient->setZero();
			return std::numeric_limits<double>::max();
		}
		phi += c * N[j];
		(*gradient)(0) += c * dN(j, 0);
		(*gradient)(1) += c * dN(j, 1);
	}
	gradient->array() *= c0.array();
	
	return phi;
}

void CubicLagrangeDiscreteGrid2D::reduceField(unsigned int field_id, Predicate pred)
{
	auto &coeffs = m_nodes[field_id];
	auto &cells = m_cells[field_id];
	auto keep = std::vector<bool>(coeffs.size());
	for (auto l = 0u; l < coeffs.size(); ++l)
	{
		auto xi = indexToNodePosition(l);
		keep[l] = pred(xi, coeffs[l]) && coeffs[l] != std::numeric_limits<double>::max();
	}

	auto &cell_map = m_cell_map[field_id];
	cell_map.resize(m_n_cells);
	std::iota(cell_map.begin(), cell_map.end(), 0u);

	const auto cells_ = cells;
	cells.clear();
	for (auto i = 0u; i < cells_.size(); ++i)
	{
		auto keep_cell = false;
		auto vals = std::vector<double>{};
		for (auto v : cells_[i])
		{
			keep_cell |= keep[v];
			vals.push_back(coeffs[v]);
		}
		if (keep_cell)
		{
			cells.push_back(cells_[i]);
			cell_map[i] = static_cast<unsigned int>(cells.size() - 1);
		}
		else
			cell_map[i] = std::numeric_limits<unsigned int>::max();
	}

	// Reduce vertices.
	auto z_values = std::vector<uint64_t>(coeffs.size());
	for (auto l = 0u; l < coeffs.size(); ++l)
	{
		auto xi = indexToNodePosition(l);
		z_values[l] = zValue(xi, 4.0 * m_inv_cell_size.minCoeff());
	}

	std::fill(keep.begin(), keep.end(), false);

	auto vertex_to_cell = std::vector<std::set<std::pair<unsigned int, unsigned int>>>(coeffs.size());
	for (auto c = 0u; c < cells.size(); ++c)
	{
		auto const &cell = cells[c];

		for (auto j = 0u; j < cell.size(); ++j)
		{
			auto v = cell[j];
			keep[v] = true;
			vertex_to_cell[v].insert({c, j});
		}
	}
	auto last_vertex = static_cast<unsigned int>(coeffs.size() - 1);
	for (auto i = static_cast<int>(coeffs.size() - 1); i >= 0; --i)
	{
		if (!keep[i])
		{
			std::swap(coeffs[i], coeffs[last_vertex]);
			std::swap(z_values[i], z_values[last_vertex]);
			std::swap(vertex_to_cell[i], vertex_to_cell[last_vertex]);
			for (auto const &kvp : vertex_to_cell[i])
			{
				cells[kvp.first][kvp.second] = i;
			}
			for (auto const &kvp : vertex_to_cell[last_vertex])
			{
				cells[kvp.first][kvp.second] = last_vertex;
			}

			last_vertex--;
		}
	}
	coeffs.resize(last_vertex + 1);
	z_values.resize(coeffs.size());

	auto sort_pattern = std::vector<unsigned int>(coeffs.size());
	std::iota(sort_pattern.begin(), sort_pattern.end(), 0u);
	std::sort(sort_pattern.begin(), sort_pattern.end(),
			  [&](unsigned int i, unsigned int j) {
				  return z_values[i] < z_values[j];
			  });

	for (auto i = 0u; i < sort_pattern.size(); ++i)
	{
		auto j = sort_pattern[i];
		for (auto const &kvp : vertex_to_cell[j])
		{
			assert(cells[kvp.first][kvp.second] == j);
			cells[kvp.first][kvp.second] = i;
		}
	}

	auto coeffs_ = coeffs;
	std::transform(sort_pattern.begin(), sort_pattern.end(), coeffs.begin(),
				   [&coeffs_](unsigned int i) { return coeffs_[i]; });
}

void CubicLagrangeDiscreteGrid2D::forEachCell(unsigned int field_id, std::function<void(unsigned int, AlignedBox2d const &, unsigned int)> const &cb) const
{
	auto n = m_resolution[0] * m_resolution[1];
	for (auto i = 0u; i < n; ++i)
	{
		auto domain = AlignedBox2d{};
		auto mi = singleToMultiIndex(i);
		domain.min() = m_domain.min() + Matrix<unsigned int, 2, 1>::Map(mi.data()).cast<double>().cwiseProduct(m_cell_size);
		domain.max() = domain.min() + m_cell_size;

		cb(i, domain, 0);
	}
}

} // namespace Discregrid
