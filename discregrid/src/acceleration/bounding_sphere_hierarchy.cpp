
#include <acceleration/bounding_sphere_hierarchy.hpp>

#include <iostream>
#include <set>

using namespace Eigen;

namespace Discregrid
{

TriangleMeshBSH::TriangleMeshBSH(
	std::vector<Vector3d> const& vertices,
	std::vector<std::array<unsigned int, 3>> const& faces)
	: super(faces.size()), m_faces(faces), m_vertices(vertices), 
		m_tri_centers(faces.size())
{
	std::transform(m_faces.begin(), m_faces.end(), m_tri_centers.begin(),
		[&](std::array<unsigned int, 3> const& f)
		{
			return 1.0 / 3.0 * (m_vertices[f[0]] + m_vertices[f[1]] +
				m_vertices[f[2]]);
		});
}

Vector3d const&
TriangleMeshBSH::entityPosition(unsigned int i) const
{
	return m_tri_centers[i];
}

struct TBSHCoordAccessor 
{
	TBSHCoordAccessor(std::vector<Vector3d> const* vertices_,
		std::set<unsigned int> const* indices_) 
		: indices(indices_), vertices(vertices_)
	{
	}

	using Pit = std::set<unsigned int>::const_iterator;
	using Cit = double const*;
	inline  Cit operator() (Pit it) const { 
		return (*vertices)[*it].data(); 
	}
	std::set<unsigned int> const* indices;
	std::vector<Vector3d> const* vertices;
};

void
TriangleMeshBSH::computeHull(unsigned int b, unsigned int n, BoundingSphere& hull) const
{
	auto vertices_subset = std::vector<Vector3d>(3 * n);
	for (unsigned int i(0); i < n; ++i)
	{
		auto const& f = m_faces[m_lst[b + i]];
		{
			vertices_subset[3 * i + 0] = m_vertices[f[0]];
			vertices_subset[3 * i + 1] = m_vertices[f[1]];
			vertices_subset[3 * i + 2] = m_vertices[f[2]];
		}
	}

	const BoundingSphere s(vertices_subset);

	hull.x() = s.x();
	hull.r() = s.r();
}

TriangleMeshBBH::TriangleMeshBBH(
	std::vector<Vector3d> const& vertices,
	std::vector<std::array<unsigned int, 3>> const& faces)
	: super(faces.size()), m_faces(faces), m_vertices(vertices), 
		m_tri_centers(faces.size())
{
	std::transform(m_faces.begin(), m_faces.end(), m_tri_centers.begin(),
		[&](std::array<unsigned int, 3> const& f)
		{
			return 1.0 / 3.0 * (m_vertices[f[0]] + m_vertices[f[1]] +
				m_vertices[f[2]]);
		});
}

Vector3d const&
TriangleMeshBBH::entityPosition(unsigned int i) const
{
	return m_tri_centers[i];
}

void
TriangleMeshBBH::computeHull(unsigned int b, unsigned int n, AlignedBox3d& hull) const
{
	for (auto i = 0u; i < n; ++i)
	{
		auto const& f = m_faces[m_lst[b + i]];
		for (auto v : f)
		{
			hull.extend(m_vertices[v]);
		}
	}
}


PointCloudBSH::PointCloudBSH()
	: super(0)
{
}

PointCloudBSH::PointCloudBSH(std::vector<Vector3d> const& vertices)
	: super(vertices.size()), m_vertices(&vertices)
{

}

Vector3d const&
PointCloudBSH::entityPosition(unsigned int i) const
{
	return (*m_vertices)[i];
}

struct PCBSHCoordAccessor
{
	PCBSHCoordAccessor(std::vector<Vector3d> const* vertices_,
		std::vector<unsigned int> const* lst_)
		: vertices(vertices_), lst(lst_)
	{
	}

	using Pit = unsigned int;
	using Cit = double const*;
	inline  Cit operator() (Pit it) const {
		return (*vertices)[(*lst)[it]].data();
	}
	std::vector<Vector3d> const* vertices;
	std::vector<unsigned int> const* lst;
};


void
PointCloudBSH::computeHull(unsigned int b, unsigned int n, BoundingSphere& hull) const
{
	auto vertices_subset = std::vector<Vector3d>(3 * n);
	for (unsigned int i(0); i < n; ++i)
		vertices_subset[3 * i + 0] = (*m_vertices)[m_lst[i]];

	const BoundingSphere s(vertices_subset);

	hull.x() = s.x();
	hull.r() = s.r();
}



}
