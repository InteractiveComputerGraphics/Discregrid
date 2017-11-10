#pragma once

#include "bounding_sphere.hpp"
#include "kd_tree.hpp"

namespace Discregrid
{

class TriangleMeshBSH : public KDTree<BoundingSphere>
{

public:

	using super = KDTree<BoundingSphere>;

	TriangleMeshBSH(std::vector<Eigen::Vector3d> const& vertices,
		std::vector<std::array<unsigned int, 3>> const& faces);

	Eigen::Vector3d const& entityPosition(unsigned int i) const final;
	void computeHull(unsigned int b, unsigned int n, BoundingSphere& hull) const final;

private:

	std::vector<Eigen::Vector3d> const& m_vertices;
	std::vector<std::array<unsigned int, 3>> const& m_faces;

	std::vector<Eigen::Vector3d> m_tri_centers;
};

class TriangleMeshBBH : public KDTree<Eigen::AlignedBox3d>
{
public:

	using super = KDTree<Eigen::AlignedBox3d>;

	TriangleMeshBBH(std::vector<Eigen::Vector3d> const& vertices,
		std::vector<std::array<unsigned int, 3>> const& faces);

	Eigen::Vector3d const& entityPosition(unsigned int i) const final;
	void computeHull(unsigned int b, unsigned int n, Eigen::AlignedBox3d& hull) const final;

private:

	std::vector<Eigen::Vector3d> const& m_vertices;
	std::vector<std::array<unsigned int, 3>> const& m_faces;

	std::vector<Eigen::Vector3d> m_tri_centers;


};

class PointCloudBSH : public KDTree<BoundingSphere>
{

public:

	using super = KDTree<BoundingSphere>;

	PointCloudBSH();
	PointCloudBSH(std::vector<Eigen::Vector3d> const& vertices);

	Eigen::Vector3d const& entityPosition(unsigned int i) const final;
	void computeHull(unsigned int b, unsigned int n, BoundingSphere& hull)
		const final;

private:

	std::vector<Eigen::Vector3d> const* m_vertices;
};

}
