#pragma once

#include <Discregrid/mesh/triangle_mesh.hpp>

namespace std {
	template <> struct hash<Eigen::Vector3d>
	{
		std::size_t operator()(Eigen::Vector3d const& x) const
		{
			std::size_t seed = 0;
			std::hash<double> hasher;
			seed ^= hasher(x[0]) + 0x9e3779b9 + (seed<<6) + (seed>>2);
			seed ^= hasher(x[1]) + 0x9e3779b9 + (seed<<6) + (seed>>2);
			seed ^= hasher(x[2]) + 0x9e3779b9 + (seed<<6) + (seed>>2);
			return seed;
		}
	};

	template <> struct less<Eigen::Vector3d>
	{
		bool operator()(Eigen::Vector3d const& left, Eigen::Vector3d const& right) const
		{
			for (auto i = 0u; i < 3u; ++i)
			{
				if (left(i) < right(i))
					return true;
				else if (left(i) > right(i))
					return false;
			}
			return false;
		}
	};
}
#include <Discregrid/utility/lru_cache.hpp>

#include <Discregrid/acceleration/bounding_sphere_hierarchy.hpp>


#include <array>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>

#include <Eigen/Dense>

namespace Discregrid
{

enum class NearestEntity;
class TriangleMesh;
class Halfedge;
class MeshDistance
{

	struct Candidate
	{
		bool operator<(Candidate const& other) const { return b < other.b; }
		unsigned int node_index;
		double b, w;
	};

public:

	MeshDistance(TriangleMesh const& mesh, bool precompute_normals = true);

	// Returns the shortest unsigned distance from a given point x to
	// the stored mesh.
	// Thread-safe function.
	double distance(Eigen::Vector3d const& x, Eigen::Vector3d* nearest_point = nullptr,
		unsigned int* nearest_face = nullptr, NearestEntity* ne = nullptr) const;

	// Requires a closed two-manifold mesh as input data.
	// Thread-safe function.
	double signedDistance(Eigen::Vector3d const& x) const;
	double signedDistanceCached(Eigen::Vector3d const& x) const;

	double unsignedDistance(Eigen::Vector3d const& x) const;
	double unsignedDistanceCached(Eigen::Vector3d const& x) const;

private:

	Eigen::Vector3d vertex_normal(unsigned int v) const;
	Eigen::Vector3d edge_normal(Halfedge const& h) const;
	Eigen::Vector3d face_normal(unsigned int f) const;

	void callback(unsigned int node_index, TriangleMeshBSH const& bsh,
		Eigen::Vector3d const& x, 
		double& dist) const;

	bool predicate(unsigned int node_index, TriangleMeshBSH const& bsh, 
		Eigen::Vector3d const& x, double& dist) const;

private:

	TriangleMesh const& m_mesh;
	TriangleMeshBSH m_bsh;

	using FunctionValueCache = LRUCache<Eigen::Vector3d, double>;
	mutable std::vector<TriangleMeshBSH::TraversalQueue> m_queues;
	mutable std::vector<unsigned int> m_nearest_face;
	mutable std::vector<FunctionValueCache> m_cache;
	mutable std::vector<FunctionValueCache> m_ucache;
	
	std::vector<Eigen::Vector3d> m_face_normals;
	std::vector<Eigen::Vector3d> m_vertex_normals;
	bool m_precomputed_normals;
};

}

