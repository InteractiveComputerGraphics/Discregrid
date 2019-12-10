
#include <geometry/mesh_distance.hpp>
#include <mesh/triangle_mesh.hpp>
#include "point_triangle_distance.hpp"

#include <limits>
#include <functional>
#include <omp.h>

using namespace Eigen;

namespace Discregrid
{

MeshDistance::MeshDistance(TriangleMesh const& mesh, bool precompute_normals)
	: m_bsh(mesh.vertex_data(), mesh.face_data()), m_mesh(mesh)
	, m_precomputed_normals(precompute_normals)
{
	auto max_threads = omp_get_max_threads();
	m_queues.resize(max_threads);
	m_nearest_face.resize(max_threads);
	m_cache.resize(max_threads, FunctionValueCache([&](Vector3d const& xi){ return signedDistance(xi);}, 10000u));
	m_ucache.resize(max_threads, FunctionValueCache([&](Vector3d const& xi){ return distance(xi);}, 10000u));

	m_bsh.construct();

	if (m_precomputed_normals)
	{
		m_face_normals.resize(m_mesh.nFaces());
		m_vertex_normals.resize(mesh.nVertices(), Vector3d::Zero());
		std::transform(m_mesh.faces().begin(), m_mesh.faces().end(),
			m_face_normals.begin(),
			[&](std::array<unsigned int, 3> const& face)
			{
				auto const& x0 = m_mesh.vertex(face[0]);
				auto const& x1 = m_mesh.vertex(face[1]);
				auto const& x2 = m_mesh.vertex(face[2]);

				auto n = (x1 - x0).cross(x2 - x0).normalized();

				auto e1 = (x1 - x0).normalized();
				auto e2 = (x2 - x1).normalized();
				auto e3 = (x0 - x2).normalized();

				auto alpha = Vector3d{
					std::acos(e1.dot(-e3)), 
					std::acos(e2.dot(-e1)),
					std::acos(e3.dot(-e2)) };
				m_vertex_normals[face[0]] += alpha[0] * n;
				m_vertex_normals[face[1]] += alpha[1] * n;
				m_vertex_normals[face[2]] += alpha[2] * n;

				return n;
			}
		);
	}
}

// Thread-safe.
double
MeshDistance::distance(Vector3d const& x, Vector3d* nearest_point, 
	unsigned int* nearest_face, NearestEntity* ne) const
{
	using namespace std::placeholders;

	auto dist_candidate = std::numeric_limits<double>::max();
	auto f = m_nearest_face[omp_get_thread_num()];
	if (f < m_mesh.nFaces())
	{
		auto t = std::array<Vector3d const*, 3>{
			&m_mesh.vertex(m_mesh.faceVertex(f, 0)),
			&m_mesh.vertex(m_mesh.faceVertex(f, 1)),
			&m_mesh.vertex(m_mesh.faceVertex(f, 2))
		};
		dist_candidate = std::sqrt(point_triangle_sqdistance(x, t));
	}

	auto pred = [&](unsigned int node_index, unsigned int)
	{
		return predicate(node_index, m_bsh, x, dist_candidate);
	};

	auto cb = [&](unsigned int node_index, unsigned int)
	{
		return callback(node_index, m_bsh, x, dist_candidate);
	};

	auto pless = [&](std::array<int, 2> const& c)
	{
		//return true;
		auto const& hull0 = m_bsh.hull(c[0]);
		auto const& hull1 = m_bsh.hull(c[1]);
		auto d0_2 = (x - hull0.x()).norm() - hull0.r();
		auto d1_2 = (x - hull1.x()).norm() - hull1.r();
		return d0_2 < d1_2;
	};

	while (!m_queues[omp_get_thread_num()].empty())
		m_queues[omp_get_thread_num()].pop();
	m_bsh.traverseDepthFirst(pred, cb, pless);

	f = m_nearest_face[omp_get_thread_num()];
	if (nearest_point)
	{
		auto t = std::array<Vector3d const*, 3>{
			&m_mesh.vertex(m_mesh.faceVertex(f, 0)),
			&m_mesh.vertex(m_mesh.faceVertex(f, 1)),
			&m_mesh.vertex(m_mesh.faceVertex(f, 2))
		};
		auto np = Vector3d{};
		auto ne_ = NearestEntity{};
		auto dist2_ = point_triangle_sqdistance(x, t, &np, &ne_);
		dist_candidate = std::sqrt(dist2_);
		if (ne)
			*ne = ne_;
		if (nearest_point)
			*nearest_point = np;
	}
	if (nearest_face)
		*nearest_face = f;
	return dist_candidate;
}

bool
MeshDistance::predicate(unsigned int node_index, 
	TriangleMeshBSH const& bsh,
	Vector3d const& x, 
	double& dist_candidate) const
{
	auto const& node = bsh.node(node_index);
	auto const& hull = bsh.hull(node_index);

	auto r = hull.r();

	auto temp = (x - hull.x()).eval();
	auto d_center2 = temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2];
	auto dmr = dist_candidate - r;
	if (dmr > 0.0)
	{
		auto temp = dist_candidate - r;
		if (temp * temp > d_center2)
			dist_candidate = std::sqrt(d_center2) + r;
	}
	else
	{
		auto d_center = std::sqrt(d_center2);
		if (dmr > d_center)
			dist_candidate = d_center + r;
	}
	auto temp_ = dist_candidate + r;
	return d_center2 <= temp_ * temp_;
}

void
MeshDistance::callback(unsigned int node_index, 
	TriangleMeshBSH const& bsh,
	Vector3d const& x, 
	double& dist_candidate) const
{
	auto const& node = m_bsh.node(node_index);
	auto const& hull = m_bsh.hull(node_index);

	if (!node.isLeaf())
		return;

	auto r = hull.r();

	auto temp = (x - hull.x()).eval();
	auto d_center2 = temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2];

	auto temp_ = dist_candidate + r;
	if (d_center2 > temp_ * temp_)
		return;

	auto dist_candidate_2 = dist_candidate * dist_candidate;
	auto changed = false;
	for (auto i = node.begin; i < node.begin + node.n; ++i)
	{
		auto f = m_bsh.entity(i);
		auto t = std::array<Vector3d const*, 3>{
			&m_mesh.vertex(m_mesh.faceVertex(f, 0)),
			&m_mesh.vertex(m_mesh.faceVertex(f, 1)),
			&m_mesh.vertex(m_mesh.faceVertex(f, 2))
		};
		auto dist2_ = point_triangle_sqdistance(x, t);
		if (dist_candidate_2 > dist2_)
		{
			dist_candidate_2 = dist2_;
			changed = true;
			m_nearest_face[omp_get_thread_num()] = f;
		}
	}
	if (changed)
	{
		dist_candidate = std::sqrt(dist_candidate_2);
	}
}

double
MeshDistance::signedDistance(Vector3d const& x) const
{
	unsigned int nf;
	auto ne = NearestEntity{};
	auto np = Vector3d{};
	auto dist = distance(x, &np, &nf, &ne);
	
	auto n = Vector3d{};
	switch (ne)
	{
	case NearestEntity::VN0:
		n = vertex_normal(m_mesh.faceVertex(nf, 0));
		break;
	case NearestEntity::VN1:
		n = vertex_normal(m_mesh.faceVertex(nf, 1));
		break;
	case NearestEntity::VN2:
		n = vertex_normal(m_mesh.faceVertex(nf, 2));
		break;
	case NearestEntity::EN0:
		n = edge_normal({nf, 0});
		break;
	case NearestEntity::EN1:
		n = edge_normal({nf, 1});
		break;
	case NearestEntity::EN2:
		n = edge_normal({nf, 2});
		break;
	case NearestEntity::FN:
		n = face_normal(nf);
		break;
	default:
		n.setZero();
		break;
	}

	if ((x - np).dot(n) < 0.0)
		dist *= -1.0;

	return dist;
}

double
MeshDistance::signedDistanceCached(Vector3d const & x) const
{
	return m_cache[omp_get_thread_num()](x);
}

double
MeshDistance::unsignedDistance(Vector3d const & x) const
{
	return distance(x);
}

double
MeshDistance::unsignedDistanceCached(Vector3d const & x) const
{
	return m_ucache[omp_get_thread_num()](x);
}

Vector3d
MeshDistance::face_normal(unsigned int f) const
{
	if (m_precomputed_normals)
		return m_face_normals[f];

	auto const& x0 = m_mesh.vertex(m_mesh.faceVertex(f, 0));
	auto const& x1 = m_mesh.vertex(m_mesh.faceVertex(f, 1));
	auto const& x2 = m_mesh.vertex(m_mesh.faceVertex(f, 2));

	return (x1 - x0).cross(x2 - x0).normalized();
}

Vector3d
MeshDistance::edge_normal(Halfedge const& h) const
{
	auto o = m_mesh.opposite(h);

	if (m_precomputed_normals)
	{
		if (o.isBoundary()) return m_face_normals[h.face()];
		return m_face_normals[h.face()] + m_face_normals[o.face()];
	}

	if (o.isBoundary()) return face_normal(h.face());
	return face_normal(h.face()) + face_normal(o.face());
}

Vector3d
MeshDistance::vertex_normal(unsigned int v) const
{
	if (m_precomputed_normals)
		return m_vertex_normals[v];

	auto const& x0 = m_mesh.vertex(v);
	auto n = Vector3d{}; n.setZero();
	for (auto h : m_mesh.incident_faces(v))
	{
		assert(m_mesh.source(h) == v);
		auto ve0 = m_mesh.target(h);
		auto e0 = (m_mesh.vertex(ve0) - x0).eval();
		e0.normalize();
		auto ve1 = m_mesh.target(h.next());
		auto e1 = (m_mesh.vertex(ve1) - x0).eval();
		e1.normalize();
		auto alpha = std::acos((e0.dot(e1)));
		n += alpha * e0.cross(e1);
	}
	return n;
}

}
