#pragma once
#include <array>
#include <limits>
#include <vector>
#include <Eigen/src/Core/Matrix.h>
#include "TriangleMeshDistance.h"
#ifdef IS_CLIPPER_ENABLED
#include "clipper2/clipper.core.h"
#endif

// This file is a derivative of TriangleMeshDistance.h providing the same functionality in 2D.
// This file can be used to get the distance of a point to a closed polygon (that can contain holes).
namespace Discregrid
{
    // Point-Edge distance definitions
    enum class NearestEntity2D { V0, V1, E };
    static double point_edge_sq_unsigned(
        NearestEntity2D& nearest_entity,
        Eigen::Vector2d& nearest_point,
        const Eigen::Vector2d& point,
        const Eigen::Vector2d& v0,
        const Eigen::Vector2d& v1);
    // -----------------------------------

    // Struct that contains the result of a 2D distance query
    struct Result2D
    {
        double distance = std::numeric_limits<double>::max();
        Eigen::Vector2d nearest_point;
        NearestEntity2D nearest_entity;
        int edge_id = -1;
    };
    // -----------------------------------

    // A class to compute signed and unsigned distances to a triangle mesh.
    class PolygonDistance
    {
    private:
        /* Definitions */
        struct BoundingSphere
        {
            Eigen::Vector2d center;
            double radius;
        };

        struct Node
        {
            BoundingSphere bv_left;
            BoundingSphere bv_right;
            int left = -1; // If left == -1, right is the edge_id
            int right = -1;
        };

        struct Edge
        {
            std::array<Eigen::Vector2d, 2> vertices;
            int id = -1;
        };

        /* Fields */
        std::vector<Eigen::Vector2d> vertices;
        std::vector<Eigen::Vector2i> edges;
        std::vector<Node> nodes;
        std::vector<Eigen::Vector2d> pseudonormals_edges;
        std::vector<Eigen::Vector2d> pseudonormals_vertices;
        BoundingSphere root_bv;
        bool is_constructed = false;

        /* Methods */
        void _construct();
        void _build_tree(const int node_id, BoundingSphere& bounding_sphere, std::vector<Edge> &edges, const int begin, const int end);
        void _query(Result2D &result, const Node &node, const Eigen::Vector2d& point) const;
        
    public:
        /* Methods */
        PolygonDistance() = default;

        /**
         * @brief Constructs a new PolygonDistance object.
         *
         * @param vertices Pointer to the vertices coordinates array in xyxy... layout.
         * @param n_vertices Number of vertices.
         * @param edges Pointer to the connectivity array in ijij... layout.
         * @param n_edges Number of edges.
        */
        template<typename FLOAT, typename INT, typename SIZE_T>
        PolygonDistance(const FLOAT* vertices, const SIZE_T n_vertices, const INT* edges, const SIZE_T n_edges);

        /**
         * @brief Constructs a new PolygonDistance object.
         *
         * @param vertices Vertices of the polygon. Y coordinate of the 3rd vertex should be accessed by `vertices[2][1]`.
         * @param edges Edges of the polygon. Index of the 2nd vertex of the 3rd edge should be accessed by `edges[2][1]`.
        */
        template<typename IndexableVector2double, typename IndexableVector2int>
        PolygonDistance(const std::vector<IndexableVector2double>& vertices, const std::vector<IndexableVector2int>& edges);

#ifdef IS_CLIPPER_ENABLED
        /**
         * @brief Constructs a new PolygonDistance object.
         *
         * @param polygon Clipper2 PathsD (= Polygon) object.
        */
        PolygonDistance(const Clipper2Lib::PathsD &polygon);
#endif

        /**
         * @brief Initializes an existing PolygonDistance object (including empty ones).
         *
         * @param vertices Pointer to the vertices coordinates array in xyxy... layout.
         * @param n_vertices Number of vertices.
         * @param edges Pointer to the conectivity array in ijij... layout.
         * @param n_edges Number of edges.
        */
        template<typename FLOAT, typename INT, typename SIZE_T>
        void construct(const FLOAT* vertices, const SIZE_T n_vertices, const INT* edges, const SIZE_T n_edges);

        /**
         * @brief Initializes an existing PolygonDistance object (including empty ones).
         *
         * @param vertices Vertices of the polygon. Y coordinate of the 3rd vertex should be access by `vertices[2][1]`.
         * @param edges Edges of the polygon. Index of the 2nd vertex of the 3rd edge should be access by `edges[2][1]`.
        */
        template<typename IndexableVector3double, typename IndexableVector3int>
        void construct(const std::vector<IndexableVector3double>& vertices, const std::vector<IndexableVector3int>& edges);

        /**
         * @brief Computes the unsigned distance from a point to the polygon. Thread safe.
         *
         * @param point to query from. Typed to `Vector2d` but can be passed as `{x, y}`.
         * 
         * @return Result containing distance, nearest point on the polygon, nearest entity and the nearest edge index.
        */
        template<typename IndexableVector2double>
        Result2D unsigned_distance(const IndexableVector2double& point) const;
        Result2D unsigned_distance(const std::array<double, 2>& point) const;

        /**
         * @brief Computes the unsigned distance from a point to the polygon. Thread safe.
         *
         * @param point to query from. Typed to `Vector2d` but can be passed as `{x, y}`.
         * 
         * @return Result containing distance, nearest point on the polygon, nearest entity and the nearest edge index.
        */
        template<typename IndexableVector2double>
        Result2D signed_distance(const IndexableVector2double& point) const;
        Result2D signed_distance(const std::array<double, 2>& point) const;
    };
}

/* ==========================================  DECLARATIONS  ========================================== */
template<typename FLOAT, typename INT, typename SIZE_T>
inline Discregrid::PolygonDistance::PolygonDistance(const FLOAT* vertices, const SIZE_T n_vertices, const INT* edges, const SIZE_T n_edges)
{
    this->construct(vertices, n_vertices, edges, n_edges);
}

template<typename IndexableVector2double, typename IndexableVector2int>
inline Discregrid::PolygonDistance::PolygonDistance(const std::vector<IndexableVector2double>& vertices, const std::vector<IndexableVector2int>& edges)
{
    this->construct(vertices, edges);
}

#ifdef IS_CLIPPER_ENABLED
inline Discregrid::PolygonDistance::PolygonDistance(const Clipper2Lib::PathsD &polygon)
{
    size_t size = 0;
    for (size_t i = 0; i < polygon.size(); i++)
    {
        size += polygon[i].size();
    }
    std::vector<Eigen::Vector2d> polygonVertices;
    polygonVertices.reserve(size);
    std::vector<Eigen::Vector2i> polygonEdges;
    polygonEdges.reserve(size);
    size_t counter = 0;
    for (size_t i = 0; i < polygon.size(); i++)
    {
        auto path = polygon[i];
        size_t subCounter = counter;
        for (size_t j = 0; j < path.size() - 1; j++)
        {
            auto point = path[j];
            polygonVertices.emplace_back(point.x, point.y);
            polygonEdges.emplace_back(counter++, counter);
        }
        auto point = path.back();
        polygonVertices.emplace_back(point.x, point.y);
        polygonEdges.emplace_back(counter++, subCounter);
    }
    this->construct(polygonVertices, polygonEdges);
}
#endif

template<typename FLOAT, typename INT, typename SIZE_T>
inline void Discregrid::PolygonDistance::construct(const FLOAT* vertices, const SIZE_T n_vertices, const INT* edges, const SIZE_T n_edges)
{
    this->vertices.resize(2 * n_vertices);
    for (size_t i = 0; i < static_cast<size_t>(n_vertices); i++)
    {
        this->vertices[i][0] = static_cast<double>(vertices[2 * i]);
        this->vertices[i][1] = static_cast<double>(vertices[2 * i + 1]);
    }

    this->edges.resize(2 * n_edges);
    for (size_t i = 0; i < static_cast<size_t>(n_edges); i++)
    {
        this->edges[i][0] = static_cast<int>(edges[2 * i]);
        this->edges[i][1] = static_cast<int>(edges[2 * i + 1]);
    }
    this->_construct();
}

template<typename IndexableVector3double, typename IndexableVector3int>
inline void Discregrid::PolygonDistance::construct(const std::vector<IndexableVector3double>& vertices, const std::vector<IndexableVector3int>& edges)
{
    this->vertices.resize(vertices.size());
    for (size_t i = 0; i < vertices.size(); i++)
    {
        this->vertices[i][0] = static_cast<double>(vertices[i][0]);
        this->vertices[i][1] = static_cast<double>(vertices[i][1]);
    }
    
    this->edges.resize(edges.size());
    for (size_t i = 0; i < edges.size(); i++)
    {
        this->edges[i][0] = static_cast<int>(edges[i][0]);
        this->edges[i][1] = static_cast<int>(edges[i][1]);
    }
    this->_construct();
}

template<typename IndexableVector2double>
inline Discregrid::Result2D Discregrid::PolygonDistance::unsigned_distance(const IndexableVector2double& point) const
{
    return this->unsigned_distance({static_cast<double>(point[0]), static_cast<double>(point[1])});
}

inline Discregrid::Result2D Discregrid::PolygonDistance::unsigned_distance(const std::array<double, 2>& point) const
{
    if (!this->is_constructed)
    {
        std::cout << "PolygonDistance error: not constructed." << std::endl;
        exit(-1);
    }

    const Eigen::Vector2d p(point[0], point[1]);
    Result2D result;
    result.distance = std::numeric_limits<double>::max();
    this->_query(result, this->nodes[0], p);
    return result;
}

template<typename IndexableVector2double>
inline Discregrid::Result2D Discregrid::PolygonDistance::signed_distance(const IndexableVector2double& point) const
{
    const Eigen::Vector2d p(point[0], point[1]);
    Result2D result = this->unsigned_distance(p);

    const Eigen::Vector2i& edge = this->edges[result.edge_id];
    Eigen::Vector2d pseudonormal;
    switch (result.nearest_entity)
    {
    case NearestEntity2D::V0:
        pseudonormal = this->pseudonormals_vertices[edge[0]];
        break;
    case NearestEntity2D::V1:
        pseudonormal = this->pseudonormals_vertices[edge[1]];
        break;
    case NearestEntity2D::E:
        pseudonormal = this->pseudonormals_edges[result.edge_id];
        break;
    default:
        break;
    }

    const Eigen::Vector2d u = p - result.nearest_point;
    result.distance *= u.dot(pseudonormal) >= 0.0 ? 1.0 : -1.0;

    return result;
}

inline Discregrid::Result2D Discregrid::PolygonDistance::signed_distance(const std::array<double, 2>& point) const
{
    return this->signed_distance({static_cast<double>(point[0]), static_cast<double>(point[1])});
}

inline void Discregrid::PolygonDistance::_construct()
{
	if (this->edges.empty())
	{
		std::cout << "PolygonDistance error: Empty edge list." << std::endl;
		exit(-1);
	}
	
	// Build the tree containing the edges
	std::vector<Edge> edges;

	edges.resize(this->edges.size());
	for (int i = 0; i < static_cast<int>(this->edges.size()); i++)
	{
		edges[i].id = i;

		const Eigen::Vector2i& edge = this->edges[i];
		edges[i].vertices[0] = this->vertices[edge[0]];
		edges[i].vertices[1] = this->vertices[edge[1]];
	}

	this->nodes.emplace_back();
	this->_build_tree(0, this->root_bv, edges, 0, static_cast<int>(edges.size()));

	// Compute
	this->pseudonormals_edges.reserve(this->edges.size());
	this->pseudonormals_vertices.resize(this->vertices.size(), { 0, 0 });

	std::unordered_map<int, std::vector<int>> vertexIndexToEdgeNormalIndices;

	// (Pseudo-)Normals edges
	for (int i = 0; i < static_cast<int>(this->edges.size()); i++) {

		// Edge
		const Eigen::Vector2i edge = this->edges[i];
		const int i0 = edge[0];
		const int i1 = edge[1];
		const Eigen::Vector2d a = this->vertices[i0];
		const Eigen::Vector2d b = this->vertices[i1];

		// Edge normal calculation routine from Clipper2 library
		if (a == b)
		{
			std::cout << "PolygonDistance::_construct error: Can not calculate edge normal because both vertices are the same!" << std::endl;
			exit(-1);
		}
		auto dx = b.x() - a.x();
		auto dy = b.y() - a.y();
		const auto inverse_hypot = 1.0 / std::sqrt(dx * dx + dy * dy);
		dx *= inverse_hypot;
		dy *= inverse_hypot;
		const Eigen::Vector2d normal(dy, -dx);

		for (const int index : edge)
		{
			auto edgeNormalIndex = static_cast<int>(this->pseudonormals_edges.size());

			if (vertexIndexToEdgeNormalIndices.find(index) == vertexIndexToEdgeNormalIndices.end())
			{
				vertexIndexToEdgeNormalIndices[index] = std::vector<int>();
			}

			vertexIndexToEdgeNormalIndices[index].push_back(edgeNormalIndex);
		}

		this->pseudonormals_edges.push_back(normal);
	}

	// Pseudonormals vertices
	for (const auto& pair : vertexIndexToEdgeNormalIndices)
	{
		int vertexIndex = pair.first;
		std::vector<int> edgeNormalIndices = pair.second;
		if (edgeNormalIndices.size() != 2)
		{
			std::cout << "PolygonDistance error: Polygon is not valid. At least one vertex does not have exactly two incident edges." << std::endl;
			exit(-1);
		}
		Eigen::Vector2i edge1 = this->edges[edgeNormalIndices[0]];
		Eigen::Vector2i edge2 = this->edges[edgeNormalIndices[1]];
		Eigen::Vector2d a = this->vertices[vertexIndex];
		Eigen::Vector2d b = this->vertices[edge1[0] == edge2[0] || edge1[0] == edge2[1] ? edge1[1] : edge1[0]];
		Eigen::Vector2d c = this->vertices[vertexIndex == edge2[0] ? edge2[1] : edge2[0]];
		double alpha0 = std::acos(std::abs((b - a).normalized().dot((c - a).normalized())));
		this->pseudonormals_vertices[vertexIndex] = (alpha0 * (this->pseudonormals_edges[edgeNormalIndices[0]] + this->pseudonormals_edges[edgeNormalIndices[1]])).normalized();
	}
	
	this->is_constructed = true;
}

inline void Discregrid::PolygonDistance::_build_tree(const int node_id, BoundingSphere& bounding_sphere, std::vector<Edge> &edges, const int begin, const int end)
{
	const int n_edges = end - begin;

	if (n_edges <= 0)
	{
		std::cout << "PolygonDistance::_construct error: Empty leaf." << std::endl;
		exit(-1);
	}
	
	if (n_edges == 1)
	{
		// Build node leaf
		this->nodes[node_id].left = -1;
		this->nodes[node_id].right = edges[begin].id;

		// Bounding sphere
		const Edge& edge = edges[begin];
		const Eigen::Vector2d helper = edge.vertices[0] + edge.vertices[1];
		const Eigen::Vector2d center(helper.x() / 2.0, helper.y() / 2.0);
		const double radius = std::max((edge.vertices[0] - center).norm(), (edge.vertices[1] - center).norm());
		bounding_sphere.center = center;
		bounding_sphere.radius = radius;
	}
	else
	{
		// Compute AABB center and largest dimension of all current edges. (AABB center is the average of all edge vertex positions)
		Eigen::Vector2d top = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
		Eigen::Vector2d bottom = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
		Eigen::Vector2d center = {0, 0};
		for (int edge_i = begin; edge_i < end; edge_i++)
		{
			for (int vertex_i = 0; vertex_i < 2; vertex_i++)
			{
				const Eigen::Vector2d& p = edges[edge_i].vertices[vertex_i];
				center += p;

				for (int coord_i = 0; coord_i < 2; coord_i++)
				{
					top[coord_i] = std::max(top[coord_i], p[coord_i]);
					bottom[coord_i] = std::min(bottom[coord_i], p[coord_i]);
				}
			}
		}
		center /= 2 * n_edges;
		const Eigen::Vector2d diagonal = top - bottom;
		const int split_dim = diagonal[0] < diagonal[1] ? 1 : 0;

		// Set node bounding sphere
		double radius_sq = 0.0;
		for (int edge_i = begin; edge_i < end; edge_i++)
		{
			for (int i = 0; i < 2; i++)
			{
				radius_sq = std::max(radius_sq, (center - edges[edge_i].vertices[i]).squaredNorm());
			}
		}
		bounding_sphere.center = center;
		bounding_sphere.radius = std::sqrt(radius_sq);

		// Sort the triangles according to their center along the split dimension
		std::sort(edges.begin() + begin, edges.begin() + end,
			[split_dim](const Edge& a, const Edge& b)
			{
				return a.vertices[0][split_dim] < b.vertices[0][split_dim];
			}
		);

		// Children
		const int mid = static_cast<int>(0.5 * (begin + end));

		this->nodes[node_id].left = static_cast<int>(this->nodes.size());
		this->nodes.emplace_back();
		this->_build_tree(this->nodes[node_id].left, this->nodes[node_id].bv_left, edges, begin, mid);

		this->nodes[node_id].right = static_cast<int>(this->nodes.size());
		this->nodes.emplace_back();
		this->_build_tree(this->nodes[node_id].right, this->nodes[node_id].bv_right, edges, mid, end);
	}
}


inline void Discregrid::PolygonDistance::_query(Result2D &result, const Node &node, const Eigen::Vector2d& point) const
{
	// End of recursion
	if (node.left == -1)
	{
		const int edge_id = node.right;
		const Eigen::Vector2i& edge = this->edges[node.right]; // If left == -1, right is the edge
		const Eigen::Vector2d& v0 = this->vertices[edge[0]];
		const Eigen::Vector2d& v1 = this->vertices[edge[1]];

		Eigen::Vector2d nearest_point;
		NearestEntity2D nearest_entity;
		const double distance_sq = point_edge_sq_unsigned(nearest_entity, nearest_point, point, v0, v1);

		if (distance_sq < result.distance * result.distance)
		{
			result.nearest_point = nearest_point;
			result.nearest_entity = nearest_entity;
			result.distance = std::sqrt(distance_sq);
			result.edge_id = edge_id;
		}
	}
	// Recursion
	else
	{
		// Find which child bounding volume is closer
		const double d_left = (point - node.bv_left.center).norm() - node.bv_left.radius;
		const double d_right = (point - node.bv_right.center).norm() - node.bv_right.radius;
		
		if (d_left < d_right)
		{
			// Overlap test
			if (d_left < result.distance)
			{
				this->_query(result, this->nodes[node.left], point);
			}

			if (d_right < result.distance)
			{
				this->_query(result, this->nodes[node.right], point);
			}
		}
		else
		{
			if (d_right < result.distance)
			{
				this->_query(result, this->nodes[node.right], point);
			}
			if (d_left < result.distance)
			{
				this->_query(result, this->nodes[node.left], point);
			}
		}
	}
}

static double Discregrid::point_edge_sq_unsigned(NearestEntity2D& nearest_entity, Eigen::Vector2d& nearest_point, const Eigen::Vector2d& point, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1)
{
	const Eigen::Vector2d e = v1 - v0;
	const double e2 = e.dot(e);
	const Eigen::Vector2d diff = point - v0;
	const double t = std::clamp(diff.dot(e) / e2, 0.0, 1.0);
	
	// Closest point on edge
	nearest_point = {v0.x() + t * e.x(), v0.y() + t * e.y()};

	if (t <= std::numeric_limits<double>::epsilon())
	{
		nearest_entity = NearestEntity2D::V0;
	}
	else if (t >= 1.0 - std::numeric_limits<double>::epsilon())
	{
		nearest_entity = NearestEntity2D::V1;
	}
	else
	{
		nearest_entity = NearestEntity2D::E;
	}

	return (point-nearest_point).squaredNorm();
}