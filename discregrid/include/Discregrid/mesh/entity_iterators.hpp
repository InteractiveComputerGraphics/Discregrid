#pragma once

#include "halfedge.hpp"

#include <iterator>
#include <array>
#include <Eigen/Core>

namespace Discregrid
{

class TriangleMesh;

class FaceContainer;
class FaceIterator
{
public:
	using iterator_category = std::random_access_iterator_tag;
    using value_type = std::array<unsigned int, 3>;
    using difference_type = unsigned int;
    using pointer = std::array<unsigned int, 3>*;
    using reference = std::array<unsigned int, 3>&;

	typedef FaceIterator _Mytype;

	FaceIterator() = delete;

	reference operator*();

	bool operator<(_Mytype const& other) const
	{
		return m_index < other.m_index;
	}
	bool operator==(_Mytype const& other) const
	{
		return m_index == other.m_index;
	}

	bool operator!=(_Mytype const& other) const
	{
		return !(*this == other);
	}

	inline _Mytype& operator++() { ++m_index; return *this; }
	inline _Mytype& operator--() { --m_index; return *this; }

	inline _Mytype operator+(_Mytype const& rhs)
	{
		return _Mytype(m_index + rhs.m_index, m_mesh);
	}
	inline difference_type operator-(_Mytype const& rhs)
	{
		return m_index - rhs.m_index;
	}
	inline _Mytype operator-(int const& rhs)
	{
		return _Mytype(m_index - rhs, m_mesh);
	}

	unsigned int vertex(unsigned int i) const;
	unsigned int& vertex(unsigned int i);

private:

	friend class FaceContainer;
	FaceIterator(unsigned int index, TriangleMesh* mesh)
		: m_index(index), m_mesh(mesh) {}

	unsigned int m_index;
	TriangleMesh* m_mesh;
};
class FaceConstIterator
{

public:
	using iterator_category = std::random_access_iterator_tag;
    using value_type = std::array<unsigned int, 3> const;
    using difference_type = unsigned int const;
    using pointer = std::array<unsigned int, 3> const*;
    using reference = std::array<unsigned int, 3> const&;

	typedef FaceConstIterator _Mytype;

	FaceConstIterator() = delete;

	reference operator*();

	bool operator<(_Mytype const& other) const
	{
		return m_index < other.m_index;
	}
	bool operator==(_Mytype const& other) const
	{
		return m_index == other.m_index;
	}

	bool operator!=(_Mytype const& other) const
	{
		return !(*this == other);
	}

	inline _Mytype& operator++() { ++m_index; return *this; }
	inline _Mytype& operator--() { --m_index; return *this; }

	inline _Mytype operator+(_Mytype const& rhs) const
	{
		return _Mytype(m_index + rhs.m_index, m_mesh);
	}
	inline difference_type operator-(_Mytype const& rhs) const
	{
		return m_index - rhs.m_index;
	}
	inline _Mytype operator-(int const& rhs) const
	{
		return _Mytype(m_index - rhs, m_mesh);
	}

	unsigned int vertex(unsigned int i) const;
	unsigned int& vertex(unsigned int i);

private:

	friend class FaceConstContainer;
	FaceConstIterator(unsigned int index, TriangleMesh const* mesh)
		: m_index(index), m_mesh(mesh) {}

	unsigned int m_index;
	TriangleMesh const* m_mesh;
};

class IncidentFaceContainer;
class IncidentFaceIterator
{

public:
	using iterator_category = std::forward_iterator_tag;
    using value_type = Halfedge;
    using difference_type = Halfedge;
    using pointer = Halfedge*;
    using reference = Halfedge&;

	typedef IncidentFaceIterator _Mytype;

	value_type operator*() { return m_h; }
	_Mytype& operator++();
	bool operator==(_Mytype const& other) const
	{
		return m_h == other.m_h;
	}

	bool operator!=(_Mytype const& other) const
	{
		return !(*this == other);
	}

private:

	friend class IncidentFaceContainer;
	IncidentFaceIterator(unsigned int v, TriangleMesh const* mesh);
	IncidentFaceIterator() : m_h(), m_begin(), m_mesh(nullptr) {}

	Halfedge m_h, m_begin;
	TriangleMesh const* m_mesh;
};


class VertexContainer;
class VertexIterator
{

public:
	using iterator_category = std::random_access_iterator_tag;
    using value_type = Eigen::Vector3d;
    using difference_type = double;
    using pointer = Eigen::Vector3d*;
    using reference = Eigen::Vector3d&;

	typedef VertexIterator _Mytype;

	VertexIterator() = delete;

	reference operator*();

	bool operator<(_Mytype const& other) const
	{
		return m_index < other.m_index;
	}
	bool operator==(_Mytype const& other) const
	{
		return m_index == other.m_index;
	}

	bool operator!=(_Mytype const& other) const
	{
		return !(*this == other);
	}

	inline _Mytype& operator++() { ++m_index; return *this; }
	inline _Mytype& operator--() { --m_index; return *this; }

	inline _Mytype operator+(_Mytype const& rhs) const
	{ 
		return _Mytype(m_index + rhs.m_index, m_mesh); 
	}
	inline difference_type operator-(_Mytype const& rhs) const
	{ 
		return m_index - rhs.m_index; 
	}
	inline _Mytype operator-(int const& rhs) const
	{ 
		return _Mytype(m_index - rhs, m_mesh); 
	}

	unsigned int index() const;

private:

	friend class VertexContainer;
	VertexIterator(unsigned int index, TriangleMesh* mesh)
		: m_index(index), m_mesh(mesh) {}

	unsigned int m_index;
	TriangleMesh* m_mesh;
};


class VertexConstContainer;
class VertexConstIterator
{

public:
	using iterator_category = std::random_access_iterator_tag;
    using value_type = Eigen::Vector3d const;
    using difference_type = double const;
    using pointer = Eigen::Vector3d* const;
    using reference = Eigen::Vector3d const&;

	typedef VertexConstIterator _Mytype;

	VertexConstIterator() = delete;

	reference operator*();

	bool operator<(_Mytype const& other) const
	{
		return m_index < other.m_index;
	}
	bool operator==(_Mytype const& other) const
	{
		return m_index == other.m_index;
	}

	bool operator!=(_Mytype const& other) const
	{
		return !(*this == other);
	}

	inline _Mytype& operator++() { ++m_index; return *this; }
	inline _Mytype& operator--() { --m_index; return *this; }

	inline _Mytype operator+(_Mytype const& rhs) const
	{ 
		return _Mytype(m_index + rhs.m_index, m_mesh); 
	}
	inline difference_type operator-(_Mytype const& rhs) const
	{ 
		return m_index - rhs.m_index; 
	}
	inline _Mytype operator-(int const& rhs) const
	{ 
		return _Mytype(m_index - rhs, m_mesh); 
	}

	unsigned int index() const;

private:

	friend class VertexConstContainer;
	VertexConstIterator(unsigned int index, TriangleMesh const* mesh)
		: m_index(index), m_mesh(mesh) {}

	unsigned int m_index;
	TriangleMesh const* m_mesh;
};
}
