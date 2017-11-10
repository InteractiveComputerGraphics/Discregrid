#pragma once

#include <Eigen/Core>

namespace Discregrid
{

class BoundingSphere
{

public:

	BoundingSphere() : m_x(Eigen::Vector3d::Zero()), m_r(0.0) {}
	BoundingSphere(Eigen::Vector3d const& x, double r) : m_x(x), m_r(r) {}

	Eigen::Vector3d const& x() const { return m_x; }
	Eigen::Vector3d& x() { return m_x; }

	double r() const { return m_r; }
	double& r() { return m_r; }

	bool overlaps(BoundingSphere const& other) const 
	{ 
		double rr = m_r + other.m_r;
		return (m_x - other.m_x).squaredNorm() < rr * rr; 
	}

	bool contains(BoundingSphere const& other) const
	{
		double rr = r() - other.r();
		return (x() - other.x()).squaredNorm() < rr * rr;
	}

	bool contains(Eigen::Vector3d const& other) const
	{
		return (x() - other).squaredNorm() < m_r * m_r;
	}

private:

	Eigen::Vector3d m_x;
	double m_r;
};

}

