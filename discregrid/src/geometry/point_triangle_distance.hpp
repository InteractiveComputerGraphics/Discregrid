
#pragma once

#include <array>
#include <Eigen/Core>

namespace Discregrid
{

enum class NearestEntity
{
	VN0, VN1, VN2, EN0, EN1, EN2, FN
};

double point_triangle_sqdistance(Eigen::Vector3d const& point, 
	std::array<Eigen::Vector3d const*, 3> const& triangle,
	Eigen::Vector3d* nearest_point = nullptr,
	NearestEntity* ne = nullptr);

}

