
#include <Discregrid/All>
#include <Eigen/Dense>
#include <cxxopts/cxxopts.hpp>

#include <string>
#include <iostream>
#include <array>

#include "sph_kernel.hpp"
#include "gauss_quadrature.hpp"

using namespace Eigen;

std::istream& operator>>(std::istream& is, std::array<unsigned int, 3>& data)
{
	is >> data[0] >> data[1] >> data[2];
	return is;
}

std::istream& operator>>(std::istream& is, AlignedBox3d& data)
{
	is	>> data.min()[0] >> data.min()[1] >> data.min()[2]
		>> data.max()[0] >> data.max()[1] >> data.max()[2];
	return is;
}


int main(int argc, char* argv[])
{
	cxxopts::Options options(argv[0], "Generates a signed distance field from a closed two-manifold triangle mesh.");
	options.positional_help("[input OBJ file]");

	options.add_options()
	("h,help", "Prints this help text")
	("r,rest_density", "Rest density rho0 of the fluid", cxxopts::value<double>()->default_value("1000.0"))
	("i,invert", "Invert field")
	("s,smoothing_length", "Kernel smoothing length", cxxopts::value<double>()->default_value("0.1"))
	("o,output", "Ouput file in cdf format", cxxopts::value<std::string>()->default_value(""))
	("no-reduction", "Disables discarding of cells for sparse layout.")
	("input", "Discrete grid file containing input SDF in field 0", cxxopts::value<std::vector<std::string>>())
	;

	try
	{
		options.parse_positional("input");
		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: GenerateSDF -r \"50 50 50\" dragon.obj" << std::endl;
			exit(0);
		}
		if (!result.count("input"))
		{
			std::cout << "ERROR: No input SDF given." << std::endl;
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: GenerateDensityMap -r \"50 50 50\" field.cdf" << std::endl;
			exit(1);
		}
		auto filename = result["input"].as<std::vector<std::string>>().front();

		if (!std::ifstream(filename).good())
		{
			std::cerr << "ERROR: Input file does not exist!" << std::endl;
			exit(1);
		}

		auto sdf = std::unique_ptr<Discregrid::DiscreteGrid>{};

		auto lastindex = filename.find_last_of(".");
		auto extension = filename.substr(lastindex + 1, filename.length() - lastindex);

		std::cout << "Load SDF...";
		if (extension == "cdf")
		{
			sdf = std::unique_ptr<Discregrid::CubicLagrangeDiscreteGrid>(
				new Discregrid::CubicLagrangeDiscreteGrid(filename));
		}
		std::cout << "DONE" << std::endl;

		auto h = result["s"].as<double>();
		auto sph_kernel = CubicKernel{};
		sph_kernel.setRadius(h);
		auto gamma = [&](Vector3d const& x)
		{
			auto ar = sph_kernel.getRadius();
			auto dist = sdf->interpolate(0u, x);
			if (dist > ar)
				return 0.0;
			return 1.0 - dist / ar;
		};
		auto int_domain = AlignedBox3d(Vector3d::Constant(-h), Vector3d::Constant(h));
		auto rho0 = result["r"].as<double>();
		auto density_func = [&](Vector3d const& x)
		{
			auto dist = sdf->interpolate(0u, x);
			if (dist > 2.0 * sph_kernel.getRadius())
			{
				return 0.0;
			}

			auto integrand = [&sph_kernel, &gamma, &x](Vector3d const& xi)
			{
				auto res = gamma(x + xi) * sph_kernel.W(xi);
				return res;
			};

			auto res = GaussQuadrature::integrate(integrand, int_domain, 30);
			return rho0 * res;
		};

		auto no_reduction = result["no-reduction"].count() > 0u;


		auto cell_diag = sdf->cellSize().norm();
		std::cout << "Generate density map..." << std::endl;
		sdf->addFunction(density_func, true, [&](Vector3d const& x_)
		{
			if (no_reduction)
			{
				return true;
			}
			auto x = x_.cwiseMax(sdf->domain().min()).cwiseMin(sdf->domain().max());
			auto dist = sdf->interpolate(0u, x);
			if (dist == std::numeric_limits<double>::max())
			{
				return false;
			}

			return -6.0 * h < dist + cell_diag && dist - cell_diag < 2.0 * h;
		});

		if (result["no-reduction"].count() == 0u)
		{
			std::cout << "Reduce discrete fields...";
			sdf->reduceField(0u, [&](const Vector3d &, double v)
			{
				return -6.0 * h < v + cell_diag && v - cell_diag < 2.0 * h;
			});
			sdf->reduceField(1u, [&](const Vector3d &, double v)
			{
				return 0.0 <= v && v <= 3.0 * rho0;
			});
			std::cout << "DONE" << std::endl;
		}

		std::cout << "Serialize discretization...";
		auto output_file = result["o"].as<std::string>();
		if (output_file == "")
		{
			output_file = filename;
			if (output_file.find(".") != std::string::npos)
			{
				auto lastindex = output_file.find_last_of(".");
				output_file = output_file.substr(0, lastindex);
			}
			output_file += ".cdm";
		}
		sdf->save(output_file);
		std::cout << "DONE" << std::endl;
	}
	catch (cxxopts::OptionException const& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
	
	return 0;
}
