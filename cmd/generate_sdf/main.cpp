
#include <Discregrid/All>
#include <Eigen/Dense>

#include "resource_path.hpp"

#include <string>
#include <iostream>
#include <array>

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

#include <cxxopts/cxxopts.hpp>

int main(int argc, char* argv[])
{
	cxxopts::Options options(argv[0], "Generates a signed distance field from a closed two-manifold triangle mesh.");
	options.positional_help("[input OBJ file]");

	options.add_options()
	("h,help", "Prints this help text")
	("r,resolution", "Grid resolution", cxxopts::value<std::array<unsigned int, 3>>()->default_value("10 10 10"))
	("d,domain", "Domain extents (bounding box), format: \"minX minY minZ maxX maxY maxZ\"", cxxopts::value<AlignedBox3d>())
	("i,invert", "Invert SDF")
	("o,output", "Ouput file in cdf format", cxxopts::value<std::string>()->default_value(""))
	("input", "OBJ file containing input triangle mesh", cxxopts::value<std::vector<std::string>>())
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
			std::cout << "ERROR: No input mesh given." << std::endl;
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: GenerateSDF -r \"50 50 50\" dragon.obj" << std::endl;
			exit(1);
		}
		auto resolution = result["r"].as<std::array<unsigned int, 3>>();
		auto filename = result["input"].as<std::vector<std::string>>().front();

		if (!std::ifstream(filename).good())
		{
			std::cerr << "ERROR: Input file does not exist!" << std::endl;
			exit(1);
		}

		std::cout << "Load mesh...";
		Discregrid::TriangleMesh mesh(filename);
		std::cout << "DONE" << std::endl;

		std::cout << "Set up data structures...";
		Discregrid::MeshDistance md(mesh);
		std::cout << "DONE" << std::endl;

		Eigen::AlignedBox3d domain;
		domain.setEmpty();
		if (result.count("d"))
		{
			domain = result["d"].as<Eigen::AlignedBox3d>();
		}
		if (domain.isEmpty())
		{
			for (auto const& x : mesh.vertices())
			{
				domain.extend(x);
			}
			domain.max() += 1.0e-3 * domain.diagonal().norm() * Vector3d::Ones();
			domain.min() -= 1.0e-3 * domain.diagonal().norm() * Vector3d::Ones();
		}

		Discregrid::CubicLagrangeDiscreteGrid sdf(domain, resolution);
		auto func = Discregrid::DiscreteGrid::ContinuousFunction{};
		if (result.count("invert"))
		{
			func = [&md](Vector3d const& xi) {return -1.0 * md.signedDistanceCached(xi); };
		}
		else
		{
			func = [&md](Vector3d const& xi) {return md.signedDistanceCached(xi); };
		}

		std::cout << "Generate discretization..." << std::endl;
		sdf.addFunction(func, true);
		std::cout << "DONE" << std::endl;

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
			output_file += ".cdf";
		}
		sdf.save(output_file);
		std::cout << "DONE" << std::endl;
	}
	catch (cxxopts::OptionException const& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
	
	return 0;
}