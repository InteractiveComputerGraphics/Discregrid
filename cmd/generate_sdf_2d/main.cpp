#include <Discregrid/All>
#include <Eigen/Dense>

#include <string>
#include <iostream>
#include <array>
#include <filesystem>

#include <clipper2/clipper.h>
#include <clipper.svg.h>
#include <clipper2/clipper.export.h>

using namespace Eigen;

std::istream& operator>>(std::istream& is, std::array<unsigned int, 2>& data)  
{  
	is >> data[0] >> data[1];  
	return is;  
}  

std::istream& operator>>(std::istream& is, AlignedBox2d& data)  
{  
	is	>> data.min()[0] >> data.min()[1]
		>> data.max()[0] >> data.max()[1];  
	return is;  
}  

#include <cxxopts/cxxopts.hpp>

int main(int argc, char* argv[])
{
	cxxopts::Options options(argv[0], "Generates a 2D signed distance field from a triangle mesh.");
	options.positional_help("[input OBJ file]");

	options.add_options()
	("h,help", "Prints this help text")
	("r,resolution", "Grid resolution", cxxopts::value<std::array<unsigned int, 2>>()->default_value("10 10"))
	("d,domain", "Domain extents (bounding box), format: \"minX minY maxX maxY\"", cxxopts::value<AlignedBox2d>())
	("i,invert", "Invert SDF")
	("o,output", "Output file in cdf format", cxxopts::value<std::string>()->default_value(""))
	("input", "OBJ file containing input triangle mesh", cxxopts::value<std::vector<std::string>>())
	;

	try
	{
		options.parse_positional("input");
		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: GenerateSDF2D -r \"50 50\" dragon.obj" << std::endl;
			exit(0);
		}
		if (!result.count("input"))
		{
			std::cout << "ERROR: No input mesh given." << std::endl;
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: GenerateSDF2D -r \"50 50\" dragon.obj" << std::endl;
			exit(1);
		}
		auto resolution = result["r"].as<std::array<unsigned int, 2>>();
		auto filename = result["input"].as<std::vector<std::string>>().front();

		if (!std::ifstream(filename).good())
		{
			std::cerr << "ERROR: Input file does not exist!" << std::endl;
			exit(1);
		}

		std::cout << "Load mesh...";
		Discregrid::TriangleMesh mesh(filename);
		std::cout << "DONE" << std::endl;

		std::cout << "Weld vertices and perform backface culling..." << std::endl;
		std::vector<Vector2d> culled_vertices;
		std::vector<unsigned int> culled_triangles;
		mesh.weldVerticesAndCullBackfaces2D(culled_vertices, culled_triangles);
		std::cout << "N Vertices: " << culled_vertices.size() << "; N Triangles: " << culled_triangles.size() / 3 << std::endl;
		std::cout << "DONE" << std::endl;

		std::cout << "Create polygon from the union of all culled triangles..." << std::endl;
		Clipper2Lib::PathsD polygon(culled_triangles.size() / 3);
		for (int i = 0; i < culled_triangles.size(); i += 3)
		{
			const auto v0 = culled_vertices[culled_triangles[i]];
			const auto v1 = culled_vertices[culled_triangles[i + 1]];
			const auto v2 = culled_vertices[culled_triangles[i + 2]];
			polygon.emplace_back(Clipper2Lib::MakePathD({v0.x(), v0.y(), v1.x(), v1.y(), v2.x(), v2.y()}));
		}
		polygon = Clipper2Lib::Union(polygon, Clipper2Lib::FillRule::NonZero, 8);
		std::cout << "DONE" << std::endl;

		std::cout << "Set up polygon distance function..." << std::endl;
		Discregrid::PolygonDistance pd(polygon);
		std::cout << "DONE" << std::endl;

		Eigen::Vector2d sample(0.0, 0.0);
		std::cout << "TEST TEST PLS REMOVE ME: Sampling at (" << sample.x() << "|" << sample.y() << ")..." << std::endl;
		const Discregrid::Result2D signed_distance = pd.signed_distance(sample);
		std::cout << "Closest point on polygon from sample: (" << signed_distance.nearest_point.x() << "|" << signed_distance.nearest_point.y() << ")" << std::endl;
		
		std::filesystem::path filePath(filename);
		std::string svg_name = filePath.filename().stem().string().append(".svg");
		std::cout << "Writing polygon as SVG to " << svg_name << "..." << std::endl;
		Clipper2Lib::SvgWriter svg;
		svg.AddPaths(polygon, false, Clipper2Lib::FillRule::NonZero, 0x4080ff9C, 0xFF003300, 1.5, false);
		Clipper2Lib::PathD signed_distance_path = Clipper2Lib::MakePathD({signed_distance.nearest_point.x(), signed_distance.nearest_point.y(), sample.x(), sample.y()});
		svg.AddPath(signed_distance_path, true, Clipper2Lib::FillRule::NonZero, 0xFFFF00FF, 0xFFFF00FF, 3, false);
		std::cout << "TEST TEST PLS REMOVE ME: Sampling from -0.5 to +0.5 in both dimensoins and save to svg..." << std::endl;
		for (double xSample = -0.5; xSample <= 0.5; xSample += 0.1)
		{
			for (double ySample = -0.5; ySample <= 0.5; ySample += 0.1)
			{
				const Discregrid::Result2D signed_distance_sample = pd.signed_distance(Vector2d{xSample, ySample});
				svg.AddPath(Clipper2Lib::MakePathD({signed_distance_sample.nearest_point.x(), signed_distance_sample.nearest_point.y(), xSample, ySample}), true, Clipper2Lib::FillRule::NonZero, 0xFFFF00FF, 0xFFFF00FF, 3, false);
			}
		}
		svg.SaveToFile(svg_name, 1000, 1000, 10);
		std::cout << "DONE" << std::endl;
		
		Eigen::AlignedBox2d domain;
		domain.setEmpty();
		if (result.count("d"))
		{
			domain = result["d"].as<Eigen::AlignedBox2d>();
		}
		if (domain.isEmpty())
		{
			for (size_t i = 0; i < polygon.size(); i++)
			{
				
				const auto &outline = polygon[i];
				for (auto const& x : outline)
				{
					domain.extend(Matrix<double, 2, 1>{x.x, x.y});
				}
			}
			
			domain.max() += 1.0e-3 * domain.diagonal().norm() * Vector2d::Ones();
			domain.min() -= 1.0e-3 * domain.diagonal().norm() * Vector2d::Ones();
		}

		Discregrid::CubicLagrangeDiscreteGrid2D sdf(domain, resolution);
		auto func = Discregrid::DiscreteGrid2D::ContinuousFunction{};
		if (result.count("invert"))
		{
			func = [&pd](Vector2d const& xi) { return -1.0 * pd.signed_distance(xi).distance; };
		}
		else
		{
			func = [&pd](Vector2d const& xi) { return pd.signed_distance(xi).distance; };
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
			output_file += ".cdf2d";
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