
#include <Discregrid/All>
#include <Eigen/Dense>
#include <cxxopts/cxxopts.hpp>

#include <string>
#include <iostream>
#include <array>

#include "bmp_file.hpp"

using namespace Eigen;

namespace
{
std::array<unsigned char, 3u> doubleToGreenBlueInverse(double v)
{
	if (v >= 0.0)
	{
		return {{0u, static_cast<unsigned char>(std::min(std::max(255.0 * (1.0 - v), 0.0), 255.0)), 0u}};
	}
	return {{0u, 0u, static_cast<unsigned char>(std::min(std::max(255.0 * (1.0 + v), 0.0), 255.0))}};
}

std::array<unsigned char, 3u> doubleToRedSequential(double v)
{
	return {{static_cast<unsigned char>(std::min(std::max(255.0 * v, 0.0), 255.0)), 0u, 0u}};
}

}

int main(int argc, char* argv[])
{
	cxxopts::Options options(argv[0], "Transforms a slice of a discrete SDF to a bitmap image.");
	options.positional_help("[input SDF file]");

	options.add_options()
	("h,help", "Prints this help text")
	("f,field_id", "ID in which the SDF to export is stored.", cxxopts::value<unsigned int>()->default_value("0"))
	("s,samples", "Number of samples in width direction", cxxopts::value<unsigned int>()->default_value("1024"))
	("p,plane", "Plane in which the image slice is extracted", cxxopts::value<std::string>()->default_value("xy"))
	("d,depth", "Relative depth value between -1 and 1 in direction of the axis orthogonal to the plane", cxxopts::value<double>()->default_value("0"))
	("o,output", "Output (in bmp format)", cxxopts::value<std::string>()->default_value(""))
	("c,colormap", "Color map options: redsequential (rs), green blue inverse diverging (gb) (suitable for visualiztion of signed distance fields)", cxxopts::value<std::string>()->default_value("gb"))
	("input", "SDF file", cxxopts::value<std::vector<std::string>>())
	;

	try
	{
		options.parse_positional("input");
		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: SDFToBitmap -p xz file.sdf" << std::endl;
			exit(0);
		}
		if (!result.count("input"))
		{
			std::cout << "ERROR: No input file given." << std::endl;
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: SDFToBitmap -p xz file.sdf" << std::endl;
			exit(1);
		}

		auto sdf = std::unique_ptr<Discregrid::DiscreteGrid>{};

		auto filename = result["input"].as<std::vector<std::string>>().front();
		auto lastindex = filename.find_last_of(".");
		auto extension = filename.substr(lastindex + 1, filename.length() - lastindex);

		std::cout << "Load SDF...";
		if (extension == "cdf" || extension == "cdm")
		{
			sdf = std::unique_ptr<Discregrid::CubicLagrangeDiscreteGrid>(
				new Discregrid::CubicLagrangeDiscreteGrid(filename));
		}
		std::cout << "DONE" << std::endl;

		auto depth = result["d"].as<double>();
		auto const& domain = sdf->domain();
		auto diag = domain.diagonal().eval();

		auto plane = result["p"].as<std::string>();
		if (plane.length() != 2 && plane[0] != plane[1])
		{
			std::cerr << "ERROR: Invalid option for plane provided. Should be one of the following options: xy, xz, yz, yx" << std::endl;
			exit(1);
		}

		auto dir = Vector3i::Zero().eval();
		if (plane[0] == 'y')
			dir(0) = 1;
		else if (plane[0] == 'z')
			dir(0) = 2;
		if (plane[1] == 'y')
			dir(1) = 1;
		else if (plane[1] == 'z')
			dir(1) = 2;
		if (dir(0) != 1 && dir(1) != 1)
			dir(2) = 1;
		if (dir(0) != 2 && dir(1) != 2)
			dir(2) = 2;

		auto xsamples = result["s"].as<unsigned int>();
		auto ysamples = static_cast<unsigned int>(std::round(diag(dir(1)) / diag(dir(0)) * static_cast<double>(xsamples)));

		auto xwidth = diag(dir(0)) / xsamples;
		auto ywidth = diag(dir(1)) / ysamples;

		auto data = std::vector<double>{};
		data.resize(xsamples * ysamples);

		auto field_id = result["f"].as<unsigned int>();

		std::cout << "Sample field...";
#pragma omp parallel for
		for (int k = 0; k < static_cast<int>(xsamples * ysamples); ++k)
		{
			auto i = k % xsamples;
			auto j = k / xsamples;

			auto xr = static_cast<double>(i) / static_cast<double>(xsamples);
			auto yr = static_cast<double>(j) / static_cast<double>(ysamples);

			auto x = domain.min()(dir(0)) + xr * diag(dir(0)) + 0.5 * xwidth;
			auto y = domain.min()(dir(1)) + yr * diag(dir(1)) + 0.5 * ywidth;

			auto sample = Vector3d{};
			sample(dir(0)) = x;
			sample(dir(1)) = y;
			sample(dir(2)) = domain.min()(dir(2)) + 0.5 * (1.0 + depth) * diag(dir(2));

			data[k] = sdf->interpolate(field_id, sample);
			if (data[k] == std::numeric_limits<double>::max())
			{
				data[k] = 0.0;
			}
		}

		std::cout << "DONE" << std::endl;

		auto min_v = *std::min_element(data.begin(), data.end());
		auto max_v = *std::max_element(data.begin(), data.end());

		auto out_file = result["o"].as<std::string>();
		if (out_file == "")
		{
			out_file = filename;
			if (out_file.find(".") != std::string::npos)
			{
				auto lastindex = out_file.find_last_of(".");
				out_file = out_file.substr(0, lastindex);
			}
			out_file += ".bmp";
		}

		std::cout << "Ouput file: " << out_file << std::endl;

		std::cout << "Export BMP...";
		std::transform(data.begin(), data.end(), data.begin(), [&max_v, &min_v](double v) {return v >= 0.0 ? v / std::abs(max_v) : v / std::abs(min_v); });

		auto pixels = std::vector<std::array<unsigned char, 3u>>(data.size());

		auto cm = result["c"].as<std::string>();
		if (cm != "gb" && cm != "rs")
		{
			std::cerr << "WARNING: Unknown color map option. Fallback to mode 'gb'." << std::endl;
		}

		if (cm == "gb")
			std::transform(data.begin(), data.end(), pixels.begin(), doubleToGreenBlueInverse);
		else if (cm == "rs")
			std::transform(data.begin(), data.end(), pixels.begin(), doubleToRedSequential);

		BmpReaderWriter::saveFile(out_file.c_str(), xsamples, ysamples, &pixels.front()[0]);
		std::cout << "DONE" << std::endl;

		std::cout << std::endl << "Statistics:" << std::endl;
		std::cout << "\tdomain         = " << domain.min().transpose() << ", " << domain.max().transpose() << std::endl;
		std::cout << "\tmin value      = " << min_v << std::endl;
		std::cout << "\tmax value      = " << max_v << std::endl;
		std::cout << "\tbmp resolution = " << xsamples << " x " << ysamples << std::endl;
	}
	catch (cxxopts::OptionException const& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
	
	return 0;
}