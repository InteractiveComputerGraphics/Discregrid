
#include <Discregrid/All>
#include <Eigen/Dense>
#include <cxxopts/cxxopts.hpp>

#include <string>
#include <iostream>
#include <array>

#include "Discregrid/utility/bmp_file.hpp"

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

std::array<unsigned char, 3u> doubleTo5ColourHeatmap(double v, double minValue, double maxValue)
{
	const double normalised = (v - minValue) / (maxValue - minValue);
	
	constexpr double h1 = 0;
	constexpr double h2 = 0.66667;
	const double h = (1.0 - normalised) * h1 + normalised * h2;
	constexpr double s = 1.0;
	constexpr double l = 0.5;

	auto hue2rgb = [](double p, double q, double t){
		if(t < 0) t += 1;
		if(t > 1) t -= 1;
		if(t < 1/6.0) return p + (q - p) * 6 * t;
		if(t < 1/2.0) return q;
		if(t < 2/3.0) return p + (q - p) * (2/3.0 - t) * 6;
		return p;
	};

	constexpr double q = l < 0.5 ? l * (1 + s) : l + s - l * s;
	constexpr double p = 2 * l - q;

	return {
		static_cast<unsigned char>(hue2rgb(p, q, h + 1 / 3.0)	* 255),
		static_cast<unsigned char>(hue2rgb(p, q, h)			* 255),
		static_cast<unsigned char>(hue2rgb(p, q, h - 1 / 3.0)	* 255)
	};
}

}

int main(int argc, char* argv[])
{
	cxxopts::Options options(argv[0], "Transforms a slice of a discrete 2D-SDF to a bitmap image.");
	options.positional_help("[input 2D-SDF file (.cdf2d)]");

	options.add_options()
	("h,help", "Prints this help text")
	("f,field_id", "ID in which the SDF to export is stored.", cxxopts::value<unsigned int>()->default_value("0"))
	("s,samples", "Number of samples in width direction", cxxopts::value<unsigned int>()->default_value("1024"))
	("d,depth", "Relative depth value between -1 and 1 in direction of the axis orthogonal to the plane", cxxopts::value<double>()->default_value("0"))
	("o,output", "Output (in bmp format)", cxxopts::value<std::string>()->default_value(""))
	("c,colormap", "Color map options: redsequential (rs), green blue inverse diverging (gb) (suitable for visualisation of signed distance fields), 5 colour heatmap (hm) (suitable for visualisation of differences/errors)", cxxopts::value<std::string>()->default_value("gb"))
	("input", "SDF file", cxxopts::value<std::vector<std::string>>())
	;

	try
	{
		options.parse_positional("input");
		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: SDFToBitmap2D -p xz file.sdf2d" << std::endl;
			exit(0);
		}
		if (!result.count("input"))
		{
			std::cout << "ERROR: No input file given." << std::endl;
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: SDFToBitmap2D -p xz file.sdf2d" << std::endl;
			exit(1);
		}

		auto sdf = std::unique_ptr<Discregrid::DiscreteGrid2D>{};

		auto filename = result["input"].as<std::vector<std::string>>().front();
		auto lastindex = filename.find_last_of(".");
		auto extension = filename.substr(lastindex + 1, filename.length() - lastindex);

		std::cout << "Load SDF...";
		if (extension == "cdf2d")
		{
			sdf = std::make_unique<Discregrid::CubicLagrangeDiscreteGrid2D>(filename);
		}
		else
		{
			std::cout << "ERROR: Input file must be a '.sdf2d' file specifically." << std::endl;
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: SDFToBitmap2D -p xz file.sdf2D" << std::endl;
			exit(1);
		}
		std::cout << "DONE" << std::endl;

		auto const& domain = sdf->domain();
		auto diag = domain.diagonal().eval();
		
		auto dir = Vector2i::Zero().eval();
		dir(1) = 1;

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

			auto sample = Vector2d{};
			sample(dir(0)) = x;
			sample(dir(1)) = y;
			
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
		if (cm != "gb" && cm != "rs" && cm != "hm")
		{
			std::cerr << "WARNING: Unknown color map option. Fallback to mode 'gb'." << std::endl;
		}

		if (cm == "gb")
			std::transform(data.begin(), data.end(), pixels.begin(), doubleToGreenBlueInverse);
		else if (cm == "rs")
			std::transform(data.begin(), data.end(), pixels.begin(), doubleToRedSequential);
		else if (cm == "hm")
		{
			const auto min_max = std::minmax_element(
				data.begin(), data.end());
			const auto min = *min_max.first;
			const auto max = *min_max.second;
			const auto& heatmap = [min, max](double v) { return doubleTo5ColourHeatmap(v, min, max); };
			std::transform(data.begin(), data.end(), pixels.begin(), heatmap);
		}
			

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