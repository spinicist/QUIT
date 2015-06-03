/*
 *  Util.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2014 Tobias Wood. All rights reserved.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_UTIL_H
#define QUIT_UTIL_H

#include <string>
#include <map>
#include <vector>
#include <random>
#include <functional>
#include <mutex>
#include <time.h>

#include <Eigen/Dense>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "Filters/VectorToImageFilter.h"
#include "Filters/ImageToVectorFilter.h"

namespace QUITK {

typedef itk::Image<float, 3> FloatImage;
typedef itk::VectorImage<float, 3> FloatVectorImage;
typedef itk::Image<float, 4> FloatTimeseries;

typedef itk::Image<std::complex<float>, 3> XFloatImage;
typedef itk::VectorImage<std::complex<float>, 3> XFloatVectorImage;
typedef itk::Image<std::complex<float>, 4> XFloatTimeseries;

typedef itk::ImageFileReader<FloatImage> ReadFloatImage;
typedef itk::ImageFileReader<FloatTimeseries> ReadFloatTimeseries;
typedef itk::ImageToVectorFilter<FloatTimeseries> FloatTimeseriesToVector;

typedef itk::ImageFileReader<XFloatImage> ReadXFloatImage;
typedef itk::ImageFileReader<XFloatTimeseries> ReadXFloatTimeseries;
typedef itk::ImageToVectorFilter<XFloatTimeseries> XFloatTimeseriesToVector;


const std::string &OutExt(); //!< Return the extension stored in $QUIT_EXT
time_t printStartTime();
time_t printElapsedTime(const time_t &start);
void printElapsedClock(const clock_t &clockStart, const int voxCount);
void printLoopTime(const clock_t &loopStart, const int voxCount);
std::mt19937_64::result_type RandomSeed(); // Thread-safe random seed

void writeResult(const typename FloatImage::Pointer img, const std::string path);
void writeResiduals(const typename FloatVectorImage::Pointer img, const std::string prefix, const bool allResids = false);

template<typename T> bool Read(const std::string &s, T &val) {
	std::istringstream stream(s);
	if (!(stream >> val)) {
		throw(std::runtime_error("Failed to parse input: " + s));
	}
	return true;
}

template<typename T> bool Read(std::istream &in, T &val) {
	std::string line;
	// Ignore comment lines. Use shell script convention
	while (in.peek() == '#') {
		if (!std::getline(in, line))
			throw(std::runtime_error("Failed to read input."));
	}
	if (!std::getline(in, line)) {
		throw(std::runtime_error("Failed to read input. Last line was: " + line));
	}
	return Read(line, val);
}

template<typename Derived> void ReadEigen(const std::string &s, const Eigen::DenseBase<Derived> &cvals) {
	std::istringstream stream(s);
	Eigen::DenseBase<Derived> &vals = const_cast<Eigen::DenseBase<Derived> &>(cvals);
	for (typename Eigen::DenseBase<Derived>::Index i = 0; i < vals.size(); i++) {
		if (!(stream >> vals[i])) {
			throw(std::runtime_error("Failed to parse input: " + s));
		}
	}
}

template<typename Derived> void ReadEigen(std::istream &in, const Eigen::DenseBase<Derived> &cvals) {
	std::string line;
	Eigen::DenseBase<Derived> &vals = const_cast<Eigen::DenseBase<Derived> &>(cvals);
	// Ignore comment lines. Use shell script convention
	while (in.peek() == '#') {
		if (!std::getline(in, line))
			throw(std::runtime_error("Failed to read input."));
	}
	if (!std::getline(in, line)) {
		throw(std::runtime_error("Failed to read input."));
	}
	ReadEigen(line, vals);
}

} // End namespace QUIT

#endif // QUIT_UTIL_H
