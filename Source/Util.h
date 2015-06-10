/*
 *  Util.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2015 Tobias Wood.
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

#include "itkVectorMagnitudeImageFilter.h"

#include "Types.h"

namespace QI {

const std::string &OutExt(); //!< Return the extension stored in $QUIT_EXT
time_t printStartTime();
time_t printElapsedTime(const time_t &start);
void printElapsedClock(const clock_t &clockStart, const int voxCount);
void printLoopTime(const clock_t &loopStart, const int voxCount);
std::mt19937_64::result_type RandomSeed(); // Thread-safe random seed

void writeResult(const typename ImageF::Pointer img, const std::string path);
void writeResiduals(const typename VectorImageF::Pointer img, const std::string prefix, const bool allResids = false);

class ProgressReport : public itk::Command
{
public:
	typedef ProgressReport           Self;
	typedef itk::Command             Superclass;
	typedef itk::SmartPointer<Self>  Pointer;
	itkNewMacro( Self );
protected:
	ProgressReport() {};
public:
	void Execute(itk::Object *caller, const itk::EventObject & event)
	{
		Execute( (const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event)
	{
    const itk::ProcessObject * filter = static_cast< const itk::ProcessObject * >( object );
    if( ! itk::ProgressEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout << "Progress: " << (filter->GetProgress()*100) << "% complete" << std::endl;
	}
};

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

template<typename Scalar>
void ReadArray(const std::string &s, Eigen::Array<Scalar, Eigen::Dynamic, 1> &array) {
	std::istringstream stream(s);
	std::vector<Scalar> vals;

	Scalar temp;
	while (stream >> temp) {
		vals.push_back(temp);
	}

	array = Eigen::Array<Scalar, Eigen::Dynamic, 1>(vals.size());
	for (int i = 0; i < vals.size(); i++) {
		array[i] = vals[i];
	}
}

template<typename Scalar>
void ReadArray(std::istream &in, Eigen::Array<Scalar, Eigen::Dynamic, 1> &array) {
	std::string line;
	// Ignore comment lines. Use shell script convention
	while (in.peek() == '#') {
		if (!std::getline(in, line))
			throw(std::runtime_error("Failed to read input."));
	}
	if (!std::getline(in, line)) {
		throw(std::runtime_error("Failed to read input."));
	}
	ReadArray(line, array);
}

} // End namespace QUIT

#endif // QUIT_UTIL_H
