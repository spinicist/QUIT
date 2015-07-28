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

#include "itkCommand.h"
#include "itkVectorMagnitudeImageFilter.h"

#include "Types.h"

namespace QI {

const std::string &OutExt(); //!< Return the extension stored in $QUIT_EXT
time_t printStartTime();
time_t printElapsedTime(const time_t &start);
void printElapsedClock(const clock_t &clockStart, const int voxCount);
void printLoopTime(const clock_t &loopStart, const int voxCount);
std::mt19937_64::result_type RandomSeed(); // Thread-safe random seed

void writeResiduals(const typename VectorImageF::Pointer img, const std::string prefix, const bool allResids = false);

template<typename TImg = QI::ImageF>
void writeResult(const typename TImg::Pointer img, const std::string path) {
    typedef itk::ImageFileWriter<TImg> TWriter;
    auto file = TWriter::New();
    file->SetFileName(path);
    file->SetInput(img);
    file->Update();
}

class GenericMonitor : public itk::Command {
public:
	typedef GenericMonitor            Self;
	typedef itk::Command            Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	itkTypeMacro(GenericMonitor, Superclass);
	itkNewMacro(Self);
protected:
	GenericMonitor() {}
public:
	void Execute(itk::Object *caller, const itk::EventObject &event) {
		Execute((const itk::Object *)caller, event);
	}
	void Execute(const itk::Object *object, const itk::EventObject & event) {
		const itk::ProcessObject *filter = static_cast<const itk::ProcessObject *>(object);
		if (typeid(event) == typeid(itk::ProgressEvent)) {
			std::cout << "Progress: " << (filter->GetProgress()*100) << "% complete" << std::endl;
		} else {
			std::cout << "Received event: " << typeid(event).name() << std::endl;
		}
	}
};

template <typename TFilter>
class SliceMonitor : public itk::Command {
public:
	typedef SliceMonitor            Self;
	typedef itk::Command            Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	itkTypeMacro(SliceMonitor, Superclass);
	itkNewMacro(Self);
protected:
	SliceMonitor() {}
public:
	void Execute(itk::Object *caller, const itk::EventObject &event) {
		Execute((const itk::Object *)caller, event);
	}
	void Execute(const itk::Object *object, const itk::EventObject & event) {
		const TFilter *filter = static_cast<const TFilter *>(object);
		if (typeid(event) == typeid(itk::ProgressEvent)) {
			std::cout << "Completed slice " << filter->GetSliceIndex() << std::endl;
		} else if (typeid(event) == typeid(itk::IterationEvent)) {
			std::cout << "Starting slice " << filter->GetSliceIndex() << std::endl;
		} else {
			std::cout << "Received event: " << typeid(event).name() << std::endl;
		}
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
