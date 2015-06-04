/*
 *  Util.cpp
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Util.h"
#include <iostream>

using namespace std;

namespace QI {

const std::string &OutExt() {
	static char *env_ext = getenv("QUIT_EXT");
	static string ext;
	static bool checked = false;
	if (!checked) {
		static map<string, string> valid_ext{
			{"NIFTI", ".nii"},
			{"NIFTI_PAIR", ".img"},
			{"NIFTI_GZ", ".nii.gz"},
			{"NIFTI_PAIR_GZ", ".img.gz"},
		};
		if (!env_ext || (valid_ext.find(env_ext) == valid_ext.end())) {
			cerr << "Environment variable QUIT_EXT is not valid, defaulting to NIFTI_GZ" << endl;
			ext = valid_ext["NIFTI_GZ"];
		} else {
			ext = valid_ext[env_ext];
		}
		checked = true;
	}
	return ext;
}

time_t printStartTime() {
	time_t theTime = time(NULL);
	char timeStr[256];
	strftime(timeStr, 256, "%H:%M:%S", localtime(&theTime));
	cout << "Started at " << timeStr << endl;
	return theTime;
}

time_t printElapsedTime(const time_t &startTime) {
    time_t theTime = time(NULL);
    char timeStr[256];
    strftime(timeStr, 512, "%H:%M:%S", localtime(&theTime));
    double elapsed = difftime(theTime, startTime);
	cout << "Finished at " << timeStr << ". Elapsed time was " << elapsed << " s." << endl;
	return theTime;
}

void printElapsedClock(const clock_t &startClock, const int voxCount) {
	clock_t endClock = clock();
	float totalMilliseconds = (endClock - startClock) * (1.e3 / CLOCKS_PER_SEC);
	cout << "Total CPU time: " << totalMilliseconds << " ms" << endl;
	cout << "Average voxel CPU time: " << totalMilliseconds / voxCount << " ms" << endl;
}

void printLoopTime(const clock_t &loopStart, const int voxCount) {
	clock_t loopEnd = clock();
	if (voxCount > 0) {
		cout << voxCount << " unmasked voxels, CPU time per voxel was "
		     << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s" << endl;
	} else {
		cout << " no voxels." << endl;
	}
}

mt19937_64::result_type RandomSeed() {
	static random_device rd;
	static mt19937_64 rng;
	static bool init = false;
	mutex seed_mtx;
	if (!init) {
		rng = mt19937_64(rd());
	}
	seed_mtx.lock();
	mt19937_64::result_type r = rng();
	seed_mtx.unlock();
	return r;
}

void writeResult(const itk::Image<float, 3>::Pointer img,
                 const string path) {
	auto file = itk::ImageFileWriter<itk::Image<float, 3>>::New();
	file->SetFileName(path);
	file->SetInput(img);
	file->Update();
}

void writeResiduals(const itk::VectorImage<float, 3>::Pointer img,
                    const string prefix,
                    const bool allResids) {
	auto magFilter = itk::VectorMagnitudeImageFilter<itk::VectorImage<float, 3>, itk::Image<float, 3>>::New();
	auto magFile = itk::ImageFileWriter<itk::Image<float, 3>>::New();
	magFilter->SetInput(img);
	magFile->SetInput(magFilter->GetOutput());
	magFile->SetFileName(prefix + "residual.nii");
	magFile->Update();
	if (allResids) {
		auto to4D = itk::VectorToImageFilter<float>::New();
		auto allFile = itk::ImageFileWriter<itk::Image<float, 4>>::New();
		to4D->SetInput(img);
		allFile->SetInput(to4D->GetOutput());
		allFile->SetFileName(prefix + "residuals.nii");
		allFile->Update();
	}
}



} // End namespace QUIT
