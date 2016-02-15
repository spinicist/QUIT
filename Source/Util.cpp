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

#include "itkVectorMagnitudeImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkDivideImageFilter.h"

using namespace std;

namespace QI {

const std::string &GetVersion() {
    // This file contains "static const std::string Version = "version";"
    #include "BuildVersion"
    return Version;
}

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

std::string StripExt(const std::string &filename) {
    std::size_t dot = filename.find_last_of(".");
    if (dot != std::string::npos) {
        /* Deal with .nii.gz files */
        if (filename.substr(dot) == "gz") {
            dot = filename.find_last_of(".", dot);
        }
        return filename.substr(0, dot);
    } else {
        return filename;
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

void WriteResiduals(const VectorImageF::Pointer img, const string prefix,
                    const bool allResids, const ImageF::Pointer scaleImage) {
	auto magFilter = itk::VectorMagnitudeImageFilter<VectorImageF, ImageF>::New();
    auto scaleFilter = itk::DivideImageFilter<ImageF, ImageF, ImageF>::New();
    auto sqrtNFilter = itk::MultiplyImageFilter<ImageF>::New();
	magFilter->SetInput(img);
    if (scaleImage != ITK_NULLPTR) {
        scaleFilter->SetInput1(magFilter->GetOutput());
        scaleFilter->SetInput2(scaleImage);
        sqrtNFilter->SetInput(scaleFilter->GetOutput());
    } else {
        sqrtNFilter->SetInput(magFilter->GetOutput());
    }
    sqrtNFilter->SetConstant(1./sqrt(img->GetNumberOfComponentsPerPixel()));
    sqrtNFilter->Update();
    WriteImage(sqrtNFilter->GetOutput(), prefix + "residual.nii");
	if (allResids) {
		auto to4D = QI::VectorToTimeseriesF::New();
		auto allFile = QI::TimeseriesWriterF::New();
		to4D->SetInput(img);
		allFile->SetInput(to4D->GetOutput());
		allFile->SetFileName(prefix + "residuals.nii");
		allFile->Update();
	}
}



} // End namespace QUIT
