/*
 *  IO.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_IO_H

#include <string>
#include <sstream>

#include <Eigen/Core>

#include "Macro.h"

namespace QI {

template<typename T> bool Read(const std::string &s, T &val) {
    std::istringstream stream(s);
    if (!(stream >> val)) {
        QI_EXCEPTION("Failed to parse input: " << s);
    }
    return true;
}

template<typename T> bool Read(std::istream &in, T &val) {
    std::string line;
    // Ignore comment lines. Use shell script convention
    while (in.peek() == '#') {
        if (!std::getline(in, line))
            QI_EXCEPTION("Failed to read input.");
    }
    if (!std::getline(in, line)) {
        QI_EXCEPTION("Failed to read input. Last line was: " << line);
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
            QI_EXCEPTION("Failed to read input.");
    }
    if (!std::getline(in, line)) {
        QI_EXCEPTION("Failed to read input.");
    }
    ReadArray(line, array);
}

template<typename Scalar>
void ReadArray(const std::string &s, Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> &array) {
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
void ReadArray(std::istream &in, Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> &array) {
    std::vector<Eigen::Array<Scalar, Eigen::Dynamic, 1>> rows;
    std::string line;
    
    Eigen::Array<Scalar, Eigen::Dynamic, 1> row;
    while (std::getline(in, line) && line != "") {
        ReadArray(line, row);
        rows.push_back(row);
    }
    
    array = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>(rows.size(), rows.at(0).size());
    for (int i = 0; i < rows.size(); i++) {
        if (rows.at(i).size() == array.cols()) {
            array.row(i) = rows.at(i);
        } else {
            QI_EXCEPTION("Inconsistent row sizes in input");
        }
    }
}

} // End namespace QUIT

#endif // QUIT_IO_H