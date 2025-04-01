/*
 *  JSON.cpp
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <fstream>

#include "JSON.h"
#include "Log.h"
#include <fstream>
#include <istream>

using namespace std::string_literals;

namespace QI {

json ReadJSON(std::istream &is) {
    return json::parse(is);
}

json ReadJSON(std::string const &path) {
    std::ifstream ifs(path);
    if (ifs) {
        QI::Info("Reading JSON from file: {}", path);
        return ReadJSON(ifs);
    } else {
        QI::Fail("Error opening file for reading: {}", path);
    }
}

std::ostream &WriteJSON(std::ostream &os, json const &doc) {
    os << doc.dump(2) << std::endl;
    return os;
}

void WriteJSON(std::string const &path, json const &doc) {
    std::ofstream ofs(path);
    if (ofs) {
        WriteJSON(ofs, doc);
    } else {
        QI::Fail("Could not open file for writing: {}", path);
    }
}

template <typename T>
Eigen::Array<T, -1, 1>
ArrayFromJSON(json const &val, std::string const &key, T const scale, Eigen::Index const sz) {
    try {
        std::vector<T> json_array;
        json_array = val[key].get<std::vector<T>>();
        if (sz > -1 && static_cast<Eigen::Index>(json_array.size()) != sz) {
            QI::Fail("JSON array {} had {} elements, expected {}", key, json_array.size(), sz);
        }
        Eigen::Array<T, -1, 1> array(json_array.size());
        for (size_t i = 0; i < json_array.size(); i++) {
            array[i] = json_array[i] * scale;
        }
        return array;
    } catch (std::exception &e) {
        QI::Fail("Error reading from JSON array {}: {}", key, e.what());
    }
}
template Eigen::ArrayXf
ArrayFromJSON(json const &val, std::string const &key, float const scale, Eigen::Index const sz);
template Eigen::ArrayXd
ArrayFromJSON(json const &val, std::string const &key, double const scale, Eigen::Index const sz);
template Eigen::ArrayXi
ArrayFromJSON(json const &val, std::string const &key, int const scale, Eigen::Index const sz);

template <typename T>
Eigen::Array<T, -1, -1> MatrixFromJSON(json const        &val,
                                       std::string const &key,
                                       T const            scale,
                                       Eigen::Index const irows,
                                       Eigen::Index const icols) {
    try {
        auto const json_rows = val[key];
        auto const nrows     = json_rows.size();
        if (irows > -1 && nrows != irows) {
            QI::Fail("JSON matrix {} had {} rows, expected {}", key, nrows, irows);
        }
        if (nrows < 1) {
            QI::Fail("JSON matrix {} must have at least 1 row", key);
        }
        auto const ncols = json_rows.front().size();
        if (icols > -1 && ncols != icols) {
            QI::Fail("JSON matrix {} had {} columns, expected {}", key, ncols, icols);
        }
        Eigen::Array<T, -1, -1> matrix(nrows, ncols);
        for (int ir = 0; ir < nrows; ir++) {
            if (json_rows[ir].size() != ncols) {
                QI::Fail("JSON matrix {} row {} had {} columns, expected {}",
                         key,
                         ir,
                         json_rows[ir].size(),
                         ncols);
            }
            std::vector<T> const row = json_rows[ir].get<std::vector<T>>();
            for (int ic = 0; ic < ncols; ic++) {
                matrix(ir, ic) = row[ic];
            }
        }
        return matrix;
    } catch (std::exception &e) {
        QI::Fail("Error reading from JSON array {}: {}", key, e.what());
    }
}
template Eigen::Array<float, -1, -1>  MatrixFromJSON(json const        &val,
                                                     std::string const &key,
                                                     float const        scale,
                                                     Eigen::Index const rows,
                                                     Eigen::Index const cols);
template Eigen::Array<double, -1, -1> MatrixFromJSON(json const        &val,
                                                     std::string const &key,
                                                     double const       scale,
                                                     Eigen::Index const rows,
                                                     Eigen::Index const cols);
template Eigen::Array<int, -1, -1>    MatrixFromJSON(json const        &val,
                                                     std::string const &key,
                                                     int const          scale,
                                                     Eigen::Index const rows,
                                                     Eigen::Index const cols);

template <typename T> void GetJSON(json const &j, std::string const &key, T &val) {
    try {
        j.at(key).get_to(val);
    } catch (std::exception &e) {
        QI::Fail("Error reading from JSON value {}: {}", key, e.what());
    }
}
template void GetJSON<float>(json const &j, std::string const &key, float &val);
template void GetJSON<double>(json const &j, std::string const &key, double &val);
template void GetJSON<int>(json const &j, std::string const &key, int &val);

} // End namespace QI
