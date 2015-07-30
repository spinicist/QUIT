/*
 *  Model.h
 *
 *  Created by Tobias Wood on 2015/03/12.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <vector>

#include <Eigen/Dense>

#include "SignalEquations.h"

using namespace std;
using namespace Eigen;

//cdbl and carrd already typedef'd in SignalEquations
typedef const VectorXd cvecd;

enum class FieldStrength { Three, Seven, User };
static const string to_string(const FieldStrength& f);

class Model {
public:

protected:
	bool m_scale_to_mean = false;
	ArrayXcd scale(const ArrayXcd &signal) const;

public:
	virtual string Name() const = 0;
	virtual size_t nParameters() const = 0;
	virtual bool ValidParameters(cvecd &params) const = 0;
	virtual const vector<string> &Names() const = 0;
	virtual ArrayXXd Bounds(const FieldStrength f, cdbl TR) const = 0;
    virtual ArrayXd Start(const FieldStrength f, cdbl T1, cdbl T2) const = 0;

	void setScaleToMean(bool s) { m_scale_to_mean = s; }

	virtual VectorXcd MultiEcho(cvecd &params, carrd &TE) const;
	virtual VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const;
    virtual VectorXcd SPGREcho(cvecd &p, carrd& a, cdbl TR, cdbl TE) const;
	virtual VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const;
	virtual VectorXcd MPRAGE(cvecd &params, cdbl a, cdbl TR, const int N, cvecd &TI, cdbl TD) const;
	virtual VectorXcd AFI(cvecd &params, cdbl a, cdbl TR1, cdbl TR2) const;
    virtual VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const;
    virtual VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const;
	virtual VectorXcd SSFPEllipse(cvecd &params, carrd &a, cdbl TR) const;
    virtual VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const;
};

#define DECLARE_INTERFACE( )\
public:\
	string Name() const override;\
	size_t nParameters() const override;\
	bool ValidParameters(cvecd &p) const override;\
	const vector<string> &Names() const override;\
	ArrayXXd Bounds(const FieldStrength f, cdbl TR) const override;\
    ArrayXd Start(const FieldStrength f, cdbl T1, cdbl T2) const override;\

class SCD : public Model {
	DECLARE_INTERFACE()

	virtual VectorXcd MultiEcho(cvecd &params, carrd &TE) const override;
    virtual VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
	virtual VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
	virtual VectorXcd MPRAGE(cvecd &params, cdbl a, cdbl TR, const int N, cvecd &TI, cdbl TD) const override;
	virtual VectorXcd AFI(cvecd &params, cdbl a, cdbl TR1, cdbl TR2) const override;
    virtual VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const override;
	virtual VectorXcd SSFPEllipse(cvecd &params, carrd &a, cdbl TR) const;
};

class MCD2 : public Model {
	DECLARE_INTERFACE()

	virtual VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    virtual VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
	virtual VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
    virtual VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const override;
};

class MCD2_NoEx : public Model {
    DECLARE_INTERFACE()

    virtual VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    virtual VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
};

class MCD3 : public Model {
	DECLARE_INTERFACE()

	virtual VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    virtual VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
	virtual VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
    virtual VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const override;
};

class MCD3_NoEx : public Model {
    DECLARE_INTERFACE()

    virtual VectorXcd SPGR(cvecd &p, carrd &a, cdbl TR) const override;
    virtual VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
    virtual VectorXcd SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const override;
    virtual VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
};

/*class MCD3_2OM : public Model {
	DECLARE_INTERFACE()

	virtual VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
	virtual VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
	virtual VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, cdbl phi) const override;
	virtual VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl phi) const override;
};*/

#undef DECLARE_INTERFACE
#endif // MODEL_H
