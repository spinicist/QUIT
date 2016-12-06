/*
 *  Option.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_OPTION_H
#define QI_OPTION_H

#include <memory>
#include <iostream>
#include <iomanip>
#include <forward_list>
#include <algorithm>
#include <string>
#include <sstream>
#include <deque>

#include "QI/Macro.h"
#include "QI/Types.h"
#include "QI/Util.h"

namespace QI {

// Forward declarations
class OptionBase;
class OptionList;
OptionList &DefaultOptions();

class OptionBase {
protected:
    const char m_short;
    const std::string m_long, m_usage;
    bool m_set = false;

public:
    OptionBase(const char s, const std::string &l, const std::string &u, OptionList &opts);

    char shortOption() const { return m_short; }
    const std::string &longOption() const { return m_long; }
    const std::string &usage() const { return m_usage; }

    virtual bool hasArgument() const = 0;
    virtual void setOption(const std::string &a) = 0;
    bool set() { return m_set; } 
};
std::ostream &operator<< (std::ostream &os, const OptionBase &o);

class OptionList {
protected:
    typedef std::forward_list<OptionBase *> TList;
    TList m_options;
    std::string m_help;

public:
    OptionList();
    OptionList(const std::string &h);
    void add(OptionBase *o);
    std::deque<std::string> parse(int argc, char *const *argv);
    void setHelp(const std::string &h);
    void print(std::ostream & os) const;
};
std::ostream &operator<< (std::ostream &os, const OptionList &l);

class Switch : public OptionBase {
protected:
    bool m_value = false;
public:
    Switch(const char s, const std::string &l, const std::string &u, OptionList &opts) :
        OptionBase(s, l, u, opts)
    {}

    bool &operator* () { return m_value; }
    virtual bool hasArgument() const override { return false; }
    virtual void setOption(const std::string &a) override {
        if (a != "") {
            QI_EXCEPTION("Switches don't have arguments");
        }
        m_value = true;
    }

};

template<typename T> class Option : public OptionBase {
protected:
    T m_value;
public:
    Option(const T &defval, const char s, const std::string &l, const std::string &u, OptionList &opts) :
        OptionBase(s, l, u, opts),
        m_value(defval)
    {}

    T &operator* () { return m_value; }
    virtual bool hasArgument() const override { return true; }
    virtual void setOption(const std::string &a) override {
        std::istringstream as(a);
        as >> m_value;
        m_set = true;
    }
};

template<typename Scalar, int ArraySize>
class ArrayOption : public OptionBase {
public:
    typedef Eigen::Array<Scalar, ArraySize, 1> TArray;
protected:
    TArray m_value;
public:
    ArrayOption(const char s, const std::string &l, const std::string &u, OptionList &opts) :
        OptionBase(s, l, u, opts),
        m_value(TArray::Zero())
    {}

    TArray &operator* () { return m_value; }
    virtual bool hasArgument() const override { return true; }
    virtual void setOption(const std::string &a) override {
        std::istringstream as(a);
        for (int i = 0; i < m_value.rows(); i++) {
            as >> m_value[i];
        }
        m_set = true;
    }
};

template<typename TRegion>
class RegionOption : public OptionBase {
protected:
    TRegion m_region;
public:
    RegionOption(const char s, const std::string &l, const std::string &u, OptionList &opts) :
        OptionBase(s, l, u, opts)
    {}

    TRegion &operator* () { return m_region; }
    virtual bool hasArgument() const override { return true; }
    virtual void setOption(const std::string &a) override {
        std::istringstream as(a);

        typename TRegion::IndexType start;
        typename TRegion::SizeType size;
        for (int i = 0; i < TRegion::ImageDimension; i++) {
            as >> start[i];
        }
        for (int i = 0; i < TRegion::ImageDimension; i++) {
            as >> size[i];
        }
        m_region.SetIndex(start);
        m_region.SetSize(size);
        m_set = true;
    }
};

template<typename TImg>
class ImageOption : public OptionBase {
public:
    typedef typename TImg::Pointer TPtr;
protected:
    TPtr m_ptr;
public:
    ImageOption(const char s, const std::string &l, const std::string &u, OptionList &opts) :
        OptionBase(s, l, u, opts),
        m_ptr(ITK_NULLPTR)
    {}

    TPtr &operator* () { return m_ptr; }
    virtual bool hasArgument() const override { return true; }
    virtual void setOption(const std::string &a) override { m_ptr = QI::ReadImage(a); m_set = true; }
};

class EnumOption : public OptionBase {
protected:
    char m_value;
    std::string m_enumValues;
public:
    EnumOption(const std::string &e, const char d,
               const char s, const std::string &l, const std::string &u, OptionList &opts) :
        OptionBase(s, l, u, opts),
        m_enumValues(e),
        m_value(d)
    {}

    char &operator* () { return m_value; }
    virtual bool hasArgument() const override { return true; }
    virtual void setOption(const std::string &a) override {
        const char v = a[0];
        if (m_enumValues.find(v) != std::string::npos) {
            m_value = v;
        } else {
            QI_EXCEPTION("Unknown enum value " + std::string{v} + " for option " + this->longOption());
        }
        m_set = true;
    }
};

class Help : public OptionBase {
private:
    const OptionList &m_opts;

public:
    Help(OptionList &opts) :
        OptionBase('h', "help", "Print the help message.", opts),
        m_opts{opts}
    {}

    virtual bool hasArgument() const override { return false; }
    virtual void setOption(const std::string &a) override {
        std::cout << m_opts << std::endl;
        exit(EXIT_FAILURE);
    }
};

} // End namespace QI

#endif // QI_OPTION_H