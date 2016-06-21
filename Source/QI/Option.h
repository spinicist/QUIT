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

#include "QI/Macro.h"
#include "QI/Types.h"
#include "QI/Util.h"

namespace QI {

class OptionBase {
protected:
    const char m_short;
    const std::string m_long, m_usage;

public:
    OptionBase(const char s, const std::string &l, const std::string &u);

    char shortOption() const { return m_short; }
    const std::string &longOption() const { return m_long; }
    const std::string &usage() const { return m_usage; }

    virtual bool hasArgument() const = 0;
    virtual void setValue() = 0;
    virtual void setArgument(const std::string &a) = 0;
};
std::ostream &operator<< (std::ostream &os, const OptionBase &o);

typedef std::forward_list<OptionBase *> OptionList;
OptionList &GetOptionList();
void ParseOptions (int argc, char *const *argv, std::vector<std::string> &nonopts);

class Switch : public OptionBase {
protected:
    bool m_value;

public:
    Switch(const char s, const std::string &l, const std::string &u) :
        OptionBase(s, l, u)
    {}

    virtual bool hasArgument() const override { return false; }
    virtual void setValue() override { m_value = true; }
    virtual void setArgument(const std::string &a) override {
        QI_EXCEPTION("Switches don't have arguments");
    }

    bool &operator* () { return m_value; }
};

template<typename T> class Option : public OptionBase {
protected:
    T m_value;

public:
    Option(const T &defval, const char s, const std::string &l, const std::string &u) :
        OptionBase(s, l, u),
        m_value(defval)
    {}

    T &operator* () { return m_value; }
    virtual bool hasArgument() const override { return true; }
    virtual void setValue() override { QI_EXCEPTION("Options must have arguments"); }
    virtual void setArgument(const std::string &a) override {
        std::istringstream as(a);
        as >> m_value;
    }
};

template<typename TImg>
class ImageOption : public OptionBase {
public:
    typedef typename TImg::Pointer TPtr;
protected:
    TPtr m_ptr;
public:
    ImageOption(const char s, const std::string &l, const std::string &u) :
        OptionBase(s, l, u),
        m_ptr(ITK_NULLPTR)
    {}

    TPtr &operator* () { return m_ptr; }
    virtual bool hasArgument() const override { return true; }
    virtual void setValue() override { QI_EXCEPTION("ImageOptions must have arguments"); }
    virtual void setArgument(const std::string &a) override {
        m_ptr = QI::ReadImage(a);
    }
};

class EnumOption : public OptionBase {
protected:
    char m_value;
    std::string m_enumValues;
public:
    EnumOption(const std::string &e, const char d,
               const char s, const std::string &l, const std::string &u) :
        OptionBase(s, l, u),
        m_enumValues(e),
        m_value(d)
    {}

    char &operator* () { return m_value; }
    virtual bool hasArgument() const override { return true; }
    virtual void setValue() override { QI_EXCEPTION("EnumOptions must have arguments"); }
    virtual void setArgument(const std::string &a) override {
        const char v = a[0];
        if (m_enumValues.find(v) != std::string::npos) {
            m_value = v;
        } else {
            QI_EXCEPTION("Unknown enum value " + std::string{v} + " for option " + this->longOption());
        }
    }
};

class Help : public OptionBase {
protected:
    std::string m_help;
public:
    Help(const std::string &help) :
        OptionBase('h', "help", "Print the help message."),
        m_help(help)
    {}

    virtual bool hasArgument() const override { return false; }
    virtual void setValue() override {
        std::cout << m_help << std::endl << std::endl << "Options: " << std::endl;
        for (OptionBase *o : GetOptionList()) {
            std::cout << *o << std::endl;
        }
        exit(EXIT_FAILURE);
    }
    virtual void setArgument(const std::string &a) override { QI_EXCEPTION("Help option does not have an argument.") };
};

} // End namespace QI

#endif // QI_OPTION_H