/*
 *  Monitor.h
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_MONITOR_H
#define QI_MONITOR_H

#include "itkCommand.h"

namespace QI {

class GenericMonitor : public itk::Command {
  public:
    typedef GenericMonitor          Self;
    typedef itk::Command            Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    itkTypeMacro(GenericMonitor, Superclass);
    itkNewMacro(Self)

  protected:
    GenericMonitor() {}

  public:
    void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE;
    void Execute(const itk::Object *object, const itk::EventObject &event) ITK_OVERRIDE;
};

} // namespace QI

#endif