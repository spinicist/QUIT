#include "Log.h"
#include "Monitor.h"
#include "itkProcessObject.h"

namespace QI {

void GenericMonitor::Execute(itk::Object *caller, const itk::EventObject &event) {
    Execute((const itk::Object *)caller, event);
}
void GenericMonitor::Execute(const itk::Object *object, const itk::EventObject &event) {
    const itk::ProcessObject *filter = static_cast<const itk::ProcessObject *>(object);
    if (typeid(event) == typeid(itk::ProgressEvent)) {
        float const p = std::round(filter->GetProgress() * 100);
        QI::Info("Progress: {}% complete", p);
    } else {
        QI::Info("Received event: {}", typeid(event).name());
    }
}

} // namespace QI
