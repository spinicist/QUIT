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
        QI::Info(true, "Progress: {}% complete", round(filter->GetProgress() * 100));
    } else {
        QI::Info(true, "Received event: {}", typeid(event).name());
    }
}

} // namespace QI
