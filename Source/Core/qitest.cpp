
#include "SequenceGroup.h"

int main(int argc, char **argv) {
    QI::SequenceGroup list;
    std::cout << "Creating temp" << std::endl;
    auto t = QI::SequenceWrapper(QI::SPGRSequence());
    std::cout << "Adding temp to list" << std::endl;
    list.addSequence(t);
    std::cout << "Direct add to list" << std::endl;
    list.addSequence(QI::SequenceWrapper(QI::SPGREchoSequence()));
    std::cout << "List sequences: " << std::endl;
    cereal::JSONOutputArchive archive(std::cout);
    archive(list);
    std::cout << "Individual sequences: " << std::endl;
    archive(QI::SPGRSequence());
    archive(QI::SequenceWrapper(QI::SPGRSequence()));
    return EXIT_SUCCESS;
}