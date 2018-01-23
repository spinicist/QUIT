
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
    {
        cereal::JSONOutputArchive archive(std::cout);
        archive(list);
    }

    std::stringstream temp;
    {
        cereal::JSONOutputArchive out_archive(temp);
        out_archive(list);
    }
    {
        cereal::JSONInputArchive  in_archive(temp);
        QI::SequenceGroup list2;
        in_archive(list2);
        cereal::JSONOutputArchive final_archive(std::cout);
        final_archive(list2);
    }
    return EXIT_SUCCESS;
}