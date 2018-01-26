
#include "SequenceGroup.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"

int main(int argc, char **argv) {
    QI::SequenceGroup list;
    list.addSequence(std::make_shared<QI::SPGREchoSequence>());
    list.addSequence(std::make_shared<QI::SSFPEchoSequence>());
    QI::SPGRSequence spgr_seq;
    std::stringstream temp;
    {
        cereal::JSONOutputArchive archive(std::cout);
        cereal::JSONOutputArchive temp_archive(temp);
        archive(cereal::make_nvp(list.name(), list));
        // archive(cereal::make_nvp("PD", std::string("PD.nii")));
        // archive(cereal::make_nvp("T2", std::string("")));
        // temp_archive(cereal::make_nvp(list.name(), list));
        // temp_archive(cereal::make_nvp("PD", std::string("PD.nii")));
        // temp_archive(cereal::make_nvp("T2", std::string("")));
    }
    std::cout << "\n" << std::endl;
    // std::cout << "\nNOW TRYING TO READ\n" << std::endl;

    // {
    //     cereal::JSONInputArchive in_archive(temp);
    //     std::string PD, T2;
    //     in_archive(cereal::make_nvp("PD", PD));
    //     in_archive(cereal::make_nvp("T2", T2));
    //     auto A = QI::ReadSequence<QI::SequenceGroup>(in_archive, true);
    //     std::cout << "\nPD = " << PD << "\nT2 = *" << T2 << "*" << std::endl;
    // }
    return EXIT_SUCCESS;
}