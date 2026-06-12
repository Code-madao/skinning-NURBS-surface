#include "../include/TestCaseSelector.h"

#include <iostream>

const TestCaseData &TestCaseSelector::Select(const TestCaseLibrary &library, int argc, char **argv)
{
    const std::string caseName = (argc > 1) ? argv[1] : library.defaultCase;
    return findTestCase(library, caseName);
}

void TestCaseSelector::PrintAvailableCases(const TestCaseLibrary &library)
{
    std::cerr << "Available cases:" << std::endl;
    for (const TestCaseData &data : library.cases)
        std::cerr << "  " << data.name << std::endl;
}
