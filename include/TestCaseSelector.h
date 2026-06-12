#pragma once

#include <string>

#include "TestCaseLoader.h"

namespace TestCaseSelector
{
    const TestCaseData &Select(const TestCaseLibrary &library, int argc, char **argv);
    void PrintAvailableCases(const TestCaseLibrary &library);
}
