#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "tinynurbs/tinynurbs.h"

struct TestCaseData
{
    std::string name;
    std::string notes;
    std::vector<tinynurbs::RationalCurve3d> curves;
};

struct TestCaseLibrary
{
    std::string defaultCase;
    std::vector<TestCaseData> cases;
};

TestCaseLibrary loadTestCaseLibrary(const std::filesystem::path &jsonPath);
const TestCaseData &findTestCase(const TestCaseLibrary &library, const std::string &caseName);
