#pragma once

#include <string>
#include <vector>
#include <cstdint>

struct ShamirSubject
{
    int_fast64_t x;
    int_fast64_t y;
};

struct ShamirParameters
{
    std::vector<ShamirSubject> subjects;
    int_fast64_t p;
};

std::vector<ShamirSubject> GetShamirSubjects(const ShamirParameters &paramaters, int_fast64_t n, std::string *steps = nullptr);
ShamirSubject DoShamirReconstruction(const ShamirParameters &paramaters, std::string *steps = nullptr);
