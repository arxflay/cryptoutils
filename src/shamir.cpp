#include "shamir.h"
#include "algos.h"

#include <fmt/core.h>
#include <algorithm>
#include <numeric>
#include <random>

ShamirSubject DoLagrangeInterpolation(const ShamirParameters &paramaters, int_fast64_t x, std::string *steps = nullptr)
{
    x = PositiveMod(x, paramaters.p);
    int_fast64_t outY = 0;

    int counter = 0;

    const std::vector<ShamirSubject> &subjects = paramaters.subjects;
    for (auto mainIt = subjects.begin(); mainIt != subjects.end(); mainIt++)
    {
        const ShamirSubject &iSubject = *mainIt;

        int_fast64_t numerator = 1;
        int_fast64_t denominator = 1;

        for (auto secondaryIt = subjects.begin(); secondaryIt != subjects.end(); secondaryIt++)
        {
            if (secondaryIt == mainIt)
                continue;

            const ShamirSubject &jSubject = *secondaryIt;
            numerator *= (x - jSubject.x);
            denominator *= (iSubject.x - jSubject.x);

            numerator %= paramaters.p;
            denominator %= paramaters.p;
            
            if (denominator == 0)
                throw std::runtime_error(fmt::format("{} - {} = 0, denominator equals zero", iSubject.x, jSubject.x));
        }

        if (steps)
            *steps += fmt::format("c{} = {} * {}^(-1) = ", counter, numerator, denominator);

        if(denominator < 0)
        {
            denominator *= -1;
            numerator *= -1;
        }

        denominator = InverseMod(denominator, paramaters.p);
        int_fast64_t ci = (numerator * denominator);
        int_fast64_t yici = (iSubject.y * ci);
        int_fast64_t ySum = outY + yici;

        if (steps)
        {
            *steps += fmt::format("{} * {} = {}\n", numerator, denominator, ci);
            *steps += fmt::format("y{}*c{} = {} * {} = {}\n", counter, counter, iSubject.y, ci, yici);
            if (mainIt == subjects.begin())
                *steps += fmt::format("y = {}\n\n", ySum);
            else
                *steps += fmt::format("y = {} + {} = {}\n\n", outY, yici, ySum);
        }

        outY = ySum;
        counter++;
    }

    int_fast64_t finalOutY = PositiveMod(outY, paramaters.p);
    
    if (steps)
        *steps += fmt::format("y = {} mod {} = {}", outY, paramaters.p, finalOutY);   

    return ShamirSubject{ x, PositiveMod(outY, paramaters.p) };
}

std::vector<ShamirSubject> GetShamirSubjects(const ShamirParameters &paramaters, int_fast64_t n, std::string *steps)
{
    std::vector<int_fast64_t> subjectsIndicies(paramaters.p - 1);
    std::iota(subjectsIndicies.begin(), subjectsIndicies.end(), 1);
    std::mt19937 generator(std::random_device{}());
    std::shuffle(subjectsIndicies.begin(), subjectsIndicies.end(), generator);

    std::vector<ShamirSubject> out;
    const std::vector<ShamirSubject> &subjects = paramaters.subjects;
    for (int_fast64_t i = 0; i < n; i++)
    {
        int_fast64_t subjectX = subjectsIndicies.at(i);

        auto it = std::find_if(subjects.begin(), subjects.end(), [&subjects, subjectX](const ShamirSubject &subject) { return subject.x == subjectX; });
        if (it != subjects.end())
        {
            out.push_back(*it);
            continue;
        }

        if (steps)
            *steps += fmt::format("<steps for x = {}>\n", subjectX);

        out.push_back(DoLagrangeInterpolation(paramaters, subjectX, steps));
        if (steps && i + 1 != n)
            *steps += "\n\n";
    }

    return out;
}

ShamirSubject DoShamirReconstruction(const ShamirParameters &paramaters, std::string *steps)
{
    return DoLagrangeInterpolation(paramaters, 0, steps);
}
