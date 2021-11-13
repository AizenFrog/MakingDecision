// lab 1

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <filesystem>

std::size_t Foo(const std::uint32_t l, const std::vector<std::uint32_t>& li)
{
    std::vector<uint32_t> snip(1, l);
    for (std::size_t i = 0; i < li.size(); ++i) {
        for (std::uint32_t j = 0; j < snip.size(); ++j) {
            if (snip[j] >= li[i]) {
                snip[j] -= li[i];
                break;
            }
            else if (snip[j] < li[i] && j == (snip.size() - 1))
                snip.push_back(l);
            else if (snip[j] < li[i])
                continue;
        }
    }
    return snip.size();
}

void OpenFile(const std::string& file_name, std::uint32_t& l, std::vector<std::uint32_t>& li)
{
    std::ifstream file(file_name);
    //std::cout << file_name << std::endl;
    if (file.is_open()) {
        file >> l;
        //std::cout << l << "\t";
        std::size_t i = 0;
        std::uint32_t tmp;
        while (!file.eof()) {
            file >> tmp;
            li.push_back(tmp);
            //std::cout << tmp << "\t";
        }
        //std::cout << std::endl;
        file.close();
    }
    else
        static_assert(true, "Error while opening file");
}

namespace base_algorithm {

void OneStep(std::vector<std::uint32_t>& li)
{
    std::sort(li.begin(), li.end(), [&](const std::uint32_t& a, const std::uint32_t& b) -> bool {
        return a < b;
    });
}

void ManySteps(std::vector<std::uint32_t>& li, const std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i) {
        std::random_device rd;
        std::mt19937 re(rd());
        std::shuffle(li.begin(), li.end(), re); // random shuffle of array
    }
}

}

namespace custom_algorithm {

void OneStep(std::vector<std::uint32_t>& li) 
{
    std::sort(li.begin(), li.end(), [&](const std::uint32_t& a, const std::uint32_t& b) -> bool {
        return a < b;
    });
    for (std::size_t i = 1; i < li.size(); i += 2)
        std::swap(li[i], li[li.size() - i]);
}

void ManySteps1(std::vector<std::uint32_t>& li, const std::size_t n)
{
    std::sort(li.begin(), li.end());
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t j = li.size() - 2;
        while (j != -1 && li[j] >= li[j + 1]) j--;
        if (j == -1)
            break;
        std::size_t k = li.size() - 1;
        while (li[j] >= li[k]) k--;
        std::swap(li[j], li[k]);
        std::size_t l = j + 1, r = li.size() - 1;
        while (l < r)
            std::swap(li[l++], li[r--]);
    }
    for (std::size_t i = 0; i < li.size(); ++i) {
        std::cout << li[i] << "\t";
        if (i == li.size() - 1)
            std::cout << std::endl << std::endl;
    }
}

void ManySteps2(std::vector<std::uint32_t>& li, std::size_t start, std::size_t end, std::size_t n)
{
    //std::cout << start << "\t" << end << "\t" << n << std::endl;
    if ((static_cast<std::int64_t>(end) - static_cast<std::int64_t>(start)) <= 1 || n-- == 0)
        return;
    std::vector<std::uint32_t>::iterator start_iter = li.begin() + start;
    std::vector<std::uint32_t>::iterator end_iter   = li.begin() + end;
    std::vector<std::uint32_t>::iterator half_iter  = li.begin() + end/2;

    if (end/2 >= start) {
        std::random_device rd;
        std::mt19937 re(rd());
        std::sort(start_iter, end_iter, [&](const std::uint32_t& a, const std::uint32_t& b) -> bool {
            return a < b;
        });
        for (std::size_t i = start + 1; i < end; i += 2)
            std::swap(li[i], li[end - i]);
        std::shuffle(start_iter, half_iter, re);
        std::shuffle(half_iter + 1, end_iter, re);
        ManySteps2(li, start, end/2, n);
        ManySteps2(li, end/2 + 1, end, n);
    }
}

}

double Optimum(std::uint32_t l, std::vector<std::uint32_t>& li)
{
    std::uint32_t sum = 0;
    for (std::size_t i = 0; i < li.size(); ++i)
        sum += li[i];
    return static_cast<double>(sum) / l;
}

void PrintResults(const std::vector<std::size_t>& res, const std::vector<double>& opt, std::vector<std::string>& file_names)
{
    size_t sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
    double sum5 = 0.0;
    std::cout << "base_os\tbase_ms\tcust_os\tcust_ms\toptimum\tfile" << std::endl;
    for (size_t i = 0; i < opt.size(); ++i) {
        std::cout << res[i * 4] << "\t"
                  << res[i * 4 + 1] << "\t"
                  << res[i * 4 + 2] << "\t"
                  << res[i * 4 + 3] << "\t"
                  << opt[i] << "\t"
                  << file_names[i] << std::endl;
        sum1 += res[i * 4];
        sum2 += res[i * 4 + 1];
        sum3 += res[i * 4 + 2];
        sum4 += res[i * 4 + 3];
        sum5 += opt[i];
    }
    std::cout << static_cast<double>(sum1) / opt.size() << "\t"
              << static_cast<double>(sum2) / opt.size() << "\t"
              << static_cast<double>(sum3) / opt.size() << "\t"
              << static_cast<double>(sum4) / opt.size() << "\t"
              << sum5 / opt.size() << std::endl;
}

int main(int argc, char** argv)
{
    std::vector<std::string> files;
    if (!std::filesystem::is_regular_file(argv[1])) {
        for (const auto& file : std::filesystem::recursive_directory_iterator(argv[1])) {
            if (std::filesystem::is_regular_file(file))
                files.push_back(std::filesystem::path(file).string());
        }
    }
    else
        static_assert(true, "Given argument is not a directory");
    
    std::vector<std::vector<std::uint32_t>> li(files.size());
    std::vector<std::uint32_t> l(files.size());
    std::vector<std::size_t> results(4 * files.size());
    std::vector<double> optimum(files.size());
    for (std::size_t i = 0; i < files.size(); ++i)
        OpenFile(files[i], l[i], li[i]);
    
    for (std::size_t i = 0; i < files.size(); ++i) {
        std::vector<std::vector<std::uint32_t>> local_li(4, li[i]);

        base_algorithm::OneStep(local_li[0]);
        results[i * 4] = Foo(l[i], local_li[0]);

        base_algorithm::ManySteps(local_li[1], 100);
        results[i * 4 + 1] = Foo(l[i], local_li[1]);

        custom_algorithm::OneStep(local_li[2]);
        results[i * 4 + 2] = Foo(l[i], local_li[2]);

        custom_algorithm::ManySteps2(local_li[3], 0, local_li[3].size() - 1, 100);
        results[i * 4 + 3] = Foo(l[i], local_li[3]);

        optimum[i] = Optimum(l[i], li[i]);
    }

    PrintResults(results, optimum, files);

    return 0;
}