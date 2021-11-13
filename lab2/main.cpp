#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>

void OpenFile(const std::string& file_name, std::size_t& A, std::vector<std::uint32_t>& ai, std::vector<std::uint32_t>& ci)
{
    std::ifstream file(file_name);
    if (file.is_open()) {
        file >> A;
        std::size_t size;
        file >> size;
        ai.resize(size);
        ci.resize(size);
        std::uint32_t tmp;
        for (std::size_t i = 0; i < size; ++i) {
            file >> tmp;
            ai[i] = tmp;
        }
        for (std::size_t i = 0; i < size; ++i) {
            file >> tmp;
            ci[i] = tmp;
        }
        file.close();
    }
    else
        static_assert(true, "Error while opening file");
}

namespace recursive {

std::size_t recursive_method(const std::vector<std::uint32_t>& ai,
                             const std::vector<std::uint32_t>& ci,
                             const std::size_t A,
                             const std::size_t k,
                             const std::size_t n,
                             std::vector<std::size_t>& cache)
{
    if (k == 1) {
        if (cache[A] == -1)
            return cache[A] = ai[0] <= A ? ci[0] : 0;
        return cache[A];
    }
    if (A >= ai[k - 1]) {
        if (cache[(k - 1) * n + A] != -1)
            return cache[(k - 1) * n + A];
        return cache[(k - 1) * n + A] = std::max(recursive_method(ai, ci, A, k - 1, n, cache), recursive_method(ai, ci, A - ai[k - 1], k - 1, n, cache) + ci[k - 1]);
    }
    else {
        if (cache[(k - 1) * n + A] != -1)
            return cache[(k - 1) * n + A];
        return cache[(k - 1) * n + A] = recursive_method(ai, ci, A, k - 1, n, cache);
    }
}

void get_orders(std::vector<bool>& xi,
                const std::vector<std::uint32_t>& ai,
                const std::vector<std::uint32_t>& ci,
                const std::size_t A,
                const std::size_t k,
                const std::size_t max_arrived)
{
    std::vector<size_t> cache(k - 1, 0);
    if (max_arrived == 0)
        return;
    std::size_t cur_arrived = max_arrived;
    for (std::size_t i = k; i-- >= 1;) {
        if (recursive_method(ai, ci, A - ai[i - 1], i - 1, i, cache) == cur_arrived)
            xi[i - 1] = false;
        else
            xi[i - 1] = true;
    }
}

}

namespace table {

std::size_t table_method(const std::vector<std::uint32_t>& ai,
                         const std::vector<std::uint32_t>& ci,
                         const std::size_t A,
                         const std::size_t k)
{
    std::size_t* first_col  = new std::size_t[A];
    std::memset(first_col, 0, sizeof(std::size_t) * A);
    std::size_t* second_col = new std::size_t[A];

    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < A; ++j) {
            if (j + 1 >= ai[i])
                second_col[j] = std::max(first_col[j], first_col[std::max(0, (int)j - (int)ai[i])] + ci[i]);
            else
                second_col[j] = first_col[j];
        }
        std::swap(first_col, second_col);
    }

    std::size_t res = first_col[A - 1];
    delete[] first_col;
    delete[] second_col;
    return res;
}

}

namespace sort {

void standart_sort(std::vector<std::uint32_t>& ai, std::vector<uint32_t>& ci)
{
    for (std::size_t i = 0; i < ai.size() - 1; ++i) {
        for (std::size_t j = i + 1; j < ai.size(); ++j)
            if (((double)ci[i] / ai[i]) < ((double)ci[j] / ai[j])) {
                std::swap(ai[i], ai[j]);
                std::swap(ci[i], ci[j]);
            }
    }
}

void custom_sort(std::vector<std::uint32_t>& ai, std::vector<uint32_t>& ci)
{
    std::size_t mean = 0;
    for (std::size_t i = 0; i < ai.size(); ++i)
        mean += ai[i];
    mean /= ai.size();
    for (std::size_t i = 0; i < ai.size() - 1; ++i)
        for (std::size_t j = i + 1; j < ai.size(); ++j)
            if (((double)ci[i] / ai[i]) < ((double)ci[j] / ai[j])) {
                std::swap(ai[i], ai[j]);
                std::swap(ci[i], ci[j]);
            }
    for (std::size_t i = 0; i < ai.size(); ++i) {
        if (ai[i] > 2 * mean) {
            std::uint32_t tmp_a = ai[i];
            std::uint32_t tmp_c = ci[i];
            ai.erase(ai.begin() + i);
            ai.push_back(tmp_a);
            ci.erase(ci.begin() + i);
            ci.push_back(tmp_c);
        }
    }
}

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

    std::vector<std::vector<std::uint32_t>> ai(files.size());
    std::vector<std::vector<std::uint32_t>> ci(files.size());
    std::vector<std::size_t> A(files.size());
    std::vector<std::size_t> results(files.size() * 4);
    for (std::size_t i = 0; i < files.size(); ++i)
        OpenFile(files[i], A[i], ai[i], ci[i]);

    for (std::size_t i = 0; i < files.size(); ++i) {
        std::vector<std::size_t> cache(ai[i].size() * (A[i] + 1), -1);
        results[4 * i]     = recursive::recursive_method(ai[i], ci[i], A[i], ai[i].size(), A[i] + 1, cache);
        results[4 * i + 1] = table::table_method(ai[i], ci[i], A[i], ai[i].size());

        sort::standart_sort(ai[i], ci[i]);
        std::vector<std::uint32_t> copy_ai_1(ai[i].begin(), ai[i].begin() + (ai[i].size() / 3));
        std::vector<std::uint32_t> copy_ci_1(ci[i].begin(), ci[i].begin() + (ci[i].size() / 3));

        sort::custom_sort(ai[i], ci[i]);
        std::vector<std::uint32_t> copy_ai_2(ai[i].begin(), ai[i].begin() + (ai[i].size() / 3));
        std::vector<std::uint32_t> copy_ci_2(ci[i].begin(), ci[i].begin() + (ci[i].size() / 3));

        std::vector<std::size_t> cache_copy_1(copy_ai_1.size() * (A[i] + 1), -1);
        std::vector<std::size_t> cache_copy_2(copy_ai_2.size() * (A[i] + 1), -1);
        results[4 * i + 2] = recursive::recursive_method(copy_ai_1, copy_ci_1, A[i], copy_ai_1.size(), A[i] + 1, cache_copy_1);
        results[4 * i + 3] = recursive::recursive_method(copy_ai_2, copy_ci_2, A[i], copy_ai_2.size(), A[i] + 1, cache_copy_2);

        std::cout << "recursive - " << results[4 * i]
                  << ",\ttable - " << results[4 * i + 1]
                  << ",\tstandart - " << results[4 * i + 2]
                  << ",\tcustom - " << results[4 * i + 3]
                  << "\t| " << files[i] << std::endl;
    }
    std::cout << "---------------------------------------------------------" << std::endl;
    double st_norm = 0.0, cust_norm = 0.0;
    for (std::size_t i = 0; i < files.size(); ++i) {
        st_norm   += (double)results[4 * i + 2] / results[4 * i];
        cust_norm += (double)results[4 * i + 3] / results[4 * i];
    }
    std::cout << "Norm standart sort - " << st_norm / files.size() << ",\nNorm custom sort - " << cust_norm / files.size() << std::endl;

    return 0;
}