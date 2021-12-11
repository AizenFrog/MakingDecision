#include <iostream>
#include <vector>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <stack>
#include <random>
#include <limits>

template <typename T>
using matrix = std::vector<std::vector<T>>;

void open_file(const std::string& file_name, std::vector<std::uint32_t>& times, matrix<std::uint32_t>& distansies)
{
    std::ifstream file(file_name);
    if (file.is_open()) {
        std::size_t size;
        file >> size;
        times.resize(size);
        distansies.resize(size + 1);
        for (std::size_t i = 0; i < size; ++i)
            file >> times[i];
        for (std::size_t i = 0; i < size + 1; ++i) {
            distansies[i].resize(size + 1);
            for (std::size_t j = 0; j < size + 1; ++j)
                file >> distansies[i][j];
        }
        file.close();
    }
    else
        static_assert(true, "Error while opening file");
}

void print_way(const std::vector<std::uint32_t>& way) {
    for (const std::uint32_t& item : way)
        std::cout << item << ", ";
    std::cout << std::endl;
}

class cell {
public:
    cell() : current_way(0), children(0) {}
    cell(const std::vector<std::uint32_t>& way) : bottom(-1), top(-1), current_way(way), free_cells(0), children(0) {}
    cell(const std::vector<std::uint32_t>& way, const std::vector<std::uint32_t>& free) : current_way(way), free_cells(free), children(0) {}
    cell* get_child(const std::size_t number) { return children[number]; }
    void set_child(cell* child) {
        children.push_back(child);
    }
    void set_element_to_current_way(const std::uint32_t elem) {
        current_way.push_back(elem);
    }
    std::uint32_t get_last_in_current_way() { return current_way[current_way.size() - 1]; }
    void set_free_cells(const std::vector<std::uint32_t>& free) {
        free_cells = free;
    }
    void delete_element_from_free_cells(const std::size_t number) {
        free_cells.erase(free_cells.begin() + number);
    }
    std::size_t get_child_count() { return children.size(); }
    std::uint32_t get_bottom() { return bottom; }
    std::uint32_t get_top() { return top; }
    void set_bottom(std::uint32_t b) { bottom = b; }
    void set_top(std::uint32_t t) { top = t; }
    static void destroy(cell* root) {
        for (std::size_t i = 0; i < root->children.size(); ++i) {
            cell* tmp = root->get_child(i);
            if (tmp != nullptr)
                destroy(tmp);
        }
        if (root != nullptr) {
            delete root;
            root = nullptr;
        }
    }
public:
    std::uint32_t bottom = -1;
    std::uint32_t top = -1;
    std::vector<std::uint32_t> current_way;
    std::vector<std::uint32_t> free_cells;
    std::vector<cell*> children;
    bool condition = true;
};

void set_free_cells(std::vector<std::uint32_t> current_way, std::vector<std::uint32_t>& free_cells) {
    std::vector<std::uint32_t> tmp(current_way.size() + free_cells.size());
    std::iota(tmp.begin(), tmp.end(), 0);
    std::sort(current_way.begin(), current_way.end());
    for (int i = current_way.size() - 1; i >= 0; --i)
        tmp.erase(tmp.begin() + current_way[i]);
    free_cells = tmp;
}

std::vector<std::uint32_t> get_finded_way(cell* root,
                                          const matrix<std::uint32_t>& distances,
                                          const std::vector<std::uint32_t>& times) {
    std::vector<std::uint32_t> current_way = root->current_way;
    std::size_t size = current_way.size();
    std::uint32_t current_sum = 0;
    std::uint32_t res = 0;
    for (std::size_t i = 0; i < size - 1; ++i) {
        current_sum += distances[current_way[i]][current_way[i + 1]];
        if (current_sum > times[current_way[i + 1] - 1])
            ++res;
    }
    std::vector<std::uint32_t> free_cells;
    for (std::size_t j = 0; size + j < distances.size(); ++j) {
        free_cells.resize(distances.size() - current_way.size());
        set_free_cells(current_way, free_cells);
        std::uint32_t tmp_sum = 0;
        std::uint32_t current_min_way = -1;
        std::uint32_t min_delta = -1;
        std::uint32_t min_way = -1;
        std::size_t min_way_index = -1;
        std::size_t index = -1;
        for (std::size_t i = 0; i < free_cells.size(); ++i) {
            tmp_sum = current_sum + distances[current_way[current_way.size() - 1]][free_cells[i]];
            if (distances[current_way[current_way.size() - 1]][free_cells[i]] < min_way) {
                min_way = distances[current_way[current_way.size() - 1]][free_cells[i]];
                min_way_index = i;
            }
            if (tmp_sum <= times[free_cells[i] - 1] && times[free_cells[i] - 1] - tmp_sum < min_delta) {
                min_delta = times[free_cells[i] - 1] - tmp_sum;
                current_min_way = distances[current_way[current_way.size() - 1]][free_cells[i]];
                index = i;
            }
        }
        if (min_delta == -1) {
            current_sum += min_way;
            current_way.push_back(free_cells[min_way_index]);
            ++res;
        }
        else {
            current_sum += current_min_way;
            current_way.push_back(free_cells[index]);
        }
    }
    return current_way;
}

namespace base {

class branch {
public:
    cell* operator()(cell* root) const {
        for (std::size_t i = 0; i < root->get_child_count(); ++i)
            if (root->children[i]->condition == true) {
                root->children[i]->condition = false;
                return root->children[i];
            }
        return nullptr;
    }
};

class bottom {
public:
    std::uint32_t operator()(cell* root,
                             const matrix<std::uint32_t>& distances,
                             const std::vector<std::uint32_t>& times) const {
        std::vector<std::uint32_t> current_way = root->current_way;
        std::size_t size = current_way.size();
        std::uint32_t current_sum = 0;
        std::uint32_t res = 0;
        for (std::size_t i = 0; i < size - 1; ++i) {
            current_sum += distances[current_way[i]][current_way[i + 1]];
            if (current_sum > times[current_way[i + 1] - 1])
                ++res;
        }
        std::vector<std::uint32_t> free_cells;
        free_cells.resize(distances.size() - current_way.size());
        set_free_cells(current_way, free_cells);
        for (std::size_t i = 0; i < free_cells.size(); ++i) {
            std::uint32_t tmp = current_sum + distances[current_way[current_way.size() - 1]][free_cells[i]];
            if (tmp >= times[free_cells[i] - 1])
                ++res;
        }
        return res;
    }
};

class top {
public:
    std::uint32_t operator()(cell* root,
                             const matrix<std::uint32_t>& distances,
                             const std::vector<std::uint32_t>& times) const {
        std::vector<std::uint32_t> current_way = root->current_way;
        std::size_t size = current_way.size();
        std::uint32_t current_sum = 0;
        std::uint32_t res = 0;
        for (std::size_t i = 0; i < size - 1; ++i) {
            current_sum += distances[current_way[i]][current_way[i + 1]];
            if (current_sum > times[current_way[i + 1] - 1])
                ++res;
        }
        std::vector<std::uint32_t> free_cells;
        for (std::size_t j = 0; size + j < distances.size(); ++j) {
            free_cells.resize(distances.size() - current_way.size());
            set_free_cells(current_way, free_cells);
            std::uint32_t tmp_sum = 0;
            std::uint32_t current_min_way = -1;
            std::uint32_t min_delta = -1;
            std::uint32_t min_way = -1;
            std::size_t min_way_index = -1;
            std::size_t index = -1;
            for (std::size_t i = 0; i < free_cells.size(); ++i) {
                tmp_sum = current_sum + distances[current_way[current_way.size() - 1]][free_cells[i]];
                if (distances[current_way[current_way.size() - 1]][free_cells[i]] < min_way) {
                    min_way = distances[current_way[current_way.size() - 1]][free_cells[i]];
                    min_way_index = i;
                }
                if (tmp_sum <= times[free_cells[i] - 1] && times[free_cells[i] - 1] - tmp_sum < min_delta) {
                    min_delta = times[free_cells[i] - 1] - tmp_sum;
                    current_min_way = distances[current_way[current_way.size() - 1]][free_cells[i]];
                    index = i;
                }
            }
            if (min_delta == -1) {
                current_sum += min_way;
                current_way.push_back(free_cells[min_way_index]);
                ++res;
            }
            else {
                current_sum += current_min_way;
                current_way.push_back(free_cells[index]);
            }
        }
        return res;
    }
};

}

namespace custom {

class branch {
public:
    cell* operator()(cell* root) const {
        std::size_t index = -1;
        std::uint32_t min = -1;
        for (std::size_t i = 0; i < root->get_child_count(); ++i) {
            if (root->children[i]->condition == true) {
                if (root->children[i]->top - root->children[i]->bottom < min) {
                    min = root->children[i]->top - root->children[i]->bottom;
                    index = i;
                }
            }
        }
        if (index != -1) {
            root->children[index]->condition = false;
            return root->children[index];
        }
        return nullptr;
    }
};

class top {
public:
    std::uint32_t operator()(cell* root,
                             const matrix<std::uint32_t>& distances,
                             const std::vector<std::uint32_t>& times) const {
        std::vector<std::uint32_t> current_way = root->current_way;
        std::size_t size = current_way.size();
        std::uint32_t current_sum = 0;
        std::uint32_t res = 0;
        for (std::size_t i = 0; i < size - 1; ++i) {
            current_sum += distances[current_way[i]][current_way[i + 1]];
            if (current_sum > times[current_way[i + 1] - 1])
                ++res;
        }
        std::vector<std::uint32_t> free_cells;
        for (std::size_t j = 0; size + j < distances.size(); ++j) {
            free_cells.resize(distances.size() - current_way.size());
            set_free_cells(current_way, free_cells);
            std::uint32_t min_way = -1;
            std::size_t min_way_index = -1;
            for (std::size_t i = 0; i < free_cells.size(); ++i) {
                if (distances[current_way[current_way.size() - 1]][free_cells[i]] < min_way) {
                    min_way = distances[current_way[current_way.size() - 1]][free_cells[i]];
                    min_way_index = i;
                }
            }
            current_sum += min_way;
            current_way.push_back(free_cells[min_way_index]);
            if (current_sum > free_cells[min_way_index])
                ++res;
        }
        return res;
    }
};

class bottom {
public:
    std::uint32_t operator()(cell* root,
                             const matrix<std::uint32_t>& distances,
                             const std::vector<std::uint32_t>& times) const {
        std::vector<std::uint32_t> current_way = root->current_way;
        std::size_t size = current_way.size();
        std::uint32_t current_sum = 0;
        std::uint32_t res = 0;
        for (std::size_t i = 0; i < size - 1; ++i) {
            current_sum += distances[current_way[i]][current_way[i + 1]];
            if (current_sum > times[current_way[i + 1] - 1])
                ++res;
        }
        std::vector<std::uint32_t> free_cells;
        free_cells.resize(distances.size() - current_way.size());
        set_free_cells(current_way, free_cells);
        std::uint32_t max = 0;
        std::size_t index = -1;
        for (std::size_t i = 0; i < free_cells.size(); ++i)
            if (max < current_sum + distances[current_way[current_way.size() - 1]][free_cells[i]]) {
                max = current_sum + distances[current_way[current_way.size() - 1]][free_cells[i]];
                index = i;
            }
        max *= 0.93;
        for (std::size_t i = 0; i < free_cells.size(); ++i)
            if (max >= times[free_cells[i] - 1])
                ++res;
        return res;
    }
};

}

bool is_in_vector(const std::size_t value, const std::vector<std::size_t>& vector) {
    for (std::size_t i = 0; i < vector.size(); ++i)
        if (vector[i] == value)
            return true;
    return false;
}

template <typename Bottom, typename Top>
void bad_solutions(cell* root,
                   const Bottom& bottom,
                   const Top& top,
                   std::vector<std::uint32_t>& free_cells,
                   const matrix<std::uint32_t>& distances,
                   const std::vector<std::uint32_t>& times,
                   std::size_t& real_cells) {
    std::vector<std::uint32_t> current_way = root->current_way;
    free_cells.resize(distances.size() - current_way.size());
    set_free_cells(current_way, free_cells);
    std::vector<std::size_t> not_perspective_step;

    for (std::size_t i = 0; i < free_cells.size(); ++i) {
        root->set_child(new cell(current_way));
        root->children[i]->set_element_to_current_way(free_cells[i]);
        root->children[i]->set_top(top(root->children[i], distances, times));
        root->children[i]->set_bottom(bottom(root->children[i], distances, times));
    }

    for (std::size_t i = 0; i < free_cells.size(); ++i) {
        for (std::size_t j = 0; j < free_cells.size(); ++j) {
            if (i == j)
                continue;
            if (root->children[i]->top == root->children[j]->bottom) {
                if (root->children[i]->bottom == root->children[j]->top &&
                    root->children[i]->bottom == root->children[j]->bottom)
                    continue;
            }
            if (root->children[i]->top < root->children[j]->bottom) {
                if (!is_in_vector(j, not_perspective_step))
                    not_perspective_step.push_back(j);
                root->children[j]->condition = false;
            }
        }
    }
    std::sort(not_perspective_step.begin(), not_perspective_step.end());
    for (int i = not_perspective_step.size() - 1; i >= 0; --i) {
        free_cells.erase(free_cells.begin() + not_perspective_step[i]);
        cell* tmp = root->children[not_perspective_step[i]];
        root->children.erase(root->children.begin() + not_perspective_step[i]);
        delete tmp;
    }
    real_cells += free_cells.size();
}

template <typename Branch, typename Bottom, typename Top>
std::vector<std::uint32_t> BAB(const Branch& branch,
                               const Bottom& bottom,
                               const Top& top,
                               cell* root,
                               const matrix<std::uint32_t>& distances,
                               const std::vector<std::uint32_t>& times,
                               std::size_t& tree_cells_count)
{
    std::size_t size = distances.size();
    cell* current_root;
    cell* current_min_cell = nullptr;
    bool not_last_flag = true;
    std::uint32_t min_max = -1;
    std::stack<cell*> min_cells;
    std::stack<cell*> parents;
    parents.push(root);
    ++tree_cells_count;
    std::vector<std::uint32_t> current_way(root->current_way);
    std::vector<std::uint32_t> free_cells(size - current_way.size());
    bad_solutions(root, bottom, top, free_cells, distances, times, tree_cells_count);
    root->set_free_cells(free_cells);
    while (true) {
        if (!parents.empty()) {
            current_root = parents.top();
            current_way = current_root->current_way;
            free_cells = current_root->free_cells;
        }
        else {
            current_root = min_cells.top();
            current_way = current_root->current_way;
            free_cells = current_root->free_cells;
            if (root->get_child_count() == root->free_cells.size()) {
                not_last_flag = false;
                return get_finded_way(current_root, distances, times);
            }
        }

        if (current_root->get_bottom() == -1) {
            current_root->set_bottom(bottom(current_root, distances, times));
            current_root->set_top(top(current_root, distances, times));
        }
        std::uint32_t max = current_root->get_top();
        std::uint32_t min = current_root->get_bottom();

        if (min > min_max && not_last_flag == false && current_root != root) {
            parents.pop();
            continue;
        }
        if (min >= min_max && not_last_flag == true && current_root != root) {
            parents.pop();
            continue;
        }
        if (max == min && not_last_flag == true && current_root != root) {
            if (min_max == max) {
                min_cells.push(current_root);
            }
            if (min_max > max) {
                current_min_cell = current_root;
                while (!min_cells.empty()) {
                    min_cells.pop();
                }
                min_cells.push(current_root);
                min_max = max;
            }
            parents.pop();
            continue;
        }

        if (current_way.size() == size - 1) {
            cell* last_item = new cell(current_way);
            tree_cells_count++;
            last_item->set_element_to_current_way(free_cells[0]);
            current_root->set_child(last_item);
            std::uint32_t t = top(last_item, distances, times);
            std::uint32_t b = bottom(last_item, distances, times);
            if (b == t && not_last_flag == false) {
                return last_item->current_way;
            }
            else
                parents.pop();
        }

        // branch
        cell* new_item = branch(current_root);
        if (new_item == nullptr) {
            parents.pop();
            continue;
        }
        parents.push(new_item);

        // cutting off unpromising solutions
        bad_solutions(new_item, bottom, top, free_cells, distances, times, tree_cells_count);
        new_item->set_free_cells(free_cells);
    }
}

std::pair<std::uint32_t, std::uint32_t> check_way(const std::vector<std::uint32_t>& way,
                                                  const matrix<std::uint32_t>& distances,
                                                  const std::vector<std::uint32_t>& times) {
    std::uint32_t way_length = 0;
    std::uint32_t violation_count = 0;

    for (std::size_t i = 0; i < way.size() - 1; ++i) {
        way_length += distances[way[i]][way[i + 1]];
        if (way_length > times[way[i + 1] - 1])
            ++violation_count;
    }

    return std::make_pair(way_length, violation_count);
}

std::size_t get_tree_cells(std::size_t n) {
    if (n > 2)
        return (n - 1) * get_tree_cells(n - 1) + 1;
    else if (n == 2)
        return 2;
    else if (n == 1)
        return 1;
}

std::size_t factorial(std::size_t n) {
    if (n > 2)
        return n * factorial(n - 1);
    if (n == 2)
        return 2;
}

double efficiency_ratio(std::size_t s_count, std::size_t r_count) {
    return (double)(s_count - r_count) / s_count;
}

int main(int argc, char** argv)
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::vector<std::string> files;
    if (!std::filesystem::is_regular_file(argv[1])) {
        for (const auto& file : std::filesystem::recursive_directory_iterator(argv[1])) {
            if (std::filesystem::is_regular_file(file))
                files.push_back(std::filesystem::path(file).string());
        }
    }
    else
        static_assert(true, "Given argument is not a directory");

    std::vector<matrix<std::uint32_t>> distansies(files.size());
    std::vector<std::vector<std::uint32_t>> times(files.size());
    std::vector<std::vector<std::uint32_t>> results(files.size() * 2);
    std::vector<double> ratio(files.size() * 2);
    for (std::size_t i = 0; i < files.size(); ++i)
        open_file(files[i], times[i], distansies[i]);

    matrix<std::uint32_t> right_way{ {0, 3, 1, 2}, // 1
                                     {0, 2, 1, 3}, // 2
                                     {0, 1, 3, 2, 4, 7, 6, 5, 8, 9, 10}, //3
                                     {0, 1, 2, 3, 4, 5, 9, 10, 6, 7, 8}, // 4
                                     {0, 1, 4, 5, 2, 3, 7, 6, 8, 9, 10}, // 5
                                     {0, 1, 2, 4, 5, 9, 7, 8, 13, 11, 10, 14, 12, 15, 3, 6}, // 6
                                     {0, 1, 2, 4, 3, 5, 12, 13, 11, 10, 15, 14, 8, 9, 6, 7}, //7
                                     {0, 15, 16, 9, 2, 1, 8, 13, 10, 7, 6, 3, 5, 11, 4, 12, 32, 17, 23, 29, 20, 21, 27, 28, 22,
                                      25, 31, 26, 19, 24, 18, 33, 30, 39, 37, 42, 40, 35, 38, 46, 49, 45, 44, 41, 36, 43, 47, 48, 50, 34, 14}, // 8
                                     {0, 3, 14, 4, 10, 11, 15, 7, 12, 5, 6, 13, 1, 2, 8, 9, 16, 22, 17, 24, 19, 29, 20, 30, 28, 23,
                                      18, 25, 26, 33, 31, 32, 27, 21, 43, 36, 35, 37, 50, 41, 34, 44, 46, 42, 48, 45, 38, 47, 39, 49, 40}, // 9
                                     {0, 14, 16, 11, 15, 13, 10, 9, 1, 5, 7, 8, 4, 3, 12, 2, 30, 17, 33, 22, 29, 24, 21, 18, 19, 31, 32, 26, 28,
                                      23, 25, 27, 20, 49, 39, 47, 36, 46, 41, 43, 50, 48, 40, 42, 37, 45, 44, 34, 38, 35, 6} }; //10

    base::branch branch_b{};
    base::bottom bottom_b{};
    base::top top_b{};

    custom::branch branch_c{};
    custom::bottom bottom_c{};
    custom::top top_c{};

    for (std::size_t i = 0; i < files.size(); ++i) {
        std::vector<std::uint32_t> tmp(distansies[i].size() - 1);
        std::iota(tmp.begin(), tmp.end(), 1);
        cell* root_b = new cell({0}, tmp);
        cell* root_c = new cell({0}, tmp);
        std::size_t real_cell_count_b = 0;
        std::size_t real_cell_count_c = 0;

        std::cout << "Number: " << i << std::endl;
        // base
        std::cout << "Base started..." << std::endl;
        results[i]                = BAB(branch_b, bottom_b, top_b, root_b, distansies[i], times[i], real_cell_count_b);
        cell::destroy(root_b);
        std::cout << "Base ended:" << std::endl;
        // custom
        std::cout << "Custom started..." << std::endl;
        results[files.size() + i] = BAB(branch_c, bottom_b, top_b, root_c, distansies[i], times[i], real_cell_count_c);
        cell::destroy(root_c);
        std::cout << "Custom ended:" << std::endl;

        std::pair<std::uint32_t, std::uint32_t> res_b = check_way(results[i], distansies[i], times[i]);
        std::pair<std::uint32_t, std::uint32_t> res_c = check_way(results[files.size() + i], distansies[i], times[i]);
        std::size_t s_count = factorial(times[i].size());
        std::cout << "Solutions count - " << s_count << std::endl;
        std::cout << "Base:" << std::endl;
        std::cout << "Real cell count - " << real_cell_count_b << std::endl;
        std::cout << "Efficiency ratio - " << (ratio[i] = efficiency_ratio(s_count, real_cell_count_b)) << std::endl;
        std::cout << "Custom:" << std::endl;
        std::cout << "Real cell count - " << real_cell_count_c << std::endl;
        std::cout << "Efficiency ratio - " << (ratio[files.size() + i] = efficiency_ratio(s_count, real_cell_count_c)) << std::endl;
        std::pair<std::uint32_t, std::uint32_t> res = check_way(right_way[i], distansies[i], times[i]);
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    double mean_b = 0.0, mean_c = 0.0;
    for (std::size_t i = 0; i < files.size(); ++i) {
        mean_b += ratio[i];
        mean_c += ratio[files.size() + i];
    }
    mean_b /= files.size();
    mean_c /= files.size();
    std::cout << "Mean base - " << mean_b << std::endl;
    std::cout << "Mean custom - " << mean_c << std::endl;

    return 0;
}
