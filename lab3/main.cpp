#include <iostream>
#include <vector>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <stack>

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

//class cell {
//public:
//    cell() : current_way(0), not_next_step(0), left(nullptr), right(nullptr), parent(nullptr) {}
//    cell(cell* parent) : current_way(0), not_next_step(0), left(nullptr), right(nullptr), parent(parent) {}
//    cell* get_right() { return right; }
//    cell* get_left() { return left; }
//    void set_right(cell* right) { this->right = right; }
//    void set_left(cell* left) { this->left = left; }
//    std::size_t get_current_way_size() { return current_way.size(); }
//    std::vector<std::uint32_t> get_current_way() { return current_way; }
//    void set_element_to_current_way(const std::uint32_t elem) {
//        current_way.push_back(elem);
//    }
//    template <class T>
//    void set_elements_to_current_way(T&& vector) {
//        current_way = vector;
//    }
//    std::size_t get_not_next_step_size() { return not_next_step.size(); }
//    std::vector<std::uint32_t> get_not_next_step() { return not_next_step; }
//    void set_element_to_not_next_step(const std::uint32_t elem) {
//        not_next_step.push_back(elem);
//    }
//    template <class T>
//    void set_elements_to_not_next_step(T&& vector) {
//        not_next_step = vector;
//    }
//    static void destroy(cell* root) {
//        cell* _left  = root->get_left();
//        cell* _right = root->get_right();
//        if (root != nullptr) {
//            delete root;
//            root = nullptr;
//        }
//        if (_left != nullptr)
//            destroy(_left);
//        if (_right != nullptr)
//            destroy(_right);
//    }
//private:
//    std::vector<std::uint32_t> current_way;
//    std::vector<std::uint32_t> not_next_step;
//    cell* left;
//    cell* right;
//    cell* parent;
//};
//
//void set_current_distances(cell* root, matrix<std::uint32_t>& distances) {
//    std::vector<std::uint32_t> current_way  = root->get_current_way();
//    std::vector<std::uint32_t> not_nex_step = root->get_not_next_step();
//    if (current_way.size() > 1) {
//        for (std::size_t i = 0; i < current_way.size() - 1; ++i) {
//            for (std::uint32_t j = 0; j < distances.size(); ++j) {
//                if (j == current_way[i + 1] || j == current_way[i]) // ?
//                    continue;
//                distances[current_way[i]][j] = -1;
//                distances[j][current_way[i + 1]] = -1;
//            }
//            distances[current_way[i + 1]][current_way[i]] = -1;
//        }
//    }
//    if (not_nex_step.size() > 0) {
//        for (std::size_t i = 0; i < not_nex_step.size(); ++i) {
//            distances[current_way[current_way.size() - 1]][not_nex_step[i]] = -1;
//        }
//    }
//}
//
//template <typename Branch, typename Bottom, typename Top>
//void BAB(const Branch& branch,
//         const Bottom& bottom,
//         const Top& top,
//         cell* root,
//         const std::uint32_t current_position,
//         std::uint32_t& record,
//         const std::vector<std::uint32_t>& time,
//         const matrix<std::uint32_t> distances)
//{
//    std::size_t size = distances.size();
//    matrix<std::uint32_t> distances_copy(distances);
//    set_current_distances(root, distances_copy);
//    std::uint32_t min_bound = bottom(distances_copy, (std::size_t)0);
//    std::uint32_t max_bound = top(distances_copy, (std::size_t)0);
//    std::cout << "record - " << record << std::endl;
//    std::cout << "size - " << size << std::endl;
//    for (std::size_t i = 0; i < size; ++i) {
//        for (std::size_t j = 0; j < size; ++j) {
//            std::cout << distances_copy[i][j] << "\t";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "min_bound - " << min_bound << std::endl;
//    std::cout << "max_bound - " << max_bound << std::endl;
//    if (max_bound < record)
//        record = max_bound;
//    if (min_bound == max_bound && (size - root->get_current_way_size()) == 1) {
//        std::vector<std::uint32_t> way = root->get_current_way();
//    }
//    else if (min_bound >= record) {
//        if (root != nullptr) {
//            delete root;
//            root = nullptr;
//        }
//        return;
//    }
//    else if (min_bound < record) {
//        std::uint32_t next = branch(root);
//        BAB(branch, bottom, top, root->get_left(), next, record, time, distances);
//        BAB(branch, bottom, top, root->get_right(), next, record, time, distances);
//    }
//}
//
//namespace base {
//
//class branch {
//private:
//    std::size_t count;
//
//    std::uint32_t get_minimum(std::vector<std::uint32_t> I, std::vector<std::uint32_t> J) const {
//        std::size_t ind1 = 0;
//        std::size_t ind2 = 0;
//        for (std::size_t i = 0; i < count; ++i) {
//            if (ind1 < I.size() && i == I[ind1]) {
//                ind1++;
//                continue;
//            }
//            if (ind2 < J.size() && i == J[ind2]) {
//                ind2++;
//                continue;
//            }
//            return i;
//        }
//        return -1;
//    }
//public:
//    branch(const std::size_t count) : count(count) {}
//    std::uint32_t operator()(cell* root) const {
//        cell* left  = new cell(root);
//        cell* right = new cell(root);
//        std::uint32_t elem = get_minimum(root->get_current_way(), root->get_not_next_step());
//        left->set_elements_to_current_way(root->get_current_way());
//        // left->set_elements_to_not_next_step(root->get_not_next_step());
//        left->set_element_to_current_way(elem);
//        right->set_elements_to_current_way(root->get_current_way());
//        right->set_elements_to_not_next_step(root->get_not_next_step());
//        right->set_element_to_not_next_step(elem);
//        root->set_left(left);
//        root->set_right(right);
//        return elem;
//    }
//};
//
//class bottom {
//public:
//    std::uint32_t operator()(const matrix<std::uint32_t>& distances, const std::size_t number) const {
//        std::vector<std::uint32_t> row_min(distances.size(), -1);
//        std::vector<std::uint32_t> col_min(distances.size(), -1);
//        // for (std::size_t j = 0; j < row_min.size(); ++j) {
//        //     std::cout << row_min[j] << "\t";
//        // }
//        // std::cout << std::endl;
//        // for (std::size_t j = 0; j < col_min.size(); ++j) {
//        //     std::cout << col_min[j] << "\t";
//        // }
//        // std::cout << std::endl;
//        for (std::size_t i = 0; i < distances.size(); ++i)
//            for (std::size_t j = 0; j < distances.size(); ++j) {
//                if (i == j)
//                    continue;
//                if (distances[i][j] < row_min[i])
//                    row_min[i] = distances[i][j];
//            }
//        // for (std::size_t j = 0; j < row_min.size(); ++j) {
//        //     std::cout << row_min[j] << "\t";
//        // }
//        // std::cout << std::endl;
//        for (std::size_t i = 0; i < distances.size(); ++i)
//            for (std::size_t j = 0; j < distances.size(); ++j) {
//                if (i == j)
//                    continue;
//                if ((distances[j][i] - row_min[j]) < col_min[i])
//                    col_min[i] = distances[j][i] - row_min[j];
//            }
//        // for (std::size_t j = 0; j < row_min.size(); ++j) {
//        //    std::cout << col_min[j] << "\t";
//        // }
//        // std::cout << std::endl;
//        std::uint32_t res = 0;
//        for (std::size_t i = 0; i < distances.size(); ++i) {
//            res += row_min[i];
//            res += col_min[i];
//        }
//        return res;
//    }
//};
//
//class top {
//private:
//    bool is_in_vector(const std::vector<uint32_t>& vec, const std::uint32_t value) const {
//        // bool res = false;
//        // for (std::size_t i = 0; i < vec.size(); ++i)
//        //     if (value == vec[i]) {
//        //         res = true;
//        //         break;
//        //     }
//        // return res;
//        return std::find(vec.begin(), vec.end(), value) != vec.end();
//    }
//public:
//    std::uint32_t operator()(const matrix<std::uint32_t>& distances, const std::size_t number) const {
//        std::uint32_t res = 0;
//        std::size_t min_index = number;
//        std::vector<std::uint32_t> way(distances.size(), -1);
//        for (std::size_t i = 0; i < distances.size() - 1; ++i) {
//            std::uint32_t cur_min = -1;
//            std::size_t cur_min_index;
//            for (std::size_t j = 0; j < distances.size(); ++j) {
//                if (min_index == j || is_in_vector(way, j))
//                    continue;
//                if (distances[min_index][j] < cur_min) {
//                    cur_min = distances[min_index][j];
//                    cur_min_index = j;
//                }
//            }
//            way[i] = min_index;
//            min_index = cur_min_index;
//            res += cur_min;
//        }
//        return res + distances[min_index][number];
//    }
//};
//
//}
//
//int main(int argc, char** argv)
//{
//    std::vector<std::string> files;
//    if (!std::filesystem::is_regular_file(argv[1])) {
//        for (const auto& file : std::filesystem::recursive_directory_iterator(argv[1])) {
//            if (std::filesystem::is_regular_file(file))
//                files.push_back(std::filesystem::path(file).string());
//        }
//    }
//    else
//        static_assert(true, "Given argument is not a directory");
//
//    std::vector<matrix<std::uint32_t>> distansies(files.size());
//    std::vector<std::vector<std::uint32_t>> times(files.size());
//    std::vector<std::vector<std::uint32_t>> results(files.size());
//    for (std::size_t i = 0; i < files.size(); ++i)
//        open_file(files[i], times[i], distansies[i]);
//
//    for (std::size_t i = 0; i < files.size(); ++i) {
//        cell* root = new cell;
//        root->set_element_to_current_way(0);
//        base::branch branch{ distansies[i].size() };
//        base::bottom bottom{};
//        base::top top{};
//        std::uint32_t record = -1;
//        BAB(branch, bottom, top, root, 0, record, times[i], distansies[i]);
//        cell::destroy(root);
//    }
//
//    return 0;
//}

class cell {
public:
    cell() : current_way(0), children(0) {}
    cell(const std::vector<std::uint32_t>& way) : current_way(way), free_cells(0), children(0) {}
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
    std::vector<std::uint32_t> get_current_way() { return current_way; }
    std::vector<std::uint32_t> get_free_cells() { return free_cells; }
    std::size_t get_child_count() { return children.size(); }
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
private:
    std::vector<std::uint32_t> current_way;
    std::vector<std::uint32_t> free_cells;
    std::vector<cell*> children;
};

void set_free_cells(const std::vector<std::uint32_t>& current_way, std::vector<std::uint32_t>& free_cells) {
    std::vector<std::uint32_t> tmp(current_way.size() + free_cells.size());
    std::iota(tmp.begin(), tmp.end(), 0);
    for (int i = current_way.size() - 1; i >= 0; --i)
        tmp.erase(tmp.begin() + current_way[i]);
    free_cells = tmp;
}

namespace base {

class branch {
public:
    std::uint32_t operator()(std::vector<std::uint32_t>& free_cells) const {
        std::sort(free_cells.begin(), free_cells.end());
        if (free_cells.size() != 0)
            return 0;
        return -1;
    }
};

class bottom {
public:
    std::uint32_t operator()(cell* root,
                             const matrix<std::uint32_t>& distances,
                             const std::vector<std::uint32_t>& times) const {
        return 0;
    }
};

class top {
public:
    std::uint32_t operator()(cell* root,
                             const matrix<std::uint32_t>& distances,
                             const std::vector<std::uint32_t>& times) const {
        std::vector<std::uint32_t> current_way = root->get_current_way();
        std::size_t size = current_way.size();
        std::uint32_t current_sum = 0;
        for (std::size_t i = 0; i < size - 1; ++i)
            current_sum += distances[current_way[i]][current_way[i + 1]];
        if (size == distances.size())
            return current_sum;
        std::vector<std::uint32_t> free_cells;
        for (std::size_t j = 0; size + j <= distances.size(); ++j) {
            free_cells.resize(distances.size() - current_way.size());
            set_free_cells(current_way, free_cells);
            std::uint32_t tmp_sum;
            std::uint32_t min_way = -1;
            std::size_t index = -1;
            for (std::size_t i = 0; i < free_cells.size(); ++i) {
                tmp_sum = current_sum + distances[current_way[current_way.size() - 1]][i];
                if (tmp_sum <= times[i] && times[i] - tmp_sum < min_way) {
                    min_way = times[i] - tmp_sum;
                    index = i;
                }
            }
            current_way.push_back(index);
        }
        return 1;
    }
};

}

template <typename Branch, typename Bottom, typename Top>
std::vector<std::uint32_t> BAB(const Branch& branch,
                               const Bottom& bottom,
                               const Top& top,
                               cell* root,
                               const matrix<std::uint32_t>& distances,
                               const std::vector<std::uint32_t>& times)
{
    std::size_t size = distances.size();
    cell* current_root;
    std::stack<cell*> parents;
    parents.push(root);
    std::vector<std::uint32_t> current_way;//(current_root->get_current_way());
    std::vector<std::uint32_t> free_cells;// (size - current_way.size());
    //set_free_cells(current_way, free_cells);
    while (true) {
        //current_way(current_root->get_current_way());
        //free_cells(size - current_way.size());
        //set_free_cells(current_way, free_cells);

        if (!parents.empty()) {
            current_root = parents.top();
            current_way = current_root->get_current_way();
            free_cells = current_root->get_free_cells();
            //if (free_cells.size() == 1)
            //    parents.pop();
        }

        // delete elements that already been here
        //std::size_t current_child_count = current_root->get_child_count();
        //free_cells.erase(free_cells->begin(), free_cells->begin() + current_child_count);

        // check stop condition
        if (free_cells.size() == 1) {
            cell* last_item = new cell(current_way);
            last_item->set_element_to_current_way(free_cells[0]);
            current_root->set_child(last_item);
            if (bottom(last_item) == top(last_item))
                return last_item->get_current_way();
            else {
                parents.pop();
                continue;
            }
        }

        // branch
        std::vector<std::uint32_t> tmp;
        for (std::size_t i = 0; i < current_root->get_child_count(); ++i) {
            tmp.push_back(current_root->get_child(i)->get_last_in_current_way());
        }
        std::vector<std::uint32_t> tmp1 = free_cells;
        for (std::size_t i = 0; i < tmp.size(); ++i) {
            for (std::size_t j = 0; j < tmp1.size(); ++j) {
                if (tmp1[j] == tmp[i]) {
                    tmp1.erase(tmp1.begin() + j);
                    break;
                }
            }
        }
        std::uint32_t next_elem_index = branch(tmp1);
        if (next_elem_index == -1) {
            parents.pop();
            continue;
        }
        current_way.push_back(tmp1[next_elem_index]);
        cell* new_item = new cell(current_way);// , free_cells);
        //new_item->set_element_to_current_way(free_cells[next_elem_index]);
        //new_item->delete_element_from_free_cells(next_elem_index);
        //free_cells.erase(free_cells.begin() + next_elem_index);
        current_root->set_child(new_item);
        parents.push(new_item);

        // cutting off unpromising solutions
        free_cells.resize(size - current_way.size());
        set_free_cells(current_way, free_cells);
        std::vector<std::size_t> not_perspective_step; // indexes in free_cells
        for (std::size_t i = 0; i < free_cells.size(); ++i)
            for (std::size_t j = 0; j < free_cells.size(); ++j) {
                if (i == j)
                    continue;
                cell a(new_item->get_current_way());
                cell b(new_item->get_current_way());
                a.set_element_to_current_way(free_cells[i]);
                b.set_element_to_current_way(free_cells[j]);
                std::uint32_t max = top(&a);
                std::uint32_t min = bottom(&b);
                if (max <= min)
                    not_perspective_step.push_back(i);
            }
        for (int i = not_perspective_step.size() - 1; i >= 0; --i)
            free_cells.erase(free_cells.begin() + not_perspective_step[i]);
        new_item->set_free_cells(free_cells);
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

    std::vector<matrix<std::uint32_t>> distansies(files.size());
    std::vector<std::vector<std::uint32_t>> times(files.size());
    std::vector<std::vector<std::uint32_t>> results(files.size());
    for (std::size_t i = 0; i < files.size(); ++i)
        open_file(files[i], times[i], distansies[i]);

    for (std::size_t i = 0; i < files.size(); ++i) {
        std::vector<std::uint32_t> tmp(distansies[i].size() - 1);
        std::iota(tmp.begin(), tmp.end(), 1);
        cell* root = new cell({0}, tmp);
        base::branch branch{};
        base::bottom bottom{};
        base::top top{};
        BAB(branch, bottom, top, root, distansies[i], times[i]);
        cell::destroy(root);
    }

    return 0;
}
