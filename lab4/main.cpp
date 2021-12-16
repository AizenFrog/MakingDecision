#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>

template <typename T>
using matrix = std::vector<std::vector<T>>;

void open_file(const std::string& file_name, std::vector<std::pair<double, double>>& distansies)
{
    std::ifstream file(file_name);
    if (file.is_open()) {
        std::size_t size;
        file >> size;
        distansies.resize(size);
        for (std::size_t i = 0; i < size; ++i) {
            std::size_t index;
            double x, y;
            file >> index >> x >> y;
            distansies[i].first = x;
            distansies[i].second = y;
        }
        file.close();
    }
    else
        assert("Error while opening file");
}

struct cluster {
    // constructors
    cluster(std::size_t size, int level) : distance(size, std::vector<double>(size)),
                                           point(size),
                                           level(level),
                                           size(size) {}

    // fields
    matrix<double> distance;
    std::vector<std::pair<double, double>> point;
    std::vector<std::size_t> way;
    std::vector<cluster*> children;
    int level;
    std::size_t size;
};

template <typename T>
inline bool is_in_vector(const std::vector<T>& vec, const T& value)
{
    for (std::size_t i = 0; i < vec.size(); ++i)
        if (vec[i] == value)
            return true;
    return false;
}

inline double euclidean_distance(const std::pair<double, double>& first, const std::pair<double, double>& second)
{
    double x = first.first - second.first;
    double y = first.second - second.second;
    return std::sqrt(x * x + y * y);
}

void set_distances(const std::vector<std::pair<double, double>>& points, matrix<double>& distances)
{
    for (std::size_t i = 0; i < points.size(); ++i)
        for (std::size_t j = 0; j < points.size(); ++j)
            distances[i][j] = euclidean_distance(points[i], points[j]);
}

template <std::uint8_t ClusterCountT>
std::vector<std::size_t> get_cluster_edges(const matrix<double>& distances)
{
    constexpr std::uint8_t count = ClusterCountT;
    std::vector<std::size_t> res;
    res.push_back(0);
    for (std::uint8_t i = 0; i < count - 1; ++i) {
        double current_max = 0;
        std::size_t max_index = -1;
        for (std::size_t j = 0; j < distances.size(); ++j) {
            if (is_in_vector(res, j))
                continue;
            double distance_sum = 0;
            for (std::size_t k = 0; k < res.size(); ++k)
                distance_sum += distances[k][j];
            if (current_max < distance_sum) {
                current_max = distance_sum;
                max_index = j;
            }
            res.push_back(max_index);
        }
    }
    return res;
}

void add_points_to_clusters(matrix<std::size_t>& clusters, const matrix<double>& distances)
{
    std::vector<std::size_t> placed_values;
    for (std::size_t i = 0; i < clusters.size(); ++i)
        placed_values.push_back(clusters[i][0]);
    for (std::size_t i = 0; i < distances.size(); ++i) {
        if (is_in_vector(placed_values, i));
        double min_distance = std::numeric_limits<double>::infinity();
        std::size_t cluster_index = -1;
        for (std::size_t j = 0; j < clusters.size(); ++j) {
            if (distances[clusters[j][0]][i] < min_distance) {
                min_distance = distances[clusters[j][0]][i];
                cluster_index = j;
            }
        }
        clusters[cluster_index].push_back(i);
    }
}

std::vector<std::pair<double, double>> get_weight_centrals(const matrix<std::size_t>& clusters,
                                                           const std::vector<std::pair<double, double>>& points)
{
    std::vector<std::pair<double, double>> res(clusters.size());
    for (std::size_t i = 0; i < clusters.size(); ++i) {
        double x = 0.0;
        double y = 0.0;
        for (std::size_t j = 0; j < clusters[i].size(); ++j) {
            x += points[clusters[i][j]].first;
            y += points[clusters[i][j]].second;
        }
        x /= clusters[i].size();
        y /= clusters[i].size();
        res[i] = std::make_pair(x, y);
    }
    return res;
}

std::size_t salesman_problem(const std::vector<std::pair<double, double>>& points, std::vector<std::size_t>& way)
{

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
        assert("Given argument is not a directory");

    std::vector<std::vector<std::pair<double, double>>> distansies(files.size());
    std::vector<std::vector<std::uint32_t>> results(files.size() * 2);
    for (std::size_t i = 0; i < files.size(); ++i)
        open_file(files[i], distansies[i]);
    return 0;
}
