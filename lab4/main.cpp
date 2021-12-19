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

template <std::uint8_t ClusterCountT>
struct cluster {
    // constructors
    cluster(std::size_t size, int level) : distance(size, std::vector<double>(size)),
                                           central_distance(ClusterCountT, std::vector<double>(ClusterCountT)),
                                           point(size),
                                           level(level),
                                           size(size),
                                           children(ClusterCountT) {}

    // fields
    matrix<double> distance;
    matrix<double> central_distance;
    std::vector<std::pair<double, double>> point;
    std::vector<std::pair<double, double>> central_point;
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
                distance_sum += distances[res[k]][j];
            if (current_max < distance_sum) {
                current_max = distance_sum;
                max_index = j;
            }
        }
        res.push_back(max_index);
    }
    return res;
}

void add_points_to_clusters(matrix<std::size_t>& clusters, const matrix<double>& distances)
{
    std::vector<std::size_t> placed_values;
    for (std::size_t i = 0; i < clusters.size(); ++i)
        placed_values.push_back(clusters[i][0]);
    for (std::size_t i = 0; i < distances.size(); ++i) {
        if (is_in_vector(placed_values, i))
            continue;
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

template <std::uint8_t ClusterCountT>
void brude_force(const matrix<double>& distances,
                 std::vector<std::uint8_t>& prev,
                 std::vector<std::uint8_t>& record_way,
                 double& record_distance,
                 const std::uint8_t level)
{
    if (level == ClusterCountT) {
        double distance = 0;
        for (std::size_t i = 0; i < prev.size() - 1; ++i)
            distance += distances[prev[i]][prev[i + 1]];
        distance += distances[prev[prev.size() - 1]][0];
        if (distance < record_distance) {
            record_distance = distance;
            record_way = prev;
        }
    }

    for (std::uint8_t i = 0; i < ClusterCountT; ++i) {
        if (is_in_vector(prev, i))
            continue;
        prev.push_back(i);
        brude_force<ClusterCountT>(distances, prev, record_way, record_distance, level + 1);
        prev.pop_back();
    }
}

void greedy_algorithm(const matrix<double>& distances, std::vector<std::size_t>& record_way)
{
    record_way.push_back(0);
    double distance = 0.0;
    for (std::size_t i = 0; i < distances.size() - 1; ++i) {
        double min_distance = std::numeric_limits<double>::infinity();
        std::size_t index = -1;
        for (std::size_t j = 0; j < distances.size(); ++j) {
            if (is_in_vector(record_way, j))
                continue;
            double tmp = distance + distances[record_way[i]][j];
            if (tmp < min_distance) {
                min_distance = tmp;
                index = j;
            }
        }
        distance = min_distance;
        record_way.push_back(index);
    }
}

template <std::uint8_t ClusterCountT>
void set_childred(const matrix<std::size_t>& clusters, cluster<ClusterCountT>* root)
{
    for (std::size_t i = 0; i < clusters.size(); ++i)
        for (std::size_t j = 0; j < clusters[i].size(); ++j)
            root->children[i]->point[j] = root->point[clusters[i][j]];
}

template <std::uint8_t ClusterCountT>
void unite_way(cluster<ClusterCountT>* root,
               const matrix<std::size_t>& clusters,
               const std::vector<std::uint8_t>& record)
{
    double min = std::numeric_limits<double>::infinity();
    std::size_t first  = -1,
                second = -1;

    for (std::size_t i = 0; i < clusters[record[0]].size(); ++i) {
        std::size_t cur_min_index = -1;
        double cur_min = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < clusters[record[1]].size(); ++j) {
            if (root->distance[clusters[record[0]][root->children[record[0]]->way[i]]][clusters[record[1]][root->children[record[1]]->way[j]]] < cur_min) {
                cur_min = root->distance[clusters[record[0]][root->children[record[0]]->way[i]]][clusters[record[1]][root->children[record[1]]->way[j]]];
                cur_min_index = root->children[record[1]]->way[j];
            }
        }
        if (cur_min < min) {
            min = cur_min;
            first = root->children[record[0]]->way[i];
            second = cur_min_index;
        }
    }

    // first and second cluster union
    for (std::size_t i = first + 1; i < clusters[record[0]].size(); ++i)
        root->way.push_back(clusters[record[0]][root->children[record[0]]->way[i]]);
    for (std::size_t i = 0; i <= first; ++i)
        root->way.push_back(clusters[record[0]][root->children[record[0]]->way[i]]);
    for (std::size_t i = second; i < clusters[record[1]].size(); ++i)
        root->way.push_back(clusters[record[1]][root->children[record[1]]->way[i]]);
    for (std::size_t i = 0; i < second; ++i)
        root->way.push_back(clusters[record[1]][root->children[record[1]]->way[i]]);

    for (std::uint8_t i = 2; i < ClusterCountT; ++i) {
        std::size_t min_index = -1;
        double min = std::numeric_limits<double>::infinity();
        // serch point with min distance
        for (std::size_t j = 0; j < clusters[record[i]].size(); ++j) {
            if (root->distance[root->way[root->way.size() - 1]][clusters[record[i]][root->children[record[i]]->way[j]]] < min) {
                min = root->distance[root->way[root->way.size() - 1]][clusters[record[i]][root->children[record[i]]->way[j]]];
                min_index = root->children[record[i]]->way[j];
            }
        }
        for (std::size_t j = min_index; j < clusters[record[i]].size(); ++j)
            root->way.push_back(clusters[record[i]][root->children[record[i]]->way[j]]);
        for (std::size_t j = 0; j < min_index; ++j)
            root->way.push_back(clusters[record[i]][root->children[record[i]]->way[j]]);
    }
}

template <std::uint8_t ClusterCountT, std::uint8_t DeepT>
std::size_t salesman(cluster<ClusterCountT>* root, const std::vector<std::pair<double, double>>& points)
{
    if (points.size() >= ClusterCountT && DeepT > 0) {
        set_distances(root->point, root->distance);
        auto cluster_edge = get_cluster_edges<ClusterCountT>(root->distance);
        matrix<std::size_t> tmp_clusters(ClusterCountT);
        for (std::uint8_t i = 0; i < ClusterCountT; ++i)
            tmp_clusters[i].push_back(cluster_edge[i]);
        add_points_to_clusters(tmp_clusters, root->distance);
        root->central_point = get_weight_centrals(tmp_clusters, root->point);

        // brude force
        set_distances(root->central_point, root->central_distance);
        std::vector<std::uint8_t> buffer, record;
        double rec = std::numeric_limits<double>::infinity();
        brude_force<ClusterCountT>(root->central_distance, buffer, record, rec, 0);

        for (std::size_t i = 0; i < ClusterCountT; ++i)
            root->children[i] = new cluster<ClusterCountT>(tmp_clusters[i].size(), root->level + 1);
        set_childred(tmp_clusters, root);

        for (std::uint8_t i = 0; i < ClusterCountT; ++i)
            salesman<ClusterCountT, DeepT - 1>(root->children[i], root->children[i]->point);
        unite_way<ClusterCountT>(root, tmp_clusters, record);

        for (std::size_t i = 0; i < ClusterCountT; ++i)
            delete root->children[i];
    }
    else {
        // greedy algorithm
        set_distances(points, root->distance);
        greedy_algorithm(root->distance, root->way);
    }
    return 0;
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

    std::vector<std::vector<std::pair<double, double>>> points(files.size());
    std::vector<std::vector<std::uint32_t>> results(files.size() * 2);
    for (std::size_t i = 0; i < files.size(); ++i)
        open_file(files[i], points[i]);

    for (std::size_t i = 0; i < files.size(); ++i) {
        cluster<8>* root = new cluster<8>(points[i].size(), 0);
        root->point = points[i];
        salesman<8, 3>(root, points[i]);
        delete root;
    }

    return 0;
}
