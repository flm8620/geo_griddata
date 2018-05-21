#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <random>
#include <vector>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <chrono>
#include <cstdio>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;


typedef CGAL::Triangulation_vertex_base_with_info_2<double, K> Vb;
typedef CGAL::Triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;

typedef CGAL::Delaunay_triangulation_2<K, Tds> DT;
typedef K::Point_2 Point;


typedef std::chrono::high_resolution_clock high_resolution_clock;
typedef std::chrono::milliseconds milliseconds;

constexpr double INVALID_VALUE = -32767;


// following code from http://cgal-discuss.949826.n4.nabble.com/Efficient-insertion-of-points-with-info-td998300.html

//spatial sort traits to use with a pair of point pointers and integer. 
struct Traits_for_spatial_sort : public DT::Geom_traits {
    typedef DT::Geom_traits Gt;
    typedef std::pair<const DT::Point*, double> Point_2;

    struct Less_x_2 {
        bool operator()(const Point_2& p, const Point_2& q) const {
            return Gt::Less_x_2()(*(p.first), *(q.first));
        }
    };

    struct Less_y_2 {
        bool operator()(const Point_2& p, const Point_2& q) const {
            return Gt::Less_y_2()(*(p.first), *(q.first));
        }
    };

    Less_x_2 less_x_2_object() const { return Less_x_2(); }
    Less_y_2 less_y_2_object() const { return Less_y_2(); }
private:
};


//function inserting points into a triangulation 
//and setting the info field to the order in the input list. 
void build_triangulation_with_value(const std::vector<Point>& pts, const std::vector<double>& values, DT& T)
{
    std::vector<std::pair<const DT::Point*, double> > points;
    points.reserve(pts.size());
    for (int i = 0; i < pts.size(); i++)
        points.emplace_back(&(pts[i]), values[i]);

    std::shuffle(points.begin(), points.end(), std::mt19937(std::random_device()()));
    spatial_sort(points.begin(), points.end(), Traits_for_spatial_sort());

    DT::Face_handle hint;
    for (auto& point : points) {
        DT::Locate_type lt;
        int li;
        DT::Face_handle f = T.locate(*(point.first), lt, li, hint);

        DT::Vertex_handle v = T.insert(*(point.first), lt, f, li);
        if (v == DT::Vertex_handle())
            hint = f;
        else {
            v->info() = point.second;
            hint = v->face();
        }
    }
}





void griddata(const std::vector<Point>& pts, const std::vector<double>& values,
    const std::vector<Point>& grid, std::vector<double>& grid_values) {
    grid_values.resize(grid.size());
    DT dt;
    auto start = high_resolution_clock::now();
    build_triangulation_with_value(pts, values, dt);

    std::cout << "\tBuild Delaunay triangulation used "
        << std::chrono::duration_cast<milliseconds>(high_resolution_clock::now() - start).count() << "ms" << std::endl;
    start = high_resolution_clock::now();
    DT::Face_handle hint = DT::Face_handle();

    for (int i = 0; i < grid.size(); i++) {
        Point query = grid[i];
        DT::Face_handle f = dt.locate(grid[i], hint);
        hint = f;
        if (dt.is_infinite(f)) {
            grid_values[i] = INVALID_VALUE;
        }
        else {

            double value[3];
            for (int k = 0; k < 3; k++) value[k] = f->vertex(k)->info();
            if (value[0] == INVALID_VALUE || value[1] == INVALID_VALUE || value[2] == INVALID_VALUE) {
                grid_values[i] = INVALID_VALUE;
                continue;
            }
            Point p[3];
            for (int k = 0; k < 3; k++) p[k] = f->vertex(k)->point();
            double area[3];
            area[0] = CGAL::area(query, p[1], p[2]);
            area[1] = CGAL::area(p[0], query, p[2]);
            area[2] = CGAL::area(p[0], p[1], query);
            double area_all = area[0] + area[1] + area[2];
            double interpo = 0.;
            for (int k = 0; k < 3; k++) interpo += area[k] / area_all*value[k];
            grid_values[i] = interpo;
        }
    }
    std::cout << "\tValue query used "
        << std::chrono::duration_cast<milliseconds>(high_resolution_clock::now() - start).count() << "ms" << std::endl;
}

void do_work(const char* FILE_NAME, const char* OUT_FILE_NAME) {
    auto start = high_resolution_clock::now();
    std::vector<double> aaigol;
    std::vector<unsigned char> mask;
    std::vector<float> latgol;
    std::vector<float> longol;
    std::cout << "\tReading file: " << FILE_NAME << std::endl;
    try {
        HighFive::File file(FILE_NAME, HighFive::File::ReadOnly);
        HighFive::DataSet aaigol_set = file.getDataSet("AAI");
        HighFive::DataSet mask_set = file.getDataSet("mask");
        HighFive::DataSet latgol_set = file.getDataSet("Latitude");
        HighFive::DataSet longol_set = file.getDataSet("Longitude");
        aaigol_set.read(aaigol);
        mask_set.read(mask);
        latgol_set.read(latgol);
        longol_set.read(longol);
    }
    catch (std::exception& e) {
        std::cout << "Error : " << e.what() << std::endl;
        return;
    }


    int N = aaigol.size();
    std::cout << "\t" << N << " points" << std::endl;
    std::vector<Point> pts(N);
    for (int i = 0; i < N; i++) {
        pts[i] = Point(latgol[i], longol[i]);
        if (mask[i] != 0) aaigol[i] = INVALID_VALUE;
    }
    std::vector<Point> grid;
    std::vector<float> new_lat, new_lon;
    for (int lat = -90; lat < 90; lat++) {
        for (int lon = -180; lon < 180; lon++) {
            grid.emplace_back(lat, lon);
            new_lat.push_back(lat);
            new_lon.push_back(lon);
        }
    }
    std::vector<double> grid_values(grid.size());
    DT dt;
    griddata(pts, aaigol, grid, grid_values);

    std::vector<unsigned char> new_mask(grid.size(), 0);
    for (int i = 0; i < grid.size(); i++) if (grid_values[i] == INVALID_VALUE) new_mask[i] = 1;

    std::cout << "\tGridding : " << std::chrono::duration_cast<milliseconds>(high_resolution_clock::now() - start).count()
        << "ms" << std::endl;

    std::vector<size_t> dims = { 180, 360 };

    std::cout << "\tWriting file: " << OUT_FILE_NAME << std::endl;
    try {
        HighFive::File file_out(OUT_FILE_NAME, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        HighFive::DataSet dataset1 = file_out.createDataSet<double>("AAI", HighFive::DataSpace(dims));
        dataset1.write<double>(grid_values.data());
        HighFive::DataSet dataset2 = file_out.createDataSet<unsigned char>("mask", HighFive::DataSpace(dims));
        dataset2.write<unsigned char>(new_mask.data());
        HighFive::DataSet dataset3 = file_out.createDataSet<float>("Latitude", HighFive::DataSpace(dims));
        dataset3.write<float>(new_lat.data());
        HighFive::DataSet dataset4 = file_out.createDataSet<float>("Longitude", HighFive::DataSpace(dims));
        dataset4.write<float>(new_lon.data());
    }
    catch (std::exception& e) {
        std::cout << "Error : " << e.what() << std::endl;
        return;
    }
}

int main() {
    for (int year = 2017; year <= 2017; year++) {
        for (int month = 1; month <= 1; month++) {
            for (int day = 1; day <= 4; day++) {
                std::cout << year << "-" << month << "-" << day << std::endl;
                char FILE_NAME[256];
                char OUT_FILE_NAME[256];
                std::sprintf(FILE_NAME, "E:/griddata/Orbit_data/OMAERO_orbit_%d-%02d-%02d.h5", year, month, day);
                std::sprintf(OUT_FILE_NAME, "E:/griddata/OMAERO_grid_%d-%02d-%02d.h5", year, month, day);
                do_work(FILE_NAME, OUT_FILE_NAME);
            }
        }
    }
    return 0;
}
