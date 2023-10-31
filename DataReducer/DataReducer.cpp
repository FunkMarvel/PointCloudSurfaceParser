// //////////////////////////////////////////////////////////////////////////
// //////////////////////////////
// //FileName: DataReducer.cpp
// //FileType: Visual C++ Source file
// //Author : Anders P. Åsbø
// //Created On : 26/10/2023
// //Last Modified On : 31/10/2023
// //Copy Rights : Anders P. Åsbø
// //Description : 
// //////////////////////////////////////////////////////////////////////////
// //////////////////////////////

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Eigen>

struct DataBounds
{
    float xmin, xmax, ymin, ymax, zmin, zmax, xExtent, yExtent, zExtent;
};

int main(int argc, char* argv[])
{
    std::string iFileName{"../RawData/merged.txt"};
    std::ifstream inFile{iFileName};

    if (!inFile.is_open())
    {
        const auto msg = "Could not read file: " + iFileName;
        std::cout << msg << std::endl;
        throw std::runtime_error(msg);
    }

    std::string oFileName{"../ProcessedData/sampledData.txt"};

    float x{}, y{}, z{};
    size_t i{};

    std::vector<Eigen::Vector3f> data;

    while (inFile >> x >> y >> z)
    {
        if (i++ % 100) continue;
        data.emplace_back(x, y, z);
    }
    inFile.close();

    DataBounds bounds{0, 0, 0, 0, 0, 0, 0, 0, 0};

    bounds.xmin = bounds.xmax = data[0].x();
    bounds.ymin = bounds.ymax = data[0].y();
    bounds.zmin = bounds.zmax = data[0].z();

    for (auto point : data)
    {
        if (bounds.xmax < point.x())
        {
            bounds.xmax = point.x();
        }
        if (bounds.xmin > point.x())
        {
            bounds.xmin = point.x();
        }

        if (bounds.ymax < point.y())
        {
            bounds.ymax = point.y();
        }
        if (bounds.ymin > point.y())
        {
            bounds.ymin = point.y();
        }

        if (bounds.zmax < point.z())
        {
            bounds.zmax = point.z();
        }
        if (bounds.zmin > point.z())
        {
            bounds.zmin = point.z();
        }
    }

    Eigen::Vector3f offset{0.5f*(bounds.xmax + bounds.xmin), 0.5f*(bounds.ymax + bounds.ymin), 0.5f*(bounds.zmax + bounds.zmin)};

    std::ofstream outFile{oFileName};
    if (!outFile.is_open())
    {
        const auto msg = "Could not open out-file: " + oFileName;
        std::cout << msg << std::endl;
        throw std::runtime_error(msg);
    }

    outFile << data.size() << "\n";

    for (auto vector : data)
    {
        vector -= offset;
        outFile << "(" << vector.x()*0.5f << ", " << vector.z()*0.5f << ", " << vector.y()*0.5f << ")\n";
    }

    outFile.close();
}
