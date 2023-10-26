
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Eigen>

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
        data.emplace_back(x,y,z);
    }
    inFile.close();

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
        vector -= data[0];
        outFile << "(" << vector.x() << ", " << vector.z() << ", " << vector.y() << ")\n";
    }

    outFile.close();

}
