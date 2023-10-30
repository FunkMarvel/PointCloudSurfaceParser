// //////////////////////////////////////////////////////////////////////////
// //////////////////////////////
// //FileName: PointCloudSurfaceParser.cpp
// //FileType: Visual C# Source file
// //Author : Anders P. Åsbø
// //Created On : 01/10/2023
// //Last Modified On : 01/10/2023
// //Copy Rights : Anders P. Åsbø
// //Description : 
// //////////////////////////////////////////////////////////////////////////
// //////////////////////////////

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Eigen>

struct DataBounds;
DataBounds findBounds(std::ifstream& file);
void writeVertexData(const std::string& filePath, const std::vector<std::vector<Eigen::Vector3f>>& dataGrid,
                     const Eigen::Vector3f& offset, float scaleX, float scaleY, float scaleZ);
void writeIndexFile(const std::string& filePath, int numX, int numY);

struct DataBounds
{
    float xmin, xmax, ymin, ymax, zmin, zmax, xExtent, yExtent, zExtent;
};

void processData(const std::string& fileName)
{
    std::ifstream inFile{fileName};

    if (!inFile.is_open())
    {
        const auto msg = "Could not read file: " + fileName;
        std::cout << msg << std::endl;
        throw std::runtime_error(msg);
    }

    const DataBounds bounds = findBounds(inFile);

    std::cout << bounds.xmin << " | " << bounds.xmax << " | " << bounds.ymin << " | " << bounds.ymax << " | " << bounds.
        xExtent << " | " << bounds.yExtent << std::endl;

    constexpr float stepLength = 10.f; // step length [m]
    const int numStepsX = static_cast<int>(ceil(bounds.xExtent / stepLength));
    const int numStepsY = static_cast<int>(ceil(bounds.yExtent / stepLength));

    std::vector<std::vector<std::vector<Eigen::Vector3f>>> buckets(
        numStepsX, std::vector<std::vector<Eigen::Vector3f>>(numStepsY, std::vector<Eigen::Vector3f>()));
    std::vector<std::vector<Eigen::Vector3f>> grid(
        numStepsX, std::vector<Eigen::Vector3f>(numStepsY, Eigen::Vector3f()));
    std::vector<std::vector<bool>> fillMask(numStepsX, std::vector<bool>(numStepsY, false));

    inFile.clear();
    inFile.seekg(0);

    float x{}, y{}, z{};

    while (inFile >> x >> y >> z)
    {
        int i = static_cast<int>((x - bounds.xmin) / stepLength);
        int j = static_cast<int>((y - bounds.ymin) / stepLength);

        if (i < 0) { i = 0; }
        else if (i >= numStepsX) { i = numStepsX - 1; }
        if (j < 0) { j = 0; }
        else if (j >= numStepsY) { j = numStepsY - 1; }

        buckets[i][j].emplace_back(x, y, z);
    }

    inFile.close();

    for (int i = 0; i < buckets.size(); ++i)
    {
        for (int j = 0; j < buckets[0].size(); ++j)
        {
            float meanPoint{};
            for (const auto& point : buckets[i][j])
            {
                meanPoint += point.z();
            }
            if(!buckets[i][j].empty())
            {
                meanPoint /= static_cast<const float&>(buckets[i][j].size());
                grid[i][j] = Eigen::Vector3f{i*stepLength, j*stepLength, meanPoint};
                fillMask[i][j] = true;
            }
        }
    }

    buckets.clear();

    for (int i = 0; i < grid.size(); ++i)
    {
        for (int j = 0; j < grid[0].size(); ++j)
        {
            if (fillMask[i][j])
            {
                continue;
            }
            
            int numPoints{};
            
            for (int xn = i-1; xn <= i+1; ++xn)
            {
                if (xn < 0 || xn >= numStepsX) continue; 
                for (int yn = j-1; yn <= j+1; ++yn)
                {
                    if (yn < 0 || yn >= numStepsY || !fillMask[xn][yn]) continue;

                    grid[i][j][2] += grid[xn][yn].z();
                    numPoints++;
                }
            }
            grid[i][j][2] /= numPoints > 0 ? numPoints: 1.f;
            fillMask[i][j] = true;
        }
        
    }

    writeVertexData("../ProcessedData/vertices.txt", grid, Eigen::Vector3f{
                        0.5f * (bounds.xExtent), 0.5f * (bounds.yExtent),
                        0.5f * (bounds.zExtent)
                    }, 0.5f, 0.5f, 0.5f);
    writeIndexFile("../ProcessedData/indices.txt", numStepsX, numStepsY);
}

void writeVertexData(const std::string& filePath, const std::vector<std::vector<Eigen::Vector3f>>& dataGrid,
                     const Eigen::Vector3f& offset, const float scaleX, const float scaleY, const float scaleZ)
{
    std::ofstream outFile{filePath};
    if (!outFile.is_open())
    {
        const auto msg = "Could not open out-file: " + filePath;
        std::cout << msg << std::endl;
        throw std::runtime_error(msg);
    }

    const size_t numVertices = dataGrid.size() * dataGrid[0].size();
    outFile << numVertices << "\n";

    for (size_t i = 0; i < dataGrid.size(); ++i)
    {
        for (size_t j = 0; j < dataGrid[0].size(); ++j)
        {
            auto point = dataGrid[i][j];
            outFile << "(" << point.x()*scaleX << ", " << point.z()*scaleZ << ", " << point.y()*scaleY << ")\n";
        }
    }

    outFile.close();
}

void writeIndexFile(const std::string& filePath, int numX, int numY)
{
    std::ofstream outFile{filePath};
    if (!outFile.is_open())
    {
        const auto msg = "Could not open out-file: " + filePath;
        std::cout << msg << std::endl;
        throw std::runtime_error(msg);
    }

    std::vector<int> indices{};

    for(unsigned int i = 0; i < numX-1; i++)       // for each row a.k.a. each strip
    {
        for (unsigned int j = 0; j < numY-1; j++) // for each column
        {
            indices.emplace_back(j+i*numY);
            indices.emplace_back((j+1)+i*numY);
            indices.emplace_back(j+(i+1)*numY);

            indices.emplace_back((j+1)+i*numY);
            indices.emplace_back((j+1)+(i+1)*numY);
            indices.emplace_back(j+(i+1)*numY);
        }
    }

    outFile << 2*(numY-1)*(numX-1) << "\n";
    std::cout << indices.size() << " " << 2*(numY-1)*(numX-1) << std::endl;

    for (size_t i = 2; i < indices.size(); i+=3)
    {
        outFile << indices[i-2] << " " << indices[i-1] << " " << indices[i] << " " << -1 << " " << -1 << " " << -1 << "\n";
    }

    outFile.close();
}

DataBounds findBounds(std::ifstream& file)
{
    file.clear();
    file.seekg(0);

    float xmin, ymin, zmin, x{}, y{}, z{};

    file >> x >> y >> z;
    float xmax = xmin = x;
    float ymax = ymin = y;
    float zmax = zmin = z;

    while (file >> x >> y >> z)
    {
        if (x < xmin)
        {
            xmin = x;
        }
        else if (x > xmax)
        {
            xmax = x;
        }
        if (y < ymin)
        {
            ymin = y;
        }
        else if (y > ymax)
        {
            ymax = y;
        }
        if (z < zmin)
        {
            zmin = z;
        }
        else if (z > zmax)
        {
            zmax = z;
        }
    }

    return DataBounds{xmin, xmax, ymin, ymax, zmin, zmax, xmax - xmin, ymax - ymin, zmax - zmin};
}

int main(int argc, char* argv[])
{
    processData("../RawData/merged.txt");
    return 0;
}
