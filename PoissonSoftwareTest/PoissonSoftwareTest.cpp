#include <fstream>
#include "geometry.h"
#include <iostream>

int main(int argc, char* argv[])
{
    // Checking for command line arguments
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    std::string filename = argv[1];
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }
    double x, y;
    std::vector<geom::Point> nodes;
    // Assuming the data is correct
    // Reading coords from a file
    while (file >> x >> y) {
        nodes.emplace_back(geom::Point{ x, y });
    }
    file.close();
    try{
        // Initialize the polygon and find its symmetry axes
        geom::ConvexPolygon poly;
        poly.SetNodes(nodes);
        auto axes = geom::FindSymmetricAxes(poly);
        if (axes.size()) {
            for (const auto& axis : axes) {
                std::cout << axis.first << " - " << axis.second << '\n';
            }
        }
        else {
            std::cout << "non-symmetric\n";
        }
    }
    catch (const geom::ConvexPolygon::IncorrectPolygonException) {
        std::cerr << "Error: Incorrectly set convex polygon\n";
    }
    return 0;
}
