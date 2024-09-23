#pragma once

#include <iostream>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <algorithm>


namespace geom {
    // 2D Point with doubles coords
    struct Point { 
        double x;
        double y;
        Point() noexcept: x(0), y(0) {}
        Point(double x, double y) noexcept: x(x), y(y) {}
        Point(int x, int y) noexcept : x(static_cast<double>(x)), y(static_cast<double>(y)) {}
        bool operator== (const Point& other) const noexcept;
    };
    bool operator< (const Point& p1, const Point& p2) noexcept;
    std::ostream& operator<<(std::ostream& os, const Point& p) noexcept;

    // Struct for hashing Point type obj.
    struct HashPoint { 
        size_t operator()(const Point& p) const noexcept;
    };

    // A class to represent polygon
    class ConvexPolygon { 
    private:
        std::vector<Point> nodes;
        // Hash set for unique nodes; Used to find items nodes
        std::unordered_set<Point, HashPoint> hash_nodes;

        //Constructing a minimal convex hull around the points of the polygon
        // This is how we get rid of unwanted dots on polyhon sides
        void Normalize();
    public:
        class IncorrectPolygonException {};
        ConvexPolygon() {}
        ConvexPolygon(const std::vector<Point>& points);
        ConvexPolygon(const std::vector<std::pair<double, double>>& points);
        const std::vector<Point>& GetNodes() const;
        const std::unordered_set<Point, HashPoint>& GetHashNodes() const;
        void SetNodes(const std::vector<Point>& nods);
    };

    // Finding the bottom left point of the polygon (For Graha'm algorithm)
    Point FindLowestPoint(const std::vector<Point>& points);

    // Finding the polar angle between two points
    double PolarAngle(const Point& origin, const Point& p);

    // Checking for convexity of a polygon (Graham's algorithm)
    bool IsConvex(const std::vector<Point>& points);

    // To find orientation of 3 ordered points
    double CrossProduct(const Point& O, const Point& A, const Point& B);

    // A convex hull on the given points set
    // Allows to get rid of unnecessary points of the polygon
    std::vector<Point> ConvexHull(const std::vector<Point>& points);

    // Calculating the remainder of division with correction for positive numbers
    constexpr size_t Pmod(size_t x, size_t y) noexcept;

    // Checking orthogonality of two vectors
    bool AreOrthoggonal(const std::pair<Point, Point>& v1, const std::pair<Point, Point>& v2) noexcept;

    // Calculates point symmetric to the point with respect to the axis
    Point ReflectPoint(const std::pair<Point, Point>& axis, const Point& p) noexcept;

    // Calculates the midpoint of a segment
    Point Midpoint(const Point& p1, const Point& p2) noexcept;

    // Checks if the given axis is the axis of symmetry for the polygon 
    bool IsSymmetricAxis(const ConvexPolygon& poly, const std::pair<Point, Point>& axis, size_t start, size_t end);

    // Finds the midpoints of all segments that the polygon includes
    std::vector<Point> FindMidpoints(const ConvexPolygon& poly);

    // Finds all symmetry axes of the polygon
    std::vector<std::pair<Point, Point>> FindSymmetricAxes(const ConvexPolygon& poly);

}
