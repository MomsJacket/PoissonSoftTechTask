#include "geometry.h"
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stack>


namespace geom{
    constexpr auto eps = 1e-6;
    bool Point::operator== (const Point& other) const noexcept {
        return std::fabs(x - other.x) < eps && std::fabs(y - other.y) < eps;
    }
    bool operator< (const Point& p1, const Point& p2) noexcept {
        return std::tie(p1.x, p1.y) < std::tie(p2.x, p2.y);
    }
    std::ostream& operator<<(std::ostream& os, const Point& p) noexcept {
        os << "(" << p.x << ", " << p.y << ")";
        return os;
    }
    size_t HashPoint::operator()(const Point& p) const noexcept {
        size_t h1 = std::hash<double>{}(std::round(p.x * 1000) / 1000);
        size_t h2 = std::hash<double>{}(std::round(p.y * 1000) / 1000);
        return h1 ^ (h2 << 1);
    }
    ConvexPolygon::ConvexPolygon(const std::vector<Point>& points) {
        SetNodes(points);
    }
    ConvexPolygon::ConvexPolygon(const std::vector<std::pair<double, double>>& points) {
        std::vector<Point> convertedPoints;
        convertedPoints.reserve(points.size());
        std::transform(points.begin(), points.end(), std::back_inserter(convertedPoints),
            [](const std::pair<double, double>& p) {
                return Point(p.first, p.second); 
            });
        SetNodes(convertedPoints);
    }
    const std::vector<Point>& ConvexPolygon::GetNodes() const {
        return nodes;
    }
    const std::unordered_set<Point, HashPoint>& ConvexPolygon::GetHashNodes() const {
        return hash_nodes;
    }
    void ConvexPolygon::SetNodes(const std::vector<Point>& nods) {
        if (!IsConvex(nods)) {
            throw IncorrectPolygonException();
        }
        nodes = nods;
        hash_nodes = std::unordered_set<Point, HashPoint>(nodes.begin(), nodes.end());
        Normalize();
    }

    double CrossProduct(const Point& O, const Point& A, const Point& B) {
        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
    }

    bool IsConvex(const std::vector<Point>& points) {
        size_t n = points.size();
        if (n < 3) return false;

        bool gotPositive = false;
        bool gotNegative = false;

        for (size_t i = 0; i < n; ++i) {
            double crossProduct = CrossProduct(points[i], points[(i + 1) % n], points[(i + 2) % n]);

            if (crossProduct > 0) gotPositive = true;
            else if (crossProduct < 0) gotNegative = true;

            if (gotPositive && gotNegative) return false;
        }

        return true;
    }

    Point FindLowestPoint(const std::vector<Point>& points) {
        Point lowest = points[0];
        for (const auto& point : points) {
            if (point.y < lowest.y || (point.y == lowest.y && point.x < lowest.x)) {
                lowest = point;
            }
        }
        return lowest;
    }

    double PolarAngle(const Point& origin, const Point& p) {
        return atan2(p.y - origin.y, p.x - origin.x);
    }

    std::vector<Point> ConvexHull(std::vector<Point>& points) {
        if (points.size() < 3) return points;

        // Find lowest point
        Point lowest = FindLowestPoint(points);

        // Sorting points by polar angle with respect of lowest
        std::sort(points.begin(), points.end(), [&lowest](const Point& a, const Point& b) {
            double angleA = PolarAngle(lowest, a);
            double angleB = PolarAngle(lowest, b);
            return angleA < angleB || (angleA == angleB && a.x < b.x);
            });

        // Using stack to store convex hull points
        std::stack<Point> hull;
        hull.push(points[0]);
        hull.push(points[1]);

        for (size_t i = 2; i < points.size(); ++i) {
            while (hull.size() > 1) {
                Point second = hull.top(); hull.pop();
                Point first = hull.top();
                // Counter-clockwise
                if (CrossProduct(first, second, points[i]) > 0) {  
                    hull.push(second);
                    break;
                }
            }
            hull.push(points[i]);
        }

        std::vector<Point> convexHull;
        while (!hull.empty()) {
            convexHull.push_back(hull.top());
            hull.pop();
        }

        return convexHull;
    }

    void ConvexPolygon::Normalize() {
        nodes = ConvexHull(nodes);
        hash_nodes.clear();
        hash_nodes.insert(nodes.begin(), nodes.end());
    }

    constexpr size_t Pmod(size_t x, size_t y) noexcept {
        return (x % y + y) % y;
    }

    bool AreOrthoggonal(const std::pair<Point, Point>& v1, const std::pair<Point, Point>& v2) noexcept {
        // Scalar product of two vectors
        return std::fabs((v1.second.x - v1.first.x) * (v2.second.x - v2.first.x)
                       + (v1.second.y - v1.first.y) * (v2.second.y - v2.first.y)) < eps;
    }

    Point ReflectPoint(const std::pair<Point, Point>& axis, const Point& p) noexcept {
        // The points a and b respectively define the axis
        // vab - Vector co-directed with the axis
        const auto& vab = std::make_pair(axis.second.x - axis.first.x, axis.second.y - axis.first.y);
        // vap - Vector from point a to point p
        const auto& vap = std::make_pair(p.x - axis.first.x, p.y - axis.first.y);
        auto norm = 1.0 / (vab.first * vab.first + vab.second * vab.second);
        auto scalar = vab.first * vap.first + vab.second * vap.second;
        // vp_ort - Projection vector of the point p on the axis 
        const auto& vp_ort = std::make_pair(norm * scalar * vab.first, norm * scalar * vab.second);
        return Point{ p.x + 2.0 * (vp_ort.first - vap.first), p.y + 2.0 * (vp_ort.second - vap.second) };
    }

    Point Midpoint(const Point& p1, const Point& p2) noexcept {
        return Point{ (p1.x + p2.x) / 2, (p1.y + p2.y) / 2 };
    }

    bool IsSymmetricAxis(const ConvexPolygon& poly, const std::pair<Point, Point>& axis, size_t start, size_t end) {
        const auto& nds = poly.GetNodes();
        const auto& point_set = poly.GetHashNodes();
        auto n = nds.size();
        for (auto i = start; Pmod(i, n) != end; ++i) {
            const auto& ref_point = ReflectPoint(axis, nds[Pmod(i, n)]);
            if (point_set.find(ref_point) == point_set.end()) return false;
        }
        return true;
    }

    std::vector<Point> FindMidpoints(const ConvexPolygon& poly) {
        std::vector<Point> midpoints;
        const auto& nodes = poly.GetNodes();
        auto n = nodes.size();
        if (n <= 2) return midpoints;
        else {
            for (size_t i = 0; i < n - 1; ++i) {
                midpoints.push_back(Midpoint(nodes[i], nodes[i + 1]));
            }
            midpoints.push_back(Midpoint(nodes[n - 1], nodes[0]));
        }
        return midpoints;
    }

    std::vector<std::pair<Point, Point>> FindSymmetricAxes(const ConvexPolygon& poly) {
        std::vector<std::pair<Point, Point>> axes;
        const auto& midpoints = FindMidpoints(poly);
        const auto& nodes = poly.GetNodes();
        const auto n = nodes.size();
        if (n < 3) return axes;
        for (size_t i = 0; i < n - 1; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                if ((j - i > 1) && (j - i < n - 1)) {
                    // Checking axes passing through two vertices of the polygon
                    const auto& axis = std::make_pair(nodes[i], nodes[j]);
                    if (IsSymmetricAxis(poly, axis, i, j)) axes.push_back(axis);
                }
                const auto& axis = std::make_pair(midpoints[i], midpoints[j]);
                const auto& midpointside1 = std::make_pair(nodes[Pmod(i, n)],
                    nodes[Pmod(i + 1, n)]);
                const auto& midpointside2 = std::make_pair(nodes[Pmod(j, n)],
                    nodes[Pmod(j + 1, n)]);
                // Checking axes passing through two midpoints of the sides
                if (AreOrthoggonal(axis, midpointside1) &&
                    AreOrthoggonal(axis, midpointside2) &&
                    IsSymmetricAxis(poly, axis, i, j)) axes.push_back(axis);
            }
        }
        for (size_t i = 0; i < n; ++i) {
            // Not checking adjacent sides
            auto f_index = Pmod(i + 1, n);
            for (auto j = f_index; j < f_index + n - 2; ++j) {
                auto modj = Pmod(j, n);
                const auto& axis = std::make_pair(nodes[i], midpoints[modj]);
                const auto& midpointside = std::make_pair(nodes[modj],
                                                          nodes[Pmod(j + 1, n)]);
                // Checking axes passing through the vertex and the midpoint of the opposite side
                if (AreOrthoggonal(midpointside, axis) &&
                    IsSymmetricAxis(poly, axis, i, modj)) axes.push_back(axis);
            }
        }
        return axes;
    }
}


