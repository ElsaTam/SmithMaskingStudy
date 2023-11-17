#pragma once

#include <set>
#include "gdt/math/vec.h"

template<int dim>
class KMeans
{
private:

    struct Point
    {
        int pointId;
        int clusterId = 0; // Initially not assigned to any cluster
        gdt::vec_t<scal, dim> values;

        Point() : pointId(0), values(0) { }
        Point(int id, gdt::vec_t<scal, dim> V) : pointId(id), values(V) { }
    };

    struct Cluster
    {
        int clusterId;
        std::vector<scal> centroid;
        std::vector<Point> points;

        Cluster(int clusterId, Point centroid);

        void addPoint(Point p);

        bool removePoint(int pointId);

        void removeAllPoints();
    };

    int K; // number of classes (user choice)
    int iters; // number of iterations (user choice)
    std::vector<Cluster> clusters;

    void clearClusters();

    int getNearestClusterId(Point point);

public:
    KMeans(int K, int iterations);

    void run(std::vector<gdt::vec_t<scal, dim>>& all_points, std::vector<int> centersID = {});

    std::set<int> getPointsIdOfCluster(int k) const;
    size_t getSizeOfCluster(int k) const;
    gdt::vec_t<scal, dim> getCentroidOfCluster(int k) const;
};