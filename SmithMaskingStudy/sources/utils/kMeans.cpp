#include <omp.h>
#include <vector>

#include "utils/console.h"
#include "utils/kMeans.h"

template<int dim>
KMeans<dim>::Cluster::Cluster(int clusterId, Point centroid)
{
    this->clusterId = clusterId;
    for (int i = 0; i < dim; i++)
    {
        this->centroid.push_back(centroid.values[i]);
    }
    this->addPoint(centroid);
}

template<int dim>
void KMeans<dim>::Cluster::addPoint(Point p)
{
    p.clusterId = this->clusterId;
    points.push_back(p);
}

template<int dim>
bool KMeans<dim>::Cluster::removePoint(int pointId)
{
    size_t size = points.size();

    for (int i = 0; i < size; i++)
    {
        if (points[i].pointId == pointId)
        {
            points.erase(points.begin() + i);
            return true;
        }
    }
    return false;
}

template<int dim>
void KMeans<dim>::Cluster::removeAllPoints() { points.clear(); }








template<int dim>
KMeans<dim>::KMeans(int K, int iterations) : K(K), iters(iterations)
{ }

template<int dim>
void KMeans<dim>::clearClusters()
{
    for (int i = 0; i < K; i++)
    {
        clusters[i].removeAllPoints();
    }
}

template<int dim>
int KMeans<dim>::getNearestClusterId(Point point)
{
    double sum = 0.0, min_dist;
    int NearestClusterId;
    if (dim == 1) {
        min_dist = abs(clusters[0].centroid[0] - point.values[0]);
    }
    else
    {
        for (int i = 0; i < dim; i++)
        {
            sum += pow(clusters[0].centroid[i] - point.values[i], 2.0);
            // sum += abs(clusters[0].getCentroidByPos(i) - point.getVal(i));
        }
        min_dist = sqrt(sum);
    }
    NearestClusterId = clusters[0].clusterId;

    for (int i = 1; i < K; i++)
    {
        double dist;
        sum = 0.0;

        if (dim == 1) {
            dist = abs(clusters[i].centroid[0] - point.values[0]);
        }
        else {
            for (int j = 0; j < dim; j++)
            {
                sum += pow(clusters[i].centroid[j] - point.values[j], 2.0);
                // sum += abs(clusters[i].getCentroidByPos(j) - point.getVal(j));
            }

            dist = sqrt(sum);
            // dist = sum;
        }
        if (dist < min_dist)
        {
            min_dist = dist;
            NearestClusterId = clusters[i].clusterId;
        }
    }

    return NearestClusterId;
}

template<int dim>
void KMeans<dim>::run(std::vector<gdt::vec_t<scal, dim>>& all_points, std::vector<int> centersID)
{
    size_t N = all_points.size();

    // Create the points
    std::vector<Point> points;
    for (int i = 0; i < N; ++i) {
        Point point;
        point.pointId = i;
        point.values = all_points[i];
        points.push_back(point);
    }

    // Initializing Clusters
    if (centersID.size() > 0) {
        for (int k = 1; k <= K && k <= centersID.size(); ++k) {
            int index = centersID[k-1]; // cluster index is id-1
            points[index].clusterId = k;
            Cluster cluster(k, points[index]);
            clusters.push_back(cluster);
        }
    }
    else {
        std::vector<int> used_pointIds;
        for (int k = 1; k <= K; ++k)
        {
            while (true)
            {
                int index = rand() % N; // randomly choose one point to initialize each cluster

                if (find(used_pointIds.begin(), used_pointIds.end(), index) == used_pointIds.end())
                {
                    used_pointIds.push_back(index);
                    points[index].clusterId = k;
                    Cluster cluster(k, points[index]);
                    clusters.push_back(cluster);
                    break;
                }
            }
        }
    }

    Console::out << Console::timePad << "Running K-Means Clustering.." << std::endl;

    int iter = 1;
    while (true)
    {
        bool done = true;

        // Add all points to their nearest cluster
#pragma omp parallel for reduction(&&: done)
        for (int i = 0; i < N; i++)
        {
            int currentClusterId = points[i].clusterId;
            int nearestClusterId = getNearestClusterId(points[i]);

            if (currentClusterId != nearestClusterId)
            {
                points[i].clusterId = nearestClusterId;
                done = false;
            }
        }

        // clear all existing clusters
        clearClusters();

        // reassign points to their new clusters
        for (int i = 0; i < N; i++)
        {
            // cluster index is ID-1
            clusters[points[i].clusterId - 1].addPoint(points[i]);
        }

        // Recalculating the center of each cluster
        for (int i = 0; i < K; i++)
        {
            size_t ClusterSize = clusters[i].points.size();

            for (int j = 0; j < dim; j++)
            {
                double sum = 0.0;
                if (ClusterSize > 0)
                {
#pragma omp parallel for reduction(+: sum)
                    for (int p = 0; p < ClusterSize; p++)
                    {
                        sum += clusters[i].points[p].values[j];
                    }
                    clusters[i].centroid[j] = sum / (float) ClusterSize;
                }
            }
        }

        if (done || iter >= iters)
        {
            Console::light << Console::timePad << "Clustering completed in iteration : " << iter << std::endl;
            break;
        }
        iter++;
    }
}


template<int dim>
std::set<int> KMeans<dim>::getPointsIdOfCluster(int k) const
{
    std::set<int> points;
    for (Point p : clusters[k].points) {
        points.insert(p.pointId);
    }
    return points;
}

template<int dim>
size_t KMeans<dim>::getSizeOfCluster(int k) const
{
    return clusters[k].points.size();
}

template<int dim>
gdt::vec_t<scal, dim> KMeans<dim>::getCentroidOfCluster(int k) const
{
    gdt::vec_t<scal, dim> centroid;
    for (int i = 0; i < dim; ++i) {
        centroid[i] = clusters[k].centroid[i];
    }
    return centroid;
}