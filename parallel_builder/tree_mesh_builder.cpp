/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  ONDŘEJ KREJČÍ <xkrejc69@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    13. 12. 2021
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"

TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{}

unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    // Suggested approach to tackle this problem is to add new method to
    // this class. This method will call itself to process the children.
    // It is also strongly suggested to first implement Octree as sequential
    // code and only when that works add OpenMP tasks to achieve parallelism.

    unsigned trianglesCount = 0;
    
    trianglesCount = octree(field, Vec3_t<float>{0, 0, 0}, mGridSize);

    return trianglesCount;
}

float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    // NOTE: This method is called from "buildCube(...)"!

    // 1. Store pointer to and number of 3D points in the field
    //    (to avoid "data()" and "size()" call in the loop).
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    // 2. Find minimum square distance from points "pos" to any point in the
    //    field.
    for(unsigned i = 0; i < count; i++)
    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        // Comparing squares instead of real distance to avoid unnecessary
        // "sqrt"s in the loop.
        value = std::min(value, distanceSquared);
    }

    // 3. Finally take square root of the minimal square distance to get the real distance
    return sqrt(value);
}

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    // NOTE: This method is called from "buildCube(...)"!

    // Store generated triangle into vector (array) of generated triangles.
    // The pointer to data in this array is return by "getTrianglesArray(...)" call
    // after "marchCubes(...)" call ends.
    mTriangles.push_back(triangle);
}

unsigned TreeMeshBuilder::octree(const ParametricScalarField &field, const Vec3_t<float>& start, unsigned len){
    if(len == 1){
        return buildCube(start, field);
    }

    unsigned totalTriangles = 0;

    const unsigned half = len / 2;

    const std::vector<Vec3_t<float>> positions{
        {start.x, start.y, start.z},
        {start.x + half, start.y, start.z},
        {start.x, start.y + half, start.z},
        {start.x + half, start.y + half, start.z},
        {start.x, start.y, start.z + half},
        {start.x + half, start.y, start.z + half},
        {start.x, start.y + half, start.z + half},
        {start.x + half, start.y + half, start.z + half}
    };

    for(unsigned i = 0; i < positions.size(); i++){
        totalTriangles += octree(field, positions[i], half);
    }

    return totalTriangles;
}

bool TreeMeshBuilder::isEmpty(const ParametricScalarField& field, const Vec3_t<float>& start, unsigned len){
    return false;
}