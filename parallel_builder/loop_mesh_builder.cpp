/**
 * @file    loop_mesh_builder.cpp
 *
 * @author  ONDŘEJ KREJČÍ <xkrejc69@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP loops
 *
 * @date    17. 12. 2021
 **/

#include <iostream>
#include <math.h>
#include <limits>
#include <omp.h>

#include "loop_mesh_builder.h"

LoopMeshBuilder::LoopMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "OpenMP Loop")
{}

unsigned LoopMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    // Loop over each coordinate in the 3D grid.
    #pragma omp parallel default(none) shared(field)
    {
        #pragma omp for schedule(static) collapse(3)
        for(size_t x = 0; x < mGridSize; x++){
            for(size_t y = 0; y < mGridSize; y++){
                for(size_t z = 0; z < mGridSize; z++){
                    buildCube(Vec3_t<float>{x, y, z}, field);
                }
            }
        }
    }

    for(const auto& vec : mTriangleVectors){
        mTriangles.insert(mTriangles.end(), vec.begin(), vec.end());
    }

    // Return total number of triangles generated.
    return mTriangles.size();
}

// NOTE: This method is called from "buildCube(...)"!
float LoopMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    // 1. Store pointer to and number of 3D points in the field
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    // 2. Find minimum square distance from points "pos" to any point in the
    //    field.
    //#pragma omp parallel for default(none) shared(value, pPoints, pos) reduction(min:value) schedule(static)
    #pragma omp simd reduction(min:value)
    for(unsigned i = 0; i < count; i++)
    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        // Comparing squares instead of real distance to avoid unnecessary "sqrt"s in the loop.
        value = std::min(value, distanceSquared);
    }

    // 3. Finally take square root of the minimal square distance to get the real distance
    return sqrt(value);
}

// NOTE: This method is called from "buildCube(...)"!
void LoopMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    // Store generated triangle into vector (array) of generated triangles.
    mTriangleVectors[omp_get_thread_num()].push_back(triangle);
}
