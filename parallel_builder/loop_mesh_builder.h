/**
 * @file    loop_mesh_builder.h
 *
 * @author  ONDŘEJ KREJČÍ <xkrejc69@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP loops
 *
 * @date    17. 12. 2021
 **/

#ifndef LOOP_MESH_BUILDER_H
#define LOOP_MESH_BUILDER_H

#include <vector>

#include "base_mesh_builder.h"

class LoopMeshBuilder : public BaseMeshBuilder
{
public:
    LoopMeshBuilder(unsigned gridEdgeSize);
    ~LoopMeshBuilder();

protected:
    unsigned marchCubes(const ParametricScalarField& field);
    float evaluateFieldAt(const Vec3_t<float>& pos, const ParametricScalarField& field);
    void emitTriangle(const Triangle_t& triangle);
    const Triangle_t* getTrianglesArray() const { return mTriangles.data(); }

    void fieldToArrays(const ParametricScalarField& field);

    const int threads = omp_get_max_threads();

    std::vector<Triangle_t>* mTriangleVectors = nullptr;
    std::vector<Triangle_t> mTriangles{}; ///< Temporary array of triangles

    unsigned fieldSize = 0;
    float* x = nullptr;
    float* y = nullptr;
    float* z = nullptr;
};

#endif // LOOP_MESH_BUILDER_H
