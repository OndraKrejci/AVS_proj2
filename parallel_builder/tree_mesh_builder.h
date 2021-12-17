/**
 * @file    tree_mesh_builder.h
 *
 * @author  ONDŘEJ KREJČÍ <xkrejc69@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    17. 12. 2021
 **/

#ifndef TREE_MESH_BUILDER_H
#define TREE_MESH_BUILDER_H

#include <math.h>
#include <omp.h>

#include "base_mesh_builder.h"

class TreeMeshBuilder : public BaseMeshBuilder
{
public:
    TreeMeshBuilder(unsigned gridEdgeSize);
    ~TreeMeshBuilder();

protected:
    unsigned marchCubes(const ParametricScalarField &field);
    float evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field);
    void emitTriangle(const Triangle_t &triangle);
    const Triangle_t *getTrianglesArray() const { return mTriangles.data(); }
    
    unsigned octree(const ParametricScalarField& field, const Vec3_t<float>& start, unsigned len);
    bool isEmpty(const ParametricScalarField& field, const Vec3_t<float>& start, unsigned len);
    void fieldToArrays(const ParametricScalarField& field);

    const int threads = omp_get_max_threads();

    std::vector<Triangle_t>* mTriangleVectors = nullptr;
    std::vector<Triangle_t> mTriangles{}; ///< Temporary array of triangles

    const float HALF_SQRT3 = sqrt(3.0f) / 2.0f;

    unsigned fieldSize = 0;
    float* x = nullptr;
    float* y = nullptr;
    float* z = nullptr;
};

#endif // TREE_MESH_BUILDER_H
