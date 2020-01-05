"""Visualization sub-module for the FireDeamon module.

This module uses very simple but fast wrappers for OpenGL to faster draw surfaces
consisting of triangles. A mesh in this sense is defined as a set of points, each of
which has 3 co-ordinates, 3 colour values (RGB) and an associated surface normal.
"""
import sys

try:
    from .cpp import *
except ImportError:
    raise ImportError("Cannot import C++ extension")


def add_mesh(data):
    """Convert a mesh for fast visualization.

    A mesh in this sense is defined as a set of points, each of which has 3
    co-ordinates, 3 colour values (RGB) and an associated surface normal. Surface
    normals are optional but they are required to use nice lighting.

    Args:
        data: (a mesh) - a list consisting of points, i.e. 3 3-element lists/tuples,
            defining the point's co-ordinates, colour values and suface normal,
            respectively. Surface normals are optional.

    Returns:
        A Python-proxy for a mesh object that can be visualized quickly using the
        function visualize_mesh
    """
    if len(data[0]) == 3:
        mesh = Trimesh(len(data), True)
        for c, p, n in data:
            mesh.AddPoint(p[0], p[1], p[2], c[0], c[1], c[2], n[0], n[1], n[2])
    elif len(data[0]) == 2:
        mesh = Trimesh(len(data), False)
        for c, p in data:
            mesh.AddPoint(p[0], p[1], p[2], c[0], c[1], c[2])
    else:
        raise ValueError(
            "Cannot convert data to Trimesh, wrong length per element, "
            "supported: 2 and 3, found: {}".format(len(data[0]))
        )
    return mesh


def visualize_mesh(mesh, scale):
    """Quickly visualize a previously converted mesh via OpenGL

    Args:
        mesh: (Trimesh object) - the mesh to be visualized
        scale: (float) - the mesh will be scaled by this factor prior to visualization
    """
    mesh.Visualize(scale)
