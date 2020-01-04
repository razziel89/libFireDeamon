import sys
from .cpp import *

def add_mesh(data):
    if len(data[0]) == 3:
        mesh = Trimesh(len(data),True)
        for c,p,n in data:
            mesh.AddPoint(p[0],p[1],p[2],c[0],c[1],c[2],n[0],n[1],n[2])
    elif len(data[0]) == 2:
        mesh = Trimesh(len(data),False)
        for c,p in data:
            mesh.AddPoint(p[0],p[1],p[2],c[0],c[1],c[2])
    else:
        raise Exception("Unhandled internal exception.")
    return mesh

def visualize_mesh(mesh,scale):
    mesh.Visualize(scale)
