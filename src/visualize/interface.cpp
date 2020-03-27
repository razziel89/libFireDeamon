#include <FireDeamon/visualize/interface.h>
#include <GL/gl.h>
#include <iostream>

void Trimesh::Visualize(float scale) {
  glBegin(GL_TRIANGLES);
  unsigned int count = 0;
  for (float *it = _mesh; count < _length; it += _spread, ++count) {
    glColor3f(*(it + 0), *(it + 1), *(it + 2));
    if (_have_normals)
      glNormal3f(*(it + 6), *(it + 7), *(it + 8));
    glVertex3f(scale * *(it + 3), scale * *(it + 4), scale * *(it + 5));
  }
  glEnd();
}

void Trimesh::AddPoint(float px, float py, float pz, float cr, float cg, float cb,
                       float nx, float ny, float nz) {
  if (!_have_normals) {
    error = true;
    return;
  }
  *(_mesh_end + 0) = cr;
  *(_mesh_end + 1) = cg;
  *(_mesh_end + 2) = cb;
  *(_mesh_end + 3) = px;
  *(_mesh_end + 4) = py;
  *(_mesh_end + 5) = pz;
  *(_mesh_end + 6) = nx;
  *(_mesh_end + 7) = ny;
  *(_mesh_end + 8) = nz;

  _mesh_end += _spread;
}
void Trimesh::AddPoint(float px, float py, float pz, float cr, float cg, float cb) {
  if (_have_normals) {
    error = true;
    return;
  }
  *(_mesh_end + 0) = cr;
  *(_mesh_end + 1) = cg;
  *(_mesh_end + 2) = cb;
  *(_mesh_end + 3) = px;
  *(_mesh_end + 4) = py;
  *(_mesh_end + 5) = pz;

  _mesh_end += _spread;
}
