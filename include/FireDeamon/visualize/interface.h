#ifndef INTERFACE_H
#define INTERFACE_H

class Trimesh {
private:
  unsigned int _length;
  unsigned int _spread;
  bool _have_normals;
  float *_mesh;
  float *_mesh_end;

public:
  bool error;
  Trimesh(int length, bool have_normals) {
    _have_normals = have_normals;
    _spread = _have_normals ? 9 : 6;
    _length = (unsigned int)length;
    _mesh = new float[_length * _spread];
    _mesh_end = _mesh;
    error = false;
  }
  ~Trimesh() { delete[] _mesh; }
  void Visualize(float scale);
  void AddPoint(float px, float py, float pz, float cr, float cg, float cb, float nx,
                float ny, float nz);
  void AddPoint(float px, float py, float pz, float cr, float cg, float cb);
};

#endif
