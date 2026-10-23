#include <Eigen/Dense>
#include <chrono>
#include <iostream>

/*
mult           0.0163621 s
solveInverse   0.0201562 s
solveFullPivLu 0.202951 s
solveColPivHh  0.547893 s
solveFullPivHh 0.621699 s
solveComplete  0.844406 s

*/

Eigen::Vector3d mult(Eigen::Matrix3d const &m, Eigen::Vector3d const &v) {
  return m * v;
}
 
Eigen::Vector3d solveInverse(Eigen::Matrix3d const &m, Eigen::Vector3d const &v) {
  return m.inverse() * v;
}
 
Eigen::Vector3d solveFullPivLu(Eigen::Matrix3d const &m, Eigen::Vector3d const &v) {
  return m.fullPivLu().solve(v);
} 

Eigen::Vector3d solveColPivHh(Eigen::Matrix3d const &m, Eigen::Vector3d const &v) {
  return m.colPivHouseholderQr().solve(v);
} 

Eigen::Vector3d solveFullPivHh(Eigen::Matrix3d const &m, Eigen::Vector3d const &v) {
  return m.fullPivHouseholderQr().solve(v);
} 

Eigen::Vector3d solveComplete(Eigen::Matrix3d const &m, Eigen::Vector3d const &v) {
  return m.completeOrthogonalDecomposition().solve(v);
} 

Eigen::Vector3d solveBdcsvd(Eigen::Matrix3d const &m, Eigen::Vector3d const &v) {
  return m.bdcSvd().solve(v);
} 

Eigen::Vector3d solveJacobiSvd(Eigen::Matrix3d const &m, Eigen::Vector3d const &v) {
  return m.jacobiSvd().solve(v);
} 

typedef Eigen::Vector3d func(Eigen::Matrix3d const &m, Eigen::Vector3d const &v);

int measure(char const *t, func &f) {
  Eigen::Matrix3d m;
  Eigen::Vector3d v,s;
  s << 0.0, 0.0, 0.0;
  int n = 0;
  auto start = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < 1000000; ++i) {
    m << ++n, ++n, ++n, ++n, ++n, ++n, ++n, ++n, ++n;
    v << ++n, ++n, ++n;
    n %= 71;
    s += f(m, v);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << t << ' ' << diff.count() << " s\n";
  return s(0);
}

int main() {
  return measure("mult", mult) +
  measure("solveInverse", solveInverse) +
  measure("solveFullPivLu", solveFullPivLu) +
  measure("solveColPivHh", solveColPivHh) +
  measure("solveFullPivHh", solveFullPivHh) +
  measure("solveComplete", solveComplete)/* +
  measure("solveBdcSvd", solveBdcsvd) +
  measure("solveJacobiSvd", solveJacobiSvd)*/;
}
