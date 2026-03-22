#include "read_file.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

SolverInput read_file(const char *filename) {

  std::fstream myfile;

  myfile.open(filename);
  std::vector<std::string> g1;
  Nodes nodes;
  Elements elems;
  std::vector<size_t> u_indices;
  std::vector<double> u;
  std::vector<double> F;
  std::unordered_map<int, int> node_id_to_index;

  if (myfile.is_open()) {
    std::string str;
    bool read = false;
    while (std::getline(myfile, str)) {
      if (str.rfind("*NODE", 0) == 0) {
        read = true;
        std::getline(myfile, str);
      }
      if (str.rfind("*", 0) == 0) {
        read = false;
      }
      if (read == true) {
        std::istringstream iss(str);
        int node_id;
        double x, y, z;
        char sep = ',';
        iss >> node_id >> sep >> x >> sep >> y >> sep >> z;
        int idx = nodes.size();
        nodes.push_back(Node{x, y});
        node_id_to_index[node_id] = idx;
      }
    }
    myfile.clear();
    myfile.seekg(0);
    read = false;
    while (std::getline(myfile, str)) {
      if (str.rfind("*ELEMENT, type=CPS6", 0) == 0) {
        read = true;
        std::getline(myfile, str);
      }
      if (str.rfind("*", 0) == 0) {
        read = false;
      }
      if (read == true) {
        std::istringstream iss(str);
        int l, n1, n2, n3, n4, n5, n6;
        char sep = ',';
        iss >> l >> sep >> n1 >> sep >> n2 >> sep >> n3 >> sep >> n4 >> sep >>
            n5 >> sep >> n6;
        elems.push_back(Element{node_id_to_index[n1], node_id_to_index[n2],
                                node_id_to_index[n3], node_id_to_index[n4],
                                node_id_to_index[n5], node_id_to_index[n6]});
      }
    }
    myfile.clear();
    myfile.seekg(0);
    read = false;
    while (std::getline(myfile, str)) {
      if (str.rfind("*BOUNDARY", 0) == 0) {
        read = true;
        std::getline(myfile, str);
      }
      if (str.rfind("*", 0) == 0) {
        read = false;
      }
      if (read == true) {
        int node, coord;
        double disp;
        std::istringstream iss(str);
        char sep = ',';
        iss >> node >> sep >> coord >> sep >> disp;
        int idx = node_id_to_index[node];
        u_indices.push_back(idx * 2 + coord);
        u.push_back(disp);
      }
    }
    myfile.clear();
    myfile.seekg(0);
    read = false;
    F.insert(F.begin(), nodes.size() * 2, 0);
    while (std::getline(myfile, str)) {
      if (str.rfind("*CLOAD", 0) == 0) {
        read = true;
        std::getline(myfile, str);
      }
      if (str == "") {
        read = false;
      }
      if (read == true) {
        int node, coord;
        double force;
        std::istringstream iss(str);
        char sep = ',';
        iss >> node >> sep >> coord >> sep >> force;
        int idx = node_id_to_index[node];
        F[idx * 2 + coord] = force;
      }
    }
    myfile.close();
  }

  SolverInput input = {nodes, elems, u_indices, u, F};
  return input;
}