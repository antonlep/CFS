#include "read_file.h"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

bool contains_letter(const std::string &s) {
  return std::any_of(s.begin(), s.end(),
                     [](unsigned char c) { return std::isalpha(c); });
}

SolverInput read_file(const char *filename) {

  std::fstream myfile;
  std::fstream myfile2;

  myfile.open(filename);
  myfile2.open(filename);
  std::vector<std::string> g1;
  Nodes nodes;
  Elements elems;
  std::vector<size_t> u_indices;
  std::vector<double> u;
  std::vector<double> F;
  std::vector<BoundaryEdge> boundary_edges;
  std::vector<TractionEdge> traction_edges;
  struct MaterialProperties material;

  std::unordered_map<int, size_t> node_id_to_index;

  if ((myfile.is_open()) && (myfile2.is_open())) {
    std::string str;
    std::string str2;
    bool read = false;
    bool read2 = false;
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
        size_t idx = nodes.size();
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
        elems.push_back(Element{static_cast<size_t>(node_id_to_index[n1]),
                                static_cast<size_t>(node_id_to_index[n2]),
                                static_cast<size_t>(node_id_to_index[n3]),
                                static_cast<size_t>(node_id_to_index[n4]),
                                static_cast<size_t>(node_id_to_index[n5]),
                                static_cast<size_t>(node_id_to_index[n6])});
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
        if (contains_letter(str)) {
          std::string nodeset, coord_str, disp_str;
          std::vector<int> nodeset_ids;
          std::istringstream iss(str);
          char sep = ',';
          std::getline(iss, nodeset, sep);
          std::getline(iss, coord_str, sep);
          std::getline(iss, disp_str, sep);
          size_t coord = static_cast<size_t>(std::stoi(coord_str));
          double disp = std::stod(disp_str);
          myfile2.clear();
          myfile2.seekg(0);
          read2 = false;
          while (std::getline(myfile2, str2)) {
            if (str2.find("NSET=" + nodeset) != std::string::npos) {
              read2 = true;
              std::getline(myfile2, str2);
            }
            if (str2.rfind("*", 0) == 0) {
              read2 = false;
            }
            if (read2 == true) {
              std::stringstream ss(str2);
              std::string item;
              while (std::getline(ss, item, ',')) {
                item.erase(0, item.find_first_not_of(" \t"));
                item.erase(item.find_last_not_of(" \t") + 1);
                if (!item.empty()) {
                  nodeset_ids.push_back(std::stoi(item));
                }
              }
            }
          }
          for (int node_id : nodeset_ids) {
            size_t idx = node_id_to_index[node_id];
            u_indices.push_back(idx * 2 + coord);
            u.push_back(disp);
          }
        } else {
          int node, coord_raw;
          double disp;
          std::istringstream iss(str);
          char sep = ',';
          iss >> node >> sep >> coord_raw >> sep >> disp;
          size_t idx = node_id_to_index[node];
          size_t coord = static_cast<size_t>(coord_raw);
          u_indices.push_back(idx * 2 + coord);
          u.push_back(disp);
        }
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
      if (str.rfind("*", 0) == 0) {
        read = false;
      }
      if (read == true) {
        if (contains_letter(str)) {
          std::string nodeset, coord_str, force_str;
          std::vector<int> nodeset_ids;
          std::istringstream iss(str);
          char sep = ',';
          std::getline(iss, nodeset, sep);
          std::getline(iss, coord_str, sep);
          std::getline(iss, force_str, sep);
          size_t coord = static_cast<size_t>(std::stoi(coord_str));
          double force = std::stod(force_str);
          myfile2.clear();
          myfile2.seekg(0);
          read2 = false;
          while (std::getline(myfile2, str2)) {
            if (str2.find("NSET=" + nodeset) != std::string::npos) {
              read2 = true;
              std::getline(myfile2, str2);
            }
            if (str2.rfind("*", 0) == 0) {
              read2 = false;
            }
            if (read2 == true) {
              std::stringstream ss(str2);
              std::string item;
              while (std::getline(ss, item, ',')) {
                item.erase(0, item.find_first_not_of(" \t"));
                item.erase(item.find_last_not_of(" \t") + 1);
                if (!item.empty()) {
                  nodeset_ids.push_back(std::stoi(item));
                }
              }
            }
          }
          for (int node_id : nodeset_ids) {
            size_t idx = node_id_to_index[node_id];
            F[idx * 2 + coord] = force;
          }
        } else {
          int node, coord_raw;
          double force;
          std::istringstream iss(str);
          char sep = ',';
          iss >> node >> sep >> coord_raw >> sep >> force;
          size_t idx = node_id_to_index[node];
          size_t coord = static_cast<size_t>(coord_raw);
          F[idx * 2 + coord] = force;
        }
      }
    }
    myfile.clear();
    myfile.seekg(0);
    read = false;
    while (std::getline(myfile, str)) {
      if (str.rfind("*SFILM", 0) == 0) {
        read = true;
        std::getline(myfile, str);
      }
      if (str.rfind("*", 0) == 0) {
        read = false;
      }
      if (read == true) {
        int n1, n2, n3;
        double h;
        double Tinf;
        std::istringstream iss(str);
        char sep = ',';
        iss >> n1 >> sep >> n2 >> sep >> n3 >> sep >> h >> sep >> Tinf;
        size_t idx1 = node_id_to_index[n1];
        size_t idx2 = node_id_to_index[n2];
        size_t idx3 = node_id_to_index[n3];
        BoundaryEdge boundary1;
        boundary1.n1 = idx1;
        boundary1.n2 = idx2;
        boundary1.n3 = idx3;
        boundary1.h = h;
        boundary1.Tinf = Tinf;
        boundary_edges.push_back(boundary1);
      }
    }
    myfile.clear();
    myfile.seekg(0);
    read = false;
    while (std::getline(myfile, str)) {
      if (str.rfind("*TRACTION", 0) == 0) {
        read = true;
        std::getline(myfile, str);
      }
      if (str.rfind("*", 0) == 0) {
        read = false;
      }
      if (read == true) {
        int n1, n2, n3;
        double tx;
        double ty;
        std::istringstream iss(str);
        char sep = ',';
        iss >> n1 >> sep >> n2 >> sep >> n3 >> sep >> tx >> sep >> ty;
        size_t idx1 = node_id_to_index[n1];
        size_t idx2 = node_id_to_index[n2];
        size_t idx3 = node_id_to_index[n3];
        TractionEdge boundary1;
        boundary1.n1 = idx1;
        boundary1.n2 = idx2;
        boundary1.n3 = idx3;
        boundary1.tx = tx;
        boundary1.ty = ty;
        traction_edges.push_back(boundary1);
      }
    }
    myfile.close();
    myfile2.close();
  }

  material.E = 210e6;
  material.nu = 0.3;
  material.t = 1.0;
  material.k = 45.0;
  material.alpha = 12e-6;
  material.T0 = 273.0;

  SolverInput input = {nodes, elems,          u_indices,      u,
                       F,     boundary_edges, traction_edges, material};

  return input;
}