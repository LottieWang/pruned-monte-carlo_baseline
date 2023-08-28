#pragma once

#include <fcntl.h>
// #include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/random.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <fstream>
#include <vector>
#include <map>

#include "utilities.h"
// #include "get_time.hpp"
using namespace std;
using namespace parlay;

struct edge {
  NodeId u;
  NodeId v;
  edge() : u(0), v(0) {}
  edge(NodeId f, NodeId s) : u(f), v(s) {}
  bool operator<(const edge& rhs) const {
    if (u != rhs.u) {
      return u < rhs.u;
    } else {
      return v < rhs.v;
    }
  }
  bool operator!=(const edge& rhs) const { return u != rhs.u || v != rhs.v; }
};

struct Graph {
  size_t n;
  size_t m;
  bool symmetric;
  sequence<EdgeId> offset;
  sequence<NodeId> E;
  sequence<EdgeId> in_offset;
  sequence<NodeId> in_E;
  sequence<float> W;
  sequence<float> in_W;
};



bool is_space(char c) {
  switch (c) {
    case '\r':
    case '\t':
    case '\n':
    case ' ':
    case 0:
      return true;
    default:
      return false;
  }
}


Graph read_binary(const char* const filename, bool enable_mmap = true) {
  Graph graph;
  if (enable_mmap) {
    struct stat sb;
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
      cerr << "Error: Cannot open file " << filename << '\n';
      abort();
    }
    if (fstat(fd, &sb) == -1) {
      cerr << "Error: Unable to acquire file stat" << filename << '\n';
      abort();
    }
    void* memory = mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    char* data = static_cast<char*>(memory);
    size_t len = sb.st_size;
    size_t n = reinterpret_cast<uint64_t*>(data)[0];
    size_t m = reinterpret_cast<uint64_t*>(data)[1];
    size_t sizes = reinterpret_cast<uint64_t*>(data)[2];
    assert(sizes == (n + 1) * 8 + m * 4 + 3 * 8);
    graph.n = n, graph.m = m;
    graph.offset = sequence<EdgeId>(n + 1);
    graph.E = sequence<NodeId>(m);
    for(size_t i = 0;i< n + 1; i++) {
      graph.offset[i] = reinterpret_cast<uint64_t*>(data + 3 * 8)[i];
    }
    for(size_t i=0; i< m;i++) {
      graph.E[i] = reinterpret_cast<uint32_t*>(data + 3 * 8 + (n + 1) * 8)[i];
    }

    if (memory != MAP_FAILED) {
      if (munmap(memory, len) != 0) {
        cerr << "Error: munmap failed" << '\n';
        abort();
      }
    }
  } else {
    ifstream ifs(filename);
    if (!ifs.is_open()) {
      cerr << "Error: Cannot open file " << filename << '\n';
      abort();
    }
    size_t n, m, sizes;
    ifs.read(reinterpret_cast<char*>(&n), sizeof(size_t));
    ifs.read(reinterpret_cast<char*>(&m), sizeof(size_t));
    ifs.read(reinterpret_cast<char*>(&sizes), sizeof(size_t));
    assert(sizes == (n + 1) * 8 + m * 4 + 3 * 8);

    graph.n = n;
    graph.m = m;
    sequence<uint64_t> offset(n + 1);
    sequence<uint32_t> edge(m);
    ifs.read(reinterpret_cast<char*>(offset.begin()), (n + 1) * 8);
    ifs.read(reinterpret_cast<char*>(edge.begin()), m * 4);
    graph.offset = sequence<EdgeId>(n + 1);
    graph.E = sequence<NodeId>(m);
    for(size_t i = 0; i< n + 1;i++) {
        graph.offset[i] = offset[i]; }
    for(size_t i=0; i< m; i++) { 
        graph.E[i] = edge[i]; }
    if (ifs.peek() != EOF) {
      cerr << "Error: Bad data\n";
      abort();
    }
    ifs.close();
  }
  return graph;
}


struct Hash_Edge {
    EdgeId hash_graph_id;
    // NodeId n;
    // NodeId id_term;
    // bool forward;
    bool operator()(NodeId u, NodeId v, float w) const {
      return _hash(hash_graph_id ^ _hash(((EdgeId)min(u,v)<< 32) + max(u,v))) < UINT_E_MAX*w;
      // if (!forward) swap(u,v);
      // return _hash(_hash(u)+v)+_hash(graph_id) < w*UINT_N_MAX;
      // return _hash( graph_id*n*n + u*n + v) < w*UINT_N_MAX;
      // return _hash(id_term + u*n + v) < w*UINT_N_MAX;

    }
};


void AssignUniWeight(Graph& graph, float w){
  graph.W = sequence<float>(graph.m);
  for(size_t i=0; i< graph.m; i++){
    graph.W[i] = w;
  }
  if (!graph.symmetric){
    graph.in_W = sequence<float>(graph.m);
    parallel_for(0, graph.m, [&](size_t i){
      graph.in_W[i]=w;
    });
  }
}

double Uniform(size_t n, size_t u, size_t v, double l, double r) {
  if (u > v) std::swap(u, v);
  double p = 1.0 * hash64(hash64(hash64(n) + u + 1) + v + 1) / std::numeric_limits<uint64_t>::max();
  return l + (r - l) * p;
}

void AssignUniformRandomWeight(Graph& graph, double l, double r) {
  graph.W = sequence<float>(graph.m);
  for(size_t u =0; u< graph.n; u++) {
    for(size_t j = graph.offset[u]; j<graph.offset[u + 1]; j++) {
      auto v = graph.E[j];
      graph.W[j] = Uniform(graph.n, u, v, l, r);
    }
  }
}

void AssignWICWeight(Graph& graph) {
  graph.W = sequence<float>(graph.m);
  for(size_t u=0; u< graph.n; u++) {
    for(size_t j = graph.offset[u]; j< graph.offset[u + 1]; j++) {
      auto v = graph.E[j];
      auto deg_u = graph.offset[u + 1] - graph.offset[u];
      auto deg_v = graph.offset[v + 1] - graph.offset[v];
      graph.W[j] = 2.0 / (deg_u + deg_v);
    }
  }
}

// void AssignIndegreeWeightSym(Graph &graph){
//   graph.W = sequence<float>(graph.m);
//   for(size_t i =0; i< graph.n; i++ ){
//     EdgeId in_degree = graph.offset[i+1] - graph.offset[i];
//     if (in_degree >0){
//       float w = 1.0/in_degree;
//       for(size_t j =graph.offset[i]; j< graph.offset[i+1]; j++){
//         graph.W[j] = w;
//       }
//     }
//   }
// }

// void AssignIndegreeWeightDir(Graph &graph){
//   graph.in_W = sequence<float>(graph.m);
//   graph.W = sequence<float>(graph.m);
//   sequence<tuple<edge,float>> edge_list = sequence<tuple<edge,float>>(graph.m);
//   parallel_for(0, graph.n, [&](size_t i){
//     EdgeId in_degree = graph.in_offset[i+1] - graph.in_offset[i];
//     if (in_degree >0){
//       float w = 1.0/in_degree;
//       parallel_for(graph.in_offset[i], graph.in_offset[i+1], [&](size_t j){
//         graph.in_W[j]=w;
//         NodeId u = graph.in_E[j];
//         edge_list[j] = make_tuple(edge(u,i),w);
//       });
//     }
//   });
//   sort_inplace(edge_list);
//   parallel_for(0,graph.m, [&](size_t i){
//     graph.W[i]=get<1>(edge_list[i]);
//   });
// }

// void AssignIndegreeWeight(Graph &graph){
//   if (graph.symmetric){
//     AssignIndegreeWeightSym(graph);
//   }else{
//     AssignIndegreeWeightDir(graph);
//   }
// }

// void WriteSampledGraph(Graph graph, const char* const OutFileName){
//   sequence<edge> edge_list = sequence<edge>(graph.m);
//   timer t;
//   t.start();
//   parallel_for(0, graph.n, [&](size_t i){
//     NodeId u = i;
//     parallel_for(graph.offset[i], graph.offset[i+1], [&](size_t j){
//       NodeId v = graph.E[j];
//       edge_list[j] = edge(u,v);
//     });
//   });
//   t.stop();
//   cout << "traverse edges costs: " << t.get_total()<<endl;

//   Hash_Edge hash_edge{0, (NodeId)graph.n ,true};
//   auto sample_edges = delayed_seq<bool>(edge_list.size(), [&](size_t i) {
//     NodeId u = edge_list[i].u;
//     NodeId v = edge_list[i].v;
//     float w = graph.W[i];
//     return hash_edge(u,v,w);
//   });
//   auto unique_edges = pack(edge_list, sample_edges);
//   edge_list.clear();
//   Graph graph_sampled = EdgeListToGraph(unique_edges, graph.n);

//   string out_name = OutFileName;
//   out_name+=".bin";
//   cout << "write to " << out_name << endl;
//   cout << "n: " << graph_sampled.n << " m: " << graph_sampled.m << endl;
//   ofstream ofs(out_name);
//   size_t sizes = (graph_sampled.n + 1) * 8 + graph_sampled.m * 4 + 3 * 8;
//   ofs.write(reinterpret_cast<char*>(&graph_sampled.n), sizeof(size_t));
//   ofs.write(reinterpret_cast<char*>(&graph_sampled.m), sizeof(size_t));
//   ofs.write(reinterpret_cast<char*>(&sizes), sizeof(size_t));
//   ofs.write(reinterpret_cast<char*>((graph_sampled.offset).data()), (graph_sampled.n + 1) * 8);
//   ofs.write(reinterpret_cast<char*>((graph_sampled.E).data()), graph_sampled.m * 4);
//   ofs.close();
// }