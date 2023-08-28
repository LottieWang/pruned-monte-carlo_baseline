#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include "../src/pmc.hpp"
#include "../src/parseCommandLine.hpp"
#include "../src/graph.hpp"

using namespace std;

int main(int argc, char **argv) {
	if (argc < 4) {
		cerr << "./pmc graph k R" << endl;
		exit(1);
	}
	char* file = argv[1];
	int k = atoi(argv[2]);
	int R = atoi(argv[3]);

	// ifstream is(file.c_str());
	// int u, v;
	// double p;
	// for (; is >> u >> v >> p;) {
	// 	if (u == v) {
	// 		continue;
	// 	}
	// 	es.push_back(make_pair(make_pair(u, v), p));
	// }
	// is.close();
	CommandLine P(argc, argv);
	double w = P.getOptionDouble("-UIC", 0.0);
	double u1 = P.getOptionDouble("-ua", 1.0);
	double u2 = P.getOptionDouble("-ub", 1.0);
	bool WIC = P.getOption("-WIC");
	float compact = P.getOptionDouble("-compact", 1.0);
	printf("w: %f ua: %f ub: %f WIC: %d k: %d R: %d compact: %f\n", w, u1, u2, WIC, k, R, compact);
	// Graph graph = read_graph(file);
	Graph graph=read_binary(file);
	if (w!=0)  AssignUniWeight(graph, w);
	else if (u2!= 1) AssignUniformRandomWeight(graph, u1 ,u2);
	else if (WIC) AssignWICWeight(graph);
	else {
		std::cout << "no weight assginment specified. -w [float] for uni weight, -u [float] for uniform (float, float+0.1) -WIC for WIC_SYM" << endl;
		return 1;
	}
	vector<pair<pair<int, int>, double> > es(graph.m);
	for (size_t i = 0; i<graph.n; i++){
		for (size_t j = graph.offset[i]; j<graph.offset[i+1]; j++){
			es[j]=make_pair(make_pair((int)i, (int)graph.E[j]),graph.W[j]);
		}
	}

	InfluenceMaximizer im;
	vector<int> seeds;
	parlay::internal::timer t("PMC", false);
	for (int i = 0; i< 4; i++){
		t.start();
		// im.run(es, k, R);
		seeds = im.run(es, k, R);
		t.next("Time");
	}
	std::cout<<"seeds: " << std::endl;
	for (int i = 0; i < k; i++) {
		// std::cout << i << "-th seed =\t" << seeds[i] << endl;
		std::cout << seeds[i] << " ";
	}
	std::cout << std::endl;

	return 0;
}
