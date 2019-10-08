// HW2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include <boost/property_map/property_map.hpp>
#include <boost/array.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/heap/priority_queue.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <limits> 
#include <stdlib.h> 
#include <cstdlib>
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include <ctime>

#define RANGE 10  // 10  100
#define COLUMNS	1000 
#define ROWS 30  // 30 60 80

using namespace std;
using namespace boost;
using namespace boost::heap;


typedef property<edge_weight_t, int> EdgeWeightProperty;

typedef adjacency_list_traits<listS, listS,
	directedS>::vertex_descriptor vertex_descriptor;

typedef property <vertex_name_t, int,
	property<vertex_index2_t, int, property<vertex_distance_t, int,
	property<vertex_predecessor_t, vertex_descriptor,
	property<vertex_color_t, int> > > > >VertexProperties;

typedef adjacency_list< listS, vecS,
	directedS,
	VertexProperties,
	EdgeWeightProperty
> Graph;
typedef graph_traits<Graph>::edge_iterator edge_iterator;

typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertex_iterator vertex_iter;
typedef pair<int, Vertex> pi;


void my_Dijksta(Graph &g,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &y,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w ) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);


	boost::heap::fibonacci_heap<pi > pq;


	std::pair<vertex_iter, vertex_iter> kp;
	kp = vertices(g);

	color[*kp.first] = 1;	//	mark s
	d[*kp.first] = 0;			//	set dist[0]
	pq.push(make_pair(d[*kp.first], *kp.first));	//	insert s into priority_queue

	pair<int, Vertex> top = pq.top();

	cout << top.first << "\t" << top.second << "\t" << color[top.second] << endl; // ektypwnei tin prwth timh tou grafhmatos

	while (!pq.empty()) {	// loop

		//	examine vertex u


		pair<int, Vertex> g1 = pq.top();
		pq.pop();
		Vertex u = g1.second;
		int dist_check = g1.first;

		cout << "g1 first: " << g1.first << "\t g1.second: " << g1.second << endl;
		cout << "top fisrt" << d[u] << endl;
		if (d[u] != dist_check)
			continue;

		graph_traits<Graph>::out_edge_iterator ei, ei_end;
		for (boost::tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
			std::cout << "( vertex: " << index[target(*ei, g)] << ")";
			std::cout << "\t weight: " << edge_w[*ei] << endl;
			cout << "distance: " << d[target(*ei, g)] << endl;
			if (edge_w[*ei] + d[u] < d[target(*ei, g)]) {

				cout << "Edge relaxed" << endl;

				d[target(*ei, g)] = edge_w[*ei] + d[u];
				cout << "update distance: " << d[target(*ei, g)] << endl;
				p[target(*ei, g)] = &u;

				if (color[target(*ei, g)] == 0) {
					color[target(*ei, g)] = 1;
					pq.push(make_pair(d[target(*ei, g)], target(*ei, g)));
				}
				else if (color[target(*ei, g)] == 1) {
					cout << endl;

					pq.push(make_pair(d[target(*ei, g)], target(*ei, g)));



				}

			}
			else {
				cout << "Edge not relaxed" << endl;
			}
		}

		color[u] = 2;

	}

}



void my_Dijksta_SP(Graph &g, Vertex s, Vertex t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &y,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);


	boost::heap::fibonacci_heap<pi > pq;


	std::pair<vertex_iter, vertex_iter> kp;
	kp = vertices(g);

	color[s] = 1;	//	mark s
	d[s] = 0;			//	set dist[0]
	pq.push(make_pair(d[s], s));	//	insert s into priority_queue

	pair<int, Vertex> top = pq.top();

	int count_k = 0;
	//cout << top.first << "\t" << top.second << "\t" << color[top.second] << endl; // ektypwnei tin prwth timh tou grafhmatos

	while (!pq.empty()) {	// loop

		//	examine vertex u


		pair<int, Vertex> g1 = pq.top();
		pq.pop();
		Vertex u = g1.second;
		int dist_check = g1.first;

		//cout << "g1 first: " << g1.first << "\t g1.second: " << g1.second << endl;
		//cout << "top fisrt" << d[u] << endl;
		if (d[u] != dist_check)
			continue;

		graph_traits<Graph>::out_edge_iterator ei, ei_end;
		for (boost::tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
			//std::cout << "( vertex: " << index[target(*ei, g)] << ")";
			//std::cout << "\t weight: " << edge_w[*ei] << endl;
			//cout << "distance: " << d[target(*ei, g)] << endl;
			if (edge_w[*ei] + d[u] < d[target(*ei, g)]) {

				//cout << "Edge relaxed" << endl;

				d[target(*ei, g)] = edge_w[*ei] + d[u];
				//cout << "update distance: " << d[target(*ei, g)] << endl;
				p[target(*ei, g)] = &u;

				if (color[target(*ei, g)] == 0) {
					color[target(*ei, g)] = 1;
					pq.push(make_pair(d[target(*ei, g)], target(*ei, g)));
				}
				else if (color[target(*ei, g)] == 1) {
					//cout << endl;

					pq.push(make_pair(d[target(*ei, g)], target(*ei, g)));



				}

			}
			else {
				//cout << "Edge not relaxed" << endl;
				count_k++;
			}
		}

		color[u] = 2;

		if (color[t] == 2) break;

	}

}

void my_A_star(Graph &g, Vertex s, Vertex t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &y,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);


	boost::heap::fibonacci_heap<pi > pq;


	std::pair<vertex_iter, vertex_iter> kp;
	kp = vertices(g);

	color[s] = 1;	//	mark s
	d[s] = 0;			//	set dist[0]
	pq.push(make_pair(d[s], s));	//	insert s into priority_queue

	pair<int, Vertex> top = pq.top();

	int count_k = 0;
	//cout << top.first << "\t" << top.second << "\t" << color[top.second] << endl; // ektypwnei tin prwth timh tou grafhmatos

	while (!pq.empty()) {	// loop

		//	examine vertex u


		pair<int, Vertex> g1 = pq.top();
		pq.pop();
		Vertex u = g1.second;
		int dist_check = g1.first;

		//cout << "g1 first: " << g1.first << "\t g1.second: " << g1.second << endl;
		//cout << "top fisrt" << d[u] << endl;
		if (d[u] != dist_check)
			continue;

		graph_traits<Graph>::out_edge_iterator ei, ei_end;
		for (boost::tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
			//std::cout << "( vertex: " << index[target(*ei, g)] << ")";
			//std::cout << "\t weight: " << edge_w[*ei] << endl;
			//cout << "distance: " << d[target(*ei, g)] << endl;

			edge_w[*ei] = edge_w[*ei] +  sqrt(pow((x[target(*ei, g)] - x[t]), 2) + pow((y[target(*ei, g)] - y[t]), 2) ) - sqrt(pow((x[u] - x[t]), 2) + pow((y[u] - y[t]), 2));

			if (edge_w[*ei] + d[u] < d[target(*ei, g)]) {

				//cout << "Edge relaxed" << endl;

				d[target(*ei, g)] = edge_w[*ei] + d[u];
				//cout << "update distance: " << d[target(*ei, g)] << endl;
				p[target(*ei, g)] = &u;

				if (color[target(*ei, g)] == 0) {
					color[target(*ei, g)] = 1;
					pq.push(make_pair(d[target(*ei, g)], target(*ei, g)));
				}
				else if (color[target(*ei, g)] == 1) {
					//cout << endl;

					pq.push(make_pair(d[target(*ei, g)], target(*ei, g)));

				}

			}
			else {
				//cout << "Edge not relaxed" << endl;
				count_k++;
			}
		}

		color[u] = 2;

		if (color[t] == 2) break;

	}

}


void create_grid_graph(Graph &g, 
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &y,
	property_map<Graph, vertex_distance_t>::type &d,
	property_map<Graph, vertex_predecessor_t>::type &p,
	property_map<Graph, vertex_color_t>::type &color,
	property_map<Graph, edge_weight_t>::type &edge_w) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);



	int num_of_nodes = ROWS * COLUMNS;

	//std::cout << "vertices(g) = ";

	std::pair<vertex_iter, vertex_iter> vp;


	//std::cout << std::endl;


	vp = vertices(g);
	int j = 0;
	for (int i = 0; i < ROWS; i++) {
		for (j = 0; j < COLUMNS - 1; j++) {
			//cout << " col |  ";
			srand(time(NULL));
			add_edge(*vp.first, *vp.first + 1, rand() % RANGE + 1, g);
			add_edge(*vp.first + 1, *vp.first, rand() % RANGE + 1, g);
			if (i < ROWS - 1) {
				add_edge(*vp.first, *vp.first + COLUMNS, rand() % RANGE + 1, g);
				add_edge(*vp.first + COLUMNS, *vp.first, rand() % RANGE + 1, g);
			}
			x[*vp.first] = i;
			y[*vp.first] = j;
			d[*vp.first] = std::numeric_limits<int>::max();
			color[*vp.first] = 0;
			//cout << " x: " << x[*vp.first] << " y: " << y[*vp.first];
			vp.first++;
		}
		//cout << " col |  ";
		x[*vp.first] = i;
		y[*vp.first] = j;
		d[*vp.first] = std::numeric_limits<int>::max();
		color[*vp.first] = 0;
		//cout << " x: " << x[*vp.first] << " y: " << y[*vp.first];
		if (i < ROWS - 1) {
			srand(time(NULL));
			add_edge(*vp.first, *vp.first + COLUMNS, rand() % RANGE + 1, g);
			add_edge(*vp.first + COLUMNS, *vp.first, rand() % RANGE + 1, g);
		}
		vp.first++;
		//cout << "\nchange row\n";
	}



	std::ofstream gout;
	gout.open("test.txt");
	boost::write_graphviz(gout, g);


}


void my_choose_vertex(Graph &g, Vertex &s, Vertex &t,
	property_map<Graph, vertex_name_t>::type &x,
	property_map<Graph, vertex_index2_t>::type &y
	) {

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);

	srand(time(NULL));
	int s_start = rand() % ROWS;
	int t_sink = rand() % ROWS;


	cout << "s_st: " << s_start << "\t t_si: " << t_sink << endl;

	std::pair<vertex_iter, vertex_iter> vp;
	for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
		if (x[*vp.first] == s_start && y[*vp.first] == 0) {
			s = *vp.first;
		}
		if (x[*vp.first] == t_sink && y[*vp.first] == (COLUMNS-1)) {
			t = *vp.first;
		}
	}



}


int main(int argc, char** argv) {

	Graph g;

	property_map<Graph, vertex_name_t>::type
		x = get(vertex_name, g);
	property_map<Graph, vertex_index2_t>::type
		y = get(vertex_index2, g);
	property_map<Graph, vertex_distance_t>::type
		d = get(vertex_distance, g);
	property_map<Graph, vertex_predecessor_t>::type
		p = get(vertex_predecessor, g);
	property_map<Graph, vertex_color_t>::type
		color = get(vertex_color, g);
	property_map<Graph, edge_weight_t>::type
		edge_w = get(edge_weight, g);

	

	cout << " Max value fot int is: " << std::numeric_limits<int>::max() << endl;

	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);
	

	create_grid_graph(g, x, y, d, p, color, edge_w);


	


	//my_Dijksta(g, x , y, d, p, color, edge_w);
	
	Vertex s, t;

	my_choose_vertex(g, s, t, x, y);

	cout << "s: " << index[s] << "\t x: " << x[s] << "\t y: " << y[s] << endl;
	cout << "t: " << index[t] << "\t x: " << x[t] << "\t y: " << y[t] << endl;
	
	std::clock_t start, stop;
	double duration;

	start = std:: clock();

	my_Dijksta_SP(g, s, t, x, y, d, p, color, edge_w);
	
	//my_A_star(g, s, t, x, y, d, p, color, edge_w);

	stop = std:: clock();

	duration = (stop - start) / (double)CLOCKS_PER_SEC;

	cout << "TIME: \t"<< duration << endl;
	

	return 0;
}

