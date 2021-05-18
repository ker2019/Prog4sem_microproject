#include <vector>
#include <iostream>
#include <string>
#include <gmsh.h>
#include <dolfin.h>
#include "geometry.hpp"

using std::vector;
using std::size_t;

int Geometry::number = 0;
int Geometry::numberOfCreated = 0;

Geometry::Geometry(int dim, const std::string &name) {
	this->dim = dim;
	this->name = name;
	if (number == 0)
		gmsh::initialize();
	number++;
}

Geometry::~Geometry() {
	if (number == 1)
		gmsh::finalize();
	number--;
}

int Geometry::addPoint(double x, double y, double z, double lc) {
	int tag = gmsh::model::geo::addPoint(x, y, z, lc);
	pointTags.push_back(tag);
	return tag;
}

int Geometry::addCircleArc(int a, int o, int b) {
	int tag = gmsh::model::geo::addCircleArc(a, o, b);
	curveTags.push_back(tag);
	return tag;
}

int Geometry::addLine(int a, int b) {
	int tag = gmsh::model::geo::addLine(a, b);
	curveTags.push_back(tag);
	return tag;
}

int Geometry::addCurveLoop(vector<int> curves) {
	return gmsh::model::geo::addCurveLoop(curves);
}

int Geometry::addPlaneSurface(vector<int> curveLoops) {
	int tag = gmsh::model::geo::addPlaneSurface(curveLoops);
	surfaceTags.push_back(tag);
	return tag;
}

void Geometry::create(double lc) {
	if (numberOfCreated > 0)
		return;
	isCreated = true;
	gmsh::model::add(name);
	createGeom(lc);
	gmsh::model::geo::synchronize();
	gmsh::model::mesh::generate(dim);
	numberOfCreated++;
}

void Geometry::destroy() {
	if (isCreated) {
		gmsh::model::remove();
		isCreated = false;
		numberOfCreated--;
	}
}

void Geometry::show() const {
	gmsh::fltk::run();
}

dolfin::Mesh Geometry::getMeshAsDolfin() const {
	int elemType = (dim == 2 ? 2 : 4);
	vector<double> nodesCoord;
	vector<size_t> nodeTags;
	vector<double> parametricCoord;
	gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

	vector<size_t> *cellTags;
	vector<int> elementTypes;
	vector<vector<size_t>> elementTags;
	vector<vector<size_t>> elementNodeTags;
	gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
	for(int i = 0; i < elementTypes.size(); i++)
		if(elementTypes[i] == elemType) {
			cellTags = &elementNodeTags[i];
			break;
		}

	dolfin::MeshEditor meshEd;
	dolfin::Mesh mesh;
	if (dim == 2) {
		meshEd.open(mesh, dolfin::CellType::Type::triangle, 2, 2, 1);
		meshEd.init_vertices(nodesCoord.size()/3);
		meshEd.init_cells(cellTags->size()/3);
		for (size_t i = 0; i < nodesCoord.size()/3; i++)
			meshEd.add_vertex(i, nodesCoord[3*i], nodesCoord[3*i + 1]);
		for (size_t i = 0; i < cellTags->size()/3; i++)
			meshEd.add_cell(i, (*cellTags)[3*i] - 1, (*cellTags)[3*i + 1] - 1, (*cellTags)[3*i + 2] - 1);
	}
	else if (dim == 3) {
		meshEd.open(mesh, dolfin::CellType::Type::tetrahedron, 3, 3, 1);
		meshEd.init_vertices(nodesCoord.size()/3);
		meshEd.init_cells(cellTags->size()/4);
		for (size_t i = 0; i < nodesCoord.size()/3; i++)
			meshEd.add_vertex(i, nodesCoord[3*i], nodesCoord[3*i + 1], nodesCoord[3*i + 2]);
		for (size_t i = 0; i < cellTags->size()/4; i++)
			meshEd.add_cell(i, (*cellTags)[4*i] - 1, (*cellTags)[4*i + 1] - 1, (*cellTags)[4*i + 2] - 1, (*cellTags)[4*i + 3] - 1);
	}
	meshEd.close();
	return mesh;
}
