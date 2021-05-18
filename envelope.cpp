#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <dolfin.h>
#include <cmath>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>

#include "envelope.h"

using namespace dolfin;
using std::make_shared;
using std::shared_ptr;
using std::vector;
using std::string;

static int k;
static double f = 1.0;
static double amplitude = 1000000000000000000;

class Left : public SubDomain {
	bool inside(const Array<double> &x, bool on_boundary) const {
		return x[0] <= DOLFIN_EPS;
	}
};

class Left_Round : public SubDomain {
	bool inside(const Array<double> &x, bool on_boundary) const {
		double x0 = 4.0, y0 = 0.0, r = 1.0;
		return x[0] <= DOLFIN_EPS || (x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0) <= r*r + DOLFIN_EPS;
	}
};

class BF_lense : public Expression {
public:
	BF_lense() : Expression(2) {}

	void eval(Array<double> &values, const Array<double> &x) const {
		values[0] = amplitude*cos(k*x[1]*x[1]/(2*f));
		values[1] = -amplitude*sin(k*x[1]*x[1]/(2*f));
	}
};

class BF_gap : public Expression {
public:
	BF_gap() : Expression(2) {}

	void eval(Array<double> &values, const Array<double> &x) const {
		if (x[1] <= 0.2 && x[1] >= -0.2)
			values[0] = amplitude;
		else
			values[0] = 0;
		values[1] = 0;
	}
};

class BF_screen : public Expression {
public:
	BF_screen() : Expression(2) {}

	void eval(Array<double> &values, const Array<double> &x) const {
		if (x[1] <= 1 && x[1] >= -1)
			values[0] = 0;
		else
			values[0] = amplitude;
		values[1] = 0;
	}
};


class BF_round : public Expression {
public:
	BF_round() : Expression(2) {}

	void eval(Array<double> &values, const Array<double> &x) const {
		if (x[0] <= DOLFIN_EPS)
			values[0] = amplitude;
		else
			values[0] = 0;
		values[1] = 0;
	}
};

class dU0 : public Expression {
	void eval(Array<double> &values, const Array<double> &r) const {
		values[0] = 0.0;
		values[1] = 0.0;
	}
};

void save(const Function &u, const std::string &filename) {
	std::fstream f(filename, f.out);
	f << "x\ty\treu\timu\tk\n";
	if (!f.is_open())
		throw 1;
	for (double x = 0.0; x <= 7.5; x += 0.1)
		for (double y = -3.0; y <= 3.0; y += 0.001)
			f << x << '\t' << y << '\t' << u[0](x, y) << '\t' << u[1](x, y) << '\t' << k << '\n';
	f.close();
}

void write(const Function &u, const std::string &filename) {
	auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto dumpPoints = vtkSmartPointer<vtkPoints>::New();

	auto v = vtkSmartPointer<vtkDoubleArray>::New();
	v->SetName("J");

	shared_ptr<const Mesh> mesh = u.function_space()->mesh();
	vector<double> coords = mesh->coordinates();
	vector<unsigned int> cells = mesh->cells();

	unsigned int number = (unsigned int)coords.size();
	if (number % 2 != 0)
		throw 1;
	number = number / 2;
	for (int i = 0; i < number; i++) {
		dumpPoints->InsertNextPoint(coords[2*i], coords[2*i + 1], 0);
		double val0 = u[0](coords[2*i], coords[2*i + 1], 0);
		double val1 = u[1](coords[2*i], coords[2*i + 1], 0);
		v->InsertNextValue(val0*val0 + val1*val1);
	}

	unstructuredGrid->SetPoints(dumpPoints);

	unsigned int cell_number = (unsigned int)cells.size();
	if (cell_number % 3 != 0) {
		std::cerr << cell_number << '\n';
		throw 1;
	}
	cell_number = cell_number / 3;
	for(unsigned int i = 0; i < cell_number; i++) {
		auto tetra = vtkSmartPointer<vtkTriangle>::New();
		tetra->GetPointIds()->SetId(0, cells[3*i]);
		tetra->GetPointIds()->SetId(1, cells[3*i + 1]);
		tetra->GetPointIds()->SetId(2, cells[3*i + 2]);
		unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
	}


	unstructuredGrid->GetPointData()->AddArray(v);

	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(unstructuredGrid);
	writer->Write();
}

void Solve(shared_ptr<Expression> bf, shared_ptr<SubDomain> boundary, const string &name) {
	auto mesh = make_shared<RectangleMesh>(Point(0.0, -3.0), Point(8.0, 3.0), 64, 512, "crossed");
	auto V = make_shared<envelope::FunctionSpace>(mesh);

	DirichletBC bc(V, bf, boundary);

	envelope::BilinearForm a(V, V);
	a.k = make_shared<Constant>(k);
	a.n = make_shared<Constant>(1);
	Function u(V);
	envelope::LinearForm L(V);
	L.zero = make_shared<Constant>(0.0);
	solve(a == L, u, bc);

	File f(name + std::to_string(k) + ".pvd");
	f << u[0];
	save(u, name + std::to_string(k) + ".tsv");
}

    
int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: wave <wave number>\n";
		return 1;
	}
	k = std::stoi(argv[1]);

	Solve(make_shared<BF_lense>(), make_shared<Left>(), string("lense"));
	Solve(make_shared<BF_gap>(), make_shared<Left>(), string("gap"));
	Solve(make_shared<BF_screen>(), make_shared<Left>(), string("screen"));
	Solve(make_shared<BF_round>(), make_shared<Left_Round>(), string("round"));

	return 0;
}
