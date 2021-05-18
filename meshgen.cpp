#include <string>
#include <cmath>
#include "geometry.hpp"

class Disk : public Geometry {
private:
	double r, R;
	int dim = 2;
	std::string name = "disk";

	void createGeom(double lc) {
		int o = addPoint(0, 0, 0, lc);
		int A = addPoint(R, 0, 0, lc);
		int B = addPoint(0, R, 0, lc);
		int C = addPoint(-R, 0, 0, lc);
		int D = addPoint(0, -R, 0, lc);
		int ABCD = addCurveLoop({addCircleArc(A, o, B), addCircleArc(B, o, C), addCircleArc(C, o, D), addCircleArc(D, o, A)});

		int a = addPoint(r, 0, 0, lc);
		int b = addPoint(0, r, 0, lc);
		int c = addPoint(-r, 0, 0, lc);
		int d = addPoint(0, -r, 0, lc);
		int abcd = addCurveLoop({addCircleArc(a, o, b), addCircleArc(b, o, c), addCircleArc(c, o, d), addCircleArc(d, o, a)});

		addPlaneSurface({ABCD, abcd});
	}
public:
	Disk(double r, double R) : Geometry(2, "disk") {
		this->r = r;
		this->R = R;
	}
};

class Rectangle : public Geometry {
private:
	double x1, y1, x2, y2;

	void createGeom(double lc) {
		int A = addPoint(x1, y1, 0, lc);
		int B = addPoint(x1, y2, 0, lc);
		int C = addPoint(x2, y2, 0, lc);
		int D = addPoint(x2, y1, 0, lc);

		int AB = addLine(A, B);
		int BC = addLine(B, C);
		int CD = addLine(C, D);
		int DA = addLine(D, A);

		addPlaneSurface({addCurveLoop({AB, BC, CD, DA})});
	}
public:
	Rectangle(double x1, double y1, double x2, double y2) : Geometry(2, "rectangle") {
		this->x1 = x1;
		this->y1 = y1;
		this->x2 = x2;
		this->y2 = y2;
	}
};

class Lense : public Geometry {
private:
	double x, r1, r2, d, w, h;

	void createGeom(double lc) {
		int O1 = addPoint(x - r1 + d/2, 0, 0, lc);
		int O2 = addPoint(x + r2 - d/2, 0, 0, lc);
		int A1 = addPoint(x - r1 + d/2 + pow(r1*r1 - h*h/4, 0.5), h/2, 0, lc);
		int B1 = addPoint(x - r1 + d/2 + pow(r1*r1 - h*h/4, 0.5), -h/2, 0, lc);
		int A2 = addPoint(x + r2 - d/2 - pow(r2*r2 - h*h/4, 0.5), h/2, 0, lc);
		int B2 = addPoint(x + r2 - d/2 - pow(r2*r2 - h*h/4, 0.5), -h/2, 0, lc);

		int A1B1 = addCircleArc(A1, O1, B1);
		int B1B2 = addLine(B1, B2);
		int B2A2 = addCircleArc(B2, O2, A2);
		int A2A1 = addLine(A2, A1);

		int A = addPoint(0, h/2, 0, lc);
		int B = addPoint(0, -h/2, 0, lc);
		int C = addPoint(w, -h/2, 0, lc);
		int D = addPoint(w, h/2, 0, lc);
		int AB = addLine(A, B);
		int BB2 = addLine(B, B2);
		int B1C = addLine(B1, C);
		int CD = addLine(C, D);
		int A2A = addLine(A2, A);
		int DA1 = addLine(D, A1);

		addPlaneSurface({addCurveLoop({AB, BB2, B2A2, A2A})});
		addPlaneSurface({addCurveLoop({-B2A2, -B1B2, -A1B1, -A2A1})});
		addPlaneSurface({addCurveLoop({A1B1, B1C, CD, DA1})});
	}
public:
	Lense(double x, double r1, double r2, double d, double w, double h) : Geometry(2, "lense") {
		this->x = x;
		this->r1 = r1;
		this->r2 = r2;
		this->d = d;
		this->w = w;
		this->h = h;
	}
};

int main(int argc, char *argv[]) {
	std::set<std::string> args(argv, argv + argc);
	Disk d(1, 3);
	Rectangle r(0, -3, 15, 3);
	Lense l(3, 6.0, 6.0, 4.0, 8, 6);

	l.create(0.02);
	if(args.count("-popup")) d.show();
	dolfin::File f1("lense.xml");
	f1 << d.getMeshAsDolfin();
	d.destroy();
	return 0;
}
