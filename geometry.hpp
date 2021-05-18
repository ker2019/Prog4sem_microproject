#include <vector>
#include <string>
#include <dolfin.h>

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

class Geometry {
private:
	static int number;
	static int numberOfCreated;
	bool isCreated = false;
	int dim;
	std::string name;
	std::vector<int> pointTags;
	std::vector<int> curveTags;
	std::vector<int> surfaceTags;
	std::vector<int> volumeTags;
	virtual void createGeom(double lc) = 0;
public:
	void create(double lc);
	void destroy();
	void show() const;
	dolfin::Mesh getMeshAsDolfin() const;
	~Geometry();
protected:
	Geometry(int dim, const std::string &name);
	int addPoint(double x, double y, double z, double lc);
	int addCircleArc(int a, int o, int b);
	int addLine(int a, int b);
	int addCurveLoop(std::vector<int> curves);
	int addPlaneSurface(std::vector<int> curveLoops);
};

#endif
