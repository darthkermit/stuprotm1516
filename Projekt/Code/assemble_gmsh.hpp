//////////////////////////////////////////////////////////////////////////////////////////////
//																							//
// !!! Diese Beschreibung enthaelt mehr Ziele als bereits Vorhandenes !!!					//
//																							//
// Diese Datei ist zZt nur dazu da die Geometrie aus einer gmsh-Datei einzulesen			//
// um daraus dann eine Steifigkeitsmatrix fuer die Waermeleitungsgleichung mit				//
// implizitem Euler und 7-Punkte-Stern zu konstruieren.										//
//																							//
// Noch zu ueberlegen: Wo werden Materialien und die zugehoerigen Konstanten hinterlegt?	//
//																							//
// Es koennen auch Startwerte(Temperatur) eingelesen werden.								//
// 																							//
//////////////////////////////////////////////////////////////////////////////////////////////

// Notiz-Zettel fuer mich:
// FreeCAD
// gmsh
// blender

#include"DIA.hpp"		// wird noch nicht benutzt
#include"Vector.hpp"	// wird noch nicht benutzt
#include<fstream>		// um Dateien einzulesen
#include<iostream>		// Ausgabe auf der Konsole
#include<string>		// fuer Dateinamen und Einlesen von Zeilen
using namespace std;	// Bequemlichkeit
// weitere includes folgen



// Liste an Parametern, die vorkommen:
//
// Fuer jedes Material:
// alpha: thermal diffusvity - Temperaturleitfaehigkeit [m/(s^2)]
// k: thermal conductivity - Waermeleitfaehigkeit [W/(mK)]
// rho: desnity - Dichte [kg/(m^3)]
// c_p: specific heat capacity - spezifische Waermekapazitaet [J/(kg K)]
//		-> alpha = k / (rho * c_p) [m^2/s]
// 


////////////////////////////////////////////////////
// Einschub: Materialklasse (verschiedene Parameter, s.o.)
template <typename datatype>
class material {
public:
	string name;
	datatype alpha;
	datatype k;
	datatype rho;
	datatype c_p;
	material(string n, datatype a, datatype setk, datatype r, datatype c): name(n), alpha(a), k(setk), rho(r), c_p(c) {};
	material(string n): name(n) {/*TODO: Bekannte Daten aus Datei einlesen*/};
	~material() {};
	void showAllData() { cout << "Data of material " /*<< name */<< ":" << endl; };	// TODO
private:
};
////////////////////////////////////////////////////

////////////////////////////////////////////////////
// Einschub: Punkt-Klasse (zZt nur Koordinaten, spaeter evtl mehr (Verbindungslinien o.ae.))
class point {
public:
	int ID;
	float x, y, z;
	point() {};
	point(float setx, float sety, float setz): x(setx), y(sety), z(setz) {};
	~point() {};
	void setID(int nr) { ID = nr; };
	void setCoords(float setx, float sety, float setz) { x = setx; y = sety; z = setz; };
	void setX(float setx) { x = setx; };
	void setY(float sety) { y = sety; };
	void setZ(float setz) { z = setz; };
};
////////////////////////////////////////////////////

////////////////////////////////////////////////////
// Einschub: Gegenstand-Klasse (Geometrie und Material)
template <typename datatype>
class object {
public:
	string name;	// zZt optional
	// evtl: unsigned int ID?
	int maxNumPoints, numPoints;
	point* points;
	material<datatype>* mat;
	object(int mnP): maxNumPoints(mnP), numPoints(0) {points = new point[maxNumPoints];};	// TODO: evtl nicht von vornherein festlegen, sondern nach und nach erweitern. Warhscheinlich aber nicht (stabil) moeglich.
	object(): numPoints(0) {};
	~object() {};
	void set(int mnP) {maxNumPoints = mnP; points = new point[maxNumPoints];};
	void addPoint(int ID) { points[numPoints].setID(ID); ++numPoints; };
	//void addPointCoord(float x, float y, float z) {points[numPoints].setCoords(x, y, z);};
	void showAllData() {
		cout << "Name: " << name << endl;
		cout << "Max points: " << maxNumPoints << ", actual num points: " << numPoints << endl;
		for (int i(0); i < numPoints; ++i)
			cout << "P" << i << ": " << points[i].x << ", " << points[i].y << ", " << points[i].z << endl;
		(*mat).showAllData();
	};
private:
};
////////////////////////////////////////////////////


void setStiffnessMatrix(/*filename*/) {
	
	/////////// Einlesen und Speichern der Eckdaten /////////////

	ifstream fin;
	fin.open("testphys.geo");
	int number(0);
	char line[100];
	object<float> objects[10];
	int objNr(0);
	do {
		fin.get(line, 3);
		if (line[0] == 'P' && line[1] == 'h') {
			fin.get(line, 14);
			fin >> number;
			fin.get(line, 6);
			objects[objNr].set(100);
			do {
				fin >> number >> line[0];
				//cout << number << ", ";
				objects[objNr].addPoint(number);
			} while (line[0] != '}');
			//cout << endl;
			fin.getline(line, 100);
			++objNr;
		}
		else {
			fin.getline(line, 100);
		}
	} while(line[0]);
	fin.close();

	fin.open("test.geo");
	float coords[3];
	do {
		fin.get(line, 3);
		if (line[0] == 'P' && line[1] == 'o') {
			fin.get(line, 5);
			fin >> number;
			fin.get(line, 6);
			fin >> coords[0];
			fin >> line[0] >> coords[1];
			fin >> line[0] >> coords[2];
			//cout << number << ": " << coords[0] << ", " << coords[1] << ", " << coords[2] << endl;
			for (int i(0); i < objNr; ++i)
				for (int j(0); j < objects[i].numPoints; ++j)
					if (objects[i].points[j].ID == number)
							objects[i].points[j].setCoords(coords[0], coords[1], coords[2]);
			fin.getline(line, 100);
		}
		else {
			fin.getline(line, 100);
		}
	} while(line[0]);
	fin.close();

	cout << endl << "Num objects: " << objNr << endl << endl;
	for (int i(0); i < objNr; ++i) {
		cout << "object " << i << ":" << endl;
		objects[i].showAllData();
		cout << endl;
	}

	//////////////// Hier kommt die eigentliche Assemblierung /////////////////


}

void setStartValues(/*filename*/) {
	// TODO
}
