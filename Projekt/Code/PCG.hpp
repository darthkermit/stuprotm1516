#pragma once

#include "Vector.hpp"
#include "DIA.hpp"
#include <iostream>
#include <math.h>
using namespace std;

// kommt bestimmt bald in Vector.hpp
template<typename data,typename data1, typename data2>
void sp(data& result, const Vector<data1>& v1, const Vector<data2>& v2)
// Anwendungsbsp: Skalarprodukt = sp(v1,v2);
{
	assert( v1.dim()==v2.dim() );

	data a(0);
	for(int i=0;i<v1._dim;i++)
	{
		a += static_cast<data>(v1[i])*static_cast<data>(v2[i]);
	}	
	result = a;
}

template <typename restype, typename mattype, typename vectype>
int CG(Vector<restype>& x, DIA<mattype>& A, Vector<vectype>& b) {
	// Fuer den Algorithmus werden einige Vektoren und Werte benoetigt. Hier die Definitionen und Initialisierungen (teilweise mit Dummi-Werten)
	Vector<restype> r(b.dim());	// r -> Residiuenvektor
	defect(r, A, b, x);	// r_0
	Vector<restype> d(r);	// Suchrichtung, d_0 := r_0
	Vector<restype> z(b.dim());	// Wird benutzt um in jedem Durchlauf A*d_k nur einmal ausrechnen zu muessen

	restype rr_old(0);	// ist im k-ten Durchlauf das Skalarprodukt von r_k*r_k
	restype rr_new(1);	// ist im k-ten Durchlauf das Skalarprodukt von r_(k+1)*r_(k+1)
	restype dz(1);		// ist im k-ten Durchlauf das Skalarprodukt von d_k*z_k
	
	sp(rr_new, r, r);	// wird fuer den ersten Schleifendurchlauf schon benoetigt
	
	Vector<restype> xold(b.dim());	// wird noch benoetigt um x_k auf x_k+1 zu aktualisieren, sollte noch aehnlich wie der Defekt neu programmiert werden um Zwischenschritte nicht speichern zu muessen
	Vector<restype> rold(b.dim());	// wird noch benoetigt um r_k auf r_k+1 zu aktualisieren, sollte noch aehnlich wie der Defekt neu programmiert werden um Zwischenschritte nicht speichern zu muessen

	double TOL(0.0000000001);	// (Vorerst eine beliebige feste) Toleranz fuer das Abbruchkriterium...

	int k(0);	// Anzahl der Iterationen. Wird zurueckgegeben, deshalb Initialisierung ausserhalb der Schleife
	for (; rr_new > TOL*TOL && k < b.dim() ; ++k) {	// laueft bis die Norm des Residuums kleiner als TOL ist oder im Notfall maximal so oft wie das System gross ist
		matvec(z, A, d);	// zwischenspeichern um es nur einmal auszurechnen
		sp(dz, d, z);		// Skalarprodukt von d_k*z_k
		
		rr_old = rr_new;	// das Skalarprodukt der letzten Iteration - r_(k+1)*r_(k+1) - ist jetzt nur noch r_k*r_k
		
		xold = x;		// x-Vektor wird aktualisiert
		x = d;			// TODO: das hier direkt programmieren um die Zwischenschritte direkt und alles ohne Zwischenspeicherungen zu rechnen
		x.scalmult(rr_old / dz);
		x.vecadd(xold);
		
		rold = r;		// Residuen-Vektor wird aktualisiert
		r = z;			// TODO: dasselbe wie beim x-Vektor, s.o.
		r.scalmult(-(rr_old / dz));
		r.vecadd(rold);
		
		sp(rr_new, r, r);
		
		d.scalmult(rr_new / rr_old);	// Suchrichtung wird aktualisiert
		d.vecadd(r);
	}
	return k;
}

template <typename restype, typename mattype, typename vectype>
int PCG_Jacobi(Vector<restype>& x, DIA<mattype>& A, Vector<vectype>& b) {
	Vector<restype> r(b.dim());				// Definiere Residuumsvektor
	defekt(r, A, b, x);						// Berechne Residuum

	Vector<restype> h(A.dim());		
	Vector<restype> diag(A.dim());
	maindiag(diag,A);						//diag ist jetzt Hauptdiagonalvektor

//***************	kann noch in Vector.hpp	***************
	for (int j=0; j<A.dim(); j++)	
	{
		h[j]=static_cast<restype>(r[j]/diag[j]);	// h=D^-1 * r   mit D=diag(A); also D^-1 Vorkonditionierungsmatrix
	}
//*********************************************************

	Vector<restype> d(h);					// d_0 = h_0
	Vector<restype> z(b.dim());					
	restype rh_old(0);
	restype rh_new(1);
	restype dz(1);
	Vector<restype> x_old(b.dim());
	Vector<restype> r_old(b.dim());

	double TOL(0.0000001);

	sp(rh_old, r, h);						// rh_old= r^T h

	int k(0);
	for (; k < b.dim() && norm(r) >= TOL; ++k) 
	{
		matvec(z, A, d);					// z=A*d_k
		sp(dz, d, z);						// dz= d^T z

		x_old = x;
		x = d;
		x.skalarmult(rh_old / dz);			
		x.vecadd(x_old);					// x_{k+1}=x_k + d*{rh_old/dz}

		r_old = r;
		r = z;
		r.skalarmult(-(rh_old / dz));
		r.vecadd(r_old);					// r_{k+1}=r_k - z*{rh_old/dz} 

		for (int j=0; j<A.dim(); j++) {
			h[j]=r[j]/diag[j];				// h_{k+1}=D^-1 * r_{k+1}
		}

		sp(rh_new, r, h);						
		d.skalarmult(rh_new / rh_old);			
		d.vecadd(h);						// d_{k+1}=h_{k+1} + {rh_new/rh_old} d_k

		rh_old=rh_new;
	}
	return k;
}
