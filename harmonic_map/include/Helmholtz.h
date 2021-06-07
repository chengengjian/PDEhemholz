#pragma once


#include <iostream>
#include <Eigen/Dense>
#include <Eigen/sparse>
#include <Eigen/QR>
#include <Eigen/LU>
#include <vector>
#include <time.h>
#include <math.h>
#include <float.h>
#include "HarmonicMap.h"
#include <complex>

double k = 10.0;
std::complex<double> iye(0, 0.1);


inline bool operator<(const Eigen::Triplet<Eigen::dcomplex>& lhs, const Eigen::Triplet<Eigen::dcomplex>& rhs)
{
	return (lhs.row() < rhs.row() || (lhs.row() == rhs.row() && lhs.col() < rhs.col()));
}

std::complex<double> JFunction(const double n, const std::complex<double>& z) {
	std::complex<double> sum(0.0, 0.0);
	std::complex<double> temp(0.0, 0.0);
	int m = 0;
	while (m <= 100)
	{
		temp = pow(-1, m) / tgamma(n + m + 1) / tgamma(m + 1) * pow(0.5 * z, 2.0 * m + n);
		sum += temp;
		m++;
	}
	return sum;
}

std::complex<double> dJFunction(const double n, const std::complex<double>& z) {
	return (-z * JFunction(n + 1, z) + n * JFunction(n, z)) / z;
}

MeshLib::CPoint Gray2Rainbox(double tmp) {
	MeshLib::CPoint p(0, 0, 0);
	if (tmp <= 51) {
		p[0] = abs(255);
		p[1] = abs(tmp * 5);
		p[2] = abs(0);
	}
	else if (tmp <= 102) {
		tmp -= 51;
		p[0] = abs(255 - tmp * 5);
		p[1] = abs(255);
		p[2] = abs(0);
	}
	else if (tmp <= 153) {
		tmp -= 102;
		p[0] = abs(0);
		p[1] = abs(255);
		p[2] = abs(tmp * 5);
	}
	else if (tmp <= 204) {
		tmp -= 153;
		p[0] = abs(0);
		p[1] = abs(unsigned char(255 - (tmp * 128.0 / 51.0)));
		p[2] = abs(255);
	}
	else {
		tmp -= 204;
		p[0] = abs(0);
		p[1] = abs(unsigned char(127 - (tmp * 128.0 / 51.0)));
		p[2] = abs(255);
	}
	return p;
}

double f( MeshLib::CPoint p) {
	double r = std::sqrt(p[0] * p[0] + p[1] * p[1]);
	return std::sin(k*r) / r;
}

std::complex<double> g(MeshLib::CPoint p, MeshLib::CHarmonicMapHalfEdge* pH) {
	double r = std::sqrt(p[0] * p[0] + p[1] * p[1]);
	std::complex<double> tmp1 = cos(k * r) / k - (std::complex<double>(cos(k), sin(k))) / (k * (JFunction(0, k) + std::complex<double>(0, 1) * JFunction(1, k))) * JFunction(0, k * r);
	tmp1 *= std::complex<double>(0, 1) * k;
	MeshLib::CPoint n = pH->normal();
	std::complex<double> tmp2 = (-sin(k * r) - (std::complex<double>(cos(k), sin(k))) / ((JFunction(0, k) + std::complex<double>(0, 1) * JFunction(1, k))) * dJFunction(0, k * r) ) * ((p / r) * n);
	return tmp1 + tmp2;
}

std::complex<double> realsolve(MeshLib::CPoint p) {
	double r = std::sqrt(p[0] * p[0] + p[1] * p[1]);
	std::complex<double> tmp1 = cos(k * r) / k - (std::complex<double>(cos(k), sin(k))) / (k * (JFunction(0, k) + std::complex<double>(0, 1) * JFunction(1, k))) * JFunction(0, k * r);
	return tmp1;
}

void assemblingMatrix(MeshLib::CHarmonicMapMesh* m_pMesh, Eigen::SparseMatrix<Eigen::dcomplex> &mat) {
	
	mat.setZero();
	std::vector<Eigen::Triplet<Eigen::dcomplex> > triplets;
	
	for (MeshLib::CHarmonicMapMesh::MeshHalfEdgeIterator hiter(m_pMesh); !hiter.end(); hiter++) {
		MeshLib::CHarmonicMapHalfEdge* pHi = *hiter;
		MeshLib::CHarmonicMapHalfEdge* pHj = (MeshLib::CHarmonicMapHalfEdge*)pHi->he_next();
		MeshLib::CHarmonicMapVertex* pVi = (MeshLib::CHarmonicMapVertex*) pHi->source();
		MeshLib::CHarmonicMapVertex* pVj = (MeshLib::CHarmonicMapVertex*) pHj->source();
		if ( ((MeshLib::CHarmonicMapVertex*)pHi->target() )->idx() != ((MeshLib::CHarmonicMapVertex*)pHj->source())->idx()) {
			std::cout << "some thing error" << std::endl;
		}
		MeshLib::CHarmonicMapFace* pF = (MeshLib::CHarmonicMapFace*)pHi->face();
		int i = pVi->idx() ;
		int j = pVj->idx() ;
		double innerProduct = 0;
		innerProduct = pHi->gredient() * pHj->gredient();
		//for (int l = 0; l < 3; l++) {
		//	innerProduct += pHi->geta(l)* pHj->geta(l);
		//}
		
		double valuereal = (innerProduct - k * k * pHi->baseFun(pF->center()) * pHj->baseFun(pF->center())) * pF->area();
		valuereal;
		double valueimage = 0.0;
		if (pHi->edge()->boundary()) {
			valueimage += k * pHi->length() * 0.25;
		}
		else {
			MeshLib::CHarmonicMapHalfEdge*  pHs = ((MeshLib::CHarmonicMapHalfEdge*)(pHi->he_sym()->he_next()));
			double a = (pHi->gredient() - pHs->gredient())* pHi->normal();
			MeshLib::CHarmonicMapHalfEdge* pHjs = ((MeshLib::CHarmonicMapHalfEdge*)(pHi->he_sym()));
			double b = (pHj->gredient() - pHjs->gredient()) * pHj->normal();
			//valueimage += 0.1 * pHi->length() * pHi->length() * a * b;
		}
		Eigen::dcomplex value(valuereal, valueimage);
		triplets.push_back(Eigen::Triplet<Eigen::dcomplex>(i, j, value));
		triplets.push_back(Eigen::Triplet<Eigen::dcomplex>(j, i, value));
		//mat(i, j) += value;
	}
	
	std::sort(triplets.begin(), triplets.end());

	std::vector<Eigen::Triplet<Eigen::dcomplex> > triplets2;
	int i = 0;
	while (i + 1 < triplets.size()) {
		if (triplets[i].row() == triplets[i + 1].row() && triplets[i].col() == triplets[i + 1].col()) {
			triplets2.push_back(Eigen::Triplet < Eigen::dcomplex>(triplets[i].row(), triplets[i].col(), triplets[i].value() + triplets[i+1].value()));
			i += 2;
		}
		else {
			triplets2.push_back(Eigen::Triplet < Eigen::dcomplex>(triplets[i].row(), triplets[i].col(), triplets[i].value() + triplets[i + 1].value()));
			i += 1;
		}
	}

	for (int i = 1; i < triplets2.size(); i++) {
		//std::cout << triplets[i].row() << " " << triplets[i].col() << " " << triplets[i].value() << std::endl;
		if (triplets2[i].row() == triplets2[i - 1].row() && triplets2[i].col() == triplets2[i - 1].col()) {
			exit(111);
			break;
		}
	}

	mat.setFromTriplets(triplets2.begin(), triplets2.end());

}

void assemblingB(MeshLib::CHarmonicMapMesh* m_pMesh, Eigen::VectorXcd& C) {
	
	C.setZero();
	
	for (MeshLib::CHarmonicMapMesh::MeshHalfEdgeIterator hiter(m_pMesh); !hiter.end(); hiter++) {
		MeshLib::CHarmonicMapHalfEdge* pH = *hiter;
		MeshLib::CHarmonicMapVertex* pV = (MeshLib::CHarmonicMapVertex* ) pH->source();
		MeshLib::CHarmonicMapVertex* pVt = (MeshLib::CHarmonicMapVertex*) pH->target();
		MeshLib::CHarmonicMapFace* pF = (MeshLib::CHarmonicMapFace* )pH->face();
		
		int id = pV->idx();
		int idt = pVt->idx();
		//assert(id >= 0);
		if (id<0 || id>=m_pMesh->numVertices()) {
			std::cout << "some thing error " << id << std::endl;
			exit(10);
		}
		C(id) += pF->area() * f(pF->center());

		if (pH->edge()->boundary()) {
			MeshLib::CPoint center = (pH->source()->point() + pH->target()->point()) / 2;
			C(id) += g(center, pH) * pH->baseFun(center) * pH->length();
			C(idt) += g(center, pH) * pH->baseFun(center) * pH->length();
		}
	}
	
}

void calcParameter(MeshLib::CHarmonicMapHalfEdge* pH) {
	MeshLib::CHarmonicMapFace* pF = (MeshLib::CHarmonicMapFace*)pH->face();
	MeshLib::CHarmonicMapVertex* pSourceV = (MeshLib::CHarmonicMapVertex*)pH->source();
	Eigen::MatrixXd A(3, 3);
	Eigen::VectorXd C(3);
	int id = 0;
	for (MeshLib::CHarmonicMapMesh::FaceHalfedgeIterator heiter(pF); !heiter.end(); heiter++) {
		MeshLib::CHarmonicMapHalfEdge* ph = *heiter;
		MeshLib::CHarmonicMapVertex* pV = (MeshLib::CHarmonicMapVertex*)ph->source();
		MeshLib::CPoint p = pV->point();
		
		//std::cout << "p: " <<p[0] << " " << p[1] << " " << p[2] << std::endl;
		
		A(id, 0) = p[0];
		A(id, 1) = p[1];
		A(id, 2) = 1.0;
		if (pV->idx() == pSourceV->idx()) {
			C(id) = 1.0;
		}
		else {
			C(id) = 0.0;
		}
		id++;
	}

	//for (int i = 0; i < 3; i++) {
	//	for (int j = 0; j < 3; j++) {
	//		std::cout << A(i, j) << " " ;
	//	}
	//	std::cout << "\n";
	//}

	Eigen::Vector3d result =  A.lu().solve(C);
	pH->geta() = MeshLib::CPoint(result[0], result[1], result[2]);
	//pH->geta(0) = result[0];
	//pH->geta(1) = result[1];
	//pH->geta(2) = result[2];

	//std::cout << "result: "<<result[0] << " " << result[1] << " " << result[2] << std::endl;
	/*
	x1 y1 z1
	x2 y2 z2   (a,b,c)^T = (1,0,0)^T
	x3 y3 z3
	*/

}

void basesFunc(MeshLib::CHarmonicMapMesh* m_pMesh) {
    
    for (MeshLib::CHarmonicMapMesh::MeshHalfEdgeIterator hiter(m_pMesh); !hiter.end(); hiter++) {
        MeshLib::CHarmonicMapHalfEdge* pH = *hiter;
        calcParameter(pH);
    }

}

void calculate_normal(MeshLib::CHarmonicMapMesh* m_pMesh)
{
	using M = MeshLib::CHarmonicMapMesh;
	for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
	{
		M::CFace* pF = *fiter;
		MeshLib::CPoint p[3];
		M::CHalfEdge* he = m_pMesh->faceHalfedge(pF);

		for (int k = 0; k < 3; k++)
		{
			M::CVertex* pv = m_pMesh->halfedgeTarget(he);
			p[k] = pv->point();
			he = m_pMesh->halfedgeNext(he);
		}

		MeshLib::CPoint fn = (p[1] - p[0]) ^ (p[2] - p[0]);
		pF->area() = fn.norm() / 2.0;
		pF->normal() = fn / fn.norm();
	}
}

void IdxInit(MeshLib::CHarmonicMapMesh* m_pMesh)
{
	using M = MeshLib::CHarmonicMapMesh;
	int i = 0;
	for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter){
		MeshLib::CHarmonicMapVertex* pV = *viter;
		pV->idx() = i;
		i++;
	}
}

void calcRealValue(MeshLib::CHarmonicMapMesh* m_pMesh) {
	using M = MeshLib::CHarmonicMapMesh;
	for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++) {
		MeshLib::CHarmonicMapVertex* pV = *viter;
		double real = realsolve(pV->point()).real();
		double image = realsolve(pV->point()).imag();
		pV->realValue() = std::sqrt(real * real + image * image);
	}
	
}

void solveHelmholz(MeshLib::CHarmonicMapMesh* m_pMesh) {

	calculate_normal(m_pMesh);
	IdxInit(m_pMesh);

	basesFunc(m_pMesh);
	

	int numVertice = m_pMesh->numVertices();

	Eigen::SparseMatrix<Eigen::dcomplex> A(numVertice,numVertice);
	//Eigen::MatrixXcd A(numVertice,numVertice);
	Eigen::VectorXcd B(numVertice);


	double max = 0.0;
	double min = 10000000.0;

	if(false){
		assemblingB(m_pMesh, B);
		assemblingMatrix(m_pMesh, A);
		Eigen::SparseLU<Eigen::SparseMatrix<Eigen::dcomplex>, Eigen::COLAMDOrdering<int>> solver;
		solver.compute(A);
		Eigen::VectorXcd X = solver.solve(B);

		for (int i = 0; i < numVertice; i++) {
			double value = std::sqrt(X(i).real() * X(i).real() + X(i).imag() * X(i).imag());
			max = value > max ? value : max;
			min = value < min ? value : min;
		}
		for (MeshLib::CHarmonicMapMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++) {
			MeshLib::CHarmonicMapVertex* pV = *viter;
			int id = pV->idx();
			double value = std::sqrt(X(id).real() * X(id).real() + X(id).imag() * X(id).imag());
			//pV->point()[2] =  100 * value;
			//pV->rgb() = MeshLib::CPoint(value * 40.0, 0, 0);
			//pV->rgb() = Gray2Rainbox(value);
			//pV->rgb() = Gray2Rainbox(value);
			pV->rgb() = Gray2Rainbox( (value-min) / (max-min) * 512);// MeshLib::CPoint(value / max * 255, 0, 0);
		}
	}
	else {
		calcRealValue(m_pMesh);
		for (MeshLib::CHarmonicMapMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++) {
			MeshLib::CHarmonicMapVertex* pV = *viter;
			//int id = pV->idx();
			double value = pV->realValue();
			max = value > max ? value : max;
			min = value < min ? value : min;
		}

		for (MeshLib::CHarmonicMapMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++) {
			MeshLib::CHarmonicMapVertex* pV = *viter;
			int id = pV->idx();
			//double value = std::sqrt(B(id).real() * B(id).real() + B(id).imag() * B(id).imag());
			pV->point()[2] = pV->realValue();
			//pV->rgb() = MeshLib::CPoint(value * 40.0, 0, 0);
			//pV->rgb() = Gray2Rainbox(pV->realValue());
			pV->rgb() = Gray2Rainbox( (pV->realValue() -min) / (max-min) * 2048);// MeshLib::CPoint(value / max * 255, 0, 0);
		}
	}
	
    //Eigen::SparseMatrix<std::complex<double>> A(10, 10);
    //std::vector<Eigen::Triplet<std::complex<double>>> A_coefficients;
    //for (int i = 0; i < 10; i++) {
        //A_coefficients.push_back(Eigen::Triplet<std::complex<double>>(i, i, std::complex<double>(1, 1)));
        //C(i) = 1;
    //}
    //A.setZero();
    //A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
    //A.makeCompressed();
    //Eigen::SparseQR< Eigen::SparseMatrix<std::complex<double>>, Eigen::COLAMDOrdering<int>> solver;
    //solver.compute(A);
    //Eigen::VectorXcd t = solver.solve(C);


	//for (MeshLib::CHarmonicMapMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++) {
	//	MeshLib::CHarmonicMapVertex* pV = *viter;
	//	int id =  pV->idx();
	//	double value = std::sqrt(B(id).real() * B(id).real() + B(id).imag() * B(id).imag());
	//	//pV->point()[2] = 100 * value;
	//	//pV->rgb() = MeshLib::CPoint(value * 40.0, 0, 0);
	//	pV->rgb() = Gray2Rainbox(pV->realValue());
	//	//pV->rgb() = Gray2Rainbox( (value-min) / (max-min) * 512);// MeshLib::CPoint(value / max * 255, 0, 0);
	//}



    //for (MeshLib::CHarmonicMapMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++) {
    //    MeshLib::CVertex* pV = *viter;
    //    MeshLib::CPoint v = pV->point();
    //    std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
    //}
   
    //exit(11);

}


