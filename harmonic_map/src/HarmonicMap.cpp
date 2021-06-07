#include <math.h>
#include <float.h>

#include <Eigen/Sparse>

#include "HarmonicMap.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279
#endif

void MeshLib::CHarmonicMap::set_mesh(CHarmonicMapMesh* pMesh)
{
    m_pMesh = pMesh;
    
    // 1. compute the weights of edges
    _calculate_edge_weight();

    // 2. map the boundary to unit circle
    _set_boundary();

    // 3. initialize the map of interior vertices to (0, 0)
    using M = CHarmonicMapMesh;
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;
		if (pV->boundary()) {
			//printf("边界点：( %f:::::%f )\n", pV->uv()[0], pV->uv()[1]);
			continue;
		}
		CPoint2 suv(0, 0);
		pV->uv() = suv;
		//printf("内点：( %f:::::%f )\n", pV->uv()[0], pV->uv()[1]);

    }
}

double MeshLib::CHarmonicMap::step_one() 
{
    if (!m_pMesh)
    {
        std::cerr << "Should set mesh first!" << std::endl;
        return DBL_MAX;
    }

    using M = CHarmonicMapMesh;

    // move each interior vertex to its weighted center of neighbors
    double max_error = -DBL_MAX;
	int i = 0;
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;
        if (pV->boundary())
            continue;

        double sw = 0;
        CPoint2 suv(0, 0);
		//printf("chen %f,:::::,%f \n", suv[0], suv[1]);
        for (M::VertexVertexIterator vviter(pV); !vviter.end(); vviter++)
        {
            M::CVertex* pW = *vviter;
			//printf("pW当前UV（ %f:::::%f） \n", pW->uv()[0], pW->uv()[1]);
            M::CEdge* pE = m_pMesh->vertexEdge(pV, pW);
			//insert your code here
			sw += pE->weight();
			suv += pW->uv()* pE->weight();
        }
				
		//printf("epoch:%d \n", i);
		//printf("SW:%f \n", sw);
		//printf("pV当前UV（ %f:::::%f） \n", pV->uv()[0], pV->uv()[1]);
		//printf("\n");
		i++;
		

        suv /= sw;
		
		//printf("chen %f,:::::,%f \n", pV->uv()[0], pV->uv()[1]);
        double error = (pV->uv() - suv).norm2();
		
        max_error = (error > max_error) ? error : max_error;
        pV->uv() = suv;
    }

    printf("Current max error is %g\n", max_error);
	//
    return max_error;
}

void MeshLib::CHarmonicMap::iterative_map(double epsilon)
{
    if (!m_pMesh)
    {
        std::cerr << "Should set mesh first!" << std::endl;
        return;
    }

    using M = CHarmonicMapMesh;

    // take steps until it converges.
    while (true)
    {
		printf("ok");
        double error = this->step_one();
		if (error < 1e-10)
			//printf("%f", error);
            break;
    }
}

void MeshLib::CHarmonicMap::map() 
{
    if (!m_pMesh)
    {
        std::cerr << "Should set mesh first!" << std::endl;
        return;
    }

    using M = CHarmonicMapMesh;

    // 1. Initialize
    int vid = 0;  // interior vertex id
    int bid = 0;  // boundary vertex id
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;

        if (pV->boundary())
            pV->idx() = bid++;
        else
            pV->idx() = vid++;
    }

    int interior_vertices = vid;
    int boundary_vertices = bid;

    // 2. Set the matrix A and B
    std::vector<Eigen::Triplet<double>> A_coefficients;
    std::vector<Eigen::Triplet<double>> B_coefficients;

    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;
        if (pV->boundary())
            continue;
        int vid = pV->idx();

        double sw = 0;
        for (M::VertexVertexIterator witer(pV); !witer.end(); ++witer)
        {
            M::CVertex* pW = *witer;
            int wid = pW->idx();

            M::CEdge* e = m_pMesh->vertexEdge(pV, pW);
            double w = e->weight();
			sw += w;
			//insert your code here
			Eigen::Triplet<double> A(vid, wid, w);

			A_coefficients.push_back(A);
			B_coefficients.push_back(A);
			//construct one element of the matrix A and B, using 
			//Eigen::Triplet<double>(i,j, val)
			//there push_back the triplet to A or B coefficients
        }
		//A_coefficients[vid, vid, 1] = -sw;
    }

    Eigen::SparseMatrix<double> A(interior_vertices, interior_vertices);
    A.setZero();
    Eigen::SparseMatrix<double> B(interior_vertices, boundary_vertices);
    B.setZero();
    A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
    B.setFromTriplets(B_coefficients.begin(), B_coefficients.end());

    // 3. Solve the equations
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    std::cerr << "Eigen Decomposition" << std::endl;
    solver.compute(A);
    std::cerr << "Eigen Decomposition Finished" << std::endl;

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Waring: Eigen decomposition failed" << std::endl;
    }

    for (int k = 0; k < 2; k++)
    {
        Eigen::VectorXd b(boundary_vertices);
        // set boundary constraints vector b
        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
        {
            M::CVertex* pV = *viter;
            if (!pV->boundary())
                continue;
            int id = pV->idx();
            b(id) = pV->uv()[k];
        }

        Eigen::VectorXd c(interior_vertices);
        c = B * b;

        Eigen::VectorXd x = solver.solve(c); // Ax=c
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "Waring: Eigen decomposition failed" << std::endl;
        }

        // set the images of the harmonic map to interior vertices
        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
        {
            M::CVertex* pV = *viter;
            if (pV->boundary())
                continue;
            int id = pV->idx();
            pV->uv()[k] = x(id);
        }
    }
}

void MeshLib::CHarmonicMap::_calculate_edge_weight() 
{
    using M = CHarmonicMapMesh;

    // 1. compute edge length
    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter)
    {
        M::CEdge* pE = *eiter;
        M::CVertex* v1 = m_pMesh->edgeVertex1(pE);
        M::CVertex* v2 = m_pMesh->edgeVertex2(pE);
        pE->length() = (v1->point() - v2->point()).norm();
    }

    // 2. compute corner angle
    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
    {
        M::CFace* pF = *fiter;
        M::CHalfEdge* pH[3];
		//insert your code here
		double Elength[3];
		int i = 0;
		for (M::FaceHalfedgeIterator fhiter(pF); !fhiter.end(); ++fhiter) {
			CHarmonicMapHalfEdge* pHe = *fhiter;
			CHarmonicMapEdge* pSymE = m_pMesh->halfedgeEdge(pHe);
			//Elength[i] = m_pMesh->halfedgeEdge(pH)->length();
			pH[i] = pHe;
			Elength[i] = pSymE->length();
			//printf("EdgeLength : %f \n", Elength[i]);
			i++;
		}
		//double l0 = Elength[0];
		//double l1 = Elength[1];
		//double l2 = Elength[2];
		//printf("l0 = %f \n", l0);
		//printf("l1 = %f \n", l1);
		//printf("l2 = %f \n", l2);
		double a = (Elength[1] * Elength[1] + Elength[2] * Elength[2] - Elength[0] * Elength[0]) / (2 * Elength[1] * Elength[2]);
		double b = (Elength[0] * Elength[0] + Elength[2] * Elength[2] - Elength[1] * Elength[1]) / (2 * Elength[0] * Elength[2]);
		double c = (Elength[1] * Elength[1] + Elength[0] * Elength[0] - Elength[2] * Elength[2]) / (2 * Elength[1] * Elength[0]);
		//printf("a = %f \n", a);
		//printf("b = %f \n", b);
		//printf("c = %f \n", c);
		//printf("\n");
		pH[0]->angle() = acos(a);
		pH[1]->angle() = acos(b);
		pH[2]->angle() = acos(c);


		//use inverse cosine law to compute the corner angles
    }

    // 3. compute edge weight
    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter)
    {
        CHarmonicMapEdge* pE = *eiter;
		CHarmonicMapHalfEdge* pH0 = (CHarmonicMapHalfEdge*)pE->halfedge(0);
		CHarmonicMapHalfEdge* pH1 = (CHarmonicMapHalfEdge*)pE->halfedge(1);

		if (pE->boundary()) {
			if (pH0 != NULL) {
				//printf("pH0->angle():%f \n", pH0->angle());
				pE->weight() = 1 / tan(pH0->angle());
				//printf("pE->weight():%f \n",pE->weight());
			}
			else {
				//printf("pH1->angle():%f \n", pH1->angle());
				pE->weight() = 1 / tan(pH1->angle());
				//printf("pE->weight():%f \n", pE->weight());
			}
			
		}
		else {
			pE->weight() = 1 / tan(pH0->angle()) + 1 / tan(pH1->angle());
		}

		//insert your code here
		//set cotangent edge weight
    }
}

void MeshLib::CHarmonicMap::_set_boundary() 
{
    using M = CHarmonicMapMesh;

    // 1. get the boundary half edge loop
    M::CBoundary boundary(m_pMesh);
    std::vector<M::CLoop*>& pLs = boundary.loops();
    if (pLs.size() != 1)
    {
        std::cerr << "Only topological disk accepted!" << std::endl;
        exit(EXIT_FAILURE);
    }
    M::CLoop* pL = pLs[0];
    std::list<M::CHalfEdge*>& pHs = pL->halfedges();
    
    // 2. compute the total length of the boundary
    double sum = 0.0;
    std::list<M::CHalfEdge*>::iterator it;
    for (it = pHs.begin(); it != pHs.end(); ++it)
    {
        M::CHalfEdge* pH = *it;
        sum += m_pMesh->halfedgeEdge(pH)->length();
    }

    // 3. parameterize the boundary using arc length parameter
    double len = 0.0;
    for (it = pHs.begin(); it != pHs.end(); ++it)
    {
        M::CHalfEdge* pH = *it;
        M::CVertex* pV = m_pMesh->halfedgeVertex(pH);

        len += m_pMesh->halfedgeEdge(pH)->length();
        double angle = len / sum * 2.0 * M_PI;
        pV->uv() = CPoint2(cos(angle), sin(angle)); 
    }
}

double MeshLib::CHarmonicMap::_inverse_cosine_law(double a, double b, double c) 
{ 
    double cs = (a * a + b * b - c * c) / (2.0 * a * b);
    assert(cs <= 1.0 && cs >= -1.0);
    return std::acos(cs);
}
