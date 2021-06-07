#ifndef _HARMONIC_MAP_MESH_H_
#define _HARMONIC_MAP_MESH_H_

#include <Eigen/Dense>
#include <Eigen/sparse>
#include <Eigen/QR>
#include <Eigen/LU>

#include "Mesh/BaseMesh.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Vertex.h"

#include "Mesh/Boundary.h"
#include "Mesh/Iterators.h"
#include "Parser/parser.h"

namespace MeshLib
{
class CHarmonicMapVertex;
class CHarmonicMapEdge;
class CHarmonicMapFace;
class CHarmonicMapHalfEdge;

/*! \brief CHarmonicMapVertex class
 *
 *   Vertex class for cut graph algoritm
 *   Trait : Vertex valence
 */
class CHarmonicMapVertex : public CVertex
{
  public:
    /*! Constructor */
    CHarmonicMapVertex() : 
        m_index(0), 
        m_rgb(229.0 / 255.0, 162.0 / 255.0, 141.0 / 255.0){};
    
    /*! Vertex index */
    int& idx() { return m_index; };


    /*! Vertex color */
    CPoint& rgb() { return m_rgb; };

    double& realValue() {
        return m_realValue;
    }
    /*!
     *	Read vertex traits to vertex string
     */
    void _from_string();



  protected:
    /*! Vertex index */
    int m_index;

    double m_realValue;
    /*! Vertex color */
    CPoint m_rgb;
};

inline void CHarmonicMapVertex::_from_string()
{
    CParser parser(m_string);
    for (std::list<CToken*>::iterator iter = parser.tokens().begin(); 
         iter != parser.tokens().end(); ++iter)
    {
        CToken* token = *iter;

        if (token->m_key == "rgb")
        {
            token->m_value >> m_rgb;
        }
    }
}


/*! \brief CHarmonicMapEdge class
 *
 *   Edge class for cut graph algorithm
 *   Trait : Edge sharp
 */
class CHarmonicMapEdge : public CEdge
{
  public:
    /*! Constructor */
    CHarmonicMapEdge() : m_length(0), m_weight(0) {};

    /*!	Edge weight */
    double& weight() { return m_weight; };
    
    /*! Edge length */
    double& length() { return m_length; };
    


  protected:
    /*!	Edge weight */
    double m_weight;



    /*! Edge length */
    double m_length;
};

/*! \brief CHarmonicMapFace class
 *
 *   Face class for cut graph algorithm
 *   Trait : Face touched flag
 */
class CHarmonicMapFace : public CFace
{
  public:
    /*! Constructor */
    CHarmonicMapFace() {};

    MeshLib::CPoint center()
    {
        MeshLib::CPoint p(0, 0, 0);
        CHalfEdge* pH = m_halfedge;
        for (int i = 0; i < 3; i++) {
            CVertex* pV =  pH->source();
            p += pV->point();
            pH = pH->he_next();
        }
        p /= 3;
        return p;
    }

	//double area() {
	//	Eigen::MatrixXd det(3, 3);
	//	
	//	CHarmonicMapVertex* pV = (CHarmonicMapVertex*)m_halfedge[0].source();
	//	CHarmonicMapVertex* pV1 = (CHarmonicMapVertex*)m_halfedge[1].source();
	//	CHarmonicMapVertex* pV2 = (CHarmonicMapVertex*)m_halfedge[2].source();
	//	
	//	//CHarmonicMapVertex* pV1 = halfedge[1].source();
	//	//CHarmonicMapVertex* pV2 = halfedge[2].source();
	//	det(0, 0) = 1; det(0, 1) = 1; det(0, 2) = 1;
	//	det(1, 0) = pV1->point()[0] - pV->point()[0]; det(1, 1) = pV1->point()[1] - pV->point()[1]; det(1, 2) = pV1->point()[2] - pV->point()[2];
	//	det(2, 0) = pV2->point()[0] - pV->point()[0]; det(2, 1) = pV2->point()[1] - pV->point()[1]; det(2, 2) = pV2->point()[2] - pV->point()[2];
	//	return 0.5 * det.determinant();
	//}
    double& area() {
        return m_area;
    }

    /*! face normal */
    CPoint& normal() { return m_normal; };
  protected:
    /*! face normal */
    CPoint m_normal;

    double m_area;
};

/*! \brief CHarmonicMapHalfEdge class
 *
 *   HalfEdge class for cut graph algorithm
 */
class CHarmonicMapHalfEdge : public CHalfEdge
{
  public:
    /*!	CHarmonicHalfEdge constructor */
    CHarmonicMapHalfEdge() : m_angle(0) {};
    
    /*!	Corner angle */
    double& angle() { return m_angle; };

    CPoint normal() {
        CHarmonicMapVertex* pV1 = (CHarmonicMapVertex * )target();
        CHarmonicMapVertex* pV2 = (CHarmonicMapVertex * )source();
        CPoint n = pV1->point() - pV2->point();
        n = n ^ CPoint(0, 0, -1);
        //std::cout << pV1->point() * (n / n.norm()) << std::endl;
        return n / n.norm();
    }

   // double& geta(int x) {
       // return a[x];
    //}
    CPoint& geta() {
        return a;
    }

    CPoint& gredient() {
        return CPoint(a(0), a(1), 0);
    }

    MeshLib::CPoint center() {
        return (this->source()->point() + this->target()->point()) / 2.0;
    }

    double length() {
        MeshLib::CHarmonicMapVertex* pV1 = (MeshLib::CHarmonicMapVertex*)this->source();
        MeshLib::CHarmonicMapVertex* pV2 = (MeshLib::CHarmonicMapVertex*)this->target();
        return (pV1->point() - pV2->point()).norm();
    }

    double baseFun(MeshLib::CPoint p) {
        return a(0) * p(0) + a(1) * p(1) + a(3);
    }

  protected:
    /*! Corner angle */
    double m_angle;

    //double a[3];
    CPoint a;
};

/*! \brief CHarmonicMapMesh class
 *
 *	Mesh class for cut graph algorithm
 *
 */
template <typename V, typename E, typename F, typename H>
class THarmonicMapMesh : public CBaseMesh<V, E, F, H>
{
  public:
    typedef V CVertex;
    typedef E CEdge;
    typedef F CFace;
    typedef H CHalfEdge;

    typedef CBoundary<V, E, F, H>                   CBoundary;
    typedef CLoop<V, E, F, H>                       CLoop;

    typedef MeshVertexIterator<V, E, F, H>          MeshVertexIterator;
    typedef MeshEdgeIterator<V, E, F, H>            MeshEdgeIterator;
    typedef MeshFaceIterator<V, E, F, H>            MeshFaceIterator;
    typedef MeshHalfEdgeIterator<V, E, F, H>        MeshHalfEdgeIterator;

    typedef VertexVertexIterator<V, E, F, H>        VertexVertexIterator;
    typedef VertexEdgeIterator<V, E, F, H>          VertexEdgeIterator;
    typedef VertexFaceIterator<V, E, F, H>          VertexFaceIterator;
    typedef VertexInHalfedgeIterator<V, E, F, H>    VertexInHalfedgeIterator;
    typedef VertexOutHalfedgeIterator<V, E, F, H>   VertexOutHalfedgeIterator;

    typedef FaceVertexIterator<V, E, F, H>          FaceVertexIterator;
    typedef FaceEdgeIterator<V, E, F, H>            FaceEdgeIterator;
    typedef FaceHalfedgeIterator<V, E, F, H>        FaceHalfedgeIterator;
};

typedef THarmonicMapMesh<CHarmonicMapVertex, CHarmonicMapEdge, CHarmonicMapFace, CHarmonicMapHalfEdge> CHarmonicMapMesh;
}

#endif // !_HARMONIC_MAP_MESH_H_
