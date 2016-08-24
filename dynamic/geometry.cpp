// This file has been generated by Py++.

#include "boost/python.hpp"

#include "indexing_suite/value_traits.hpp"

#include "indexing_suite/container_suite.hpp"

#include "indexing_suite/vector.hpp"

#include "geometry_headers.hpp"

namespace bp = boost::python;

namespace boost { namespace python { namespace indexing {

template<>
struct value_traits< boost::numeric::ublas::c_vector< double, 3 > >{

    static bool const equality_comparable = false;
    

    static bool const less_than_comparable = false;
    

    template<typename PythonClass, typename Policy>
    static void visit_container_class(PythonClass &, Policy const &){
        
    }

};

}/*indexing*/ } /*python*/ } /*boost*/

namespace boost { namespace python { namespace indexing {

template<>
struct value_traits< boost::shared_ptr< Facet > >{

    static bool const equality_comparable = false;
    

    static bool const less_than_comparable = false;
    

    template<typename PythonClass, typename Policy>
    static void visit_container_class(PythonClass &, Policy const &){
        
    }

};

}/*indexing*/ } /*python*/ } /*boost*/

namespace boost { namespace python { namespace indexing {

template<>
struct value_traits< boost::shared_ptr< Polygon > >{

    static bool const equality_comparable = false;
    

    static bool const less_than_comparable = false;
    

    template<typename PythonClass, typename Policy>
    static void visit_container_class(PythonClass &, Policy const &){
        
    }

};

}/*indexing*/ } /*python*/ } /*boost*/

namespace boost { namespace python { namespace indexing {

template<>
struct value_traits< boost::shared_ptr< Vertex > >{

    static bool const equality_comparable = false;
    

    static bool const less_than_comparable = false;
    

    template<typename PythonClass, typename Policy>
    static void visit_container_class(PythonClass &, Policy const &){
        
    }

};

}/*indexing*/ } /*python*/ } /*boost*/

BOOST_PYTHON_MODULE(_chaste_project_Angiogenesis_geometry){
    bp::class_< std::vector< unsigned int > >("vector_less__unsigned_int__greater_")    
        .def( bp::indexing::vector_suite< std::vector< unsigned int > >() );

    bp::class_< std::vector< std::pair<unsigned int, unsigned int> > >("vector_less__std_scope_pair_less_unsigned_int_comma__unsigned_int_greater___greater_")    
        .def( bp::indexing::vector_suite< std::vector< std::pair<unsigned int, unsigned int> > >() );

    { //::std::vector< boost::shared_ptr<Vertex> >
        typedef bp::class_< std::vector< boost::shared_ptr<Vertex> > > vector_less__boost_scope_shared_ptr_less_Vertex_greater___greater__exposer_t;
        vector_less__boost_scope_shared_ptr_less_Vertex_greater___greater__exposer_t vector_less__boost_scope_shared_ptr_less_Vertex_greater___greater__exposer = vector_less__boost_scope_shared_ptr_less_Vertex_greater___greater__exposer_t( "vector_less__boost_scope_shared_ptr_less_Vertex_greater___greater_" );
        bp::scope vector_less__boost_scope_shared_ptr_less_Vertex_greater___greater__scope( vector_less__boost_scope_shared_ptr_less_Vertex_greater___greater__exposer );
        vector_less__boost_scope_shared_ptr_less_Vertex_greater___greater__exposer.def( bp::indexing::vector_suite< std::vector< boost::shared_ptr<Vertex> > >() );
    }

    { //::std::vector< boost::shared_ptr<Polygon> >
        typedef bp::class_< std::vector< boost::shared_ptr<Polygon> > > vector_less__boost_scope_shared_ptr_less_Polygon_greater___greater__exposer_t;
        vector_less__boost_scope_shared_ptr_less_Polygon_greater___greater__exposer_t vector_less__boost_scope_shared_ptr_less_Polygon_greater___greater__exposer = vector_less__boost_scope_shared_ptr_less_Polygon_greater___greater__exposer_t( "vector_less__boost_scope_shared_ptr_less_Polygon_greater___greater_" );
        bp::scope vector_less__boost_scope_shared_ptr_less_Polygon_greater___greater__scope( vector_less__boost_scope_shared_ptr_less_Polygon_greater___greater__exposer );
        vector_less__boost_scope_shared_ptr_less_Polygon_greater___greater__exposer.def( bp::indexing::vector_suite< std::vector< boost::shared_ptr<Polygon> > >() );
    }

    { //::std::vector< boost::shared_ptr<Facet> >
        typedef bp::class_< std::vector< boost::shared_ptr<Facet> > > vector_less__boost_scope_shared_ptr_less_Facet_greater___greater__exposer_t;
        vector_less__boost_scope_shared_ptr_less_Facet_greater___greater__exposer_t vector_less__boost_scope_shared_ptr_less_Facet_greater___greater__exposer = vector_less__boost_scope_shared_ptr_less_Facet_greater___greater__exposer_t( "vector_less__boost_scope_shared_ptr_less_Facet_greater___greater_" );
        bp::scope vector_less__boost_scope_shared_ptr_less_Facet_greater___greater__scope( vector_less__boost_scope_shared_ptr_less_Facet_greater___greater__exposer );
        vector_less__boost_scope_shared_ptr_less_Facet_greater___greater__exposer.def( bp::indexing::vector_suite< std::vector< boost::shared_ptr<Facet> > >() );
    }

    { //::std::vector< boost::numeric::ublas::c_vector<double, 3> >
        typedef bp::class_< std::vector< boost::numeric::ublas::c_vector<double, 3> > > __type_exposer_t;
        __type_exposer_t __type_exposer = __type_exposer_t( "__type" );
        bp::scope __type_scope( __type_exposer );
        __type_exposer.def( bp::indexing::vector_suite< std::vector< boost::numeric::ublas::c_vector<double, 3> > >() );
    }

    { //::Facet
        typedef bp::class_< Facet > Facet_exposer_t;
        Facet_exposer_t Facet_exposer = Facet_exposer_t( "Facet", bp::init< std::vector< boost::shared_ptr<Polygon> > >(( bp::arg("polygons") )) );
        bp::scope Facet_scope( Facet_exposer );
        bp::implicitly_convertible< std::vector< boost::shared_ptr<Polygon> >, Facet >();
        Facet_exposer.def( bp::init< boost::shared_ptr< Polygon > >(( bp::arg("pPolygon") )) );
        bp::implicitly_convertible< boost::shared_ptr< Polygon >, Facet >();
        { //::Facet::AddPolygon
        
            typedef void ( ::Facet::*AddPolygon_function_type)( ::boost::shared_ptr< Polygon > ) ;
            
            Facet_exposer.def( 
                "AddPolygon"
                , AddPolygon_function_type( &::Facet::AddPolygon )
                , ( bp::arg("pPolygon") ) );
        
        }
        { //::Facet::AddPolygons
        
            typedef void ( ::Facet::*AddPolygons_function_type)( ::std::vector< boost::shared_ptr<Polygon> > ) ;
            
            Facet_exposer.def( 
                "AddPolygons"
                , AddPolygons_function_type( &::Facet::AddPolygons )
                , ( bp::arg("polygons") ) );
        
        }
        { //::Facet::ContainsPoint
        
            typedef bool ( ::Facet::*ContainsPoint_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Facet_exposer.def( 
                "ContainsPoint"
                , ContainsPoint_function_type( &::Facet::ContainsPoint )
                , ( bp::arg("location") ) );
        
        }
        { //::Facet::Create
        
            typedef ::boost::shared_ptr< Facet > ( *Create_function_type )( ::std::vector< boost::shared_ptr<Polygon> > );
            
            Facet_exposer.def( 
                "Create"
                , Create_function_type( &::Facet::Create )
                , ( bp::arg("polygons") ) );
        
        }
        { //::Facet::Create
        
            typedef ::boost::shared_ptr< Facet > ( *Create_function_type )( ::boost::shared_ptr< Polygon > );
            
            Facet_exposer.def( 
                "Create"
                , Create_function_type( &::Facet::Create )
                , ( bp::arg("pPolygon") ) );
        
        }
        { //::Facet::GetBoundingBox
        
            typedef ::boost::numeric::ublas::c_vector< double, 6 > ( ::Facet::*GetBoundingBox_function_type)(  ) ;
            
            Facet_exposer.def( 
                "GetBoundingBox"
                , GetBoundingBox_function_type( &::Facet::GetBoundingBox ) );
        
        }
        { //::Facet::GetCentroid
        
            typedef ::boost::numeric::ublas::c_vector< double, 3 > ( ::Facet::*GetCentroid_function_type)(  ) ;
            
            Facet_exposer.def( 
                "GetCentroid"
                , GetCentroid_function_type( &::Facet::GetCentroid ) );
        
        }
        { //::Facet::GetData
        
            typedef double ( ::Facet::*GetData_function_type)( ::std::string const & ) ;
            
            Facet_exposer.def( 
                "GetData"
                , GetData_function_type( &::Facet::GetData )
                , ( bp::arg("label") ) );
        
        }
        { //::Facet::GetDistance
        
            typedef double ( ::Facet::*GetDistance_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Facet_exposer.def( 
                "GetDistance"
                , GetDistance_function_type( &::Facet::GetDistance )
                , ( bp::arg("location") ) );
        
        }
        { //::Facet::GetNormal
        
            typedef ::boost::numeric::ublas::c_vector< double, 3 > ( ::Facet::*GetNormal_function_type)(  ) ;
            
            Facet_exposer.def( 
                "GetNormal"
                , GetNormal_function_type( &::Facet::GetNormal ) );
        
        }
        { //::Facet::GetPlane
        
            typedef ::vtkSmartPointer< vtkPlane > ( ::Facet::*GetPlane_function_type)(  ) ;
            
            Facet_exposer.def( 
                "GetPlane"
                , GetPlane_function_type( &::Facet::GetPlane ) );
        
        }
        { //::Facet::GetPolygons
        
            typedef ::std::vector< boost::shared_ptr<Polygon> > ( ::Facet::*GetPolygons_function_type)(  ) ;
            
            Facet_exposer.def( 
                "GetPolygons"
                , GetPolygons_function_type( &::Facet::GetPolygons ) );
        
        }
        { //::Facet::GetVertices
        
            typedef ::std::vector< boost::shared_ptr<Vertex> > ( ::Facet::*GetVertices_function_type)(  ) ;
            
            Facet_exposer.def( 
                "GetVertices"
                , GetVertices_function_type( &::Facet::GetVertices ) );
        
        }
        { //::Facet::GetVtkVertices
        
            typedef ::std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkIdTypeArray > > ( ::Facet::*GetVtkVertices_function_type)(  ) ;
            
            Facet_exposer.def( 
                "GetVtkVertices"
                , GetVtkVertices_function_type( &::Facet::GetVtkVertices ) );
        
        }
        { //::Facet::RotateAboutAxis
        
            typedef void ( ::Facet::*RotateAboutAxis_function_type)( ::boost::numeric::ublas::c_vector< double, 3 >,double ) ;
            
            Facet_exposer.def( 
                "RotateAboutAxis"
                , RotateAboutAxis_function_type( &::Facet::RotateAboutAxis )
                , ( bp::arg("axis"), bp::arg("angle") ) );
        
        }
        { //::Facet::SetData
        
            typedef void ( ::Facet::*SetData_function_type)( ::std::string const &,double ) ;
            
            Facet_exposer.def( 
                "SetData"
                , SetData_function_type( &::Facet::SetData )
                , ( bp::arg("label"), bp::arg("value") ) );
        
        }
        { //::Facet::Translate
        
            typedef void ( ::Facet::*Translate_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Facet_exposer.def( 
                "Translate"
                , Translate_function_type( &::Facet::Translate )
                , ( bp::arg("translationVector") ) );
        
        }
        { //::Facet::UpdateVertices
        
            typedef void ( ::Facet::*UpdateVertices_function_type)(  ) ;
            
            Facet_exposer.def( 
                "UpdateVertices"
                , UpdateVertices_function_type( &::Facet::UpdateVertices ) );
        
        }
        Facet_exposer.staticmethod( "Create" );
        bp::register_ptr_to_python< boost::shared_ptr< Facet > >();
    }

    { //::Part< 3 >
        typedef bp::class_< Part< 3 > > Part3_exposer_t;
        Part3_exposer_t Part3_exposer = Part3_exposer_t( "Part3", bp::init< >() );
        bp::scope Part3_scope( Part3_exposer );
        { //::Part< 3 >::AddCircle
        
            typedef Part< 3 > exported_class_t;
            typedef ::boost::shared_ptr< Polygon > ( exported_class_t::*AddCircle_function_type)( double,::boost::numeric::ublas::c_vector< double, 3 >,unsigned int ) ;
            
            Part3_exposer.def( 
                "AddCircle"
                , AddCircle_function_type( &::Part< 3 >::AddCircle )
                , ( bp::arg("radius"), bp::arg("centre"), bp::arg("numSegments") ) );
        
        }
        { //::Part< 3 >::AddCuboid
        
            typedef Part< 3 > exported_class_t;
            typedef void ( exported_class_t::*AddCuboid_function_type)( double,double,double,::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Part3_exposer.def( 
                "AddCuboid"
                , AddCuboid_function_type( &::Part< 3 >::AddCuboid )
                , ( bp::arg("sizeX"), bp::arg("sizeY"), bp::arg("sizeZ"), bp::arg("origin") ) );
        
        }
        { //::Part< 3 >::AddCylinder
        
            typedef Part< 3 > exported_class_t;
            typedef void ( exported_class_t::*AddCylinder_function_type)( double,double,::boost::numeric::ublas::c_vector< double, 3 >,unsigned int ) ;
            
            Part3_exposer.def( 
                "AddCylinder"
                , AddCylinder_function_type( &::Part< 3 >::AddCylinder )
                , ( bp::arg("radius"), bp::arg("depth"), bp::arg("centre"), bp::arg("numSegments") ) );
        
        }
        { //::Part< 3 >::AddHoleMarker
        
            typedef Part< 3 > exported_class_t;
            typedef void ( exported_class_t::*AddHoleMarker_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Part3_exposer.def( 
                "AddHoleMarker"
                , AddHoleMarker_function_type( &::Part< 3 >::AddHoleMarker )
                , ( bp::arg("location") ) );
        
        }
        { //::Part< 3 >::AddPolygon
        
            typedef Part< 3 > exported_class_t;
            typedef ::boost::shared_ptr< Polygon > ( exported_class_t::*AddPolygon_function_type)( ::std::vector< boost::shared_ptr<Vertex> >,bool,::boost::shared_ptr< Facet > ) ;
            
            Part3_exposer.def( 
                "AddPolygon"
                , AddPolygon_function_type( &::Part< 3 >::AddPolygon )
                , ( bp::arg("vertices"), bp::arg("newFacet"), bp::arg("pFacet") ) );
        
        }
        { //::Part< 3 >::AddPolygon
        
            typedef Part< 3 > exported_class_t;
            typedef ::boost::shared_ptr< Polygon > ( exported_class_t::*AddPolygon_function_type)( ::boost::shared_ptr< Polygon >,bool,::boost::shared_ptr< Facet > ) ;
            
            Part3_exposer.def( 
                "AddPolygon"
                , AddPolygon_function_type( &::Part< 3 >::AddPolygon )
                , ( bp::arg("pPolygon"), bp::arg("newFacet"), bp::arg("pFacet") ) );
        
        }
        { //::Part< 3 >::AddRectangle
        
            typedef Part< 3 > exported_class_t;
            typedef ::boost::shared_ptr< Polygon > ( exported_class_t::*AddRectangle_function_type)( double,double,::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Part3_exposer.def( 
                "AddRectangle"
                , AddRectangle_function_type( &::Part< 3 >::AddRectangle )
                , ( bp::arg("sizeX"), bp::arg("sizeY"), bp::arg("origin") ) );
        
        }
        { //::Part< 3 >::AddVesselNetwork
        
            typedef Part< 3 > exported_class_t;
            typedef void ( exported_class_t::*AddVesselNetwork_function_type)( ::boost::shared_ptr< VesselNetwork< 3 > >,bool ) ;
            
            Part3_exposer.def( 
                "AddVesselNetwork"
                , AddVesselNetwork_function_type( &::Part< 3 >::AddVesselNetwork )
                , ( bp::arg("pVesselNetwork"), bp::arg("surface") ) );
        
        }
        { //::Part< 3 >::Create
        
            typedef Part< 3 > exported_class_t;
            typedef ::boost::shared_ptr< Part< 3 > > ( *Create_function_type )(  );
            
            Part3_exposer.def( 
                "Create"
                , Create_function_type( &::Part< 3 >::Create ) );
        
        }
        { //::Part< 3 >::Extrude
        
            typedef Part< 3 > exported_class_t;
            typedef void ( exported_class_t::*Extrude_function_type)( ::boost::shared_ptr< Polygon >,double ) ;
            
            Part3_exposer.def( 
                "Extrude"
                , Extrude_function_type( &::Part< 3 >::Extrude )
                , ( bp::arg("pPolygon"), bp::arg("distance") ) );
        
        }
        { //::Part< 3 >::GetBoundingBox
        
            typedef Part< 3 > exported_class_t;
            typedef ::boost::numeric::ublas::c_vector< double, 6 > ( exported_class_t::*GetBoundingBox_function_type)(  ) ;
            
            Part3_exposer.def( 
                "GetBoundingBox"
                , GetBoundingBox_function_type( &::Part< 3 >::GetBoundingBox ) );
        
        }
        { //::Part< 3 >::GetContainingGridIndices
        
            typedef Part< 3 > exported_class_t;
            typedef ::std::vector< unsigned int > ( exported_class_t::*GetContainingGridIndices_function_type)( unsigned int,unsigned int,unsigned int,double ) ;
            
            Part3_exposer.def( 
                "GetContainingGridIndices"
                , GetContainingGridIndices_function_type( &::Part< 3 >::GetContainingGridIndices )
                , ( bp::arg("num_x"), bp::arg("num_y"), bp::arg("num_z"), bp::arg("spacing") ) );
        
        }
        { //::Part< 3 >::GetFacets
        
            typedef Part< 3 > exported_class_t;
            typedef ::std::vector< boost::shared_ptr<Facet> > ( exported_class_t::*GetFacets_function_type)(  ) ;
            
            Part3_exposer.def( 
                "GetFacets"
                , GetFacets_function_type( &::Part< 3 >::GetFacets ) );
        
        }
        { //::Part< 3 >::GetHoleMarkers
        
            typedef Part< 3 > exported_class_t;
            typedef ::std::vector< boost::numeric::ublas::c_vector<double, 3> > ( exported_class_t::*GetHoleMarkers_function_type)(  ) ;
            
            Part3_exposer.def( 
                "GetHoleMarkers"
                , GetHoleMarkers_function_type( &::Part< 3 >::GetHoleMarkers ) );
        
        }
        { //::Part< 3 >::GetPolygons
        
            typedef Part< 3 > exported_class_t;
            typedef ::std::vector< boost::shared_ptr<Polygon> > ( exported_class_t::*GetPolygons_function_type)(  ) ;
            
            Part3_exposer.def( 
                "GetPolygons"
                , GetPolygons_function_type( &::Part< 3 >::GetPolygons ) );
        
        }
        { //::Part< 3 >::GetSegmentIndices
        
            typedef Part< 3 > exported_class_t;
            typedef ::std::vector< std::pair<unsigned int, unsigned int> > ( exported_class_t::*GetSegmentIndices_function_type)(  ) ;
            
            Part3_exposer.def( 
                "GetSegmentIndices"
                , GetSegmentIndices_function_type( &::Part< 3 >::GetSegmentIndices ) );
        
        }
        { //::Part< 3 >::GetVertexLocations
        
            typedef Part< 3 > exported_class_t;
            typedef ::std::vector< boost::numeric::ublas::c_vector<double, 3> > ( exported_class_t::*GetVertexLocations_function_type)(  ) ;
            
            Part3_exposer.def( 
                "GetVertexLocations"
                , GetVertexLocations_function_type( &::Part< 3 >::GetVertexLocations ) );
        
        }
        { //::Part< 3 >::GetVertices
        
            typedef Part< 3 > exported_class_t;
            typedef ::std::vector< boost::shared_ptr<Vertex> > ( exported_class_t::*GetVertices_function_type)(  ) ;
            
            Part3_exposer.def( 
                "GetVertices"
                , GetVertices_function_type( &::Part< 3 >::GetVertices ) );
        
        }
        { //::Part< 3 >::GetVtk
        
            typedef Part< 3 > exported_class_t;
            typedef ::vtkSmartPointer< vtkPolyData > ( exported_class_t::*GetVtk_function_type)( bool ) ;
            
            Part3_exposer.def( 
                "GetVtk"
                , GetVtk_function_type( &::Part< 3 >::GetVtk )
                , ( bp::arg("update") ) );
        
        }
        { //::Part< 3 >::IsPointInPart
        
            typedef Part< 3 > exported_class_t;
            typedef bool ( exported_class_t::*IsPointInPart_function_type)( ::boost::numeric::ublas::c_vector< double, 3 >,bool ) ;
            
            Part3_exposer.def( 
                "IsPointInPart"
                , IsPointInPart_function_type( &::Part< 3 >::IsPointInPart )
                , ( bp::arg("location"), bp::arg("update") ) );
        
        }
        { //::Part< 3 >::MergeCoincidentVertices
        
            typedef Part< 3 > exported_class_t;
            typedef void ( exported_class_t::*MergeCoincidentVertices_function_type)(  ) ;
            
            Part3_exposer.def( 
                "MergeCoincidentVertices"
                , MergeCoincidentVertices_function_type( &::Part< 3 >::MergeCoincidentVertices ) );
        
        }
        { //::Part< 3 >::Translate
        
            typedef Part< 3 > exported_class_t;
            typedef void ( exported_class_t::*Translate_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Part3_exposer.def( 
                "Translate"
                , Translate_function_type( &::Part< 3 >::Translate )
                , ( bp::arg("vector") ) );
        
        }
        { //::Part< 3 >::Write
        
            typedef Part< 3 > exported_class_t;
            typedef void ( exported_class_t::*Write_function_type)( ::std::string const & ) ;
            
            Part3_exposer.def( 
                "Write"
                , Write_function_type( &::Part< 3 >::Write )
                , ( bp::arg("rFilename") ) );
        
        }
        { //::Part< 3 >::WriteStl
        
            typedef Part< 3 > exported_class_t;
            typedef void ( exported_class_t::*WriteStl_function_type)( ::std::string const & ) ;
            
            Part3_exposer.def( 
                "WriteStl"
                , WriteStl_function_type( &::Part< 3 >::WriteStl )
                , ( bp::arg("rFilename") ) );
        
        }
        Part3_exposer.staticmethod( "Create" );
        bp::register_ptr_to_python< boost::shared_ptr< Part<3> > >();
    }

    { //::Polygon
        typedef bp::class_< Polygon > Polygon_exposer_t;
        Polygon_exposer_t Polygon_exposer = Polygon_exposer_t( "Polygon", bp::init< std::vector< boost::shared_ptr<Vertex> > >(( bp::arg("vertices") )) );
        bp::scope Polygon_scope( Polygon_exposer );
        bp::implicitly_convertible< std::vector< boost::shared_ptr<Vertex> >, Polygon >();
        Polygon_exposer.def( bp::init< boost::shared_ptr< Vertex > >(( bp::arg("pVertex") )) );
        bp::implicitly_convertible< boost::shared_ptr< Vertex >, Polygon >();
        { //::Polygon::AddVertex
        
            typedef void ( ::Polygon::*AddVertex_function_type)( ::boost::shared_ptr< Vertex > ) ;
            
            Polygon_exposer.def( 
                "AddVertex"
                , AddVertex_function_type( &::Polygon::AddVertex )
                , ( bp::arg("pVertex") ) );
        
        }
        { //::Polygon::AddVertices
        
            typedef void ( ::Polygon::*AddVertices_function_type)( ::std::vector< boost::shared_ptr<Vertex> > ) ;
            
            Polygon_exposer.def( 
                "AddVertices"
                , AddVertices_function_type( &::Polygon::AddVertices )
                , ( bp::arg("vertices") ) );
        
        }
        { //::Polygon::ContainsPoint
        
            typedef bool ( ::Polygon::*ContainsPoint_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Polygon_exposer.def( 
                "ContainsPoint"
                , ContainsPoint_function_type( &::Polygon::ContainsPoint )
                , ( bp::arg("location") ) );
        
        }
        { //::Polygon::Create
        
            typedef ::boost::shared_ptr< Polygon > ( *Create_function_type )( ::std::vector< boost::shared_ptr<Vertex> > );
            
            Polygon_exposer.def( 
                "Create"
                , Create_function_type( &::Polygon::Create )
                , ( bp::arg("vertices") ) );
        
        }
        { //::Polygon::Create
        
            typedef ::boost::shared_ptr< Polygon > ( *Create_function_type )( ::boost::shared_ptr< Vertex > );
            
            Polygon_exposer.def( 
                "Create"
                , Create_function_type( &::Polygon::Create )
                , ( bp::arg("pVertex") ) );
        
        }
        { //::Polygon::GetBoundingBox
        
            typedef ::boost::numeric::ublas::c_vector< double, 6 > ( ::Polygon::*GetBoundingBox_function_type)(  ) ;
            
            Polygon_exposer.def( 
                "GetBoundingBox"
                , GetBoundingBox_function_type( &::Polygon::GetBoundingBox ) );
        
        }
        { //::Polygon::GetCentroid
        
            typedef ::boost::numeric::ublas::c_vector< double, 3 > ( ::Polygon::*GetCentroid_function_type)(  ) ;
            
            Polygon_exposer.def( 
                "GetCentroid"
                , GetCentroid_function_type( &::Polygon::GetCentroid ) );
        
        }
        { //::Polygon::GetDistance
        
            typedef double ( ::Polygon::*GetDistance_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Polygon_exposer.def( 
                "GetDistance"
                , GetDistance_function_type( &::Polygon::GetDistance )
                , ( bp::arg("location") ) );
        
        }
        { //::Polygon::GetDistanceToEdges
        
            typedef double ( ::Polygon::*GetDistanceToEdges_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Polygon_exposer.def( 
                "GetDistanceToEdges"
                , GetDistanceToEdges_function_type( &::Polygon::GetDistanceToEdges )
                , ( bp::arg("location") ) );
        
        }
        { //::Polygon::GetNormal
        
            typedef ::boost::numeric::ublas::c_vector< double, 3 > ( ::Polygon::*GetNormal_function_type)(  ) ;
            
            Polygon_exposer.def( 
                "GetNormal"
                , GetNormal_function_type( &::Polygon::GetNormal ) );
        
        }
        { //::Polygon::GetPlane
        
            typedef ::vtkSmartPointer< vtkPlane > ( ::Polygon::*GetPlane_function_type)(  ) ;
            
            Polygon_exposer.def( 
                "GetPlane"
                , GetPlane_function_type( &::Polygon::GetPlane ) );
        
        }
        { //::Polygon::GetVertex
        
            typedef ::boost::shared_ptr< Vertex > ( ::Polygon::*GetVertex_function_type)( unsigned int ) ;
            
            Polygon_exposer.def( 
                "GetVertex"
                , GetVertex_function_type( &::Polygon::GetVertex )
                , ( bp::arg("idx") ) );
        
        }
        { //::Polygon::GetVertices
        
            typedef ::std::vector< boost::shared_ptr<Vertex> > ( ::Polygon::*GetVertices_function_type)(  ) ;
            
            Polygon_exposer.def( 
                "GetVertices"
                , GetVertices_function_type( &::Polygon::GetVertices ) );
        
        }
        { //::Polygon::GetVtkPolygon
        
            typedef ::vtkSmartPointer< vtkPolygon > ( ::Polygon::*GetVtkPolygon_function_type)(  ) ;
            
            Polygon_exposer.def( 
                "GetVtkPolygon"
                , GetVtkPolygon_function_type( &::Polygon::GetVtkPolygon ) );
        
        }
        { //::Polygon::GetVtkVertices
        
            typedef ::std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkIdTypeArray > > ( ::Polygon::*GetVtkVertices_function_type)(  ) ;
            
            Polygon_exposer.def( 
                "GetVtkVertices"
                , GetVtkVertices_function_type( &::Polygon::GetVtkVertices ) );
        
        }
        { //::Polygon::ReplaceVertex
        
            typedef void ( ::Polygon::*ReplaceVertex_function_type)( unsigned int,::boost::shared_ptr< Vertex > ) ;
            
            Polygon_exposer.def( 
                "ReplaceVertex"
                , ReplaceVertex_function_type( &::Polygon::ReplaceVertex )
                , ( bp::arg("idx"), bp::arg("pVertex") ) );
        
        }
        { //::Polygon::RotateAboutAxis
        
            typedef void ( ::Polygon::*RotateAboutAxis_function_type)( ::boost::numeric::ublas::c_vector< double, 3 >,double ) ;
            
            Polygon_exposer.def( 
                "RotateAboutAxis"
                , RotateAboutAxis_function_type( &::Polygon::RotateAboutAxis )
                , ( bp::arg("axis"), bp::arg("angle") ) );
        
        }
        { //::Polygon::Translate
        
            typedef void ( ::Polygon::*Translate_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Polygon_exposer.def( 
                "Translate"
                , Translate_function_type( &::Polygon::Translate )
                , ( bp::arg("translationVector") ) );
        
        }
        Polygon_exposer.staticmethod( "Create" );
        bp::register_ptr_to_python< boost::shared_ptr< Polygon > >();
    }

    { //::Vertex
        typedef bp::class_< Vertex > Vertex_exposer_t;
        Vertex_exposer_t Vertex_exposer = Vertex_exposer_t( "Vertex", bp::init< bp::optional< double, double, double > >(( bp::arg("x")=0., bp::arg("y")=0., bp::arg("z")=0. )) );
        bp::scope Vertex_scope( Vertex_exposer );
        bp::implicitly_convertible< double, Vertex >();
        Vertex_exposer.def( bp::init< boost::numeric::ublas::c_vector< double, 3 > >(( bp::arg("coords") )) );
        bp::implicitly_convertible< boost::numeric::ublas::c_vector< double, 3 >, Vertex >();
        { //::Vertex::Create
        
            typedef ::boost::shared_ptr< Vertex > ( *Create_function_type )( double,double,double );
            
            Vertex_exposer.def( 
                "Create"
                , Create_function_type( &::Vertex::Create )
                , ( bp::arg("x")=0., bp::arg("y")=0., bp::arg("z")=0. ) );
        
        }
        { //::Vertex::Create
        
            typedef ::boost::shared_ptr< Vertex > ( *Create_function_type )( ::boost::numeric::ublas::c_vector< double, 3 > );
            
            Vertex_exposer.def( 
                "Create"
                , Create_function_type( &::Vertex::Create )
                , ( bp::arg("coords") ) );
        
        }
        { //::Vertex::GetIndex
        
            typedef unsigned int ( ::Vertex::*GetIndex_function_type)(  ) ;
            
            Vertex_exposer.def( 
                "GetIndex"
                , GetIndex_function_type( &::Vertex::GetIndex ) );
        
        }
        { //::Vertex::RotateAboutAxis
        
            typedef void ( ::Vertex::*RotateAboutAxis_function_type)( ::boost::numeric::ublas::c_vector< double, 3 >,double ) ;
            
            Vertex_exposer.def( 
                "RotateAboutAxis"
                , RotateAboutAxis_function_type( &::Vertex::RotateAboutAxis )
                , ( bp::arg("axis"), bp::arg("angle") ) );
        
        }
        { //::Vertex::SetIndex
        
            typedef void ( ::Vertex::*SetIndex_function_type)( unsigned int ) ;
            
            Vertex_exposer.def( 
                "SetIndex"
                , SetIndex_function_type( &::Vertex::SetIndex )
                , ( bp::arg("index") ) );
        
        }
        { //::Vertex::Translate
        
            typedef void ( ::Vertex::*Translate_function_type)( ::boost::numeric::ublas::c_vector< double, 3 > ) ;
            
            Vertex_exposer.def( 
                "Translate"
                , Translate_function_type( &::Vertex::Translate )
                , ( bp::arg("translationVector") ) );
        
        }
        Vertex_exposer.staticmethod( "Create" );
        bp::register_ptr_to_python< boost::shared_ptr< Vertex > >();
        bp::implicitly_convertible< boost::shared_ptr< Vertex >, boost::shared_ptr< boost::enable_shared_from_this< Vertex > > >();
        bp::implicitly_convertible< boost::shared_ptr< Vertex >, boost::shared_ptr< ChastePoint< 3 > > >();
    }
}
