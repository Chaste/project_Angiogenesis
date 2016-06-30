/*

 Copyright (c) 2005-2015, University of Oxford.
 All rights reserved.

 University of Oxford means the Chancellor, Masters and Scholars of the
 University of Oxford, having an administrative office at Wellington
 Square, Oxford OX1 2JD, UK.

 This file is part of Chaste.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifdef CHASTE_ANGIOGENESIS_PYTHON
#ifndef BP_CONVERTERS_HPP_
#define BP_CONVERTERS_HPP_

#include <iostream>
#include <sstream>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/stl_iterator.hpp>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION //Ignore numpy depracated warnings
#include <numpy/ndarrayobject.h>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Ignore vtk deprecated warnings
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkObjectBase.h>
#include "Debug.hpp"
//#include <vtkPythonUtil.h>

#include "UblasIncludes.hpp"

using namespace boost::python;

/* Collection of converters to and from Python objects
 */

/******************************
 * To Python converters
 */

/* VTK pointers to Python VTK objects
 */
template<class T>
struct VtkSmartPointerToPython
{
    static PyObject* convert(const vtkSmartPointer<T> &p)
    {
        std::ostringstream oss;
        oss << (void*) p.GetPointer();
        std::string address_str = oss.str();

        boost::python::object obj = import("vtk").attr("vtkObjectBase")(address_str);
        return incref(obj.ptr());
    }
};

//#define VTK_PYTHON_CONVERSION(type) converter::registry::insert(&extract_vtk_wrapped_pointer, type_id<type>());

//vtkObjectBase* GetImageDataPtr(PyObject *  obj)
//{
//    return vtkPythonUtil::GetPointerFromObject(PyObject *  obj, const char *    classname)
//}

//vtkObjectBase*

void* extract_vtk_wrapped_pointer(PyObject* obj)
{
    char thisStr[] = "__this__";
    //first we need to get the __this__ attribute from the Python Object
    if (!PyObject_HasAttrString(obj, thisStr))
    return NULL;

    PyObject* thisAttr = PyObject_GetAttrString(obj, thisStr);
    if (thisAttr == NULL)
    return NULL;

    const char* str = PyString_AsString(thisAttr);
    if(str == 0 || strlen(str) < 1)
    return NULL;

    char hex_address[32], *pEnd;
    const char *_p_ = strstr(str, "_p_vtk");
    if(_p_ == NULL) return NULL;
    const char *class_name = strstr(_p_, "vtk");
    if(class_name == NULL) return NULL;
    strcpy(hex_address, str+1);
    hex_address[_p_-str-1] = '\0';
    long address = strtol(hex_address, &pEnd, 16);

    vtkObjectBase* vtk_object = (vtkObjectBase*)((void*)address);
    if(vtk_object->IsA(class_name))
    {
        vtk_object->Register(NULL);
        return vtk_object;
    }
    return NULL;
}

/* c_vector to numpy array
 */
template<class T>
struct CVectorToNumpyArray
{
    static PyObject* convert(T const& vec)
    {
        npy_intp size = vec.size();
        double * data = size ? const_cast<double *>(&vec[0]) : static_cast<double *>(NULL);
        PyObject * pyObj = PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE, data);
        boost::python::handle<> handle( pyObj );
        boost::python::numeric::array arr( handle );
        return incref(arr.ptr());
    }
};

template<class T>
struct StdVectorDoubleToNumpyArray
{
    static PyObject* convert(T const& vec)
    {
        npy_intp size = vec.size();
        double * data = size ? const_cast<double *>(&vec[0]) : static_cast<double *>(NULL);
        PyObject * pyObj = PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE, data);
        boost::python::handle<> handle( pyObj );
        boost::python::numeric::array arr( handle );
        return incref(arr.ptr());
    }
};

/* STL Iterators to Python Iterators
 */
template <class Container>
class vector_ptr_indexing_suite : public vector_indexing_suite<Container, true, vector_ptr_indexing_suite<Container> >
{
public:

    template <class Class>
    static void extension_def(Class & cl)
    {
        vector_indexing_suite<Container, true, vector_ptr_indexing_suite<Container> >::extension_def(cl);
        cl.def("__iter__", iterator<Container, return_value_policy<copy_non_const_reference> >());
    }
};

/*****************************
 * From Python converters
 */

struct TupleToCVector
{
    template <unsigned DIM>
    TupleToCVector& from_python()
    {
        boost::python::converter::registry::push_back(&TupleToCVector::convertible,
                &TupleToCVector::construct<DIM>,
                boost::python::type_id<c_vector<double, DIM> >());
        return *this;
    }

    static void* convertible(PyObject* obj_ptr)
    {
        if (PyTuple_Check(obj_ptr))
        {
            return obj_ptr;
        }
        return NULL;
    }

    // Do the conversion
    template <unsigned DIM>
    static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        typedef boost::python::converter::rvalue_from_python_storage<c_vector<double, DIM> > storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;
        new (storage) c_vector<double, DIM>;
        c_vector<double, DIM>* vec = (c_vector<double, DIM>*) storage;

        // Populate the vector
        int tuple_size = PyTuple_Size(obj_ptr);
        for(int idx = 0; idx<DIM; idx++)
        {
            if(idx < tuple_size)
            {
                (*vec)[idx] = extract<double>(PyTuple_GetItem(obj_ptr, idx));;
            }
            else
            {
                (*vec)[idx] = 0.0;
            }
        }
        data->convertible = storage;
    }
};

struct NumpyArrayToCVector
{
    template <unsigned DIM>
    NumpyArrayToCVector& from_python()
    {
        boost::python::converter::registry::push_back(&NumpyArrayToCVector::convertible,
                &NumpyArrayToCVector::construct<DIM>,
                boost::python::type_id<c_vector<double, DIM> >());
        return *this;
    }

    // Determine if a c_vector can be generated
    static void* convertible(PyObject* obj_ptr)
    {
        if (PyArray_Check(obj_ptr))
        {
            PyArrayObject* arrayObjPtr = reinterpret_cast<PyArrayObject*>(obj_ptr);
            if(PyArray_NDIM(arrayObjPtr) == 1)
            {
                return obj_ptr;
            }
        }
        return NULL;
    }

    // Do the conversion
    template <unsigned DIM>
    static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        typedef boost::python::converter::rvalue_from_python_storage<c_vector<double, DIM> > storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;
        new (storage) c_vector<double, DIM>;
        c_vector<double, DIM>* vec = (c_vector<double, DIM>*) storage;

        // Populate the vector
        PyArrayObject* arrayObjPtr = reinterpret_cast<PyArrayObject*>(obj_ptr);
        int size_array = PyArray_DIM(arrayObjPtr, 0);
        for(int idx = 0; idx<DIM; idx++)
        {
            if(idx < size_array)
            {
                (*vec)[idx] = extract<double>(PyArray_GETITEM(arrayObjPtr, (const char*)PyArray_GETPTR1(arrayObjPtr, idx)));
            }
            else
            {
                (*vec)[idx] = 0.0;
            }
        }
        data->convertible = storage;
    }
};

/* Python Iterables to C++
 *
 */
struct PythonIterableToStl
{

    template <typename Container>
    PythonIterableToStl&
    from_python()
    {
        boost::python::converter::registry::push_back(
                &PythonIterableToStl::convertible,
                &PythonIterableToStl::construct<Container>,
                boost::python::type_id<Container>());

        // Support chaining.
        return *this;
    }

    // Check if PyObject is iterable.
    static void* convertible(PyObject* object)
    {
        return PyObject_GetIter(object) ? object : NULL;
    }

    template <typename Container>
    static void construct(PyObject* object, boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        namespace python = boost::python;
        python::handle<> handle(python::borrowed(object));

        typedef python::converter::rvalue_from_python_storage<Container> storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;
        typedef python::stl_input_iterator<typename Container::value_type>iterator;

        new (storage) Container(
                iterator(python::object(handle)), // begin
                iterator());// end
        data->convertible = storage;
    }
};

#endif /* BP_CONVERTERS_HPP_ */
#endif // CHASTE_ANGIOGENESIS_PYTHON
