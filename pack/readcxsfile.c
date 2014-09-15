#include <Python.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "cio.h"
#define PY_ARRAY_UNIQUE_SYMBOL cio_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>


#define VERTICES 0
#define FACE 1
#define ATOMS 2
#define DISTANCE 3
#define INDICES 4
#define MOMENTS 5
#define COEFFICIENTS 6
#define INVARIANTS 7


#ifndef __APPLE__
  #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#pragma GCC diagnostic ignored "-Wuninitialized"

int line_startswith(char * pre, char * line)
{
    return strncmp(pre, line, sizeof(char) * strlen(pre));
}

char * get_formula(char * buf)
{
    char * s = calloc(BUFSIZ, sizeof(char));
    if(s == NULL) return NULL;
    strncpy(s,buf,BUFSIZ-1);
    return s;
}

int readvals(FILE * f, int kind, int count, void * s)
{
    char buf[BUFSIZ];
    char * res = NULL;
    int n = 0;
    switch (kind) {
    case VERTICES: {
        float * v = s;
        for(int i = 0; i < count; i++) {
            res = fgets(buf,BUFSIZ,f);
            //Treat arr as a pointer to an Nx3 array
            float pts[3];
            sscanf(buf,"%e %e %e\n",&pts[0],&pts[1],&pts[2]);
            memcpy(&v[3*i], pts, sizeof(float) * 3);
            n++;
        }
        return n;
        break;
    }
    case MOMENTS: {
        float * v = s;
        for(int i = 0; i < count; i++) {
            res = fgets(buf,BUFSIZ,f);
            sscanf(buf,"%*d %e",&v[i]);
            n++;
        }
        return n;
        break;
    }
    case INDICES: {
        int * x = s;
        for(int i = 0; i < count; i++) {
            res = fgets(buf,BUFSIZ,f);
            int ind[3];
            sscanf(buf,"%d %d %d\n",&ind[0],&ind[1],&ind[2]);
            memcpy(&x[3*i], ind, sizeof(int) *3);
            n++;
        }
        return n;
        break;
    }
    case FACE: {
        int * face = s;
        for(int i = 0; i < count; i++) {
            res = fgets(buf, BUFSIZ, f);
            sscanf(buf,"%d %*s",&face[i]);
            n++;
        }
        return n;
        break;
    }
    case ATOMS: {
        char * atom = s;
        for(int i = 0; i < count; i++) {
            res = fgets(buf, BUFSIZ, f);
            char tmp[3] = {'\0','\0','\0'};
            sscanf(buf,"%*s %2s ",(char *) &tmp);
            memcpy(&atom[2*i],tmp,2*sizeof(char));
            n++;
        }
        return n;
        break;
    }
    case DISTANCE: {
        float * d = s;
        for(int i = 0; i < count; i++) {
            res = fgets(buf, BUFSIZ, f);
            sscanf(buf,"%e",&d[i]);
            n++;
        }
        return n;
        break;
    }
    case COEFFICIENTS: {
        float * d = s;
        for(int i = 0; i < count; i++) {
            res = fgets(buf, BUFSIZ, f);
            float pts[2];
            sscanf(buf,"%e %e\n",&pts[0],&pts[1]);
            memcpy(&d[2*i], pts, sizeof(float) * 3);
            n++;
        }
        return n;
        break;
    }
    case INVARIANTS: {
        float * d = s;
        for(int i = 0; i < count; i++) {
            res = fgets(buf, BUFSIZ, f);
            sscanf(buf,"%e \n",&d[i]);
            n++;
        }
    }
    if(res == NULL && feof(f) == 0) {
        fprintf(stderr,"problem reading file ");
        perror("error:");
    }
    }
    return count;
}





PyObject * readcxsfile(char * fname)
{
    //File processing vars
    FILE * inputFile = fopen(fname, "r");
    if(inputFile == NULL) {
        fprintf(stderr, "Error opening file:  %s",fname);
        perror("");
        goto FAIL;
    }
    char buf[BUFSIZ];
    char * res = NULL;
    int n = 0;
    //RETURN VALUES
    float * devals = NULL, * divals = NULL, * vertices = NULL;
    float * dnorm_vals = NULL, * dnorm_evals = NULL, * dnorm_ivals = NULL;
    int * indices = NULL;
    char * external = NULL, *internal = NULL;
    // Temporary variables
    char * atoms = NULL, * formula = NULL;
    int * de_face_atoms = NULL, * di_face_atoms = NULL;
    int * atoms_outside = NULL, * atoms_inside = NULL;
    float *coefficients = NULL;
    float *invariants = NULL;
    int count = 0;
    int nfaces = 0;
    int nvertices = 0;
    int ncoefficients = 0, ninvariants = 0;

    while (!feof(inputFile)) {
        int r = 0;
        res = fgets(buf, BUFSIZ, inputFile);
        // ERRORS READING FILE
        if (res == NULL && feof(inputFile) == 0) {
            fprintf(stderr,"error reading file %s:", fname);
            perror("error:");
            goto FAIL;
        }

        if (line_startswith("begin unit_cell",buf) == 0) {
            sscanf(buf,"begin unit_cell %d\n",&count);
            atoms = malloc(sizeof(char) * count * 2);
            if(atoms == NULL) goto FAIL;
            r = readvals(inputFile, ATOMS, count, (void *) atoms);
        }
        else if(line_startswith("   formula = ",buf) == 0) {
            formula = get_formula(buf);
            if(formula == NULL) goto FAIL;
        }
        // SURFACE
        else if( line_startswith("begin vertices",buf) == 0) {
            sscanf(buf,"begin vertices %d\n",&count);
            vertices = malloc(sizeof(float) * count * 3);
            if(vertices == NULL) goto FAIL;
            r = readvals(inputFile, VERTICES, count, (void *) vertices);
            nvertices = count;
        }
        else if (line_startswith("begin indices",buf) == 0) {
            sscanf(buf,"begin indices %d\n",&count);
            indices = malloc(sizeof(int) * count * 3);
            if(indices == NULL) goto FAIL;
            r = readvals(inputFile, INDICES, count, (void *) indices);
        }
        // EXTERNAL AND INTERNAL ATOM STUFF

        else if (line_startswith("begin atoms_inside",buf) == 0) {
            sscanf(buf,"begin atoms_inside_surface %d\n",&count);
            atoms_inside = malloc(sizeof(int) * count);
            if(atoms_inside == NULL) goto FAIL;
            r = readvals(inputFile, FACE, count, (void *) atoms_inside);
        }
        else if (line_startswith("begin atoms_outside",buf) == 0) {
            sscanf(buf,"begin atoms_outside_surface %d\n",&count);
            atoms_outside = malloc(sizeof(int) * count);
            if(atoms_outside == NULL) goto FAIL;
            r = readvals(inputFile, FACE, count, (void *) atoms_outside);
        }
        else if (line_startswith("begin d_i_face_atoms",buf) == 0) {
            sscanf(buf,"begin d_i_face_atoms %d\n",&count);
            di_face_atoms = malloc(sizeof(*di_face_atoms) * count);
            if(di_face_atoms == NULL) goto FAIL;
            r = readvals(inputFile, FACE, count, (void *) di_face_atoms);
            internal = malloc(sizeof(char) * count * 2);
            if(internal == NULL) goto FAIL;
        }
        else if (line_startswith("begin d_e_face_atoms",buf) == 0) {
            sscanf(buf,"begin d_e_face_atoms %d\n",&count);
            de_face_atoms = malloc(sizeof(*de_face_atoms) * count);
            if(de_face_atoms == NULL) goto FAIL;
            r = readvals(inputFile, FACE, count, (void *) de_face_atoms);
            external = malloc(sizeof(char) * count * 2);
            if(external == NULL) goto FAIL;
            nfaces = count;
        }
        // SURFACE PROPERTIES

        else if (line_startswith("begin d_i ",buf) == 0) {
            sscanf(buf,"begin d_i %d\n",&count);
            divals = malloc(sizeof(float) * count);
            if(divals == NULL) goto FAIL;
            r = readvals(inputFile, DISTANCE , count, (void *) divals);
        }
        else if (line_startswith("begin d_e ",buf) == 0) {
            sscanf(buf,"begin d_e %d\n",&count);
            devals = malloc(sizeof(float) * count);
            if(devals == NULL) goto FAIL;
            r = readvals(inputFile, DISTANCE , count, (void *) devals);
        }
        else if (line_startswith("begin coefficients", buf) == 0) {
            sscanf(buf, "begin coefficients %d\n", &count);
            coefficients = malloc(sizeof(*coefficients) * count * 2);
            if(coefficients == NULL) goto FAIL;
            r = readvals(inputFile, COEFFICIENTS, count, (void *) coefficients);
            ncoefficients = count;
        }
        else if (line_startswith("begin invariants", buf) == 0) {
            sscanf(buf, "begin invariants %d\n", &count);
            invariants = malloc(sizeof(*invariants) * count);
            if(invariants == NULL) goto FAIL;
            r = readvals(inputFile, INVARIANTS, count, (void *) invariants);
            ninvariants = count;
        }

        if(r > 0 && r < count) {
            error_string = "CRITICAL: read less values than count";
            goto FAIL;
        }

        n++;
    }
    if(nfaces == 0 || nvertices == 0) {
        error_string = "CRITICAL: nfaces or nvertices == 0";
        goto FAIL;
    }

    fclose(inputFile);
    //NOW THAT WE HAVE DATA IN ARRAYS, MANIPULATE AND RETURN;
    if(external == NULL || internal == NULL
            || atoms == NULL || atoms_outside == NULL
            || atoms_inside == NULL || de_face_atoms == NULL
            || di_face_atoms == NULL || divals == NULL
            || devals == NULL || invariants == NULL
            || coefficients == NULL) {
        goto FAIL;
    }
    // DEREFERENCE ALL THE ATOMS TO THEIR CHEMICAL SYMBOLS
    for(int i = 0; i < nfaces; i++) {
        int ai_index = di_face_atoms[i] - 1;
        int ao_index = de_face_atoms[i] - 1;
        if (ao_index < 0 || ai_index < 0) {
            error_string = "critical: indexing below 0 (face_atoms)";
            goto FAIL;
        }
        int a1_index = atoms_outside[ao_index] - 1;
        int a2_index = atoms_inside[ai_index] -1;
        if(a1_index < 0 || a2_index < 0) {
            error_string = "critical: indexing below 0 (atoms outside)";
            goto FAIL;
        }
        memcpy(&external[2*i],
               &atoms[(atoms_outside[de_face_atoms[i] - 1] -1) * 2 ],
               2*sizeof(char));
        memcpy(&internal[2*i],
               &atoms[(atoms_inside[di_face_atoms[i] - 1] -1) * 2 ],
               2*sizeof(char));
    }
    // free unused stuff
    free(atoms_inside);
    free(atoms_outside);
    free(de_face_atoms);
    free(di_face_atoms);
    free(atoms);

    //build result and then return it

    // PYTHON STUFF
    PyObject * formulaobj =  PyString_FromString(formula);
    //dimensions etc
    npy_intp didims[1] = {nvertices};
    npy_intp vdims[2] = {nvertices, 3};
    npy_intp idims[2] = {nfaces, 3};
    npy_intp cdims[2] = {ncoefficients, 2};
    npy_intp exdims[1] = {nfaces};
    npy_intp stride[1] = {2*sizeof(char)};
    npy_intp invdims[1] = {ninvariants};

    //Only way to create an array of c strings with length 2 like we have
    PyArray_Descr * desc = PyArray_DescrNewFromType(NPY_STRING);
    desc->elsize = 2;
    const int FLAGS = NPY_CARRAY | NPY_OWNDATA;

    //Construct our numpy arrays
    PyObject * divalsobj = PyArray_SimpleNewFromData(1, didims, NPY_FLOAT, divals);
    PyObject * devalsobj = PyArray_SimpleNewFromData(1, didims, NPY_FLOAT, devals);
    PyObject * verticesobj = PyArray_SimpleNewFromData(2, vdims, NPY_FLOAT, vertices);
    PyObject * indicesobj = PyArray_SimpleNewFromData(2, idims, NPY_INT, indices);
    PyObject * coefficientsobj = PyArray_SimpleNewFromData(2, cdims, NPY_FLOAT, coefficients);
    PyObject * invariantsobj = PyArray_SimpleNewFromData(1, invdims, NPY_FLOAT, invariants);

    PyObject * internalobj = PyArray_NewFromDescr(&PyArray_Type, desc,
                          1, exdims, stride, internal, FLAGS, NULL);
    PyObject * externalobj = PyArray_NewFromDescr(&PyArray_Type, desc,
                          1, exdims, stride, external, FLAGS, NULL);

    //UPDATE THE FLAGS ON SIMPLE ARRAYS
    PyArray_UpdateFlags((PyArrayObject * ) verticesobj, FLAGS);
    PyArray_UpdateFlags((PyArrayObject * ) indicesobj, FLAGS);
    PyArray_UpdateFlags((PyArrayObject * ) divalsobj, FLAGS);
    PyArray_UpdateFlags((PyArrayObject * ) devalsobj, FLAGS);
    PyArray_UpdateFlags((PyArrayObject * ) coefficientsobj, FLAGS);

    PyObject * h = Py_BuildValue("OO", coefficientsobj, invariantsobj);
    PyObject * x = Py_BuildValue("OOOOO",formulaobj, verticesobj, indicesobj,
                                 internalobj, externalobj);
    PyObject * rslt = Py_BuildValue("OOOO", divalsobj, devalsobj, x, h);

    return rslt;

    //failure point, free memory and return null
FAIL:
    free(atoms_inside);
    free(atoms_outside);
    free(de_face_atoms);
    free(di_face_atoms);
    free(atoms);
    free(internal);
    free(external);

    free(coefficients);
    free(invariants);


    free(divals);
    free(devals);
    free(dnorm_evals);
    free(dnorm_ivals);
    free(dnorm_vals);

    free(formula);
    free(vertices);
    free(indices);
    error_string = "Problem reading file, not enough data?";
    printf("Failure in readcxsfile c module");
    return NULL;
}
#ifndef __APPLE__
  #pragma GCC diagnostic pop
#endif
