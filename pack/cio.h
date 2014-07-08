#ifndef HEADER_INCLUDED
#define HEADER_INCLUDED
#ifdef MAIN_FILE
char * error_string;
#else
extern char * error_string;
#endif
#endif

typedef struct {
    float * devals, *divals; //arrays of nvertices floats
    float *vertices; //array of nvertices *3 floats, with a stride of 3
    int * indices; //array of nfaces indices * 3 ints, with a stride of 3
    char * internal, * external;// array of nfaces *2 chars with a stride of 2
    int nfaces, nvertices;
    char * formula;

} CXS_DATA;

CXS_DATA * readcxsfile(char * fname);
void cross3D(double v1[], double v2[], double r[]);
double normal3D(double v[]);
