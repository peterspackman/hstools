#ifndef HEADER_INCLUDED
#define HEADER_INCLUDED
#ifdef MAIN_FILE
char * error_string;
#else
extern char * error_string;
#endif
#endif
PyObject * readcxsfile(char * fname);
void cross3D(double v1[], double v2[], double r[]);
double normal3D(double v[]);
