typedef struct
{
  float * devals, *divals; //arrays of nvertices floats
  float *vertices; //array of nvertices *3 floats, with a stride of 3
  int * indices; //array of nfaces indices * 3 ints, with a stride of 3
  char * internal, * external;// array of nfaces *2 chars with a stride of 2
  int nfaces, nvertices;
  char * formula;

} CXS_DATA;

CXS_DATA * readcxsfile(char * fname);
