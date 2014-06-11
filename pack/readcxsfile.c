#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "readcxsfile.h"
#define VERTICES 0
#define FACE 1
#define ATOMS 2
#define DISTANCE 3
#define INDICES 4


int line_startswith(char * pre, char * line) {
  return strncmp(pre, line, sizeof(char) * strlen(pre));
}
char * get_formula(char * buf) {
  char * s = malloc(sizeof(char)*BUFSIZ);
  strncpy(s,buf,BUFSIZ);
  return s;
}

int readvals(FILE * f, int kind, int count, void * s) {
  char buf[BUFSIZ];
  int n = 0;
  switch (kind)
  {
    case VERTICES: {
      float * v = s;
      for(int i = 0; i < count; i++) {
        fgets(buf,BUFSIZ,f);
        //Treat arr as a pointer to an Nx3 array
        float pts[3];
        sscanf(buf,"%e %e %e\n",&pts[0],&pts[1],&pts[2]);
        memcpy(&v[3*i], pts, sizeof(float) * 3);
        n++;
      }
      return n;
      break;
    }
    case INDICES: {
      int * x = s;
      for(int i = 0; i < count; i++) {
        fgets(buf,BUFSIZ,f);
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
        fgets(buf, BUFSIZ, f);
        sscanf(buf,"%d %*s",&face[i]);
        n++;
      }
      return n;
      break;
    }
    case ATOMS: {
      char * atom = s;
      for(int i = 0; i < count; i++) {
        fgets(buf, BUFSIZ, f);
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
        fgets(buf, BUFSIZ, f);
        sscanf(buf,"%e",&d[i]);
        n++;
      }
      return n;
      break;
    }
  }
  return count;
}

CXS_DATA * readcxsfile(char * fname) {
  //File processing vars
  FILE * inputFile;
  char buf[BUFSIZ];
  int n = 0;
  //RETURN VALUES
  float * devals, * divals, * vertices;
  int * indices;
  char * external, *internal;
  // Temporary variables
  char * atoms;
  int * de_face_atoms;
  int * di_face_atoms;
  int * atoms_outside;
  int * atoms_inside;
  char * formula;
  inputFile = fopen(fname, "r");
  int count = 0;
  int nfaces;
  int nvertices;

  while (!feof(inputFile)) {
    fgets(buf, BUFSIZ, inputFile);
    int r = 0;
    if (line_startswith("begin unit_cell",buf) == 0) {
      sscanf(buf,"begin unit_cell %d\n",&count);
      atoms = malloc(sizeof(char) * count * 2);
      r =readvals(inputFile, ATOMS, count, (void *) atoms);
    }
    else if( line_startswith("begin vertices",buf) == 0){
      sscanf(buf,"begin vertices %d\n",&count);
      vertices = malloc(sizeof(float) * count * 3);
      r = readvals(inputFile, VERTICES, count, (void *) vertices);
      nvertices = count;
    }

    else if (line_startswith("begin indices",buf) == 0) {
      sscanf(buf,"begin indices %d\n",&count);
      indices = malloc(sizeof(int) * count * 3);
      r = readvals(inputFile, INDICES, count, (void *) indices);


    }
    else if (line_startswith("begin atoms_inside",buf) == 0) {
      sscanf(buf,"begin atoms_inside_surface %d\n",&count);
      atoms_inside = malloc(sizeof(int) * count);
      r = readvals(inputFile, FACE, count, (void *) atoms_inside);
    }

    else if (line_startswith("begin atoms_outside",buf) == 0) {
      sscanf(buf,"begin atoms_outside_surface %d\n",&count);
      atoms_outside = malloc(sizeof(int) * count);
      r = readvals(inputFile, FACE, count, (void *) atoms_outside);
    }

    else if (line_startswith("begin d_i_face_atoms",buf) == 0) {
      sscanf(buf,"begin d_i_face_atoms %d\n",&count);
      di_face_atoms = malloc(sizeof(*di_face_atoms) * count);
      r = readvals(inputFile, FACE, count, (void *) di_face_atoms);
      internal = malloc(sizeof(char) * count * 2);
    }

    else if (line_startswith("begin d_e_face_atoms",buf) == 0) {
      sscanf(buf,"begin d_e_face_atoms %d\n",&count);
      de_face_atoms = malloc(sizeof(*de_face_atoms) * count);
      r = readvals(inputFile, FACE, count, (void *) de_face_atoms);
      external = malloc(sizeof(char) * count * 2);
      nfaces = count;
    }

    else if (line_startswith("begin d_i ",buf) == 0) {
      sscanf(buf,"begin d_i %d\n",&count);
      divals = malloc(sizeof(float) * count);
      r = readvals(inputFile, DISTANCE , count, (void *) divals);
    }

    else if (line_startswith("begin d_e ",buf) == 0) {
      sscanf(buf,"begin d_e %d\n",&count);
      devals = malloc(sizeof(float) * count);
      r = readvals(inputFile, DISTANCE , count, (void *) devals);
    }
    else if(line_startswith("   formula = ",buf) == 0) {
      formula = get_formula(buf);
    }
    if(r > 0 && r < count) {
      printf("read less values than count!!\n");
      printf("r = %d, count = %d\n",r,count);
      return NULL;
    }
    n++;
  }
  if(nfaces == 0 || nvertices == 0) {
    printf("nfaces or nvertices == 0");
    return NULL;
  }
  //NOW THAT WE HAVE DATA IN ARRAYS, MANIPULATE AND RETURN;
  for(int i = 0; i < nfaces; i++) {
    memcpy(&external[2*i],&atoms[(atoms_outside[de_face_atoms[i] - 1] -1) * 2 ],2*sizeof(char));
    memcpy(&internal[2*i],&atoms[(atoms_inside[di_face_atoms[i] - 1] -1) * 2 ],2*sizeof(char));
  }
  free(atoms_outside);
  free(atoms_inside);
  free(de_face_atoms);
  free(di_face_atoms);
  free(atoms);
  //divals, devals, internal, external, vertices, indices
  // turn these into lists
  CXS_DATA * rslt = malloc(sizeof(CXS_DATA));
  rslt->internal = internal;
  rslt->external = external;
  rslt->divals = divals;
  rslt->devals = devals;
  rslt->vertices = vertices;
  rslt->indices = indices;
  rslt->nfaces = nfaces;
  rslt->nvertices =nvertices;
  rslt->formula = formula;
  return rslt;
  //NEED TO REMEMBER TO free memory, and decrement ref counts after packed into tuple
}
/*
int main(int argc, char ** argv) {
  clock_t start, end;
  start= clock();
  CXS_DATA * cxs = readcxsfile("test.cxs");
  float * vertices = cxs->vertices;
  int * indices = cxs->indices;
  char * external = cxs->external;
  printf("%d vertices, %d faces\n",cxs->nvertices,cxs->nfaces);
  end = clock();
  double time_spent = (double)(end - start) /CLOCKS_PER_SEC;
  printf("Time spent: %e seconds\n",time_spent);
  return 0;
} */
