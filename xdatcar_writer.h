#ifndef XDATCAR_WRITER_H
#define XDATCAR_WRITER_H

#include "poscar_reader.h"
#include "structure.h"

void start_writing_xdatcar(const char *filename, Topology *structure, atomic_site *sites, Atom *atoms, int num_atoms);
void write_config(const char *filename, int* counter_xdatcar, atomic_site *sites, Atom *atoms, int num_atoms );

#endif // 
