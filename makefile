CC=gcc
CFLAGS= -O4
LDFLAGS=
VERSION=54

all: plantri plantri_hamiltonian plantri_longpath

plantri: plantri.c
	${CC} -o plantri ${CFLAGS} plantri_modified.c ${LDFLAGS}

plantri_hamiltonian: plantri_modified.c hamiltonian_filter.c
	${CC} -o plantri_hamiltonian ${CFLAGS} '-DPLUGIN="hamiltonian_filter.c"' plantri_modified.c ${LDFLAGS}

plantri_longpath: plantri.c long_path_filter.c
	${CC} -o plantri_longpath ${CFLAGS} '-DPLUGIN="long_path_filter.c"' plantri.c ${LDFLAGS}
