CC=gcc
CFLAGS= -O4
LDFLAGS=
VERSION=1

all: plantri plantri_filter stellation
stellation: all_ham_cycles all_longest_paths

plantri: plantri.c
	${CC} -o plantri ${CFLAGS} plantri.c ${LDFLAGS}

plantri_filter: plantri.c filter.c
	${CC} -o plantri_filter ${CFLAGS} '-DPLUGIN="filter.c"' plantri.c ${LDFLAGS}

all_ham_cycles: stellations/bin/all_ham_cycles.c
	${CC} -o stellations/bin/all_ham_cycles stellations/bin/all_ham_cycles.c ${LDFLAGS}
all_longest_paths: stellations/bin/all_longest_paths.c
	${CC} -o stellations/bin/all_longest_paths stellations/bin/all_longest_paths.c ${LDFLAGS}