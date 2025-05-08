CC=gcc
CFLAGS= -O4
LDFLAGS=
VERSION=54

all: plantri plantri_filter

plantri: plantri.c
	${CC} -o plantri ${CFLAGS} plantri.c ${LDFLAGS}

plantri_filter: plantri.c filter.c
	${CC} -o plantri_filter ${CFLAGS} '-DPLUGIN="filter.c"' plantri.c ${LDFLAGS}

