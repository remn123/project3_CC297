#
#	Makefile
#
#	Project 3 - CC297
#

# compiler to use
CC = clang
#CC = gcc

# Flags (-Wall, -Werror nao)
UNUSEDARGS=-Wunused-value
#UNUSEDARGS=-Ounused-arguments

CFLAGS = -ggdb3 -O0 $(UNUSEDARGS) -std=c11 -Wall #-Werror

# nome do executavel
EXE = project3

# lista de arquivos 'Headers', separados por 'espacos'.
HDRS = structs.h settings.h functions.h

# lista de 'Libraries', separados por espacos.
# cada qual com acompanhada pelo prefixo -l
LIBS = -lm

# lista de arquivos fonte, separados por espaco.
SRCS = project3.c functions.c

# lista automatica de arquivos .o
OBJS = $(SRCS:.c=.o)

# target principal (que vai ser executado em default pelo comando 'make')
$(EXE): $(OBJS) $(HDRS) Makefile
		$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)
		
# dependencias
$(OBJS): $(HDRS) Makefile

# apagar tudo!
clean:
		rm -f core $(EXE) *.o
