#
# Makefile: BSM algorithm for the set covering problem
# Programmed by Takahiro Kan <u123826e@ecs.osaka-u.ac.jp>
# Date: 2016/11/8
#

CC	= gcc
CFLAGS	= -O3
#CFLAGS  = -Wall
#CFLAGS	= -O2 -pg
#CFLAGS	= -O3 -fomit-frame-pointer
#CFLAGS	= -O4 -Wall -g
#CFLAGS	= 
RM	= rm -f
LIBS	= -lm
TARGET	= s
HDR	= 

#SRCS	= SM_with_lb.c cpu_time.c mt19937ar.c LB_SM.c

#SRCS	= u_SM.c cpu_time.c mt19937ar.c LB_class.c
#SRCS	= u_SSM.c cpu_time.c mt19937ar.c LB_class.c
#SRCS	= u_PSM.c cpu_time.c mt19937ar.c LB_class.c
#SRCS	= u_PSSM.c cpu_time.c mt19937ar.c LB_class.c

#SRCS	= SM.c cpu_time.c mt19937ar.c LB_class.c
#SRCS	= SSM.c cpu_time.c mt19937ar.c LB_class.c
#SRCS	= PSM.c cpu_time.c mt19937ar.c LB_class.c
#SRCS	= PSSM.c cpu_time.c mt19937ar.c LB_class.c
SRCS	= PriceS.c cpu_time.c mt19937ar.c LB_class.c

OBJS	= $(SRCS:.c=.o)
PROGS	= $(SRCS:.c=)
CORE	= $(TARGET).core $(EXE).stackdump
EXE	= $(TARGET).exe

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $? -o $(TARGET) $(LIBS)
	$(RM) $(OBJS)

$(OBJS):

clean:
	$(RM) *~ $(PROGS) $(OBJS) $(TARGET) $(CORE) $(EXE)

#
# End of file
#
