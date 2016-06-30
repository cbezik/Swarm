CMPLR = mpic++
CFLGS = -O3 -w -g -std=c++11

EXECU = ./swarm

SRCES = mueller.cpp\
		swarm.cpp \
		main.cpp \

OBJCS := $(addsuffix .o, $(basename $(SRCES)))

all: ${EXECU}

%.o : %.cpp
	${CMPLR} ${CFLGS} -c $< -o $@

%.o : %.c
	${CMPLR} ${CFLGS} -c $< -o $@

$(EXECU): ${OBJCS}
	${CMPLR} ${CFLGS} ${OBJCS} -o $@

clean:
	rm -f ${OBJCS} *.out swarm
